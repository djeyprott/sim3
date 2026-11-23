"""Builds sim3's world data file from public datasets.

Runs on GitHub's machines (see .github/workflows/data.yml), which can reach
the data hosts. Output: one binary file with terrain, sea and ice for every
map cell on 1 January 10,000 BCE and today, plus a JSON report.

Inputs (downloaded by the workflow):
  grid.bin       cell centres and areas from `gridgen 8` (f64 lat, lon, m2)
  surface.tif    ETOPO 2022 60 arc-second ice-surface elevation (NOAA NCEI)
  bed.tif        ETOPO 2022 60 arc-second bedrock elevation
  ice_0.nc, ice_11.5.nc, ice_12.nc
                 ICE-6G_C (VM5a) 10 arc-minute topography and ice mask
                 (Peltier, Argus & Drummond 2015), via PMIP4
  chelsa_tas_MM.tif (MM = 01..12)
                 CHELSA V2.1 monthly mean near-surface air temperature,
                 1981-2010, 30 arc-seconds (Karger et al. 2017; v2.1 2021)
  orbit91        Berger & Loutre (1991) orbital parameters, 1 kyr steps
  co2_composite.txt
                 Antarctic ice-core CO2 composite (Bereiter et al. 2015)
  co2_mlo.txt    Mauna Loa annual mean CO2 (NOAA GML)

Method
------
Today: each cell gets the latitude-weighted mean, minimum and maximum of the
~230 ETOPO points that fall in it, the share of those points above sea level,
and the mean ice thickness (ice surface minus bedrock; for floating ice
shelves, the thickness at which the ice floats).

10,000 BCE (= 11,950 years before 1950): today's ETOPO terrain plus the
ICE-6G change in topography between then and now, the standard way to put a
coarse reconstruction on detailed modern terrain. ICE-6G is interpolated
linearly between its 12,000 and 11,500 year slices.

Ice thickness in 10,000 BCE is not in the files we can reach, so it is
derived from isostasy (the crust floating on the mantle): a load of mass M
per square metre pushes the ground down by M / rho_mantle. Then, the ground
carried grounded ice of thickness H (load rho_ice * H). Today it carries
either grounded ice H0 (load rho_ice * H0) or, under open sea and floating
ice shelves, the water column of depth d (load rho_seawater * d; floating
ice weighs exactly what it displaces). With surface s (then, relative to then
sea level), today's bedrock b0 and the global sea-level drop dSL:

    s = b0 + dSL - (rho_ice * H - load_now) / rho_mantle + H
    =>  H = (s - b0 - dSL - load_now / rho_mantle) / (1 - rho_ice / rho_mantle) This assumes the ground had fully adjusted to
the ice; in reality it lags, so H may be somewhat underestimated. The report
checks the result: the extra ice volume should match the sea-level drop.
"""

import json
import sys

import netCDF4
import numpy as np
import rasterio
from rasterio.windows import Window
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

# Physical constants, each with its source.
RHO_ICE = 917.0  # kg/m3, glacier ice (Cuffey & Paterson 2010, The Physics of Glaciers, 4th ed.)
RHO_SEAWATER = 1028.0  # kg/m3, mean surface seawater
RHO_MANTLE = 3300.0  # kg/m3, upper mantle (Turcotte & Schubert 2014, Geodynamics, 3rd ed.)
OCEAN_AREA_M2 = 3.618e14  # present ocean area (Eakins & Sharman 2010, NOAA NGDC)

# 1 January 10,000 BCE is 11,950 years before 1950 (BP), between the ICE-6G
# slices at 12,000 and 11,500 BP.
YEARS_BP = 11_950
W_115 = (12_000 - YEARS_BP) / 500.0  # weight on the 11.5 ka slice

FORMAT_VERSION = 2
CHUNK_ROWS = 270


def unit_vectors(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    c = np.cos(lat)
    return np.stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)], axis=-1)


def load_grid(path):
    g = np.fromfile(path, dtype="<f8").reshape(-1, 3)
    return g[:, 0], g[:, 1], g[:, 2]


def ice6g_interpolator(path_a, path_b, w_b, var):
    """Bilinear interpolator over (lat, lon) for (1 - w_b) * a + w_b * b."""
    with netCDF4.Dataset(path_a) as a, netCDF4.Dataset(path_b) as b:
        lat = np.array(a["lat"][:], dtype=np.float64)
        lon = np.array(a["lon"][:], dtype=np.float64)
        field = (1 - w_b) * np.array(a[var][:], dtype=np.float64) + w_b * np.array(b[var][:], dtype=np.float64)
    return periodic_interpolator(lat, lon, field)


def ice6g_field(path, var):
    with netCDF4.Dataset(path) as d:
        return (
            np.array(d["lat"][:], dtype=np.float64),
            np.array(d["lon"][:], dtype=np.float64),
            np.array(d[var][:], dtype=np.float64),
        )


def periodic_interpolator(lat, lon, field):
    if lat[0] > lat[-1]:
        lat, field = lat[::-1], field[::-1]
    # Wrap longitude so cells near 0/360 interpolate across the seam.
    lon_ext = np.concatenate([lon[-1:] - 360.0, lon, lon[:1] + 360.0])
    f_ext = np.concatenate([field[:, -1:], field, field[:, :1]], axis=1)
    return RegularGridInterpolator((lat, lon_ext), f_ext, bounds_error=False, fill_value=None)


def sample(interp, lat, lon):
    return interp(np.stack([lat, np.mod(lon, 360.0)], axis=-1))


def main(grid_path, surface_path, bed_path, ice0, ice115, ice12, chelsa_pattern, orbit_path, co2_path, mlo_path,
         out_path, report_path):
    lat_c, lon_c, area = load_grid(grid_path)
    n = lat_c.size
    print(f"{n} cells", flush=True)
    tree = cKDTree(unit_vectors(lat_c, lon_c))

    # ICE-6G change in topography (then minus now) and ice mask then, at
    # each cell centre.
    lat6, lon6, topo0 = ice6g_field(ice0, "Topo")
    _, _, topo12 = ice6g_field(ice12, "Topo")
    _, _, topo115 = ice6g_field(ice115, "Topo")
    delta_grid = (1 - W_115) * topo12 + W_115 * topo115 - topo0
    delta = sample(periodic_interpolator(lat6, lon6, delta_grid), lat_c, lon_c)
    ice_frac_then = np.clip(sample(ice6g_interpolator(ice12, ice115, W_115, "sftgif"), lat_c, lon_c) / 100.0, 0, 1)

    # Global sea-level drop: in the deep open ocean far from ice, the
    # change in topography is just the change in sea level.
    lat_g, lon_g = np.meshgrid(lat6, lon6, indexing="ij")
    deep = (topo0 < -4000) & (np.abs(lat_g) < 40)
    dsl = float(np.median(delta_grid[deep]))
    print(f"global sea level then: {-dsl:.1f} m relative to today", flush=True)

    # Aggregate ETOPO points into cells.
    sum_w = np.zeros(n)
    sum_s = np.zeros(n)
    sum_b = np.zeros(n)
    sum_ice = np.zeros(n)
    land_now = np.zeros(n)
    land_then = np.zeros(n)
    s_min = np.full(n, np.inf)
    s_max = np.full(n, -np.inf)

    with rasterio.open(surface_path) as sf, rasterio.open(bed_path) as bf:
        rows, cols = sf.height, sf.width
        t = sf.transform
        lons = t.c + t.a * (np.arange(cols) + 0.5)
        for r0 in range(0, rows, CHUNK_ROWS):
            nr = min(CHUNK_ROWS, rows - r0)
            win = Window(0, r0, cols, nr)
            s = sf.read(1, window=win).astype(np.float64).ravel()
            b = bf.read(1, window=win).astype(np.float64).ravel()
            lats = t.f + t.e * (np.arange(r0, r0 + nr) + 0.5)
            lat2 = np.repeat(lats, cols)
            lon2 = np.tile(lons, nr)
            _, idx = tree.query(unit_vectors(lat2, lon2), workers=-1)
            w = np.cos(np.radians(lat2))  # pixel area is proportional to cos(latitude)

            # Ice thickness today: surface minus bed, but floating ice is only
            # as thick as flotation allows.
            h = np.maximum(s - b, 0.0)
            floating_h = np.where(b < 0, np.maximum(s, 0) * RHO_SEAWATER / (RHO_SEAWATER - RHO_ICE), np.inf)
            h = np.minimum(h, floating_h)

            sum_w += np.bincount(idx, w, n)
            sum_s += np.bincount(idx, w * s, n)
            sum_b += np.bincount(idx, w * b, n)
            sum_ice += np.bincount(idx, w * h, n)
            land_now += np.bincount(idx, w * (s > 0), n)
            land_then += np.bincount(idx, w * ((s + delta[idx]) > 0), n)
            np.minimum.at(s_min, idx, s)
            np.maximum.at(s_max, idx, s)
            print(f"rows {r0 + nr}/{rows}", flush=True)

    assert (sum_w > 0).all(), "every cell should contain ETOPO points"
    s_now = sum_s / sum_w
    b_now = sum_b / sum_w
    ice_now = sum_ice / sum_w
    s_then = s_now + delta

    # Load on the ground today (kg/m2): grounded ice, or else the water
    # column (open sea, or ice floating on it).
    depth_now = np.clip(-b_now, 0, None)
    grounded_now = RHO_ICE * ice_now > RHO_SEAWATER * depth_now
    load_now = np.where(grounded_now, RHO_ICE * ice_now, RHO_SEAWATER * depth_now)
    r = RHO_ICE / RHO_MANTLE
    h_iso = (s_then - b_now - dsl - load_now / RHO_MANTLE) / (1 - r)
    # Only where ICE-6G has ice: elsewhere small differences between the
    # datasets would otherwise turn into phantom ice.
    ice_then = np.where(ice_frac_then > 0, np.clip(h_iso, 0, None), 0.0)
    extra_volume = float(np.sum((ice_then - ice_now).clip(0) * area))
    sle = extra_volume * RHO_ICE / RHO_SEAWATER / OCEAN_AREA_M2
    land_share_now = float(np.sum(land_now) / np.sum(sum_w))
    report = {
        "cells": n,
        "years_bp": YEARS_BP,
        "sea_level_then_m": -dsl,
        "extra_ice_volume_km3": extra_volume / 1e9,
        "extra_ice_sea_level_equivalent_m": sle,
        "sea_level_check_ratio": sle / dsl if dsl > 0 else None,
        "land_share_today": land_share_now,
        "land_share_then": float(np.sum(land_then) / np.sum(sum_w)),
        "max_ice_then_m": float(ice_then.max()),
    }

    def i16(x):
        return np.clip(np.rint(x), -32768, 32767).astype("<i2")

    def frac_u8(num):
        return np.clip(np.rint(num / sum_w * 255), 0, 255).astype("u1")

    def u16(x):
        return np.clip(np.rint(x), 0, 65535).astype("<u2")

    with open(out_path, "wb") as f:
        f.write(b"SIM3WRLD")
        grid_level = int(round(np.log((n - 2) / 10) / np.log(4)))  # cells = 10 * 4^level + 2
        assert 10 * 4**grid_level + 2 == n
        np.array([FORMAT_VERSION, grid_level, n, 2], dtype="<u4").tofile(f)
        for year, mean, lo, hi, land, ice in [
            (-9999, s_then, s_min + delta, s_max + delta, land_then, ice_then),
            (2022, s_now, s_min, s_max, land_now, ice_now),
        ]:
            np.array([year], dtype="<i4").tofile(f)
            i16(mean).tofile(f)
            i16(lo).tofile(f)
            i16(hi).tofile(f)
            frac_u8(land).tofile(f)
            u16(ice).tofile(f)

        # Observations are compared on land only (cells at least half land
        # today), so only those are kept; this also keeps the file small.
        obs_cells, obs_t = observed_temperature(chelsa_pattern, tree, n, land_now / sum_w >= 0.5)
        report["observed_temperature_cells"] = int(obs_cells.size)
        section(f, b"OBST", [np.array([obs_cells.size], "<u4"), obs_cells.astype("<u4"),
                             np.clip(np.rint(obs_t * 100), -32768, 32767).astype("<i2")])
        co2 = read_co2(co2_path)
        section(f, b"CO2 ", [np.array([len(co2)], "<u4"), np.array(co2, "<f4")])
        co2_now = read_mlo_mean(mlo_path, 1981, 2010)
        report["co2_1981_2010_ppm"] = co2_now
        section(f, b"CO2P", [np.array([co2_now], "<f4")])
        orbit = read_orbit(orbit_path)
        section(f, b"ORBT", [np.array([len(orbit)], "<u4"), np.array(orbit, "<f4")])
    print(json.dumps(report, indent=2), flush=True)
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)


def section(f, tag, arrays):
    """A tagged block: 4-byte tag, u32 byte length, payload."""
    payload = b"".join(a.tobytes() for a in arrays)
    f.write(tag)
    np.array([len(payload)], "<u4").tofile(f)
    f.write(payload)


def observed_temperature(pattern, tree, n, keep_mask):
    """CHELSA 1981-2010 monthly mean temperature (deg C) averaged over each
    map cell, for cells in `keep_mask` that CHELSA covers at least half of.

    CHELSA is 30 arc-seconds; each 2x2 block of pixels is assigned to the
    cell containing the block's centre (blocks are ~1.8 km, cells ~28 km).
    """
    sums = np.zeros((12, n))
    wsum = np.zeros(n)
    wall = np.zeros(n)
    files = [rasterio.open(pattern.replace("MM", f"{m:02d}")) for m in range(1, 13)]
    f0 = files[0]
    t = f0.transform
    scale, offset = f0.scales[0], f0.offsets[0]
    rows, cols = f0.height, f0.width
    rows -= rows % 2
    block_lons = t.c + t.a * (np.arange(0, cols, 2) + 1.0)
    step = 400  # full-resolution rows per chunk (even)
    for r0 in range(0, rows, step):
        nr = min(step, rows - r0)
        block_lats = t.f + t.e * (np.arange(r0, r0 + nr, 2) + 1.0)
        lat2 = np.repeat(block_lats, block_lons.size)
        lon2 = np.tile(block_lons, block_lats.size)
        _, bidx = tree.query(unit_vectors(lat2, lon2), workers=-1)
        bidx = bidx.reshape(block_lats.size, block_lons.size)
        idx = np.repeat(np.repeat(bidx, 2, axis=0), 2, axis=1)[:nr, :cols].ravel()
        lat_px = np.repeat(t.f + t.e * (np.arange(r0, r0 + nr) + 0.5), cols)
        w = np.cos(np.radians(lat_px))
        wall += np.bincount(idx, w, n)
        valid = None
        for m, fm in enumerate(files):
            raw = fm.read(1, window=Window(0, r0, cols, nr)).ravel()
            ok = (raw > 0) & (raw < 65535)  # 0 and 65535 mark no data (sea)
            if valid is None:
                valid = ok
                wsum += np.bincount(idx[ok], w[ok], n)
            temp = raw[ok].astype(np.float64) * scale + offset
            sums[m] += np.bincount(idx[ok], w[ok] * temp, n)
        print(f"CHELSA rows {r0 + nr}/{rows}", flush=True)
    for fm in files:
        fm.close()
    # Cells mostly covered by CHELSA; for cells partly outside its coverage
    # (north of 84 N there is none), wall undercounts, which is fine.
    keep = np.flatnonzero((wall > 0) & (wsum >= 0.5 * wall) & keep_mask)
    temps = (sums[:, keep] / wsum[keep]).T  # cells x 12
    return keep, temps.ravel()


def read_co2(path):
    """(age in years before 1950, CO2 ppm) pairs for the last 15,000 years."""
    rows = []
    header = None
    for line in open(path, encoding="utf-8-sig", errors="replace"):
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if header is None:
            header = [h.lower() for h in parts]
            continue
        try:
            age = float(parts[header.index(next(h for h in header if h.startswith("age")))])
            ppm = float(parts[header.index(next(h for h in header if h.startswith("co2") and "1s" not in h and "sigma" not in h))])
        except (ValueError, StopIteration, IndexError):
            continue
        if age <= 15000:
            rows.append((age, ppm))
    rows.sort()
    assert len(rows) > 50, f"CO2 record parse failed ({len(rows)} rows); header {header}"
    return rows


def read_mlo_mean(path, y0, y1):
    vals = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if y0 <= int(float(p[0])) <= y1:
            vals.append(float(p[1]))
    assert len(vals) == y1 - y0 + 1, "Mauna Loa record incomplete"
    return float(np.mean(vals))


def read_orbit(path):
    """(kyr relative to 1950, eccentricity, longitude of perihelion measured
    from the northern spring equinox in degrees, obliquity in degrees) for
    the last 15,000 years. The file's OMEGA is measured from the autumn
    equinox, so 180 degrees are added (Berger 1978, appendix; as in climlab).
    """
    rows = []
    for line in open(path):
        p = line.split()
        if len(p) < 4:
            continue
        try:
            k, ecc, omega, obl = int(p[0]), float(p[1]), float(p[2]), float(p[3])
        except ValueError:
            continue
        if -15 <= k <= 0:
            rows.append((k, ecc, (omega + 180.0) % 360.0, obl))
    rows.sort()
    assert len(rows) == 16, f"orbit table parse failed: {rows}"
    return rows


if __name__ == "__main__":
    main(*sys.argv[1:])
