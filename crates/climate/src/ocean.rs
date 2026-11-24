//! Wind-driven ocean currents.
//!
//! The winds' stress on the sea surface drives the great loops (gyres) in
//! each ocean basin. In the interior of a basin the depth-integrated flow
//! follows Sverdrup's balance (Sverdrup 1947, Proc. Natl. Acad. Sci. 33,
//! 318-326): β V = curl τ / ρ₀, with V the north-south transport per metre
//! of width and β the northward change of the Coriolis parameter. Starting
//! from each basin's eastern shore and adding up westward gives the
//! transport streamfunction ψ; at the western shore the loop closes in a
//! narrow, fast western boundary current (Stommel 1948): the Gulf Stream,
//! the Kuroshio, the Agulhas, the Brazil and East Australian currents.
//!
//! Only the winds and the coastlines go in, so any map gets its own
//! currents (project rule). The basins are found from the land/sea pattern,
//! row by row, on a 1° latitude-longitude raster.
//!
//! Recomputed once a year from the year's average winds: the gyres take
//! years to adjust to the winds anyway.
//!
//! The streamfunction is sampled at the corners of the 440 km cells, so
//! the water crossing each cell edge is the difference of ψ at its two
//! ends: what leaves one cell enters the next exactly, and the currents
//! neither make nor lose water.

use world::grid::{from_lat_lon, parent_map, to_lat_lon, Grid, Vec3, EARTH_RADIUS_M};

use crate::atmosphere::{AIR_DENSITY, DRAG_COEF, GUSTINESS, OMEGA_EARTH};

/// Seawater density (kg/m³) and specific heat (J/kg/K) (as in model.rs).
pub const SEA_DENSITY: f64 = 1025.0;
pub const SEA_SPECIFIC_HEAT: f64 = 3985.0;
/// Ekman's balance (the wind-driven surface layer moves at right angles
/// to the wind, transport τ × ẑ / (ρ₀ f)) breaks down near the equator;
/// closer than 5° the Coriolis parameter at 5° is used.
pub const EKMAN_MIN_LAT_DEG: f64 = 5.0;
/// Overturning (Stommel 1961, Tellus 13, 224-230): a warm and a cold box
/// per basin, flow q = k (α ΔT - β ΔS). Thermal and haline expansion of
/// seawater (Gill 1982, Atmosphere-Ocean Dynamics, app. 3).
pub const THERMAL_EXPANSION: f64 = 2.0e-4; // 1/K
pub const HALINE_CONTRACTION: f64 = 7.6e-4; // 1/psu
/// Mean salinity of seawater, psu.
pub const MEAN_SALINITY: f64 = 34.7;
/// Depth of the water each overturning box stands for: the mean ocean
/// depth, 3,700 m (Charette & Smith 2010, Oceanography 23, 112-114).
pub const OVERTURNING_DEPTH_M: f64 = 3700.0;
/// Overturning strength per unit density difference, m³/s. Tuned once so
/// the basin reaching the North Atlantic carries today's 17 Sv (RAPID
/// array 2004-2012; McCarthy et al. 2015, Prog. Oceanogr. 130, 91-111).
pub const OVERTURNING_K: f64 = 6.0e9;
/// Water rising into the surface layer comes from about 150 m down (the
/// source depth of coastal upwelling, 50-300 m; Huyer 1983, Prog.
/// Oceanogr. 12, 259-284). Below the surface layer temperature falls
/// toward the deep-water value with an e-folding depth of about 500 m, the
/// scale of the main thermocline (Pedlosky 1996). So rising water is at
/// T_deep + (T_surface - T_deep) exp(-150/500).
pub const UPWELLING_SOURCE_DEPTH_M: f64 = 150.0;
pub const THERMOCLINE_SCALE_M: f64 = 500.0;
/// Response time of the overturning, years: it changes over decades
/// (e.g. Delworth, Manabe & Stouffer 1993, J. Climate 6, 1993-2011, find
/// multidecadal variability); 20 years is used.
pub const OVERTURNING_ADJUSTMENT_YEARS: f64 = 20.0;
/// Time for the gyres to settle after the winds change, years: the
/// crossing time of long planetary (Rossby) waves over a subtropical
/// basin, about 5 years (Anderson & Gill 1975).
pub const GYRE_ADJUSTMENT_YEARS: f64 = 5.0;
/// Depth of the wind-driven upper ocean (the main thermocline), m: about
/// 1 km in the subtropical gyres (Pedlosky 1996, Ocean Circulation Theory).
/// Used only to show a typical current speed.
pub const GYRE_DEPTH_M: f64 = 1000.0;

const N_LON: usize = 360;
const N_LAT: usize = 180;

pub struct Ocean {
    /// Level-4 cells: per neighbour entry, water flowing out across that
    /// edge, m³/s (negative = flowing in).
    pub outflow: Vec<f64>,
    /// Per level-4 cell: upwind heat-carrying coefficients for the
    /// temperature solve, W/m²/K per neighbour entry (inflow only).
    pub heat_adv: Vec<f64>,
    /// Per level-4 cell: typical current (east, north) in the upper ocean,
    /// m/s, for drawing.
    pub current: Vec<(f64, f64)>,
    offsets: Vec<u32>,
    nbrs: Vec<u32>,
    corner_ids: Vec<u32>,
    corner_ll: Vec<(f64, f64)>,
    area: Vec<f64>,
    lat_deg: Vec<f64>,
    lon_deg: Vec<f64>,
    /// Raster point → climate cell.
    raster_fine: Vec<u32>,
    fine_lat: Vec<f64>,
    fine_lon: Vec<f64>,
    fine_area: Vec<f64>,
    /// Climate cell → level-4 cell.
    parent: Vec<u32>,
    /// Year's sums of wind stress (east, north) per level-4 cell, N/m², and
    /// the months summed.
    stress_sum: Vec<(f64, f64)>,
    months: u32,
    /// Streamfunction on the raster, m³/s.
    pub psi: Vec<f64>,
    /// Upwelling at the base of the mixed layer per level-4 cell, m/s (0
    /// where water sinks).
    pub upwelling: Vec<f64>,
    edge_dir: Vec<(f64, f64)>,
    edge_len: Vec<f64>,
    /// Overturning basins and their state.
    pub basins: Vec<Basin>,
    /// Temperature of the water that rises from below, °C: that of the
    /// water sinking in the overturning, in its coldest month.
    pub deep_temperature: f64,
    raster_sea_init: bool,
    psi_ready: bool,
}

/// One overturning basin: a cold box (its sea at 50-70° latitude, one
/// stretch of longitude) and a warm box (the sea between 40° S and 40° N
/// in the same longitudes, which supplies the warm water: the Atlantic
/// overturning draws on the South Atlantic too).
pub struct Basin {
    pub name_lat_lon: (f64, f64),
    pub cold: Vec<u32>,
    pub warm: Vec<u32>,
    pub cold_area: f64,
    pub warm_area: f64,
    /// Salinity of the warm box minus the cold box, psu.
    pub salinity_diff: f64,
    /// Overturning flow, m³/s.
    pub flow: f64,
    /// Heat carried from the warm to the cold box, W.
    pub heat_w: f64,
    /// Rain minus evaporation over the warm and cold boxes, m per year,
    /// and their mean temperatures, °C.
    pub fresh_m_per_year: (f64, f64),
    pub temperatures: (f64, f64),
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl Ocean {
    /// `fine` is the climate grid; `coarse` the level-4 grid with
    /// `coarse_area` its cells' areas.
    pub fn new(fine: &Grid, coarse: &Grid, coarse_area: &[f64]) -> Ocean {
        let n = coarse.len();
        let mut raster_fine = Vec::with_capacity(N_LON * N_LAT);
        let mut last = 0;
        for r in 0..N_LAT {
            let lat = 89.5 - r as f64;
            for c in 0..N_LON {
                let lon = -179.5 + c as f64;
                last = fine.nearest_from(from_lat_lon(lat, lon), last);
                raster_fine.push(last as u32);
            }
        }
        let corner_ll = coarse.corners.iter().map(|&c| to_lat_lon(c)).collect();
        let lat_lon: Vec<(f64, f64)> = coarse.centers.iter().map(|&c| to_lat_lon(c)).collect();
        // Per neighbour entry: unit direction to the neighbour in the
        // cell's (east, north) frame and the shared edge's length.
        let mut edge_dir = Vec::with_capacity(coarse.nbrs.len());
        let mut edge_len = Vec::with_capacity(coarse.nbrs.len());
        for p in 0..n {
            let c = coarse.centers[p];
            let e = {
                let v = [-c[1], c[0], 0.0];
                let l = (v[0] * v[0] + v[1] * v[1]).sqrt().max(1e-12);
                [v[0] / l, v[1] / l, 0.0]
            };
            let nn = [
                c[1] * e[2] - c[2] * e[1],
                c[2] * e[0] - c[0] * e[2],
                c[0] * e[1] - c[1] * e[0],
            ];
            let (s, t) = (coarse.offsets[p] as usize, coarse.offsets[p + 1] as usize);
            for j in s..t {
                let q = coarse.centers[coarse.nbrs[j] as usize];
                let d = [q[0] - c[0], q[1] - c[1], q[2] - c[2]];
                let (x, y) = (dot(d, e), dot(d, nn));
                let l = (x * x + y * y).sqrt();
                edge_dir.push((x / l, y / l));
                let prev = if j == s { t - 1 } else { j - 1 };
                let a = coarse.corners[coarse.corner_ids[prev] as usize];
                let b = coarse.corners[coarse.corner_ids[j] as usize];
                edge_len.push(dot(a, b).clamp(-1.0, 1.0).acos() * EARTH_RADIUS_M);
            }
        }
        Ocean {
            outflow: vec![0.0; coarse.nbrs.len()],
            heat_adv: vec![0.0; coarse.nbrs.len()],
            current: vec![(0.0, 0.0); n],
            offsets: coarse.offsets.clone(),
            nbrs: coarse.nbrs.clone(),
            corner_ids: coarse.corner_ids.clone(),
            corner_ll,
            area: coarse_area.to_vec(),
            lat_deg: lat_lon.iter().map(|v| v.0).collect(),
            lon_deg: lat_lon.iter().map(|v| v.1).collect(),
            raster_fine,
            fine_lat: fine.centers.iter().map(|&c| to_lat_lon(c).0).collect(),
            fine_lon: fine.centers.iter().map(|&c| to_lat_lon(c).1).collect(),
            fine_area: fine.areas_m2.clone(),
            parent: parent_map(fine, coarse),
            stress_sum: vec![(0.0, 0.0); n],
            months: 0,
            psi: vec![0.0; N_LON * N_LAT],
            upwelling: vec![0.0; n],
            edge_dir,
            edge_len,
            basins: Vec::new(),
            deep_temperature: 2.0,
            raster_sea_init: false,
            psi_ready: false,
        }
    }

    /// Adds a month's surface winds (east, north; m/s per level-4 cell) to
    /// the year's wind stress: τ = ρ_air C_D (|U| + gusts) U (bulk formula).
    pub fn add_month(&mut self, wind: &[(f64, f64)]) {
        for (s, &(u, v)) in self.stress_sum.iter_mut().zip(wind) {
            let k = AIR_DENSITY * DRAG_COEF * ((u * u + v * v).sqrt() + GUSTINESS);
            s.0 += k * u;
            s.1 += k * v;
        }
        self.months += 1;
    }

    /// Recomputes the currents from the year's average wind stress.
    /// `sea` is each climate cell's share of open or ice-covered sea.
    pub fn update(&mut self, sea: &[f64]) {
        if self.months == 0 {
            return;
        }
        let m = self.months as f64;
        let stress: Vec<(f64, f64)> = self.stress_sum.iter().map(|s| (s.0 / m, s.1 / m)).collect();
        self.stress_sum.iter_mut().for_each(|s| *s = (0.0, 0.0));
        self.months = 0;
        // Raster of sea and stress.
        let wet: Vec<bool> = self
            .raster_fine
            .iter()
            .map(|&f| sea[f as usize] >= 0.5)
            .collect();
        let tau: Vec<(f64, f64)> = self
            .raster_fine
            .iter()
            .map(|&f| stress[self.parent[f as usize] as usize])
            .collect();
        let dlam = 1f64.to_radians();
        let at = |r: usize, c: usize| r * N_LON + c;
        let mut fresh = vec![0.0; N_LON * N_LAT];
        for r in 1..N_LAT - 1 {
            let phi = (89.5 - r as f64).to_radians();
            let cosp = phi.cos();
            let beta = 2.0 * OMEGA_EARTH * cosp / EARTH_RADIUS_M;
            let dx = EARTH_RADIUS_M * cosp * dlam;
            // A row that is sea all the way round has no eastern shore
            // (the Southern Ocean): no gyre there.
            if (0..N_LON).all(|c| wet[at(r, c)]) {
                continue;
            }
            // curl τ = (1/(R cos φ)) [∂τy/∂λ - ∂(τx cos φ)/∂φ].
            let curl = |c: usize| {
                let (e, w) = ((c + 1) % N_LON, (c + N_LON - 1) % N_LON);
                let dty = (tau[at(r, e)].1 - tau[at(r, w)].1) / (2.0 * dlam);
                let cn = ((89.5 - (r - 1) as f64).to_radians()).cos();
                let cs = ((89.5 - (r + 1) as f64).to_radians()).cos();
                let dtx = (tau[at(r - 1, c)].0 * cn - tau[at(r + 1, c)].0 * cs) / (2.0 * dlam);
                (dty - dtx) / (EARTH_RADIUS_M * cosp)
            };
            // Start just west of some land and walk west once round, so
            // every stretch of sea is reached from its eastern shore.
            let start = (0..N_LON).find(|&c| !wet[at(r, c)]).unwrap();
            let mut acc = 0.0;
            for k in 1..=N_LON {
                let c = (start + N_LON * 2 - k) % N_LON;
                if !wet[at(r, c)] {
                    acc = 0.0;
                    continue;
                }
                let v = curl(c) / (SEA_DENSITY * beta);
                acc -= v * dx;
                fresh[at(r, c)] = acc;
            }
        }
        // The gyres adjust to new winds over years, as planetary waves
        // cross the basin (Anderson & Gill 1975, Deep-Sea Res. 22, 583-596):
        // ψ relaxes toward this year's balance with a time constant of
        // `GYRE_ADJUSTMENT_YEARS`. Land is set at once.
        let keep = if self.psi_ready {
            (-1.0 / GYRE_ADJUSTMENT_YEARS).exp()
        } else {
            0.0
        };
        self.psi_ready = true;
        for (k, v) in self.psi.iter_mut().enumerate() {
            *v = if wet[k] {
                keep * *v + (1.0 - keep) * fresh[k]
            } else {
                0.0
            };
        }
        // Water across each level-4 edge from ψ at its two corners.
        let psi_at = |lat: f64, lon: f64| {
            let y = (89.5 - lat).clamp(0.0, (N_LAT - 1) as f64);
            let x = (lon + 179.5).rem_euclid(N_LON as f64);
            let (r0, c0) = (y.floor() as usize, x.floor() as usize);
            let (r1, c1) = ((r0 + 1).min(N_LAT - 1), (c0 + 1) % N_LON);
            let (fy, fx) = (y - r0 as f64, x - x.floor());
            let p = &self.psi;
            (p[at(r0, c0)] * (1.0 - fx) + p[at(r0, c1)] * fx) * (1.0 - fy)
                + (p[at(r1, c0)] * (1.0 - fx) + p[at(r1, c1)] * fx) * fy
        };
        let corner_psi: Vec<f64> = self.corner_ll.iter().map(|&(a, b)| psi_at(a, b)).collect();
        let n = self.area.len();
        for p in 0..n {
            let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
            for j in s..e {
                let prev = if j == s { e - 1 } else { j - 1 };
                let a = corner_psi[self.corner_ids[prev] as usize];
                let b = corner_psi[self.corner_ids[j] as usize];
                self.outflow[j] = OUTFLOW_SIGN * (b - a);
            }
            // Heat carried in from upstream: ρ c_p × inflow / area, times
            // the share of the surface temperature difference the moving
            // water carries on average over the gyre's depth (temperature
            // falls off below the surface with the thermocline scale z:
            // (z / H)(1 - exp(-H / z)), about 0.43).
            let depth_share = THERMOCLINE_SCALE_M / GYRE_DEPTH_M
                * (1.0 - (-GYRE_DEPTH_M / THERMOCLINE_SCALE_M).exp());
            for j in s..e {
                self.heat_adv[j] =
                    depth_share * SEA_DENSITY * SEA_SPECIFIC_HEAT * (-self.outflow[j]).max(0.0)
                        / self.area[p];
            }
            // Typical current from the water crossing the cell's edges: for
            // a uniform flow U, Σ (outflow / edge length) d = (n/2) U H.
            let (mut ux, mut uy) = (0.0, 0.0);
            for j in s..e {
                let k = self.outflow[j] / self.edge_len[j];
                ux += k * self.edge_dir[j].0;
                uy += k * self.edge_dir[j].1;
            }
            let scale = 2.0 / ((e - s) as f64 * GYRE_DEPTH_M);
            self.current[p] = (ux * scale, uy * scale);
        }
    }

    /// The month's upwelling from the winds: the Ekman transport of the
    /// surface layer (τ × ẑ / (ρ₀ f)) spreading apart pulls water up from
    /// below, along coasts where winds blow along the shore and wherever
    /// winds diverge. `sea` is each level-4 cell's share of sea.
    pub fn update_upwelling(&mut self, wind: &[(f64, f64)], sea: &[f64]) {
        let n = self.area.len();
        let f_min = 2.0 * OMEGA_EARTH * EKMAN_MIN_LAT_DEG.to_radians().sin();
        let ekman: Vec<(f64, f64)> = (0..n)
            .map(|p| {
                if sea[p] < 0.5 {
                    return (0.0, 0.0);
                }
                let (u, v) = wind[p];
                let k = AIR_DENSITY * DRAG_COEF * ((u * u + v * v).sqrt() + GUSTINESS);
                let (tx, ty) = (k * u, k * v);
                let f = 2.0 * OMEGA_EARTH * self.lat_deg[p].to_radians().sin();
                let f = if f.abs() < f_min {
                    f_min.copysign(f)
                } else {
                    f
                };
                (ty / (SEA_DENSITY * f), -tx / (SEA_DENSITY * f))
            })
            .collect();
        for p in 0..n {
            if sea[p] < 0.5 {
                self.upwelling[p] = 0.0;
                continue;
            }
            let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
            let mut out = 0.0;
            for j in s..e {
                let q = self.nbrs[j] as usize;
                let (dx, dy) = self.edge_dir[j];
                let m = 0.5 * ((ekman[p].0 + ekman[q].0) * dx + (ekman[p].1 + ekman[q].1) * dy);
                out += m * self.edge_len[j];
            }
            self.upwelling[p] = (out / self.area[p]).max(0.0);
        }
    }

    /// Finds the overturning basins from the sea pattern (each stretch of
    /// sea along 60° latitude reaching from 50° to 70°).
    fn find_basins(&mut self, sea: &[f64]) {
        let wet = |lat: f64, lon: f64| {
            let r = ((89.5 - lat).round() as usize).min(N_LAT - 1);
            let c = ((lon + 179.5).round().rem_euclid(N_LON as f64)) as usize % N_LON;
            sea[self.raster_fine[r * N_LON + c] as usize] >= 0.5
        };
        let mut basins = Vec::new();
        for hemi in [1.0f64, -1.0] {
            let row = 60.0 * hemi;
            let w: Vec<bool> = (0..N_LON).map(|c| wet(row, -179.5 + c as f64)).collect();
            let mut segs: Vec<(usize, usize)> = Vec::new();
            if w.iter().all(|&x| x) {
                segs.push((0, N_LON));
            } else if let Some(start) = (0..N_LON).find(|&c| !w[c]) {
                let mut c = start;
                let mut k = 0;
                while k < N_LON {
                    if w[c] {
                        let s0 = c;
                        let mut len = 0;
                        while w[c] && len < N_LON {
                            c = (c + 1) % N_LON;
                            len += 1;
                            k += 1;
                        }
                        if len >= 10 {
                            segs.push((s0, len));
                        }
                    } else {
                        c = (c + 1) % N_LON;
                        k += 1;
                    }
                }
            }
            for (s0, len) in segs {
                // The stretch must connect by sea, within its own
                // longitudes, down to 40°: a closed bay doesn't overturn.
                let mut reach: Vec<bool> =
                    (0..N_LON).map(|c| (c + N_LON - s0) % N_LON < len).collect();
                for step in 1..=20 {
                    let lat = (60.0 - step as f64) * hemi;
                    let row: Vec<bool> = (0..N_LON).map(|c| wet(lat, -179.5 + c as f64)).collect();
                    let mut next: Vec<bool> = (0..N_LON)
                        .map(|c| reach[c] && row[c] && (c + N_LON - s0) % N_LON < len)
                        .collect();
                    // Spread sideways along the row within the stretch.
                    for _ in 0..len {
                        let mut changed = false;
                        for c in 0..N_LON {
                            if next[c] || !row[c] || (c + N_LON - s0) % N_LON >= len {
                                continue;
                            }
                            if next[(c + 1) % N_LON] || next[(c + N_LON - 1) % N_LON] {
                                next[c] = true;
                                changed = true;
                            }
                        }
                        if !changed {
                            break;
                        }
                    }
                    reach = next;
                }
                if !reach.iter().any(|&x| x) {
                    continue;
                }
                let in_lon = |lon: f64| {
                    let c = ((lon + 179.5).round().rem_euclid(N_LON as f64)) as usize % N_LON;
                    (c + N_LON - s0) % N_LON < len
                };
                let mut b = Basin {
                    name_lat_lon: (row, -179.5 + ((s0 + len / 2) % N_LON) as f64),
                    cold: Vec::new(),
                    warm: Vec::new(),
                    cold_area: 0.0,
                    warm_area: 0.0,
                    salinity_diff: 0.0,
                    flow: 0.0,
                    heat_w: 0.0,
                    fresh_m_per_year: (0.0, 0.0),
                    temperatures: (0.0, 0.0),
                };
                for (i, (&lat, &lon)) in self.fine_lat.iter().zip(&self.fine_lon).enumerate() {
                    if sea[i] < 0.5 || !in_lon(lon) {
                        continue;
                    }
                    let a = lat * hemi;
                    if (50.0..70.0).contains(&a) {
                        b.cold.push(i as u32);
                        b.cold_area += sea[i] * self.fine_area[i];
                    } else if (-40.0..40.0).contains(&a) {
                        b.warm.push(i as u32);
                        b.warm_area += sea[i] * self.fine_area[i];
                    }
                }
                if b.cold_area > 0.0 && b.warm_area > 0.0 {
                    basins.push(b);
                }
            }
        }
        self.basins = basins;
    }

    /// One year of the overturning. `t_annual` and `t_coldest` are each
    /// climate cell's yearly mean and coldest-month temperature (°C);
    /// `fresh` its yearly rain minus evaporation (m of water per second);
    /// `sea` its sea share.
    pub fn update_overturning(
        &mut self,
        t_annual: &[f64],
        t_coldest: &[f64],
        fresh: &[f64],
        sea: &[f64],
    ) {
        if !self.raster_sea_init {
            self.find_basins(sea);
            self.raster_sea_init = true;
        }
        let year = 365.25 * 86_400.0;
        let (mut tq, mut wq) = (0.0, 0.0);
        let first_year = self.basins.iter().all(|b| b.flow == 0.0);
        for b in &mut self.basins {
            let mean = |cells: &[u32], f: &[f64], area: f64| {
                cells
                    .iter()
                    .map(|&i| f[i as usize] * sea[i as usize] * self.fine_area[i as usize])
                    .sum::<f64>()
                    / area
            };
            let (tw, tc) = (
                mean(&b.warm, t_annual, b.warm_area),
                mean(&b.cold, t_annual, b.cold_area),
            );
            let (fw, fc) = (
                mean(&b.warm, fresh, b.warm_area),
                mean(&b.cold, fresh, b.cold_area),
            );
            let (vw, vc) = (
                b.warm_area * OVERTURNING_DEPTH_M,
                b.cold_area * OVERTURNING_DEPTH_M,
            );
            // Salinity: rain freshens, evaporation salts; the flow mixes
            // the boxes. From a standing start, settle the first year.
            let steps = if first_year { 2000 } else { 1 };
            for _ in 0..steps {
                let q = OVERTURNING_K
                    * (THERMAL_EXPANSION * (tw - tc) - HALINE_CONTRACTION * b.salinity_diff);
                // Implicit in the mixing term, so any box size is stable.
                let forcing = -MEAN_SALINITY * (fw * b.warm_area / vw - fc * b.cold_area / vc);
                let mixing = q.abs() * (1.0 / vw + 1.0 / vc);
                b.salinity_diff = (b.salinity_diff + forcing * year) / (1.0 + mixing * year);
                // The overturning takes decades to respond (its flow
                // relaxes toward the box balance with
                // `OVERTURNING_ADJUSTMENT_YEARS`); from a standing start it
                // settles at once.
                let keep = if first_year {
                    0.0
                } else {
                    (-1.0 / OVERTURNING_ADJUSTMENT_YEARS).exp()
                };
                b.flow = keep * b.flow + (1.0 - keep) * q;
            }
            b.heat_w = SEA_DENSITY * SEA_SPECIFIC_HEAT * b.flow.abs() * (tw - tc);
            b.fresh_m_per_year = (fw * year, fc * year);
            b.temperatures = (tw, tc);
            tq += b.flow.abs() * mean(&b.cold, t_coldest, b.cold_area);
            wq += b.flow.abs();
        }
        if wq > 0.0 {
            self.deep_temperature = (tq / wq).max(crate::model::SEA_FREEZING_C);
        }
    }

    /// Heat added by the overturning per climate cell, W/m² of the cell:
    /// taken evenly from each basin's warm box and given to its cold box.
    pub fn overturning_heating(&self, out: &mut [f64], sea: &[f64]) {
        out.iter_mut().for_each(|v| *v = 0.0);
        for b in &self.basins {
            for &i in &b.cold {
                out[i as usize] += sea[i as usize] * b.heat_w / b.cold_area;
            }
            for &i in &b.warm {
                out[i as usize] -= sea[i as usize] * b.heat_w / b.warm_area;
            }
        }
    }

    /// Every level-4 cell's position and typical upper-ocean current
    /// (east, north; m/s), for drawing. Land cells have none.
    pub fn currents(&self) -> Vec<(f64, f64, f64, f64)> {
        (0..self.area.len())
            .map(|p| {
                (
                    self.lat_deg[p],
                    self.lon_deg[p],
                    self.current[p].0,
                    self.current[p].1,
                )
            })
            .collect()
    }

    /// Largest streamfunction value in a latitude band, Sv (10⁶ m³/s), for
    /// checks: the strength of the gyres there.
    pub fn gyre_strength_sv(&self, lat_from: f64, lat_to: f64) -> f64 {
        let mut best = 0.0f64;
        for r in 0..N_LAT {
            let lat = 89.5 - r as f64;
            if lat < lat_from || lat > lat_to {
                continue;
            }
            for c in 0..N_LON {
                best = best.max(self.psi[r * N_LON + c].abs());
            }
        }
        best / 1e6
    }
}

/// Sign that makes ψ(b) - ψ(a) the outward flow across the edge from corner
/// a to corner b (corners go round each cell in the grid's order; checked
/// by the `eastward_flow_leaves_through_the_east` test).
const OUTFLOW_SIGN: f64 = -1.0;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eastward_flow_leaves_through_the_east() {
        // ψ = -U y gives u = U eastward: a cell's eastern edges must show
        // water leaving.
        let fine = Grid::new(5);
        let coarse = Grid::new(4);
        let area = coarse.areas_m2.clone();
        let mut o = Ocean::new(&fine, &coarse, &area);
        for r in 0..N_LAT {
            let lat = (89.5 - r as f64).to_radians();
            for c in 0..N_LON {
                o.psi[r * N_LON + c] = -EARTH_RADIUS_M * lat;
            }
        }
        // Reuse update's corner sampling by calling it with a fake year:
        // compute outflows directly instead.
        let p = coarse.nearest_from(from_lat_lon(20.0, 30.0), 0);
        let (s, e) = (coarse.offsets[p] as usize, coarse.offsets[p + 1] as usize);
        let c = coarse.centers[p];
        let east = [-c[1], c[0], 0.0];
        for j in s..e {
            let prev = if j == s { e - 1 } else { j - 1 };
            let la = to_lat_lon(coarse.corners[coarse.corner_ids[prev] as usize])
                .0
                .to_radians();
            let lb = to_lat_lon(coarse.corners[coarse.corner_ids[j] as usize])
                .0
                .to_radians();
            let out = OUTFLOW_SIGN * (-EARTH_RADIUS_M * lb + EARTH_RADIUS_M * la);
            let q = coarse.centers[coarse.nbrs[j] as usize];
            let dir = dot([q[0] - c[0], q[1] - c[1], q[2] - c[2]], east);
            if dir.abs() > 0.02 {
                assert_eq!(out > 0.0, dir > 0.0, "edge {j}: out {out} dir {dir}");
            }
        }
        let _ = &mut o;
    }
}
