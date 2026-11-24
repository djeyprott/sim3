//! Winds, moisture and rain on the 440 km grid (the same coarse grid that
//! carries heat; see solver.rs). Runs once a month.
//!
//! Everything here is driven by the temperature pattern and the terrain,
//! never by knowledge of where Earth's continents are, so it works on any
//! map (project rule).
//!
//! **Pressure.** Sea-level pressure has two parts:
//! - the planet-wide belts of a rotating Earth-like planet (a low along the
//!   tropical rain belt, highs near the edges of the tropical overturning
//!   cells, lows in the storm belts, highs at the poles). Their positions
//!   follow the "thermal equator" (the warmest latitude, found per 10° of
//!   longitude so it swings north and south over continents);
//! - thermal lows and highs: a column warmer than its latitude's average
//!   expands, lowering surface pressure by ρ g H ΔT / T (hydrostatic), with
//!   H the depth of the heated layer.
//!
//! Pressure differences are spread over the distance gravity waves carry
//! them (the Rossby radius).
//!
//! **Winds** follow the pressure gradient, turned by the Earth's rotation
//! and slowed by surface friction (bulk drag law).
//!
//! **Moisture.** Each month the column water vapour W (kg/m²) is solved as a
//! balance: evaporation in (from the sea, and from wet ground), carried by
//! the winds and by large-scale mixing, and rained out. Rain comes from
//! - humid columns: rain falls more readily the closer the column is to
//!   saturation (the satellite-observed relationship of Bretherton et al.
//!   2004);
//! - converging winds, where air rises (Kuo 1974);
//! - winds blowing up slopes.
//!
//! The last two only once the rising air passes its condensation level.

use crate::model::LAPSE_RATE;
use world::grid::{parent_map, to_lat_lon, Grid, Vec3, EARTH_RADIUS_M};

/// Angular speed of the Earth's rotation, rad/s (sidereal day 86,164 s).
pub const OMEGA_EARTH: f64 = 7.292_115e-5;
/// Near-surface air density, kg/m³ (standard atmosphere, sea level).
pub const AIR_DENSITY: f64 = 1.225;
/// Standard gravity, m/s².
pub const GRAVITY: f64 = 9.806_65;
/// Drag coefficient of the sea surface for momentum, 1.2e-3 (Large & Pond
/// 1981, J. Phys. Oceanogr. 11, 324-336; land is rougher but its winds are
/// measured lower; one value is used).
pub const DRAG_COEF: f64 = 1.2e-3;
/// Depth of the well-mixed layer next to the surface whose winds are
/// modelled: about 500 m, the cloud-base height of the trade-wind layer
/// (Stevens 2005, Annu. Rev. Earth Planet. Sci. 33, 605-643).
pub const MIXED_LAYER_DEPTH_M: f64 = 500.0;
/// Depth of the layer whose heating makes thermal lows: monsoon heat lows
/// reach about 700 hPa, roughly 3 km; 2.5 km is used (an Asian summer
/// anomaly of +10 K then gives about -10 hPa, as observed).
pub const THERMAL_LOW_DEPTH_M: f64 = 2_500.0;
/// Temperature range over which latitudes count toward the thermal
/// equator, K (see `update_winds`).
pub const THERMAL_EQUATOR_WIDTH_K: f64 = 1.0;
/// Buoyancy frequency × depth of the heated layer, m/s: the speed of the
/// gravity waves that spread pressure changes (N ≈ 0.01 /s in the
/// troposphere; Gill 1982, Atmosphere-Ocean Dynamics). Pressure differences
/// made by heating are spread over the Rossby radius NH/f, or the
/// equatorial radius sqrt(NH/β) near the equator (Gill 1980, Q. J. R.
/// Meteorol. Soc. 106, 447-462).
pub const GRAVITY_WAVE_SPEED: f64 = 0.01 * THERMAL_LOW_DEPTH_M;
/// Scale height of water vapour, m: column vapour ≈ surface saturation
/// vapour density × 2 km (Hartmann, Global Physical Climatology, 2016).
pub const VAPOUR_SCALE_HEIGHT_M: f64 = 2_000.0;
/// Heat capacity of the air column, cp·ps/g, J/m²/K.
pub const COLUMN_HEAT_CAPACITY: f64 = 1004.0 * 101_325.0 / 9.81;
/// Bulk transfer coefficient for evaporation over the sea (Large & Pond
/// 1982, J. Phys. Oceanogr. 12, 464-482): 1.2e-3.
pub const EVAPORATION_COEF: f64 = 1.2e-3;
/// Relative humidity of air just above the surface, sea and moist land
/// (observed marine boundary layer, ~80%).
pub const SURFACE_RELATIVE_HUMIDITY: f64 = 0.8;
/// Extra wind speed for evaporation from gusts and weather that monthly-mean
/// winds miss (bulk-flux "gustiness", commonly 3 m/s; Godfrey & Beljaars
/// 1991, J. Geophys. Res. 96, 22043).
pub const GUSTINESS: f64 = 3.0;
/// Bretherton, Peters & Back (2004, J. Climate 17, 1517-1528): monthly
/// precipitation over the tropical oceans P = exp(a (r - b)) mm/day, r the
/// column relative humidity. Applied everywhere here (an extrapolation
/// outside the tropics; validated against observed rainfall).
pub const BRETHERTON_A: f64 = 15.6;
pub const BRETHERTON_B: f64 = 0.603;

/// Values a model this simple can't derive, set once from today's world.
#[derive(Clone, Copy, Debug)]
pub struct AtmosphereTuning {
    /// Belt strengths, Pa: tropical trough, subtropical highs, storm-belt
    /// lows, polar highs. Set to roughly match the observed yearly
    /// zonal-mean sea-level pressure (trough ~3 hPa below the global mean,
    /// subtropical highs ~6 above, subpolar lows ~9 below averaged over
    /// both hemispheres, polar highs ~4 above; Peixoto & Oort 1992,
    /// Physics of Climate, ch. 7). Not yet fitted to pressure data.
    pub trough_pa: f64,
    pub subtropical_pa: f64,
    pub storm_pa: f64,
    pub polar_pa: f64,
    /// Sideways mixing of moisture by weather systems, m²/s. Tuned to
    /// today's rainfall: of 2, 3, 4.5, 5, 8 and 10 × 10⁶ m²/s, 2 × 10⁶
    /// matched the observed land pattern best (correlation 0.60 vs 0.59,
    /// 0.58, ...). Typical of eddy diffusivities in moist energy balance
    /// models (~10⁶ m²/s).
    pub moisture_mixing: f64,
}

pub const ATMOSPHERE_TUNED: AtmosphereTuning = AtmosphereTuning {
    trough_pa: 300.0,
    subtropical_pa: 600.0,
    storm_pa: 900.0,
    polar_pa: 400.0,
    moisture_mixing: 2.0e6,
};

/// Saturation vapour pressure over water, Pa (Bolton 1980, Mon. Wea. Rev.
/// 108, 1046-1053), T in °C.
pub fn saturation_vapour_pressure(t_c: f64) -> f64 {
    611.2 * (17.67 * t_c / (t_c + 243.5)).exp()
}

/// Saturation vapour density, kg/m³, T in °C.
pub fn saturation_vapour_density(t_c: f64) -> f64 {
    // e / (R_v T), R_v = 461.5 J/kg/K.
    saturation_vapour_pressure(t_c) / (461.5 * (t_c + 273.15))
}

/// A function tabulated in 0.05 °C steps over -90..+60 °C.
pub struct Table {
    v: Vec<f64>,
}

impl Table {
    const MIN: f64 = -90.0;
    const STEP: f64 = 0.05;
    pub fn new(f: impl Fn(f64) -> f64) -> Table {
        let n = ((60.0 - Self::MIN) / Self::STEP) as usize + 2;
        Table {
            v: (0..n)
                .map(|k| f(Self::MIN + k as f64 * Self::STEP))
                .collect(),
        }
    }
    #[inline]
    pub fn at(&self, t: f64) -> f64 {
        let x = ((t - Self::MIN) / Self::STEP).clamp(0.0, (self.v.len() - 2) as f64);
        let k = x as usize;
        let f = x - k as f64;
        self.v[k] * (1.0 - f) + self.v[k + 1] * f
    }
}

pub struct Atmosphere {
    pub tuning: AtmosphereTuning,
    n: usize,
    lat_deg: Vec<f64>,
    offsets: Vec<u32>,
    nbrs: Vec<u32>,
    /// Edge length / cell area (1/m) and edge length / (centre distance ×
    /// cell area) (1/m²) per neighbour entry.
    l_over_a: Vec<f64>,
    g: Vec<f64>,
    /// Unit vector from each cell towards each neighbour, in the cell's
    /// (east, north) frame.
    dir: Vec<(f64, f64)>,
    /// Least-squares gradient weights per neighbour entry: ∇f ≈ Σ w (f_j - f_i).
    grad_w: Vec<(f64, f64)>,
    pub area: Vec<f64>,
    /// Mean surface height of each coarse cell (sea = 0), m.
    coarse_height: Vec<f64>,
    /// Fine (climate) cell → coarse cell.
    pub parent: Vec<u32>,
    pub fine_area: Vec<f64>,
    /// Climate-cell terrain slope (east, north), m/m.
    fine_slope: Vec<(f64, f64)>,
    fine_height: Vec<f64>,
    fine_lift: Vec<f64>,
    // State and outputs per coarse cell.
    /// Pressure anomaly, Pa.
    pub pressure: Vec<f64>,
    /// Heat carried by the air from upwind neighbours, W/m²/K per
    /// neighbour entry: C_air × inflow speed × edge / area.
    pub heat_adv: Vec<f64>,
    /// Near-surface wind (east, north), m/s.
    pub wind: Vec<(f64, f64)>,
    /// Vapour-weighted wind of the whole column (east, north), m/s.
    pub carry_wind: Vec<(f64, f64)>,
    /// Wind above the surface layer (geostrophic + thermal), m/s.
    pub free_wind: Vec<(f64, f64)>,
    /// Column water vapour, kg/m².
    pub vapour: Vec<f64>,
    /// Rain rate (kg/m²/s) per climate cell, from the last solve.
    pub precip: Vec<f64>,
    /// Evaporation (kg/m²/s) per coarse cell, from the last solve.
    pub evaporation: Vec<f64>,
    /// Evaporation from wet ground or open water (kg/m²/s) per coarse cell.
    pub potential_evaporation: Vec<f64>,
    sea_evaporation: Vec<f64>,
    /// Largest evaporation wet ground can sustain from its energy supply
    /// per coarse cell, kg/m²/s (Priestley-Taylor; set each month).
    pub land_energy_limit: Vec<f64>,
    sat_density: Table,
    // Work arrays.
    flow: Vec<f64>,
    surface_flow: Vec<f64>,
    wsat: Vec<f64>,
    uplift: Vec<f64>,
    pub thermal_equator_deg: f64,
    /// Thermal equator per 10° of longitude from 180°W.
    pub sector_equator: [f64; 36],
    lon_deg: Vec<f64>,
    pub moisture_sweeps: usize,
    vapour_by_month: Vec<Vec<f64>>,
    prev_vapour: Vec<f64>,
    pressure_by_month: Vec<Vec<f64>>,
    coef: Vec<f64>,
    lcl: Vec<f64>,
    /// Surface convergence per coarse cell, 1/s (0 where winds diverge).
    pub convergence: Vec<f64>,
    diag: Vec<f64>,
    inv_ws: Vec<f64>,
    /// Squared spreading radius per cell (m²) and the raw pressure.
    radius2: Vec<f64>,
    raw_pressure: Vec<f64>,
    /// For each neighbour entry p→q, the index of the entry q→p.
    reverse: Vec<u32>,
}

fn east_north(c: Vec3) -> (Vec3, Vec3) {
    let e = [-c[1], c[0], 0.0];
    let len = (e[0] * e[0] + e[1] * e[1]).sqrt();
    let e = if len < 1e-12 {
        [1.0, 0.0, 0.0]
    } else {
        [e[0] / len, e[1] / len, 0.0]
    };
    let n = [
        c[1] * e[2] - c[2] * e[1],
        c[2] * e[0] - c[0] * e[2],
        c[0] * e[1] - c[1] * e[0],
    ];
    (e, n)
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn arc(a: Vec3, b: Vec3) -> f64 {
    dot(a, b).clamp(-1.0, 1.0).acos() * EARTH_RADIUS_M
}

/// Least-squares gradient weights for a cell from its neighbours' offsets.
fn gradient_weights(offs: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for &(x, y) in offs {
        sxx += x * x;
        sxy += x * y;
        syy += y * y;
    }
    let det = sxx * syy - sxy * sxy;
    offs.iter()
        .map(|&(x, y)| {
            if det.abs() < 1e-30 {
                (0.0, 0.0)
            } else {
                ((syy * x - sxy * y) / det, (sxx * y - sxy * x) / det)
            }
        })
        .collect()
}

impl Atmosphere {
    /// `fine` is the climate grid; `fine_height` its cells' mean surface
    /// height (m, sea = 0).
    pub fn new(
        fine: &Grid,
        fine_height: &[f64],
        level: u32,
        tuning: AtmosphereTuning,
    ) -> Atmosphere {
        let coarse = Grid::new(level);
        let n = coarse.len();
        let parent = parent_map(fine, &coarse);
        let mut area = vec![0.0; n];
        let mut coarse_height = vec![0.0; n];
        for (i, &p) in parent.iter().enumerate() {
            area[p as usize] += fine.areas_m2[i];
            coarse_height[p as usize] += fine.areas_m2[i] * fine_height[i];
        }
        for p in 0..n {
            coarse_height[p] /= area[p];
        }
        let mut l_over_a = Vec::new();
        let mut g = Vec::new();
        let mut dir = Vec::new();
        let mut grad_w = Vec::new();
        for i in 0..n {
            let c = coarse.centers[i];
            let (e, nn) = east_north(c);
            let (s, t) = (coarse.offsets[i] as usize, coarse.offsets[i + 1] as usize);
            let mut offs = Vec::new();
            for j in s..t {
                let prev = if j == s { t - 1 } else { j - 1 };
                let edge = arc(
                    coarse.corners[coarse.corner_ids[prev] as usize],
                    coarse.corners[coarse.corner_ids[j] as usize],
                );
                let q = coarse.centers[coarse.nbrs[j] as usize];
                let dist = arc(c, q);
                l_over_a.push(edge / area[i]);
                g.push(edge / (dist * area[i]));
                let d = [q[0] - c[0], q[1] - c[1], q[2] - c[2]];
                let (x, y) = (dot(d, e), dot(d, nn));
                let len = (x * x + y * y).sqrt();
                dir.push((x / len, y / len));
                offs.push((x * EARTH_RADIUS_M, y * EARTH_RADIUS_M));
            }
            grad_w.extend(gradient_weights(&offs));
        }
        // Fine-cell slopes from fine-cell heights.
        let mut fine_slope = Vec::with_capacity(fine.len());
        for i in 0..fine.len() {
            let c = fine.centers[i];
            let (e, nn) = east_north(c);
            let nb = fine.neighbors(i);
            let offs: Vec<(f64, f64)> = nb
                .iter()
                .map(|&j| {
                    let q = fine.centers[j as usize];
                    let d = [q[0] - c[0], q[1] - c[1], q[2] - c[2]];
                    (dot(d, e) * EARTH_RADIUS_M, dot(d, nn) * EARTH_RADIUS_M)
                })
                .collect();
            let w = gradient_weights(&offs);
            let (mut sx, mut sy) = (0.0, 0.0);
            for (k, &j) in nb.iter().enumerate() {
                let dh = fine_height[j as usize] - fine_height[i];
                sx += w[k].0 * dh;
                sy += w[k].1 * dh;
            }
            fine_slope.push((sx, sy));
        }
        let lat_deg: Vec<f64> = coarse.centers.iter().map(|&c| to_lat_lon(c).0).collect();
        let radius2 = lat_deg
            .iter()
            .map(|&l| {
                let phi = l.to_radians();
                let f = 2.0 * OMEGA_EARTH * phi.sin();
                let beta = 2.0 * OMEGA_EARTH * phi.cos() / EARTH_RADIUS_M;
                GRAVITY_WAVE_SPEED * GRAVITY_WAVE_SPEED / (f * f + beta * GRAVITY_WAVE_SPEED)
            })
            .collect();
        let mut reverse = vec![0u32; coarse.nbrs.len()];
        for p in 0..n {
            for j in coarse.offsets[p] as usize..coarse.offsets[p + 1] as usize {
                let q = coarse.nbrs[j] as usize;
                let (s, e) = (coarse.offsets[q] as usize, coarse.offsets[q + 1] as usize);
                reverse[j] = (s..e)
                    .find(|&k| coarse.nbrs[k] as usize == p)
                    .expect("symmetric neighbours") as u32;
            }
        }
        let m = coarse.nbrs.len();
        Atmosphere {
            tuning,
            n,
            lat_deg,
            offsets: coarse.offsets.clone(),
            nbrs: coarse.nbrs.clone(),
            l_over_a,
            g,
            dir,
            grad_w,
            area,
            coarse_height,
            parent,
            fine_area: fine.areas_m2.clone(),
            fine_slope,
            fine_height: fine_height.to_vec(),
            fine_lift: vec![0.0; fine.len()],
            pressure: vec![0.0; n],
            wind: vec![(0.0, 0.0); n],
            carry_wind: vec![(0.0, 0.0); n],
            free_wind: vec![(0.0, 0.0); n],
            heat_adv: vec![0.0; m],
            surface_flow: vec![0.0; m],
            vapour: vec![10.0; n],
            precip: vec![0.0; fine.len()],
            evaporation: vec![0.0; n],
            potential_evaporation: vec![0.0; n],
            sea_evaporation: vec![0.0; n],
            land_energy_limit: vec![f64::MAX; n],
            sat_density: Table::new(saturation_vapour_density),
            flow: vec![0.0; m],
            wsat: vec![0.0; n],
            uplift: vec![0.0; n],
            thermal_equator_deg: 0.0,
            sector_equator: [0.0; 36],
            lon_deg: coarse.centers.iter().map(|&c| to_lat_lon(c).1).collect(),
            moisture_sweeps: 0,
            vapour_by_month: vec![Vec::new(); 12],
            prev_vapour: Vec::new(),
            pressure_by_month: vec![Vec::new(); 12],
            coef: vec![0.0; m],
            lcl: vec![0.0; n],
            convergence: vec![0.0; n],
            diag: vec![0.0; n],
            inv_ws: vec![0.0; n],
            reverse,
            radius2,
            raw_pressure: vec![0.0; n],
        }
    }

    pub fn lat(&self, p: usize) -> f64 {
        self.lat_deg[p]
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Pressure and winds from the coarse sea-level-equivalent temperatures.
    pub fn update_winds(&mut self, month: usize, t_sea_level: &[f64]) {
        let n = self.n;

        // Zonal-mean temperature in 1° bands, smoothed over ±5° (the grid's
        // rows of cells don't fill every band evenly).
        let bands = 180;
        let mut sum = vec![0.0; bands];
        let mut w = vec![0.0; bands];
        for p in 0..n {
            let b = ((self.lat_deg[p] + 90.0) as usize).min(bands - 1);
            sum[b] += t_sea_level[p] * self.area[p];
            w[b] += self.area[p];
        }
        let zonal: Vec<f64> = (0..bands as i64)
            .map(|b| {
                let (mut s, mut ww) = (0.0, 0.0);
                for d in -5i64..=5 {
                    let k = b + d;
                    if (0..bands as i64).contains(&k) {
                        let f = 6.0 - d.abs() as f64;
                        s += f * sum[k as usize];
                        ww += f * w[k as usize];
                    }
                }
                s / ww
            })
            .collect();
        // The thermal equator (the warmest latitude, where the rain belt
        // sits) is found per 10° of longitude, from the temperatures within
        // 20° of longitude either side, so the belt can swing far over
        // continents (monsoons) and little over oceans. Kept within ±25°.
        const SECTORS: usize = 36;
        const HALF_WINDOW_DEG: f64 = 20.0;
        let tropic_bands = 25;
        let mut st = vec![[0.0f64; 25]; SECTORS];
        let mut sw = vec![[0.0f64; 25]; SECTORS];
        for p in 0..n {
            let lat = self.lat_deg[p];
            if lat.abs() >= 25.0 {
                continue;
            }
            let b = (((lat + 25.0) / 2.0) as usize).min(tropic_bands - 1);
            let first = ((self.lon_deg[p] - HALF_WINDOW_DEG + 180.0) / 10.0 - 0.5).ceil() as i64;
            let last = ((self.lon_deg[p] + HALF_WINDOW_DEG + 180.0) / 10.0 - 0.5).floor() as i64;
            for k in first..=last {
                let k = k.rem_euclid(SECTORS as i64) as usize;
                st[k][b] += t_sea_level[p] * self.area[p];
                sw[k][b] += self.area[p];
            }
        }
        // Not just the single warmest band, which would make the belt jump
        // between hemispheres: the average latitude of the bands, each
        // weighted by exp(-(T_max - T) / W) with W = THERMAL_EQUATOR_WIDTH_K,
        // so it glides as the seasons change.
        let scale = THERMAL_EQUATOR_WIDTH_K;
        let mut eqs = [0.0f64; SECTORS];
        for k in 0..SECTORS {
            let mut tmax = f64::MIN;
            for b in 0..tropic_bands {
                if sw[k][b] > 0.0 {
                    tmax = tmax.max(st[k][b] / sw[k][b]);
                }
            }
            let (mut num, mut den) = (0.0, 0.0);
            for b in 0..tropic_bands {
                if sw[k][b] > 0.0 {
                    let wgt = (-(tmax - st[k][b] / sw[k][b]) / scale).exp();
                    num += wgt * (-25.0 + 2.0 * b as f64 + 1.0);
                    den += wgt;
                }
            }
            eqs[k] = if den > 0.0 { num / den } else { 0.0 };
        }
        self.thermal_equator_deg = eqs.iter().sum::<f64>() / SECTORS as f64;
        self.sector_equator = eqs;
        let t = &self.tuning;
        for p in 0..n {
            let lat = self.lat_deg[p];
            let x_sec = (self.lon_deg[p] + 180.0) / 10.0 - 0.5;
            let k0 = x_sec.floor();
            let f = x_sec - k0;
            let k0 = (k0 as i64).rem_euclid(SECTORS as i64) as usize;
            let k1 = (k0 + 1) % SECTORS;
            let eq = eqs[k0] * (1.0 - f) + eqs[k1] * f;
            // Latitude measured from the thermal equator, stretched so the
            // poles stay at ±90.
            let x = if lat >= eq {
                (lat - eq) / (90.0 - eq) * 90.0
            } else {
                (lat - eq) / (90.0 + eq) * 90.0
            };
            let ax = x.abs();
            let gauss = |c: f64, wd: f64| (-((ax - c) / wd).powi(2)).exp();
            let belts = -t.trough_pa * gauss(0.0, 8.0) + t.subtropical_pa * gauss(30.0, 10.0)
                - t.storm_pa * gauss(60.0, 10.0)
                + t.polar_pa * gauss(90.0, 10.0);
            let b = ((lat + 90.0) as usize).min(bands - 1);
            let anomaly = t_sea_level[p] - zonal[b];
            let thermal =
                -AIR_DENSITY * GRAVITY * THERMAL_LOW_DEPTH_M * anomaly / (zonal[b] + 273.15);
            self.raw_pressure[p] = belts + thermal;
        }
        // Spread: solve (1 - L² ∇²) p = p_raw by over-relaxed sweeps,
        // starting from the same month last year.
        if self.pressure_by_month[month].len() == n {
            self.pressure
                .copy_from_slice(&self.pressure_by_month[month]);
        }
        for _ in 0..300 {
            let mut max_change = 0.0f64;
            for p in 0..n {
                let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
                let (mut sg, mut sgp) = (0.0, 0.0);
                for j in s..e {
                    let a = self.radius2[p] * self.g[j];
                    sg += a;
                    sgp += a * self.pressure[self.nbrs[j] as usize];
                }
                let gs = (self.raw_pressure[p] + sgp) / (1.0 + sg);
                let new = self.pressure[p] + 1.6 * (gs - self.pressure[p]);
                max_change = max_change.max((new - self.pressure[p]).abs());
                self.pressure[p] = new;
            }
            if max_change < 1.0 {
                break;
            }
        }
        self.pressure_by_month[month].clone_from(&self.pressure);
        for p in 0..n {
            let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
            let (mut px, mut py) = (0.0, 0.0);
            for j in s..e {
                let dp = self.pressure[self.nbrs[j] as usize] - self.pressure[p];
                px += self.grad_w[j].0 * dp;
                py += self.grad_w[j].1 * dp;
            }
            let f = 2.0 * OMEGA_EARTH * self.lat_deg[p].to_radians().sin();
            // Friction from the bulk drag law: the surface takes momentum
            // out of the mixed layer at C_D (|U| + gusts) / h per second,
            // so stronger winds meet more drag. Solved by fixed-point
            // iteration (converges in a few steps).
            let g2 = px * px + py * py;
            let mut speed = 5.0;
            let mut r = 0.0;
            for _ in 0..6 {
                r = DRAG_COEF * (speed + GUSTINESS) / MIXED_LAYER_DEPTH_M;
                speed = (g2 / (r * r + f * f)).sqrt() / AIR_DENSITY;
            }
            let den = AIR_DENSITY * (r * r + f * f);
            self.wind[p] = ((-r * px - f * py) / den, (f * px - r * py) / den);
        }
        // The wind that carries vapour. Only the share of the vapour inside
        // the mixed layer, 1 - exp(-h/H_v), moves with the surface wind; the
        // rest sits higher, where there is no surface friction and the wind
        // is in geostrophic balance with the pressure, plus the thermal wind
        // that grows with height over a temperature contrast (storm tracks,
        // monsoons; Holton, An Introduction to Dynamic Meteorology, ch. 3).
        // That vapour's average height is h + H_v. Near the equator, where
        // the balance breaks down, f² is replaced by f² + βc (the same
        // equatorial scale as the pressure spreading above), so the
        // free-air wind fades out there instead of growing without bound.
        //
        // Where surface air spreads out, air sinks from above to replace it
        // at h × divergence, pushing the vapour above down into the mixed
        // layer at the rate the column turns over, h / H_v of it; the
        // surface wind then carries that too. This is what keeps the
        // subtropical deserts dry.
        let low_share = 1.0 - (-MIXED_LAYER_DEPTH_M / VAPOUR_SCALE_HEIGHT_M).exp();
        let surface_share =
            low_share + (1.0 - low_share) * MIXED_LAYER_DEPTH_M / VAPOUR_SCALE_HEIGHT_M;
        let z = MIXED_LAYER_DEPTH_M + VAPOUR_SCALE_HEIGHT_M;
        for p in 0..n {
            let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
            let (mut px, mut py, mut tx, mut ty) = (0.0, 0.0, 0.0, 0.0);
            for j in s..e {
                let q = self.nbrs[j] as usize;
                let dp = self.pressure[q] - self.pressure[p];
                let dt = t_sea_level[q] - t_sea_level[p];
                px += self.grad_w[j].0 * dp;
                py += self.grad_w[j].1 * dp;
                tx += self.grad_w[j].0 * dt;
                ty += self.grad_w[j].1 * dt;
            }
            let phi = self.lat_deg[p].to_radians();
            let f = 2.0 * OMEGA_EARTH * phi.sin();
            let beta = 2.0 * OMEGA_EARTH * phi.cos() / EARTH_RADIUS_M;
            let g = f / (f * f + beta * GRAVITY_WAVE_SPEED);
            let k = GRAVITY * z / (t_sea_level[p] + 273.15);
            let free = (
                -g * (py / AIR_DENSITY + k * ty),
                g * (px / AIR_DENSITY + k * tx),
            );
            self.carry_wind[p] = (
                surface_share * self.wind[p].0 + (1.0 - low_share) * free.0,
                surface_share * self.wind[p].1 + (1.0 - low_share) * free.1,
            );
            self.free_wind[p] = free;
        }
        // Velocity across each edge (positive outward), for the surface wind
        // (convergence) and the vapour-carrying wind. Each edge's speed is
        // computed once, from the lower-numbered cell, and mirrored, so what
        // leaves one cell enters the other exactly.
        for p in 0..n {
            let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
            for j in s..e {
                let q = self.nbrs[j] as usize;
                if q < p {
                    let r = self.reverse[j] as usize;
                    self.flow[j] = -self.flow[r];
                    self.surface_flow[j] = -self.surface_flow[r];
                    continue;
                }
                let (dx, dy) = self.dir[j];
                // Average the two cells' winds projected on the direction p→q
                // (the neighbour's wind is expressed in p's frame approximately:
                // over ~440 km the frames differ by a few degrees).
                let across = |w: &[(f64, f64)]| {
                    0.5 * (w[p].0 * dx + w[p].1 * dy + w[q].0 * dx + w[q].1 * dy)
                };
                self.surface_flow[j] = across(&self.wind);
                self.flow[j] = across(&self.carry_wind);
            }
        }
        // The air column's heat moves with the wind above the surface
        // layer (geostrophic plus thermal wind), the part that circles the
        // globe. The near-surface inflow to the rain belts goes back out
        // aloft, carrying little heat on balance; its effect is in the
        // large-scale mixing.
        for p in 0..n {
            let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
            for j in s..e {
                let q = self.nbrs[j] as usize;
                let (dx, dy) = self.dir[j];
                let w = &self.free_wind;
                let un = 0.5 * (w[p].0 * dx + w[p].1 * dy + w[q].0 * dx + w[q].1 * dy);
                self.heat_adv[j] = COLUMN_HEAT_CAPACITY * (-un).max(0.0) * self.l_over_a[j];
            }
        }
    }

    /// Solves the month's moisture balance.
    /// - `t_surface`: coarse surface temperature, °C.
    /// - `sea`: coarse open-water share (0-1), `soil_wetness`: coarse land
    ///   evaporation factor β (0-1, soil moisture / capacity), `land`: land share.
    /// - `rain_rate` for fine cells is written to `self.precip` (kg/m²/s).
    pub fn solve_moisture(
        &mut self,
        month: usize,
        t_surface: &[f64],
        sea: &[f64],
        land: &[f64],
        soil_wetness: &[f64],
    ) {
        let n = self.n;
        // Evaporation (kg/m²/s): ρ C_E (|U| + gust) q_sat (1 - RH), with
        // q_sat = ρ_v,sat / ρ_air.
        for p in 0..n {
            let speed = (self.wind[p].0.powi(2) + self.wind[p].1.powi(2)).sqrt() + GUSTINESS;
            let rho_v = self.sat_density.at(t_surface[p]);
            // Over the sea the air just above the water stays near 80%
            // humid; over land it is as humid as the air column above
            // (well mixed by daytime turbulence), from last month.
            let sea_rate = EVAPORATION_COEF * speed * rho_v * (1.0 - SURFACE_RELATIVE_HUMIDITY);
            // Vapour is carried as a sea-level-equivalent column (what the
            // column would hold if its air reached down to sea level), so
            // moist lowland air spreading onto a plateau doesn't count the
            // part of its column that lies below the plateau. Its
            // saturation value uses the sea-level-equivalent temperature.
            self.wsat[p] = self
                .sat_density
                .at(t_surface[p] + LAPSE_RATE * self.coarse_height[p])
                * VAPOUR_SCALE_HEIGHT_M;
            let column_rh = (self.vapour[p] / self.wsat[p]).clamp(0.0, SURFACE_RELATIVE_HUMIDITY);
            // Wet ground can't evaporate faster than its energy supply
            // allows (`land_energy_limit`, set by the temperature model).
            let land_rate = (EVAPORATION_COEF * speed * rho_v * (1.0 - column_rh))
                .min(self.land_energy_limit[p]);
            self.potential_evaporation[p] = land_rate;
            self.sea_evaporation[p] = sea_rate;
            self.evaporation[p] = sea_rate * sea[p] + land_rate * land[p] * soil_wetness[p];
        }
        // Mean upslope wind of each coarse cell's fine cells (m/s per m of
        // column scale → used as w / H).
        // Lifted air only rains once it has risen past its condensation
        // level, about 125 m per °C of dew-point depression, which is about
        // (100 - RH%)/5 °C (Lawrence 2005, Bull. Amer. Meteor. Soc. 86,
        // 225-233). Only the part of a slope above that level counts, using
        // last month's humidity.
        for p in 0..n {
            let rh = (self.vapour[p] / self.wsat[p].max(1e-6)).min(1.0);
            self.lcl[p] = 125.0 * (100.0 - 100.0 * rh) / 5.0;
            self.uplift[p] = 0.0;
        }
        for i in 0..self.parent.len() {
            let p = self.parent[i] as usize;
            let h = self.fine_height[i];
            let lcl = self.lcl[p];
            if h <= lcl {
                self.fine_lift[i] = 0.0;
                continue;
            }
            let (u, v) = self.wind[p];
            let (sx, sy) = self.fine_slope[i];
            let above = 1.0 - lcl / h;
            self.fine_lift[i] = (u * sx + v * sy).max(0.0) * above;
            self.uplift[p] += self.fine_area[i] * self.fine_lift[i];
        }
        for p in 0..n {
            self.uplift[p] /= self.area[p];
        }
        // Steady state per cell p, per unit area (kg/m²/s):
        //   out_p W_p - Σ in_pq W_q + K Σ g (W_p - W_q) + P(W_p) = E_p
        // upwind fluxes; P(W) = exp(a (W/Wsat - b)) mm/day + W uplift / H.
        let k = self.tuning.moisture_mixing;
        let mm_day = 1.0 / 86_400.0; // kg/m²/s per mm/day
        self.moisture_sweeps = 0;
        // The month's coefficients: what each cell loses per kg/m² it holds
        // (`diag`: outflow and mixing) and gains per kg/m² its neighbours
        // hold (`self.coef`: inflow and mixing).
        for p in 0..n {
            let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
            // Air converging near the surface rises and its vapour condenses
            // (Kuo 1974, J. Atmos. Sci. 31, 1232-1240): rain at the rate the
            // convergence brings vapour in, W × convergence.
            let conv = -(s..e)
                .map(|j| self.surface_flow[j] * self.l_over_a[j])
                .sum::<f64>();
            // As with slopes, only air lifted past its condensation level
            // rains; converging air rises through the depth of the heated
            // layer, so shallow heat lows over dry land stay dry.
            let rh = (self.vapour[p] / self.wsat[p].max(1e-6)).min(1.0);
            let lcl = 125.0 * (100.0 - 100.0 * rh) / 5.0;
            let above = (1.0 - lcl / THERMAL_LOW_DEPTH_M).max(0.0);
            self.convergence[p] = conv.max(0.0) * above;
            let mut out_c = (self.uplift[p] / VAPOUR_SCALE_HEIGHT_M) + self.convergence[p];
            for j in s..e {
                let un = self.flow[j];
                let mix = k * self.g[j];
                out_c += un.max(0.0) * self.l_over_a[j] + mix;
                self.coef[j] = (-un).max(0.0) * self.l_over_a[j] + mix;
            }
            self.diag[p] = out_c;
            self.inv_ws[p] = BRETHERTON_A / self.wsat[p].max(1e-6);
        }
        let rain0 = (-BRETHERTON_A * BRETHERTON_B).exp() * mm_day;
        // Start from the same month last year, the closest guess.
        // Starting guess: last month's answer plus last year's change
        // between the same two months.
        if self.vapour_by_month[month].len() == n && self.prev_vapour.len() == n {
            for p in 0..n {
                self.vapour[p] = (self.vapour[p] + self.vapour_by_month[month][p]
                    - self.prev_vapour[p])
                    .max(1e-6);
            }
        }
        // Passes alternate direction (symmetric Gauss-Seidel), so moisture
        // carried either way along the grid order spreads quickly. Each pass
        // takes one Newton step per cell on
        //   diag W + rain(W) = E + Σ coef W_neighbour,
        // rain(W) = exp(a (W/Wsat - b)) mm/day.
        for sweep in 0..400 {
            self.moisture_sweeps += 1;
            let mut max_change = 0.0f64;
            for idx in 0..n {
                let p = if sweep % 2 == 0 { idx } else { n - 1 - idx };
                let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
                let mut rhs = self.evaporation[p];
                for j in s..e {
                    rhs += self.coef[j] * self.vapour[self.nbrs[j] as usize];
                }
                let a = self.inv_ws[p];
                let wmax = 1.5 * BRETHERTON_A / a;
                let wv = self.vapour[p].min(wmax);
                let pb = rain0 * (a * wv).exp();
                let f = self.diag[p] * wv + pb - rhs;
                let df = self.diag[p] + pb * a;
                // Never wetter than 1.5 × saturation (the rain law grows
                // exponentially beyond saturation; real columns stay below).
                let next = (wv - f / df).clamp(1e-6, wmax);
                max_change = max_change.max(f.abs());
                self.vapour[p] = next;
            }
            // Converged when every column's balance is out by less than
            // 0.3 mm/day (rain is 0-30 mm/day; the global balance of
            // evaporation and rain then holds to about 1%).
            if max_change < 0.3 * mm_day {
                break;
            }
        }
        std::mem::swap(&mut self.prev_vapour, &mut self.vapour_by_month[month]);
        self.vapour_by_month[month].clone_from(&self.vapour);
        // Rain per fine cell: humidity-driven part of its coarse cell plus
        // its own upslope part (which averages to the coarse cell's).
        for p in 0..n {
            let ws = self.wsat[p].max(1e-6);
            self.uplift[p] = (BRETHERTON_A * (self.vapour[p] / ws - BRETHERTON_B)).exp() * mm_day;
        }
        for i in 0..self.parent.len() {
            let p = self.parent[i] as usize;
            let wv = self.vapour[p];
            let pb = self.uplift[p];
            let lift = self.fine_lift[i];
            self.precip[i] = pb + wv * (lift / VAPOUR_SCALE_HEIGHT_M + self.convergence[p]);
        }
    }

    /// Evaporation from open sea in a coarse cell, kg/m²/s (mm/s).
    pub fn potential_evaporation_sea(&self, p: usize) -> f64 {
        self.sea_evaporation[p]
    }

    /// Column relative humidity, evaporation (mm/day) and convergence
    /// (1/day) of the coarse cell containing a fine cell.
    pub fn moisture_at_fine(&self, fine_cell: usize) -> (f64, f64, f64) {
        let p = self.parent[fine_cell] as usize;
        (
            self.vapour[p] / self.wsat[p].max(1e-6),
            self.evaporation[p] * 86_400.0,
            self.convergence[p] * 86_400.0,
        )
    }

    /// Column relative humidity of a coarse cell (0-1), from the last
    /// moisture solve.
    pub fn relative_humidity(&self, p: usize) -> f64 {
        (self.vapour[p] / self.wsat[p].max(1e-6)).clamp(0.0, 1.0)
    }

    /// Rain per coarse cell from the last moisture solve, kg/m²/s.
    pub fn coarse_rain(&self) -> Vec<f64> {
        let mut out = vec![0.0; self.n];
        for (i, &p) in self.parent.iter().enumerate() {
            out[p as usize] += self.precip[i] * self.fine_area[i] / self.area[p as usize];
        }
        out
    }

    /// Pressure anomaly (Pa) of the coarse cell containing a fine cell.
    pub fn pressure_at_fine(&self, fine_cell: usize) -> f64 {
        self.pressure[self.parent[fine_cell] as usize]
    }

    /// Every coarse cell's position (lat, lon in degrees) and wind (east,
    /// north; m/s), for drawing.
    pub fn winds(&self) -> Vec<(f64, f64, f64, f64)> {
        (0..self.n)
            .map(|p| {
                (
                    self.lat_deg[p],
                    self.lon_deg[p],
                    self.wind[p].0,
                    self.wind[p].1,
                )
            })
            .collect()
    }

    /// Wind of the coarse cell containing a fine cell.
    pub fn wind_at_fine(&self, fine_cell: usize) -> (f64, f64) {
        self.wind[self.parent[fine_cell] as usize]
    }

    /// Area-weighted global means: (evaporation, precipitation) in mm/day.
    pub fn global_balance(&self) -> (f64, f64) {
        let total: f64 = self.area.iter().sum();
        let e: f64 = (0..self.n)
            .map(|p| self.evaporation[p] * self.area[p])
            .sum::<f64>()
            / total;
        let p: f64 = (0..self.parent.len())
            .map(|i| self.precip[i] * self.fine_area[i])
            .sum::<f64>()
            / total;
        (e * 86_400.0, p * 86_400.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn saturation_vapour_values() {
        // Bolton: 611.2 Pa at 0 °C, about 3,170 Pa at 25 °C.
        assert!((saturation_vapour_pressure(0.0) - 611.2).abs() < 0.1);
        assert!((saturation_vapour_pressure(25.0) - 3168.0).abs() < 15.0);
        // Tropical column: ~50 kg/m² at 27 °C.
        let w = saturation_vapour_density(27.0) * VAPOUR_SCALE_HEIGHT_M;
        assert!((45.0..60.0).contains(&w), "{w}");
        let t = Table::new(saturation_vapour_density);
        assert!((t.at(12.34) - saturation_vapour_density(12.34)).abs() < 1e-5);
    }

    #[test]
    fn moisture_is_conserved_and_winds_are_plausible() {
        let fine = Grid::new(5);
        let heights = vec![0.0; fine.len()];
        let mut a = Atmosphere::new(&fine, &heights, 3, ATMOSPHERE_TUNED);
        let n = a.len();
        // Idealized: warm tropics, cold poles, all sea.
        let t: Vec<f64> = a
            .lat_deg
            .iter()
            .map(|l| 28.0 - 45.0 * (l.to_radians().sin()).powi(2))
            .collect();
        a.update_winds(0, &t);
        let ones = vec![1.0; n];
        let zeros = vec![0.0; n];
        a.solve_moisture(0, &t, &ones, &zeros, &zeros);
        let (e, p) = a.global_balance();
        // Steady state: global rain equals global evaporation (within the
        // solver tolerance), and is of Earth's order (~2.7 mm/day).
        assert!((p - e).abs() / e < 0.05, "E {e:.2} P {p:.2}");
        assert!((1.0..6.0).contains(&e), "E {e:.2}");
        // Trade winds blow from the east near 15°.
        let i = (0..n)
            .min_by(|&x, &y| {
                (a.lat_deg[x] - 15.0)
                    .abs()
                    .total_cmp(&(a.lat_deg[y] - 15.0).abs())
            })
            .unwrap();
        assert!(a.wind[i].0 < 0.0, "trade wind {:?}", a.wind[i]);
        // Westerlies near 45°.
        let i = (0..n)
            .min_by(|&x, &y| {
                (a.lat_deg[x] - 45.0)
                    .abs()
                    .total_cmp(&(a.lat_deg[y] - 45.0).abs())
            })
            .unwrap();
        assert!(a.wind[i].0 > 0.0, "westerly {:?}", a.wind[i]);
    }
}
