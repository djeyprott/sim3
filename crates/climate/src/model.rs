//! The monthly temperature model (Phase 1, Step 2): a two-dimensional energy
//! balance on the 40,962-cell climate grid.
//!
//! Each month every cell warms by the sunlight it absorbs and cools by the
//! heat it radiates to space; heat stored in land (thin) or sea (the mixed
//! top layer) carries over between months; large-scale transport by the
//! atmosphere and ocean moves heat from warm to cold regions (see
//! `solver.rs`).
//!
//! The balance is solved for sea-level-equivalent temperature: heat leaves
//! for space from the atmosphere above, whose temperature follows the
//! column, not the ground. A place's surface temperature is that value minus
//! 6.5 °C per km of its height. (This is the usual treatment where energy
//! balance models meet high ice sheets.) Map cells (28 km) use their own
//! height.
//!
//! Not yet included, and added in later steps: rain and seasonal snow
//! (Step 3), ocean currents and proper sea ice (Step 4), random variability
//! (Step 5), changing ice sheets and sea level (Step 6).

use crate::orbit::{month_days, MonthlySunshine, Orbit, SOLAR_CONSTANT};
use crate::solver::Solver;
use world::data::Epoch;
use world::grid::{to_lat_lon, Grid, EARTH_RADIUS_M};

// ---- Physical constants, each with its source ----

/// Outgoing longwave radiation at the top of the atmosphere as a linear
/// function of surface temperature, OLR = A + B·T (T in °C), fitted to
/// satellite data (North, Cahalan & Coakley 1981, Rev. Geophys. 19, 91-121).
pub const OLR_A: f64 = 203.3; // W/m²
pub const OLR_B: f64 = 2.09; // W/m²/K
/// CO2 level of the satellite period those values describe (Mauna Loa,
/// mid-1970s, NOAA GML).
pub const OLR_CO2_PPM: f64 = 330.0;
/// Radiative forcing of CO2: 5.35 ln(C/C0) W/m² (Myhre et al. 1998,
/// Geophys. Res. Lett. 25, 2715-2718).
pub const CO2_FORCING_COEF: f64 = 5.35;
/// Heat capacity of the whole atmospheric column, cp·ps/g
/// = 1004 J/kg/K × 101,325 Pa / 9.81 m/s² (standard values). Land stores
/// little heat beyond the air above it over a month.
pub const HEAT_CAPACITY_AIR: f64 = 1004.0 * 101_325.0 / 9.81; // J/m²/K
/// Heat capacity of the ocean's mixed surface layer: seawater density
/// 1025 kg/m³ × specific heat 3985 J/kg/K × 50 m, a typical global mixed-
/// layer depth (de Boyer Montégut et al. 2004, J. Geophys. Res. 109, C12003).
pub const HEAT_CAPACITY_SEA: f64 = 1025.0 * 3985.0 * 50.0; // J/m²/K
/// Until Step 4 adds sea ice properly, sea is treated as ice-covered
/// (bright, storing little heat) where the annual mean temperature is below
/// -10 °C: the classic ice-edge rule of energy balance models (Budyko 1969,
/// Tellus 21, 611-619; North 1975, J. Atmos. Sci. 32, 2033-2043). Using the
/// annual mean avoids flipping between ice and open water month by month.
pub const SEA_ICE_ANNUAL_MEAN_C: f64 = -10.0;
/// Temperature falls 6.5 °C per km of height (ICAO standard atmosphere).
pub const LAPSE_RATE: f64 = 6.5e-3; // K/m

/// Values that a model this simple cannot derive and that are tuned once to
/// today's observations (see `tools/tuning` notes in check.rs). Never tuned
/// to the past.
#[derive(Clone, Copy, Debug)]
pub struct Tuning {
    /// Planetary albedo (share of sunlight reflected, mostly by clouds)
    /// without ice: a0 + a2·P2(sin latitude).
    pub albedo_a0: f64,
    pub albedo_a2: f64,
    /// Planetary albedo over ice and snow.
    pub albedo_ice: f64,
    /// Heat transport strength D on the unit sphere, W/m²/K.
    pub diffusion: f64,
}

/// Tuned against today's world: see check.rs (`tune` test) for the
/// procedure and results.
pub const TUNED: Tuning = Tuning {
    albedo_a0: 0.30,
    albedo_a2: 0.10,
    albedo_ice: 0.70,
    diffusion: 0.5,
};

/// Length of the averaging used for "recent climate": an exponentially
/// weighted average with a 30-year time constant (the standard length of a
/// climate normal).
const CLIMATE_NORMAL_YEARS: f64 = 30.0;

pub struct Climate {
    tuning: Tuning,
    n: usize,
    /// Ice-free planetary albedo per cell (depends only on latitude).
    free_albedo: Vec<f64>,
    /// Row of the sunshine table for each cell.
    sun_row: Vec<usize>,
    /// Whether each cell's sea is ice-covered this year.
    frozen: Vec<bool>,
    buf_c: Vec<f64>,
    buf_b: Vec<f64>,
    area: Vec<f64>,
    land: Vec<f64>,
    ice: Vec<f64>,
    /// Current temperature per climate cell, °C (this month's mean).
    t: Vec<f64>,
    transport: Vec<f64>,
    solver: Solver,
    sunshine: MonthlySunshine,
    /// Monthly means of the current year, °C.
    this_year: Vec<[f32; 12]>,
    /// Recent-climate monthly means (30-year weighted), °C.
    normal: Vec<[f32; 12]>,
    years: u32,
    /// Area-weighted global mean temperature of each completed year.
    pub global_history: Vec<(i64, f32)>,
    /// Optional averaging window (first year, last year) and its sums.
    window: Option<(i64, i64)>,
    window_sum: Vec<[f64; 12]>,
    window_years: u32,
    // Map downscaling.
    map_parent: Vec<u32>,
    /// Surface height of each map cell (sea counts as 0), m.
    map_height: Vec<f32>,
    /// Mean surface height of each climate cell, m.
    height: Vec<f64>,
}

impl Climate {
    /// `map_parent[i]` is the climate cell containing map cell i.
    pub fn new(
        climate_grid: &Grid,
        map_grid: &Grid,
        map_parent: Vec<u32>,
        epoch: &Epoch,
        orbit: &Orbit,
        tuning: Tuning,
    ) -> Climate {
        let n = climate_grid.len();
        let mut w = vec![0.0; n];
        let mut land = vec![0.0; n];
        let mut ice = vec![0.0; n];
        let mut elev = vec![0.0; n];
        for (i, &p) in map_parent.iter().enumerate() {
            let p = p as usize;
            let a = map_grid.areas_m2[i];
            w[p] += a;
            land[p] += a * epoch.land_share[i] as f64 / 255.0;
            ice[p] += if epoch.ice_thickness_m[i] > 0 { a } else { 0.0 };
            elev[p] += a * (epoch.elev_mean_m[i] as f64).max(0.0);
        }
        for p in 0..n {
            land[p] /= w[p];
            ice[p] /= w[p];
            elev[p] /= w[p];
        }
        let map_height = epoch
            .elev_mean_m
            .iter()
            .map(|&e| (e as f32).max(0.0))
            .collect();
        let lat_deg: Vec<f64> = climate_grid
            .centers
            .iter()
            .map(|&c| to_lat_lon(c).0)
            .collect();
        // Start from a rough equilibrium so spin-up is short.
        let sunshine = MonthlySunshine::new(orbit, SOLAR_CONSTANT, &month_days());
        let t_start: Vec<f64> = lat_deg
            .iter()
            .map(|&l| {
                let q: f64 = (0..12).map(|m| sunshine.at(l, m) as f64).sum::<f64>() / 12.0;
                ((1.0 - 0.3) * q - OLR_A) / OLR_B
            })
            .collect();
        let free_albedo = lat_deg
            .iter()
            .map(|&l| {
                let x = l.to_radians().sin();
                (tuning.albedo_a0 + tuning.albedo_a2 * 0.5 * (3.0 * x * x - 1.0)).clamp(0.0, 0.9)
            })
            .collect();
        let sun_row = lat_deg.iter().map(|&l| MonthlySunshine::row(l)).collect();
        let frozen = t_start.iter().map(|&v| v < SEA_ICE_ANNUAL_MEAN_C).collect();
        Climate {
            tuning,
            n,
            free_albedo,
            sun_row,
            frozen,
            buf_c: vec![0.0; n],
            buf_b: vec![0.0; n],
            area: climate_grid.areas_m2.clone(),
            land,
            ice,
            t: t_start,
            transport: vec![0.0; n],
            solver: Solver::new(climate_grid),
            sunshine,
            this_year: vec![[0.0; 12]; n],
            normal: vec![[0.0; 12]; n],
            years: 0,
            global_history: Vec::new(),
            window: None,
            window_sum: vec![[0.0; 12]; n],
            window_years: 0,
            map_parent,
            map_height,
            height: elev,
        }
    }

    /// Updates the sunshine table for a new orbit (call every century or so;
    /// the orbit changes by well under 1% per century).
    pub fn set_orbit(&mut self, orbit: &Orbit) {
        self.sunshine = MonthlySunshine::new(orbit, SOLAR_CONSTANT, &month_days());
    }

    /// Averages each month over the years `first..=last` (inclusive) for
    /// comparison with observations.
    pub fn set_window(&mut self, first: i64, last: i64) {
        self.window = Some((first, last));
    }

    fn albedo(&self, i: usize) -> f64 {
        let sea_ice = if self.frozen[i] {
            1.0 - self.land[i]
        } else {
            0.0
        };
        let icy = (self.ice[i] + sea_ice * (1.0 - self.ice[i])).min(1.0);
        icy * self.tuning.albedo_ice + (1.0 - icy) * self.free_albedo[i]
    }

    fn heat_capacity(&self, i: usize) -> f64 {
        let sea = (1.0 - self.land[i]) * (1.0 - self.ice[i]);
        let sea_c = if self.frozen[i] {
            HEAT_CAPACITY_AIR
        } else {
            HEAT_CAPACITY_SEA
        };
        sea * sea_c + (1.0 - sea) * HEAT_CAPACITY_AIR
    }

    /// Advances one month (0 = January) of `days` days at CO2 `co2_ppm`.
    pub fn step_month(&mut self, year: i64, month: usize, days: u32, co2_ppm: f64) {
        let dt = days as f64 * 86_400.0;
        let a_eff = OLR_A - CO2_FORCING_COEF * (co2_ppm / OLR_CO2_PPM).ln();
        let mut c = std::mem::take(&mut self.buf_c);
        let mut b = std::mem::take(&mut self.buf_b);
        let this = &*self;
        c.iter_mut()
            .zip(b.iter_mut())
            .enumerate()
            .for_each(|(i, (ci, bi))| {
                let cap = this.heat_capacity(i) / dt;
                let absorbed =
                    (1.0 - this.albedo(i)) * this.sunshine.table[this.sun_row[i]][month] as f64;
                *ci = cap + OLR_B;
                *bi = cap * this.t[i] + absorbed - a_eff;
            });
        let k = self.tuning.diffusion * EARTH_RADIUS_M * EARTH_RADIUS_M;
        let mut t = std::mem::take(&mut self.t);
        self.solver
            .solve(month, &c, &b, k, &mut t, &mut self.transport);
        self.t = t;
        self.buf_c = c;
        self.buf_b = b;

        for i in 0..self.n {
            self.this_year[i][month] = self.t[i] as f32;
        }
        if let Some((first, last)) = self.window {
            if (first..=last).contains(&year) {
                for i in 0..self.n {
                    self.window_sum[i][month] += self.t[i];
                }
                if month == 11 {
                    self.window_years += 1;
                }
            }
        }
        if month == 11 {
            self.finish_year(year);
        }
    }

    fn finish_year(&mut self, year: i64) {
        let weight = if self.years == 0 {
            1.0
        } else {
            1.0 / CLIMATE_NORMAL_YEARS
        } as f32;
        let mut sum = 0.0;
        let mut area = 0.0;
        for i in 0..self.n {
            for m in 0..12 {
                self.normal[i][m] += weight * (self.this_year[i][m] - self.normal[i][m]);
            }
            let annual: f32 = self.this_year[i].iter().sum::<f32>() / 12.0;
            sum += (annual as f64 - LAPSE_RATE * self.height[i]) * self.area[i];
            area += self.area[i];
        }
        self.global_history.push((year, (sum / area) as f32));
        self.years += 1;
        for i in 0..self.n {
            let annual = self.normal[i].iter().map(|&v| v as f64).sum::<f64>() / 12.0;
            self.frozen[i] = annual < SEA_ICE_ANNUAL_MEAN_C;
        }
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Surface temperature (°C) of a map cell from its climate cell's
    /// sea-level-equivalent value.
    fn at_map(&self, map_cell: usize, sea_level_value: f32) -> f32 {
        sea_level_value - (LAPSE_RATE as f32) * self.map_height[map_cell]
    }

    /// Per map cell: this month's temperature, °C.
    pub fn map_current(&self) -> Vec<f32> {
        (0..self.map_parent.len())
            .map(|i| self.at_map(i, self.t[self.map_parent[i] as usize] as f32))
            .collect()
    }

    /// Per map cell: recent-climate temperature for `month`, or the annual
    /// mean when `month` is None, °C.
    pub fn map_normal(&self, month: Option<usize>) -> Vec<f32> {
        (0..self.map_parent.len())
            .map(|i| {
                let row = &self.normal[self.map_parent[i] as usize];
                let v = match month {
                    Some(m) => row[m],
                    None => row.iter().sum::<f32>() / 12.0,
                };
                self.at_map(i, v)
            })
            .collect()
    }

    /// Recent-climate monthly temperatures of one map cell, °C.
    pub fn cell_normal(&self, map_cell: usize) -> [f32; 12] {
        let row = self.normal[self.map_parent[map_cell] as usize];
        row.map(|v| self.at_map(map_cell, v))
    }

    /// Monthly means over the averaging window for one map cell, if the
    /// window has completed at least one year.
    pub fn cell_window(&self, map_cell: usize) -> Option<[f32; 12]> {
        if self.window_years == 0 {
            return None;
        }
        let row = self.window_sum[self.map_parent[map_cell] as usize];
        Some(row.map(|s| self.at_map(map_cell, (s / self.window_years as f64) as f32)))
    }

    /// Years averaged so far in the window.
    pub fn window_years(&self) -> u32 {
        self.window_years
    }

    /// Global area-weighted mean over the window, °C.
    pub fn window_global_mean(&self) -> Option<f64> {
        if self.window_years == 0 {
            return None;
        }
        let total: f64 = self.area.iter().sum();
        Some(
            (0..self.n)
                .map(|i| {
                    let sea_level =
                        self.window_sum[i].iter().sum::<f64>() / 12.0 / self.window_years as f64;
                    (sea_level - LAPSE_RATE * self.height[i]) * self.area[i]
                })
                .sum::<f64>()
                / total,
        )
    }

    /// Passes the transport solver needed last month (for profiling).
    pub fn solver_passes(&self) -> usize {
        self.solver.last_sweeps
    }

    pub fn years_simulated(&self) -> u32 {
        self.years
    }
}
