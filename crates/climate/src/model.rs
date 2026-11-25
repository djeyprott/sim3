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
//! Winds and rain (Step 3) come from `atmosphere.rs`; soil water and
//! lying snow (kept as water until it melts) are tracked here.
//!
//! Ocean currents and sea ice (Step 4) are in `ocean.rs` and below; eruptions and the El Niño-like
//! cycle (Step 5) in `variability.rs`. Not yet included: changing ice sheets and sea level (Step 6).

use crate::atmosphere::{Atmosphere, AtmosphereTuning, Table};
use crate::orbit::{month_days, MonthlySunshine, Orbit, SOLAR_CONSTANT};
use crate::solver::{Solver, TRANSPORT_LEVEL};
use crate::variability::{
    Enso, EnsoEvent, Eruption, EruptionStats, Rng, Volcanism, ENSO_LAT_WIDTH_DEG,
};
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
/// Sea ice, as a single slab (the "zero-layer" model of Semtner 1976, J.
/// Phys. Oceanogr. 6, 379-389): open sea that cools to the freezing point
/// freezes instead of cooling further; ice thickens as heat is conducted
/// up through it from the water below (which stays at freezing) and melts
/// when its surface warms past 0 °C. While ice covers the sea it reflects
/// sunlight like other ice and only the air above stores heat.
///
/// Freezing point of seawater at salinity 35, °C.
pub const SEA_FREEZING_C: f64 = -1.8;
/// Thermal conductivity of sea ice, W/m/K (Maykut & Untersteiner 1971,
/// J. Geophys. Res. 76, 1550-1575).
pub const ICE_CONDUCTIVITY: f64 = 2.03;
/// Density of ice (kg/m³) and latent heat of fusion (J/kg).
pub const ICE_DENSITY: f64 = 917.0;
pub const LATENT_HEAT_FUSION: f64 = 3.34e5;
/// Thinnest ice used for conduction, m: new ice is a few cm to tens of cm
/// thick within days, and conduction through thinner ice is limited by
/// the air above it.
pub const MIN_CONDUCTING_ICE_M: f64 = 0.1;
// ---- Two layers: the ground or sea surface, and the air above ----
//
// Each place has a surface temperature (ground, ice, or the sea's mixed
// layer) and the air above it has its own (solved on the 440 km grid, where
// the winds and large-scale mixing carry it). They trade heat by turbulence
// and by longwave radiation; some of the surface's radiation escapes
// straight to space through the "window" where air doesn't absorb. On
// clear, dry winter nights that loss lets frozen ground sit far below the
// air above it; over the sea the two stay close.

/// Effective heat capacity of ground over a seasonal cycle, J/m²/K: soil
/// volumetric heat capacity ~2.0 MJ/m³/K times the depth the yearly
/// temperature wave reaches (√(2κ/ω) ≈ 2.6 m for κ = 0.7 mm²/s), over √2
/// for a sinusoid (Hartmann 2016, Global Physical Climatology, ch. 4).
pub const HEAT_CAPACITY_GROUND: f64 = 2.0e6 * 2.6 / std::f64::consts::SQRT_2;
/// Share of absorbed sunlight taken up by the air itself (water vapour,
/// clouds, ozone): 78 of 239 W/m² (Trenberth, Fasullo & Kiehl 2009, Bull.
/// Amer. Meteor. Soc. 90, 311-323).
pub const AIR_SUNLIGHT_SHARE: f64 = 78.0 / 239.0;
/// Surface radiation escaping through the atmospheric window in today's
/// world on average, W/m² (Trenberth et al. 2009). The air's own outgoing
/// radiation is the satellite fit minus this, so the total matches.
pub const WINDOW_REFERENCE: f64 = 40.0;
/// Stefan-Boltzmann constant, W/m²/K⁴.
pub const STEFAN_BOLTZMANN: f64 = 5.670_374e-8;
/// Bulk transfer coefficient for heat between surface and air, 1.2e-3
/// (Large & Pond 1982), and the air's density and specific heat.
pub const HEAT_TRANSFER_COEF: f64 = 1.2e-3;
pub const AIR_DENSITY_SURFACE: f64 = 1.225;
pub const AIR_SPECIFIC_HEAT: f64 = 1004.0;
/// Latent heat of vaporisation, J/kg.
pub const LATENT_HEAT_VAPOUR: f64 = 2.5e6;
/// Priestley & Taylor (1972, Mon. Wea. Rev. 100, 81-92): wet ground
/// evaporates at most 1.26 × Δ/(Δ + γ) of its net radiation as latent
/// heat, Δ the slope of saturation vapour pressure and γ = 66 Pa/K the
/// psychrometric constant.
pub const PRIESTLEY_TAYLOR: f64 = 1.26;
pub const PSYCHROMETRIC_PA_PER_K: f64 = 66.0;

/// Priestley-Taylor share of net radiation that wet ground at `t_c` (°C)
/// can turn into evaporation: 1.26 Δ/(Δ + γ).
fn priestley_taylor_factor(t_c: f64) -> f64 {
    let es = crate::atmosphere::saturation_vapour_pressure(t_c);
    let slope = es * 17.67 * 243.5 / (t_c + 243.5).powi(2);
    PRIESTLEY_TAYLOR * slope / (slope + PSYCHROMETRIC_PA_PER_K)
}

/// Height at which the surface layer's stability is judged, m (the
/// standard 10 m of bulk formulas).
pub const STABILITY_HEIGHT_M: f64 = 10.0;

/// Clear-sky emissivity of the air seen from the ground, from the
/// near-surface vapour pressure e (hPa): 0.52 + 0.065 √e (Brunt 1932, Q. J.
/// R. Meteorol. Soc. 58, 389-420). Dry winter air lets much more of the
/// ground's heat escape than humid tropical air.
///
/// Clouds close the window almost completely, so only the clear part of
/// the sky lets radiation out.
fn air_emissivity(vapour_pressure_hpa: f64) -> f64 {
    let clear = (0.52 + 0.065 * vapour_pressure_hpa.max(0.0).sqrt()).min(1.0);
    clear + (1.0 - clear) * CLOUD_COVER
}

/// Global mean cloud cover, 0.67 (ISCCP; Rossow & Schiffer 1999, Bull.
/// Amer. Meteor. Soc. 80, 2261-2287). With it the world's mean window
/// comes to ~35 W/m², close to the 40 W/m² observed.
pub const CLOUD_COVER: f64 = 0.67;

/// Reduction of turbulent exchange when air is warmer than the ground
/// below (a stable layer damps turbulence): Louis (1979, Boundary-Layer
/// Meteorol. 17, 187-202), 1 / (1 + 10 Ri / √(1 + 5 Ri)) for bulk
/// Richardson number Ri > 0; no change otherwise.
fn stability_factor(air_minus_surface: f64, speed: f64, t_kelvin: f64) -> f64 {
    let ri = 9.81 * STABILITY_HEIGHT_M * air_minus_surface / (t_kelvin * speed * speed);
    if ri <= 0.0 {
        1.0
    } else {
        1.0 / (1.0 + 10.0 * ri / (1.0 + 5.0 * ri).sqrt())
    }
}

/// Temperature falls 6.5 °C per km of height (ICAO standard atmosphere).
pub const LAPSE_RATE: f64 = 6.5e-3; // K/m
/// Water the soil can hold for plants and evaporation: 150 mm, the "bucket"
/// of Manabe (1969, Mon. Wea. Rev. 97, 739-774). Evaporation runs at the
/// full rate while the bucket is at least 75% full and in proportion below.
pub const SOIL_CAPACITY_MM: f64 = 150.0;
pub const SOIL_FULL_RATE_SHARE: f64 = 0.75;
/// Snowmelt per degree-day above 0 °C, mm of water per °C per day
/// (Braithwaite 1995, J. Glaciol. 41, 153-160: 3-5 for snow; low end).
pub const DEGREE_DAY_FACTOR: f64 = 3.0;
/// Spread (standard deviation) of daily temperatures around the monthly
/// mean, used for melt and for splitting rain from snow: 5 °C (Reeh 1991,
/// Polarforschung 59, 113-128; the usual value in ice-sheet models).
pub const DAILY_SPREAD_C: f64 = 5.0;

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
    albedo_a0: 0.33,
    albedo_a2: 0.075,
    albedo_ice: 0.84,
    diffusion: 0.2,
};

/// Length of the averaging used for "recent climate": an exponentially
/// weighted average with a 30-year time constant (the standard length of a
/// climate normal).
const CLIMATE_NORMAL_YEARS: f64 = 30.0;

/// Complementary error function (Abramowitz & Stegun 1964, 7.1.26;
/// error below 1.5e-7).
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * z);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    let r = poly * (-z * z).exp();
    if x >= 0.0 {
        r
    } else {
        2.0 - r
    }
}

/// Mean degrees above 0 °C per day in a month of mean temperature `t`
/// (°C), with daily temperatures spread normally by `DAILY_SPREAD_C`.
fn positive_degrees(t: f64) -> f64 {
    let s = DAILY_SPREAD_C;
    s / (2.0 * std::f64::consts::PI).sqrt() * (-t * t / (2.0 * s * s)).exp()
        + 0.5 * t * erfc(-t / (s * std::f64::consts::SQRT_2))
}

/// A place and its rain as a share of normal: (lat, lon, %).
pub type RainPlace = (f64, f64, f64);

/// What the monthly temperature step needs to know about a climate cell
/// that doesn't change during a run, packed together so the step reads one
/// record per cell instead of many separate lists.
#[derive(Clone, Copy, Default)]
struct CellFacts {
    /// The 440 km cell it belongs to, and its share of that cell's area.
    parent: u32,
    weight: f64,
    sun_row: u32,
    height: f64,
    free_albedo: f64,
    /// Ground, ice sheet and sea shares.
    shares: [f64; 3],
}

#[derive(Clone)]
pub struct Climate {
    tuning: Tuning,
    n: usize,
    facts: Vec<CellFacts>,
    /// Share of each 440 km cell that is sea (open or iced).
    coarse_sea_share: Vec<f64>,
    /// Ice-free planetary albedo per cell (depends only on latitude).
    free_albedo: Vec<f64>,
    /// Row of the sunshine table for each cell.
    sun_row: Vec<usize>,
    /// Sea ice thickness over each cell's sea, m (0 = open water).
    sea_ice_m: Vec<f64>,
    /// Centre of each cell is north of the equator.
    northern: Vec<bool>,
    /// Sea ice extent (km², north and south) per month, summed over the
    /// averaging window's years.
    window_ice_extent: [[f64; 2]; 12],
    buf_c: Vec<f64>,
    buf_b: Vec<f64>,
    area: Vec<f64>,
    land: Vec<f64>,
    ice: Vec<f64>,
    /// Current temperature per climate cell, °C (this month's mean).
    t: Vec<f64>,
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
    pub atmosphere: Atmosphere,
    pub ocean: crate::ocean::Ocean,
    /// Year-to-year variability (Step 5): eruptions and the El Niño-like
    /// cycle, their random numbers, and what they do this month.
    pub volcanism: Volcanism,
    pub enso: Enso,
    rng: Rng,
    /// Share of sunlight the volcanic veil takes away, per 440 km cell.
    dimming: Vec<f64>,
    /// El Niño-like sea surface anomaly, per 440 km cell, °C.
    enso_offset: Vec<f64>,
    /// What happened this month, for the event log.
    pub eruptions: Vec<Eruption>,
    pub enso_events: Vec<EnsoEvent>,
    /// The widest basin's anomaly each month since the start, °C.
    pub enso_history: Vec<f32>,
    /// Each climate cell's rain this month and the two before, and the
    /// normal for those months, mm (where strong El Niño years bring
    /// drought and floods).
    rain_recent: [Vec<f32>; 3],
    normal_recent: [Vec<f32>; 3],
    /// Climate cell centres (lat, lon), degrees.
    cell_lat_lon: Vec<(f64, f64)>,
    /// Largest height step from each climate cell to a neighbour, m.
    relief: Vec<f64>,
    /// The month (0-11) most recently stepped.
    pub last_month: usize,
    /// During the comparison window: the widest basin's monthly anomaly,
    /// and land rain and its normal at the basin's western and eastern
    /// ends in El Niño months (anomaly at least +0.5 °C), mm.
    pub window_enso: Vec<f32>,
    enso_ends: Option<[Vec<usize>; 2]>,
    window_enso_rain: [[f64; 2]; 2],
    /// Air temperature per 440 km cell (sea-level equivalent), °C, and the
    /// open sea's mixed-layer temperature there, °C.
    pub air: Vec<f64>,
    pub sea_t: Vec<f64>,
    /// Last solution for each calendar month (starting guesses).
    air_by_month: Vec<Vec<f64>>,
    sea_by_month: Vec<Vec<f64>>,
    prev_air: Vec<f64>,
    prev_sea: Vec<f64>,
    /// Surface temperatures of each climate cell's ground, ice sheet and
    /// sea ice (sea-level equivalent), °C.
    t_ground: Vec<f64>,
    t_glacier: Vec<f64>,
    t_seaice: Vec<f64>,
    /// Passes of the last temperature solve.
    last_passes: usize,
    /// Heat brought by the overturning, W/m² of each cell.
    overturning_heat: Vec<f64>,
    /// This year's rain minus evaporation over each cell's sea so far, m.
    fresh_m: Vec<f64>,
    /// Tabulated share of a month's precipitation falling as snow, and
    /// mean daily degrees above freezing, by monthly mean temperature.
    snow_share: Table,
    pt_factor: Table,
    positive_degrees: Table,
    /// Seconds spent so far in winds, temperature and water (profiling).
    pub timing: [f64; 3],
    /// Water in the soil (mm) and lying snow (mm of water) per climate cell.
    soil_mm: Vec<f64>,
    snow_mm: Vec<f64>,
    /// Rain and snowfall this month, mm, per climate cell.
    precip_mm: Vec<f64>,
    this_year_p: Vec<[f32; 12]>,
    normal_p: Vec<[f32; 12]>,
    normal_snow: Vec<[f32; 12]>,
    this_year_snow: Vec<[f32; 12]>,
    window_sum_p: Vec<[f64; 12]>,
    // Coarse-grid work arrays.
    co_t: Vec<f64>,
    co_sea: Vec<f64>,
    co_land: Vec<f64>,
    co_beta: Vec<f64>,
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
        atmosphere_tuning: AtmosphereTuning,
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
        let sea_ice_m = t_start
            .iter()
            .map(|&v| if v < SEA_FREEZING_C { 1.0 } else { 0.0 })
            .collect();
        let atmosphere = Atmosphere::new(climate_grid, &elev, TRANSPORT_LEVEL, atmosphere_tuning);
        let co_start = {
            let mut v = vec![0.0; atmosphere.len()];
            for (i, &p) in atmosphere.parent.iter().enumerate() {
                v[p as usize] +=
                    t_start[i] * climate_grid.areas_m2[i] / atmosphere.area[p as usize];
            }
            v
        };
        let ocean =
            crate::ocean::Ocean::new(climate_grid, &Grid::new(TRANSPORT_LEVEL), &atmosphere.area);
        let nc = atmosphere.len();
        let mut c = Climate {
            atmosphere,
            air: co_start.clone(),
            sea_t: co_start,
            air_by_month: vec![Vec::new(); 12],
            sea_by_month: vec![Vec::new(); 12],
            prev_air: Vec::new(),
            prev_sea: Vec::new(),
            t_ground: t_start.clone(),
            t_glacier: t_start.clone(),
            t_seaice: t_start.clone(),
            last_passes: 0,
            ocean,
            volcanism: Volcanism::new(None, Vec::new()),
            enso: Enso::default(),
            rng: Rng::new(0),
            dimming: vec![0.0; nc],
            enso_offset: vec![0.0; nc],
            eruptions: Vec::new(),
            enso_events: Vec::new(),
            enso_history: Vec::new(),
            rain_recent: [vec![0.0; n], vec![0.0; n], vec![0.0; n]],
            normal_recent: [vec![0.0; n], vec![0.0; n], vec![0.0; n]],
            cell_lat_lon: climate_grid
                .centers
                .iter()
                .map(|&c| to_lat_lon(c))
                .collect(),
            window_enso: Vec::new(),
            enso_ends: None,
            window_enso_rain: [[0.0; 2]; 2],
            last_month: 0,
            relief: (0..n)
                .map(|i| {
                    climate_grid
                        .neighbors(i)
                        .iter()
                        .map(|&q| (elev[i] - elev[q as usize]).abs())
                        .fold(0.0, f64::max)
                })
                .collect(),
            overturning_heat: vec![0.0; n],
            fresh_m: vec![0.0; n],
            pt_factor: Table::new(priestley_taylor_factor),
            snow_share: Table::new(|t| 0.5 * erfc(t / (DAILY_SPREAD_C * std::f64::consts::SQRT_2))),
            positive_degrees: Table::new(positive_degrees),
            timing: [0.0; 3],
            soil_mm: vec![SOIL_CAPACITY_MM * SOIL_FULL_RATE_SHARE; n],
            snow_mm: vec![0.0; n],
            precip_mm: vec![0.0; n],
            this_year_p: vec![[0.0; 12]; n],
            normal_p: vec![[0.0; 12]; n],
            normal_snow: vec![[0.0; 12]; n],
            this_year_snow: vec![[0.0; 12]; n],
            window_sum_p: vec![[0.0; 12]; n],
            co_t: vec![0.0; nc],
            co_sea: vec![0.0; nc],
            co_land: vec![0.0; nc],
            co_beta: vec![0.0; nc],
            tuning,
            n,
            facts: Vec::new(),
            coarse_sea_share: Vec::new(),
            free_albedo,
            sun_row,
            sea_ice_m,
            northern: lat_deg.iter().map(|&l| l >= 0.0).collect(),
            window_ice_extent: [[0.0; 2]; 12],
            buf_c: vec![0.0; n],
            buf_b: vec![0.0; n],
            area: climate_grid.areas_m2.clone(),
            land,
            ice,
            t: t_start,
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
        };
        c.rebuild_facts();
        c
    }

    /// Repacks the per-cell facts (call again if land or ice sheets change).
    fn rebuild_facts(&mut self) {
        self.facts = (0..self.n)
            .map(|i| {
                let p = self.atmosphere.parent[i];
                CellFacts {
                    parent: p,
                    weight: self.area[i] / self.atmosphere.area[p as usize],
                    sun_row: self.sun_row[i] as u32,
                    height: self.height[i],
                    free_albedo: self.free_albedo[i],
                    shares: self.tiles(i),
                }
            })
            .collect();
        let mut share = vec![0.0; self.atmosphere.len()];
        for f in &self.facts {
            share[f.parent as usize] += f.weight * f.shares[2];
        }
        self.coarse_sea_share = share;
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

    /// Share of the cell that is open sea or sea ice (not land, not an ice
    /// sheet).
    /// This month's eruptions and El Niño-like cycle: the veil's dimming and
    /// the sea surface anomaly, laid on before the temperatures are solved.
    fn variability_month(&mut self) {
        let lat_lon = &self.cell_lat_lon;
        self.eruptions = self.volcanism.draw(&mut self.rng, |c| lat_lon[c as usize]);
        self.volcanism.step();
        self.enso_events = self.enso.step(&mut self.rng);
        if let Some(b) = self.enso.main_basin() {
            self.enso_history.push(self.enso.basins[b].t as f32);
        }
        for p in 0..self.dimming.len() {
            let (lat, lon) = self.atmosphere.lat_lon(p);
            self.dimming[p] = self.volcanism.dimming(lat);
            self.enso_offset[p] = self.enso.anomaly(lat, lon);
            self.sea_t[p] += self.enso_offset[p];
        }
    }

    /// Comparison-window record of the El Niño-like cycle and its rain.
    fn window_enso_month(&mut self, month: usize) {
        let Some(b) = self.enso.main_basin() else {
            return;
        };
        let t = self.enso.basins[b].t;
        self.window_enso.push(t as f32);
        if self.enso_ends.is_none() {
            // Land within one anomaly width (latitude and longitude) of the
            // basin's western and eastern ends on the equator.
            let (w, e) = (self.enso.basins[b].west_lon, self.enso.basins[b].east_lon);
            let near = |lon: f64, end: f64, west: bool| {
                let d = (lon - end + 540.0).rem_euclid(360.0) - 180.0;
                if west {
                    (-ENSO_LAT_WIDTH_DEG..=0.0).contains(&d)
                } else {
                    (0.0..=ENSO_LAT_WIDTH_DEG).contains(&d)
                }
            };
            let pick = |west: bool, end: f64| -> Vec<usize> {
                (0..self.n)
                    .filter(|&i| {
                        let (lat, lon) = self.cell_lat_lon[i];
                        self.land[i] >= 0.5
                            && lat.abs() <= ENSO_LAT_WIDTH_DEG
                            && near(lon, end, west)
                    })
                    .collect()
            };
            self.enso_ends = Some([pick(true, w), pick(false, e)]);
        }
        if t >= 0.5 {
            let ends = self.enso_ends.as_ref().unwrap();
            for k in 0..2 {
                for &i in &ends[k] {
                    self.window_enso_rain[k][0] += self.area[i] * self.precip_mm[i];
                    self.window_enso_rain[k][1] += self.area[i] * self.normal_p[i][month] as f64;
                }
            }
        }
    }

    /// Rain on land at the widest basin's western and eastern ends in El
    /// Niño months, % of normal (None without data).
    pub fn window_enso_rain_percent(&self) -> Option<[f64; 2]> {
        let r = self.window_enso_rain;
        if r[0][1] <= 0.0 || r[1][1] <= 0.0 {
            return None;
        }
        Some([100.0 * r[0][0] / r[0][1], 100.0 * r[1][0] / r[1][1]])
    }

    /// The world's average surface temperature right now, °C.
    pub fn global_mean_now(&self) -> f64 {
        let total: f64 = self.area.iter().sum();
        (0..self.n)
            .map(|i| (self.t[i] - LAPSE_RATE * self.height[i]) * self.area[i])
            .sum::<f64>()
            / total
    }

    /// Where the last three months' rain on land fell furthest below and
    /// above normal, in mm: ((lat, lon, % of normal), (lat, lon, % of
    /// normal)). Only cells that are mostly land count.
    pub fn rain_extremes(&self) -> Option<(RainPlace, RainPlace)> {
        let mut dry: Option<(usize, f32)> = None;
        let mut wet: Option<(usize, f32)> = None;
        for i in 0..self.n {
            if self.land[i] < 0.5 {
                continue;
            }
            let r: f32 = (0..3).map(|k| self.rain_recent[k][i]).sum();
            let nrm: f32 = (0..3).map(|k| self.normal_recent[k][i]).sum();
            let d = r - nrm;
            if dry.is_none_or(|(_, v)| d < v) {
                dry = Some((i, d));
            }
            if wet.is_none_or(|(_, v)| d > v) {
                wet = Some((i, d));
            }
        }
        let place = |i: usize| {
            let r: f32 = (0..3).map(|k| self.rain_recent[k][i]).sum();
            let nrm: f32 = (0..3).map(|k| self.normal_recent[k][i]).sum();
            let (lat, lon) = self.cell_lat_lon[i];
            (lat, lon, 100.0 * r as f64 / (nrm as f64).max(1.0))
        };
        Some((place(dry?.0), place(wet?.0)))
    }

    /// Turns on year-to-year variability. `eruptions` gives how often and
    /// how big eruptions are (None: no eruptions); `seed` picks this run's
    /// random history.
    pub fn set_variability(&mut self, eruptions: Option<EruptionStats>, seed: u64) {
        // Eruptions can happen on land, more likely where the ground is
        // steep: weight = land share × the largest height step to a
        // neighbouring cell (m). Real volcanoes follow plate edges, which
        // the terrain alone doesn't give; mountains are the closest proxy.
        let sites = (0..self.n)
            .map(|i| {
                (
                    i as u32,
                    self.land[i] * (1.0 - self.ice[i]) * self.relief[i],
                )
            })
            .collect();
        self.volcanism = Volcanism::new(eruptions, sites);
        let sea: Vec<f64> = (0..self.n).map(|i| self.sea_share(i)).collect();
        self.enso = Enso::new(&self.ocean.equator_wet(&sea));
        self.rng = Rng::new(seed);
    }

    fn sea_share(&self, i: usize) -> f64 {
        (1.0 - self.land[i]) * (1.0 - self.ice[i])
    }

    /// The month's temperatures: air on the 440 km grid, carried by winds
    /// and mixing; the open sea's mixed layer there, carried by currents;
    /// and each climate cell's ground, ice sheet and sea ice, which follow
    /// the air above them. Then sea ice freezes, grows and melts.
    fn temperature_month(&mut self, month: usize, dt: f64, co2_ppm: f64) {
        use crate::ocean::{
            SEA_DENSITY, SEA_SPECIFIC_HEAT, THERMOCLINE_SCALE_M, UPWELLING_SOURCE_DEPTH_M,
        };
        let nc = self.atmosphere.len();
        let a_air = OLR_A - CO2_FORCING_COEF * (co2_ppm / OLR_CO2_PPM).ln() - WINDOW_REFERENCE;
        let c_air = HEAT_CAPACITY_AIR / dt;
        let c_sea = HEAT_CAPACITY_SEA / dt;
        let c_ground = HEAT_CAPACITY_GROUND / dt;
        // Air: c_a T_a - Σ a (T_q - T_a) - open_γ S = b_a; sea: c_s S - γ_s T_a - Σ a_s (S_q - S) = b_s.
        let mut ca = vec![c_air + OLR_B; nc];
        let mut ba: Vec<f64> = self.air.iter().map(|&t| c_air * t - a_air).collect();
        // Latent heat: what each surface loses by evaporation, the air
        // above gains. (Heat released by rain in rising air is mostly
        // carried off again by the rising air itself; Neelin & Held 1987,
        // Mon. Wea. Rev. 115, 3-12. Putting it into the air where the rain
        // falls makes the rain belts swing wildly.)
        // Open sea per coarse cell: share, sunlight, coupling, window.
        let mut open = vec![0.0; nc];
        let mut sea_sun = vec![0.0; nc];
        let mut sea_ot = vec![0.0; nc];
        let mut eps = vec![0.0; nc];
        let mut speed = vec![0.0; nc];
        for p in 0..nc {
            let (u, v) = self.atmosphere.wind[p];
            speed[p] = (u * u + v * v).sqrt() + crate::atmosphere::GUSTINESS;
            // Near-surface vapour pressure from the column humidity and
            // the air temperature at sea level.
            let e = self.atmosphere.relative_humidity(p)
                * crate::atmosphere::saturation_vapour_pressure(self.air[p])
                / 100.0;
            eps[p] = air_emissivity(e);
        }
        let mut land_limit = vec![0.0; nc];
        let mut land_w = vec![0.0; nc];
        // Per climate cell: each surface tile ends up as T = α + β T_air.
        let mut alpha = std::mem::take(&mut self.buf_c);
        let mut beta = std::mem::take(&mut self.buf_b);
        alpha.resize(3 * self.n, 0.0);
        beta.resize(3 * self.n, 0.0);
        let ice_albedo = self.tuning.albedo_ice;
        let sun_month: Vec<f64> = self
            .sunshine
            .table
            .iter()
            .map(|r| r[month] as f64)
            .collect();
        let pe = &self.atmosphere.potential_evaporation;
        for i in 0..self.n {
            let f = self.facts[i];
            let p = f.parent as usize;
            let w = f.weight;
            let sun = sun_month[f.sun_row as usize] * (1.0 - self.dimming[f.parent as usize]);
            let iced = self.sea_ice_m[i] > 0.0;
            let (air_p, eps_p, speed_p) = (self.air[p], eps[p], speed[p]);
            for k in 0..3 {
                let share = f.shares[k];
                if share <= 0.0 {
                    continue;
                }
                // k: 0 ground, 1 ice sheet, 2 sea (open or iced).
                let albedo = if k == 1 || (k == 2 && iced) {
                    ice_albedo
                } else {
                    f.free_albedo
                };
                let absorbed = (1.0 - albedo) * sun;
                ba[p] += w * share * AIR_SUNLIGHT_SHARE * absorbed;
                let surface_sun = (1.0 - AIR_SUNLIGHT_SHARE) * absorbed;
                if k == 2 && !iced {
                    open[p] += w * share;
                    sea_sun[p] += w * share * surface_sun;
                    sea_ot[p] += w * self.overturning_heat[i];
                    continue;
                }
                let (t_old, height) = match k {
                    0 => (self.t_ground[i], f.height),
                    1 => (self.t_glacier[i], f.height),
                    _ => (self.t_seaice[i], 0.0),
                };
                let tk = t_old - LAPSE_RATE * height + 273.15;
                let tk3 = tk * tk * tk;
                let window = (1.0 - eps_p) * STEFAN_BOLTZMANN * tk3 * tk;
                let window_slope = 4.0 * (1.0 - eps_p) * STEFAN_BOLTZMANN * tk3;
                let turb = AIR_DENSITY_SURFACE
                    * AIR_SPECIFIC_HEAT
                    * HEAT_TRANSFER_COEF
                    * speed_p
                    * stability_factor(air_p - t_old, speed_p, tk);
                let g = turb + 4.0 * eps_p * STEFAN_BOLTZMANN * tk3;
                // Sea ice: heat conducted up from the water below at freezing.
                let cond = if k == 2 {
                    ICE_CONDUCTIVITY / self.sea_ice_m[i].max(MIN_CONDUCTING_ICE_M)
                } else {
                    0.0
                };
                // Evaporation from wet ground cools it: last month's rate
                // from the winds, capped by what the energy supply allows.
                let evap = if k == 0 {
                    let limit = (self.pt_factor.at(tk - 273.15) * (surface_sun - window)).max(0.0);
                    land_limit[p] += w * share * limit * (1.0 / LATENT_HEAT_VAPOUR);
                    land_w[p] += w * share;
                    (LATENT_HEAT_VAPOUR * pe[p]).min(limit) * self.beta(i)
                } else {
                    0.0
                };
                let inv_den = 1.0 / (c_ground + g + window_slope + cond);
                let a = (c_ground * t_old + surface_sun - window
                    + window_slope * t_old
                    + cond * SEA_FREEZING_C
                    - evap)
                    * inv_den;
                let bt = g * inv_den;
                alpha[3 * i + k] = a;
                beta[3 * i + k] = bt;
                ca[p] += w * share * g * (1.0 - bt);
                ba[p] += w * share * (g * a + evap);
            }
        }
        for p in 0..nc {
            self.atmosphere.land_energy_limit[p] = if land_w[p] > 0.0 {
                land_limit[p] / land_w[p]
            } else {
                f64::MAX
            };
        }
        // Sea coefficients per coarse cell (per m² of open sea).
        let k_deep = (-UPWELLING_SOURCE_DEPTH_M / THERMOCLINE_SCALE_M).exp();
        let mut cs = vec![1.0; nc];
        let mut bs = vec![0.0; nc];
        let mut gs = vec![0.0; nc];
        for p in 0..nc {
            if open[p] < 1e-3 {
                bs[p] = 0.0;
                continue;
            }
            let s_old = self.sea_t[p];
            let tk = s_old + 273.15;
            let window = (1.0 - eps[p]) * STEFAN_BOLTZMANN * tk.powi(4);
            let window_slope = 4.0 * window / tk;
            let turb = AIR_DENSITY_SURFACE
                * AIR_SPECIFIC_HEAT
                * HEAT_TRANSFER_COEF
                * speed[p]
                * stability_factor(self.air[p] - s_old, speed[p], tk);
            let g = turb + 4.0 * eps[p] * STEFAN_BOLTZMANN * tk.powi(3);
            let up = SEA_DENSITY * SEA_SPECIFIC_HEAT * self.ocean.upwelling[p];
            let s_year = if self.sea_by_month.iter().all(|v| v.len() == nc) {
                self.sea_by_month.iter().map(|v| v[p]).sum::<f64>() / 12.0
            } else {
                s_old
            };
            let deep = self.ocean.deep_temperature;
            let t_up = deep + (s_year.max(SEA_FREEZING_C) - deep) * k_deep;
            let evap = LATENT_HEAT_VAPOUR * self.atmosphere.potential_evaporation_sea(p);
            cs[p] = c_sea + g + window_slope + up;
            bs[p] = c_sea * s_old + sea_sun[p] / open[p] - window
                + window_slope * s_old
                + up * t_up
                + sea_ot[p] / open[p]
                - evap;
            gs[p] = g;
            ca[p] += open[p] * g;
            ba[p] += open[p] * evap;
        }
        // Solve both together, cell by cell (a 2 × 2 system each), with
        // over-relaxed passes until the largest change is under 0.01 °C.
        let (offsets, nbrs, geo) = self.solver.coarse();
        let k_mix = self.tuning.diffusion * EARTH_RADIUS_M * EARTH_RADIUS_M;
        let air_adv = &self.atmosphere.heat_adv;
        let sea_adv = &self.ocean.heat_adv;
        // Starting guess: last month's answer plus the change between the
        // same two months last year (follows both the seasons and any slow
        // drift of the oceans and ice).
        if self.air_by_month[month].len() == nc && self.prev_air.len() == nc {
            for p in 0..nc {
                self.air[p] += self.air_by_month[month][p] - self.prev_air[p];
                self.sea_t[p] += self.sea_by_month[month][p] - self.prev_sea[p];
            }
        }
        // The month's couplings, once: air (mixing + winds) and sea
        // (currents, per m² of open sea) per neighbour entry, and their sums.
        let m = nbrs.len();
        let mut ka = vec![0.0; m];
        let mut ks = vec![0.0; m];
        let mut sum_a = vec![0.0; nc];
        let mut sum_s = vec![0.0; nc];
        for p in 0..nc {
            for j in offsets[p] as usize..offsets[p + 1] as usize {
                let q = nbrs[j] as usize;
                ka[j] = k_mix * geo[j] + air_adv[j];
                sum_a[p] += ka[j];
                if open[p] >= 1e-3 && open[q] >= 1e-3 {
                    ks[j] = sea_adv[j] / open[p];
                    sum_s[p] += ks[j];
                }
            }
        }
        // Plain Gauss-Seidel: the starting guesses are close, and with
        // currents and winds the system isn't symmetric, so over-relaxing
        // tends to overshoot.
        // Each cell's 2 × 2 matrix is the same in every pass: invert once.
        let inv: Vec<[f64; 4]> = (0..nc)
            .map(|p| {
                if open[p] >= 1e-3 {
                    let (m11, m12) = (ca[p] + sum_a[p], open[p] * gs[p]);
                    let (m21, m22) = (gs[p], cs[p] + sum_s[p]);
                    let inv_det = 1.0 / (m11 * m22 - m12 * m21);
                    [m22 * inv_det, m12 * inv_det, m21, 1.0 / m22]
                } else {
                    [1.0 / (ca[p] + sum_a[p]), 0.0, 0.0, 0.0]
                }
            })
            .collect();
        let omega: f64 = 1.0;
        let mut passes = 0;
        for _ in 0..500 {
            passes += 1;
            let mut max_change = 0.0f64;
            for p in 0..nc {
                let (s0, e0) = (offsets[p] as usize, offsets[p + 1] as usize);
                let (mut ra, mut rs) = (ba[p], bs[p]);
                for j in s0..e0 {
                    let q = nbrs[j] as usize;
                    ra += ka[j] * self.air[q];
                    rs += ks[j] * self.sea_t[q];
                }
                let [c1, c2, c3, c4] = inv[p];
                let (ta, ts) = if open[p] >= 1e-3 {
                    let ta = ra * c1 + rs * c2;
                    (ta, (rs + c3 * ta) * c4)
                } else {
                    // No open sea here: its sea temperature is unused.
                    (ra * c1, self.sea_t[p])
                };
                let na = self.air[p] + omega * (ta - self.air[p]);
                let ns = self.sea_t[p] + omega * (ts - self.sea_t[p]);
                max_change = max_change
                    .max((na - self.air[p]).abs())
                    .max((ns - self.sea_t[p]).abs());
                self.air[p] = na;
                self.sea_t[p] = ns;
            }
            if max_change < 1e-2 {
                break;
            }
        }
        self.last_passes = passes;
        // Each cell's surfaces from the air above, smoothed between 440 km
        // cells so maps don't show the coarse grid.
        let per_m = ICE_DENSITY * LATENT_HEAT_FUSION;
        // Sea temperature handed down to cells with sea: where a 440 km
        // cell has no open sea of its own, use the mean of its neighbours
        // that do (or its air, if none do).
        let mut sea_fill = vec![0.0; nc];
        for p in 0..nc {
            sea_fill[p] = if open[p] >= 1e-3 {
                self.sea_t[p]
            } else {
                let (mut n_open, mut sum) = (0.0, 0.0);
                for &q in &nbrs[offsets[p] as usize..offsets[p + 1] as usize] {
                    if open[q as usize] >= 1e-3 {
                        n_open += 1.0;
                        sum += self.sea_t[q as usize];
                    }
                }
                if n_open > 0.0 {
                    sum / n_open
                } else {
                    self.air[p]
                }
            };
        }
        for i in 0..self.n {
            let (idx, wts) = self.solver.blend(i);
            let air_i: f64 = (0..7).map(|k| wts[k] * self.air[idx[k] as usize]).sum();
            let shares = self.facts[i].shares;
            let sea_i = if shares[2] > 0.0 {
                (0..7).map(|k| wts[k] * sea_fill[idx[k] as usize]).sum()
            } else {
                air_i
            };
            // (A cell's ground or ice-sheet temperature is only kept where
            // it has some; elsewhere it is never used.)
            if shares[0] > 0.0 {
                self.t_ground[i] = alpha[3 * i] + beta[3 * i] * air_i;
            }
            if shares[1] > 0.0 {
                self.t_glacier[i] = alpha[3 * i + 1] + beta[3 * i + 1] * air_i;
            }
            let p = self.facts[i].parent as usize;
            if self.sea_ice_m[i] > 0.0 && shares[2] > 0.0 {
                let mut t_ice = alpha[3 * i + 2] + beta[3 * i + 2] * air_i;
                let h = self.sea_ice_m[i].max(MIN_CONDUCTING_ICE_M);
                let up_flux = ICE_CONDUCTIVITY * (SEA_FREEZING_C - t_ice) / h;
                let ocean = SEA_DENSITY
                    * SEA_SPECIFIC_HEAT
                    * self.ocean.upwelling[p]
                    * (self.ocean.deep_temperature - SEA_FREEZING_C).max(0.0)
                    + self.overturning_heat[i] / shares[2];
                self.sea_ice_m[i] += (up_flux - ocean) * dt / per_m;
                if t_ice > 0.0 {
                    // The surface can't warm past melting: the extra heat
                    // melts ice instead.
                    let den = c_ground + 4.0 * STEFAN_BOLTZMANN * 273.15f64.powi(3);
                    self.sea_ice_m[i] -= den * t_ice * dt / per_m;
                    t_ice = 0.0;
                }
                self.t_seaice[i] = t_ice;
                if self.sea_ice_m[i] <= 0.0 {
                    self.sea_ice_m[i] = 0.0;
                    self.t_seaice[i] = SEA_FREEZING_C;
                }
            } else if shares[2] > 0.0 && sea_i < SEA_FREEZING_C {
                // Open water cooled below freezing: that heat came from
                // freezing water instead.
                self.sea_ice_m[i] = HEAT_CAPACITY_SEA * (SEA_FREEZING_C - sea_i) / per_m;
                self.t_seaice[i] = SEA_FREEZING_C;
            }
            let sea_surface = if self.sea_ice_m[i] > 0.0 {
                self.t_seaice[i]
            } else {
                sea_i.max(SEA_FREEZING_C)
            };
            self.t[i] = shares[0] * self.t_ground[i]
                + shares[1] * self.t_glacier[i]
                + shares[2] * sea_surface;
        }
        for p in 0..nc {
            if self.sea_t[p] < SEA_FREEZING_C {
                self.sea_t[p] = SEA_FREEZING_C;
            }
        }
        // Keep last year's value for this month (next month's guess uses
        // it), then store this year's.
        std::mem::swap(&mut self.prev_air, &mut self.air_by_month[month]);
        std::mem::swap(&mut self.prev_sea, &mut self.sea_by_month[month]);
        self.air_by_month[month].clone_from(&self.air);
        self.sea_by_month[month].clone_from(&self.sea_t);
        self.buf_c = alpha;
        self.buf_b = beta;
    }

    /// Shares of a climate cell that are ground, ice sheet and sea.
    fn tiles(&self, i: usize) -> [f64; 3] {
        let ground = self.land[i] * (1.0 - self.ice[i]);
        [ground, self.ice[i], self.sea_share(i)]
    }

    /// Advances one month (0 = January) of `days` days at CO2 `co2_ppm`.
    pub fn step_month(&mut self, year: i64, month: usize, days: u32, co2_ppm: f64) {
        let dt = days as f64 * 86_400.0;
        let clock = std::time::Instant::now();
        self.last_month = month;
        self.variability_month();
        // Winds from last month's air temperatures.
        let air = std::mem::take(&mut self.air);
        self.atmosphere.update_winds(month, &air);
        self.air = air;
        self.ocean.add_month(&self.atmosphere.wind);
        self.ocean
            .update_upwelling(&self.atmosphere.wind, &self.coarse_sea_share);
        self.timing[0] += clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        self.temperature_month(month, dt, co2_ppm);
        self.timing[1] += clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        self.water_month(month, dt, days);
        self.timing[2] += clock.elapsed().as_secs_f64();
        // The El Niño-like anomaly was laid on the sea for this month only.
        for p in 0..self.sea_t.len() {
            self.sea_t[p] -= self.enso_offset[p];
        }
        self.rain_recent.rotate_left(1);
        self.normal_recent.rotate_left(1);
        for i in 0..self.n {
            self.rain_recent[2][i] = self.precip_mm[i] as f32;
            self.normal_recent[2][i] = self.normal_p[i][month];
        }
        if self.window.is_some_and(|(a, b)| (a..=b).contains(&year)) {
            self.window_enso_month(month);
        }

        for i in 0..self.n {
            self.this_year[i][month] = self.t[i] as f32;
            self.this_year_p[i][month] = self.precip_mm[i] as f32;
            self.this_year_snow[i][month] = self.snow_mm[i] as f32;
        }
        if let Some((first, last)) = self.window {
            if (first..=last).contains(&year) {
                let ext = self.sea_ice_extent_km2();
                self.window_ice_extent[month][0] += ext.0;
                self.window_ice_extent[month][1] += ext.1;
                for i in 0..self.n {
                    self.window_sum[i][month] += self.t[i];
                    self.window_sum_p[i][month] += self.precip_mm[i];
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

    /// The month's moisture: rain and snow, soil water, snowpack.
    fn water_month(&mut self, month: usize, dt: f64, days: u32) {
        let nc = self.atmosphere.len();
        self.co_t.iter_mut().for_each(|v| *v = 0.0);
        self.co_sea.iter_mut().for_each(|v| *v = 0.0);
        self.co_land.iter_mut().for_each(|v| *v = 0.0);
        self.co_beta.iter_mut().for_each(|v| *v = 0.0);
        for i in 0..self.n {
            let p = self.atmosphere.parent[i] as usize;
            let a = self.area[i];
            self.co_t[p] += a * (self.t[i] - LAPSE_RATE * self.height[i]);
            let open = if self.sea_ice_m[i] > 0.0 { 0.0 } else { 1.0 };
            self.co_sea[p] += a * (1.0 - self.land[i]) * (1.0 - self.ice[i]) * open;
            let land = self.land[i] * (1.0 - self.ice[i]);
            self.co_land[p] += a * land;
            self.co_beta[p] += a * land * self.beta(i);
        }
        for p in 0..nc {
            let area = self.atmosphere.area[p];
            self.co_t[p] /= area;
            if self.co_land[p] > 0.0 {
                self.co_beta[p] /= self.co_land[p];
            }
            self.co_sea[p] /= area;
            self.co_land[p] /= area;
        }
        let (co_t, sea, land, beta) = (
            std::mem::take(&mut self.co_t),
            std::mem::take(&mut self.co_sea),
            std::mem::take(&mut self.co_land),
            std::mem::take(&mut self.co_beta),
        );
        self.atmosphere
            .solve_moisture(month, &co_t, &sea, &land, &beta);
        (self.co_t, self.co_sea, self.co_land, self.co_beta) = (co_t, sea, land, beta);
        for i in 0..self.n {
            let p = self.atmosphere.parent[i] as usize;
            let fall = self.atmosphere.precip[i] * dt;
            self.precip_mm[i] = fall;
            // Rain minus evaporation over the sea, m of water.
            self.fresh_m[i] += (fall - self.atmosphere.potential_evaporation_sea(p) * dt) / 1000.0;
            if self.land[i] <= 0.0 {
                continue;
            }
            let t_surface = self.t[i] - LAPSE_RATE * self.height[i];
            // Days within a month vary around its mean: the share of
            // snowfall and the degree-days of melt follow from a normal
            // spread of daily temperatures (Reeh 1991).
            let snowfall = fall * self.snow_share.at(t_surface);
            let rain = fall - snowfall;
            let melt = (DEGREE_DAY_FACTOR * self.positive_degrees.at(t_surface) * days as f64)
                .min(self.snow_mm[i] + snowfall);
            self.snow_mm[i] += snowfall - melt;
            let evap = (self.atmosphere.potential_evaporation[p] * dt * self.beta(i))
                .min(self.soil_mm[i] + rain + melt);
            self.soil_mm[i] = (self.soil_mm[i] + rain + melt - evap).min(SOIL_CAPACITY_MM);
        }
    }

    /// Land evaporation as a share of the full rate (Manabe 1969).
    fn beta(&self, i: usize) -> f64 {
        (self.soil_mm[i] * (1.0 / (SOIL_CAPACITY_MM * SOIL_FULL_RATE_SHARE))).min(1.0)
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
                self.normal_p[i][m] += weight * (self.this_year_p[i][m] - self.normal_p[i][m]);
                self.normal_snow[i][m] +=
                    weight * (self.this_year_snow[i][m] - self.normal_snow[i][m]);
            }
            let annual: f32 = self.this_year[i].iter().sum::<f32>() / 12.0;
            sum += (annual as f64 - LAPSE_RATE * self.height[i]) * self.area[i];
            area += self.area[i];
        }
        self.global_history.push((year, (sum / area) as f32));
        self.years += 1;
        let sea: Vec<f64> = (0..self.n).map(|i| self.sea_share(i)).collect();
        self.ocean.update(&sea);
        let year_s = 365.25 * 86_400.0;
        let t_annual: Vec<f64> = self
            .this_year
            .iter()
            .map(|r| r.iter().map(|&v| v as f64).sum::<f64>() / 12.0)
            .collect();
        let t_coldest: Vec<f64> = self
            .this_year
            .iter()
            .map(|r| r.iter().cloned().fold(f32::MAX, f32::min) as f64)
            .collect();
        let fresh: Vec<f64> = self.fresh_m.iter().map(|v| v / year_s).collect();
        self.fresh_m.iter_mut().for_each(|v| *v = 0.0);
        self.ocean
            .update_overturning(&t_annual, &t_coldest, &fresh, &sea);
        let mut heat = std::mem::take(&mut self.overturning_heat);
        self.ocean.overturning_heating(&mut heat, &sea);
        self.overturning_heat = heat;
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

    /// Per map cell: recent-climate rainfall (mm) for `month`, or the annual
    /// total when None.
    pub fn map_precip(&self, month: Option<usize>) -> Vec<f32> {
        self.map_parent
            .iter()
            .map(|&p| {
                let row = &self.normal_p[p as usize];
                match month {
                    Some(m) => row[m],
                    None => row.iter().sum(),
                }
            })
            .collect()
    }

    /// Per map cell: recent-climate lying snow (mm of water) at the end of
    /// `month`, or the yearly maximum when None.
    pub fn map_snow(&self, month: Option<usize>) -> Vec<f32> {
        self.map_parent
            .iter()
            .map(|&p| {
                let row = &self.normal_snow[p as usize];
                match month {
                    Some(m) => row[m],
                    None => row.iter().copied().fold(0.0, f32::max),
                }
            })
            .collect()
    }

    /// Area of sea covered by ice now, km² (northern, southern hemisphere).
    pub fn sea_ice_extent_km2(&self) -> (f64, f64) {
        let (mut n, mut s) = (0.0, 0.0);
        for i in 0..self.n {
            if self.sea_ice_m[i] > 0.0 {
                let a = self.sea_share(i) * self.area[i] / 1e6;
                if self.northern[i] {
                    n += a;
                } else {
                    s += a;
                }
            }
        }
        (n, s)
    }

    /// Average sea ice extent for `month` over the averaging window, km²
    /// (north, south).
    pub fn window_sea_ice_extent_km2(&self, month: usize) -> Option<(f64, f64)> {
        if self.window_years == 0 {
            return None;
        }
        let y = self.window_years as f64;
        Some((
            self.window_ice_extent[month][0] / y,
            self.window_ice_extent[month][1] / y,
        ))
    }

    /// Per map cell: sea ice thickness now, m (0 on open sea and land).
    pub fn map_sea_ice(&self) -> Vec<f32> {
        self.map_parent
            .iter()
            .map(|&p| self.sea_ice_m[p as usize] as f32)
            .collect()
    }

    /// Per map cell: rain and snowfall this month, mm.
    pub fn map_precip_current(&self) -> Vec<f32> {
        self.map_parent
            .iter()
            .map(|&p| self.precip_mm[p as usize] as f32)
            .collect()
    }

    /// Per map cell: lying snow now, mm of water.
    pub fn map_snow_current(&self) -> Vec<f32> {
        self.map_parent
            .iter()
            .map(|&p| self.snow_mm[p as usize] as f32)
            .collect()
    }

    /// Per map cell: climate zone of the recent climate, or of the
    /// averaging window when `window` is set (None before it has data).
    pub fn map_zones(&self, window: bool, northern: &[bool]) -> Vec<Option<crate::koppen::Zone>> {
        (0..self.map_parent.len())
            .map(|i| {
                let (t, p) = if window {
                    (self.cell_window(i)?, self.cell_precip_window(i)?)
                } else {
                    if self.years == 0 {
                        return None;
                    }
                    (self.cell_normal(i), self.cell_precip_normal(i))
                };
                Some(crate::koppen::classify(&t, &p, northern[i]))
            })
            .collect()
    }

    /// Recent-climate lying snow of one map cell by month, mm of water.
    pub fn cell_snow_normal(&self, map_cell: usize) -> [f32; 12] {
        self.normal_snow[self.map_parent[map_cell] as usize]
    }

    /// Recent-climate monthly rainfall of one map cell, mm.
    pub fn cell_precip_normal(&self, map_cell: usize) -> [f32; 12] {
        self.normal_p[self.map_parent[map_cell] as usize]
    }

    /// Monthly rainfall over the averaging window for one map cell, mm.
    pub fn cell_precip_window(&self, map_cell: usize) -> Option<[f32; 12]> {
        if self.window_years == 0 {
            return None;
        }
        let row = self.window_sum_p[self.map_parent[map_cell] as usize];
        Some(row.map(|s| (s / self.window_years as f64) as f32))
    }

    /// Wind (east, north; m/s) over a map cell this month.
    pub fn map_wind(&self, map_cell: usize) -> (f64, f64) {
        self.atmosphere
            .wind_at_fine(self.map_parent[map_cell] as usize)
    }

    /// Passes the transport solver needed last month (for profiling).
    pub fn solver_passes(&self) -> usize {
        self.last_passes
    }

    pub fn years_simulated(&self) -> u32 {
        self.years
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degree_days() {
        // Far above freezing every day is warm: degrees = the mean.
        assert!((positive_degrees(20.0) - 20.0).abs() < 1e-3);
        // Far below: no melt.
        assert!(positive_degrees(-20.0) < 1e-3);
        // At a 0 °C mean, half the days thaw: σ/√(2π) ≈ 2.0 °C.
        assert!((positive_degrees(0.0) - 1.9947).abs() < 1e-3);
        assert!((erfc(0.0) - 1.0).abs() < 1e-6 && (erfc(1.0) - 0.157_299).abs() < 1e-6);
    }
}
