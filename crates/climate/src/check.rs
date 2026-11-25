//! The present-day check: the model run in today's world, compared with
//! observed 1981-2010 temperatures (CHELSA V2.1, land) and with the
//! observed global mean.

use crate::model::Climate;
use serde::Serialize;
use world::data::{ObservedPrecipitation, ObservedTemperature};
use world::grid::{to_lat_lon, Grid};

/// Observed global mean near-surface air temperature 1981-2010:
/// 14.3 ± 0.5 °C (Jones et al. 1999 climatology adjusted to 1981-2010;
/// reanalyses give 14.05-14.55 °C. See Jones & Harpham 2013, J. Geophys.
/// Res. 118, 3213-3217, and the comparison at RealClimate, Aug 2017).
pub const OBSERVED_GLOBAL_MEAN_C: f64 = 14.3;
pub const OBSERVED_GLOBAL_MEAN_UNCERTAINTY_C: f64 = 0.5;

/// Agreed pass marks (Phase 1 plan, "How we'll know it's realistic").
pub const PASS_GLOBAL_MEAN_C: f64 = 1.0;
pub const PASS_ANNUAL_TYPICAL_ERROR_C: f64 = 3.0;
pub const PASS_SEASONAL_SWING_ERROR_C: f64 = 4.0;

#[derive(Serialize, Clone, Debug)]
pub struct Band {
    pub lat_from: f64,
    pub lat_to: f64,
    pub model_c: f64,
    pub observed_c: f64,
    /// Warmest minus coldest month, land average.
    pub model_swing_c: f64,
    pub observed_swing_c: f64,
}

#[derive(Serialize, Clone, Debug)]
pub struct Scores {
    pub years_averaged: u32,
    pub global_mean_model_c: f64,
    pub global_mean_observed_c: f64,
    /// Average of model minus observed over land cells (area-weighted).
    pub land_bias_c: f64,
    /// Root-mean-square error of the annual mean over land cells.
    pub annual_typical_error_c: f64,
    /// Root-mean-square error of warmest-minus-coldest month over land.
    pub seasonal_swing_error_c: f64,
    pub land_cells_compared: usize,
    pub bands: Vec<Band>,
    pub pass_global_mean: bool,
    pub pass_annual: bool,
    pub pass_seasonal: bool,
}

/// Compares land cells only (at least half land today, `land_share` >= 128),
/// as agreed: people live on land, and observations there are densest.
pub fn scores(
    climate: &Climate,
    obs: &ObservedTemperature,
    map_grid: &Grid,
    land_share: &[u8],
) -> Option<Scores> {
    let global = climate.window_global_mean()?;
    let (mut w_sum, mut bias, mut sq_annual, mut sq_swing) = (0.0, 0.0, 0.0, 0.0);
    let band_deg = 10.0;
    let nb = (180.0 / band_deg) as usize;
    let mut bands = vec![(0.0, 0.0, 0.0, 0.0, 0.0); nb];
    let mut compared = 0;
    for (k, &cell) in obs.cells.iter().enumerate() {
        let cell = cell as usize;
        if land_share[cell] < 128 {
            continue;
        }
        compared += 1;
        let model = climate.cell_window(cell)?;
        let o = &obs.monthly_c[k * 12..k * 12 + 12];
        let w = map_grid.areas_m2[cell];
        let m_ann = model.iter().map(|&v| v as f64).sum::<f64>() / 12.0;
        let o_ann = o.iter().map(|&v| v as f64).sum::<f64>() / 12.0;
        let swing = |v: &mut dyn Iterator<Item = f64>| {
            let (lo, hi) = v.fold((f64::MAX, f64::MIN), |(a, b), x| (a.min(x), b.max(x)));
            hi - lo
        };
        let m_sw = swing(&mut model.iter().map(|&v| v as f64));
        let o_sw = swing(&mut o.iter().map(|&v| v as f64));
        w_sum += w;
        bias += w * (m_ann - o_ann);
        sq_annual += w * (m_ann - o_ann).powi(2);
        sq_swing += w * (m_sw - o_sw).powi(2);
        let lat = to_lat_lon(map_grid.centers[cell]).0;
        let b = (((lat + 90.0) / band_deg) as usize).min(nb - 1);
        bands[b].0 += w * m_ann;
        bands[b].1 += w * o_ann;
        bands[b].2 += w;
        bands[b].3 += w * m_sw;
        bands[b].4 += w * o_sw;
    }
    let annual = (sq_annual / w_sum).sqrt();
    let seasonal = (sq_swing / w_sum).sqrt();
    Some(Scores {
        years_averaged: climate.window_years(),
        global_mean_model_c: global,
        global_mean_observed_c: OBSERVED_GLOBAL_MEAN_C,
        land_bias_c: bias / w_sum,
        annual_typical_error_c: annual,
        seasonal_swing_error_c: seasonal,
        land_cells_compared: compared,
        bands: bands
            .iter()
            .enumerate()
            .filter(|(_, b)| b.2 > 0.0)
            .map(|(i, b)| Band {
                lat_from: -90.0 + i as f64 * band_deg,
                lat_to: -90.0 + (i + 1) as f64 * band_deg,
                model_c: b.0 / b.2,
                observed_c: b.1 / b.2,
                model_swing_c: b.3 / b.2,
                observed_swing_c: b.4 / b.2,
            })
            .collect(),
        pass_global_mean: (global - OBSERVED_GLOBAL_MEAN_C).abs() <= PASS_GLOBAL_MEAN_C,
        pass_annual: annual <= PASS_ANNUAL_TYPICAL_ERROR_C,
        pass_seasonal: seasonal <= PASS_SEASONAL_SWING_ERROR_C,
    })
}

/// Agreed Step 3 pass marks (Phase 1 plan).
pub const PASS_RAIN_CORRELATION: f64 = 0.7;
pub const PASS_ZONE_AGREEMENT: f64 = 0.7;

/// A place the model must get right in kind (dry or wet).
#[derive(Serialize, Clone, Debug)]
pub struct Landmark {
    pub name: &'static str,
    pub lat: f64,
    pub lon: f64,
    pub expect: &'static str,
    /// Model and observed annual rainfall (mm); for the monsoon, the share
    /// of the year's rain falling June-September.
    pub model: f64,
    pub observed: f64,
    pub pass: bool,
}

#[derive(Serialize, Clone, Debug)]
pub struct RainScores {
    /// Area-weighted correlation of annual rainfall over land cells.
    pub correlation: f64,
    pub land_mean_model_mm: f64,
    pub land_mean_observed_mm: f64,
    /// Share of land where model and observed climate zones have the same
    /// main class (A, B, C, D, E).
    pub zone_agreement: f64,
    /// Area share of each zone (order of `koppen::Zone::ALL`), model and
    /// observed.
    pub zone_share_model: Vec<f64>,
    pub zone_share_observed: Vec<f64>,
    pub landmarks: Vec<Landmark>,
    pub pass_correlation: bool,
    pub pass_zones: bool,
    pub pass_landmarks: bool,
}

/// Places from the Phase 1 plan: deserts must be dry (under 250 mm a year,
/// the usual desert limit), rainforests wet (over 1,500 mm), and central
/// India must get most of its rain in the summer monsoon (June-September)
/// and enough of it to farm without irrigation (at least 500 mm a year).
const LANDMARKS: [(&str, f64, f64, &str); 8] = [
    ("Sahara", 23.0, 10.0, "dry"),
    ("Arabian Desert", 22.0, 47.0, "dry"),
    ("Atacama", -23.0, -69.5, "dry"),
    ("Australian interior", -25.0, 133.0, "dry"),
    ("Amazon", -3.0, -60.0, "wet"),
    ("Congo basin", 0.0, 22.0, "wet"),
    ("Borneo", 1.0, 114.0, "wet"),
    ("Central India", 22.0, 78.0, "monsoon"),
];

/// Rainfall and climate-zone scores over land cells (at least half land).
pub fn rain_scores(
    climate: &Climate,
    obs_t: &ObservedTemperature,
    obs_p: &ObservedPrecipitation,
    map_grid: &Grid,
    land_share: &[u8],
) -> Option<RainScores> {
    use crate::koppen::classify;
    // Observed temperature and rain are on the same cells (same mask).
    let index: std::collections::HashMap<u32, usize> = obs_t
        .cells
        .iter()
        .enumerate()
        .map(|(k, &c)| (c, k))
        .collect();
    let (mut sw, mut sm, mut so, mut smm, mut soo, mut smo) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut agree, mut zone_m, mut zone_o) = (0.0, vec![0.0; 6], vec![0.0; 6]);
    for (k, &cell) in obs_p.cells.iter().enumerate() {
        let c = cell as usize;
        if land_share[c] < 128 {
            continue;
        }
        let Some(&kt) = index.get(&cell) else {
            continue;
        };
        let mp = climate.cell_precip_window(c)?;
        let mt = climate.cell_window(c)?;
        let op: [f32; 12] = obs_p.monthly_mm[k * 12..k * 12 + 12].try_into().unwrap();
        let ot: [f32; 12] = obs_t.monthly_c[kt * 12..kt * 12 + 12].try_into().unwrap();
        let w = map_grid.areas_m2[c];
        let (m, o) = (mp.iter().sum::<f32>() as f64, op.iter().sum::<f32>() as f64);
        sw += w;
        sm += w * m;
        so += w * o;
        smm += w * m * m;
        soo += w * o * o;
        smo += w * m * o;
        let north = to_lat_lon(map_grid.centers[c]).0 >= 0.0;
        let (zm, zo) = (classify(&mt, &mp, north), classify(&ot, &op, north));
        zone_m[zm.index() as usize] += w;
        zone_o[zo.index() as usize] += w;
        if zm.main_class() == zo.main_class() {
            agree += w;
        }
    }
    if sw == 0.0 {
        return None;
    }
    let (mm, mo) = (sm / sw, so / sw);
    let cov = smo / sw - mm * mo;
    let correlation = cov / ((smm / sw - mm * mm).sqrt() * (soo / sw - mo * mo).sqrt());
    let mut landmarks = Vec::new();
    for (name, lat, lon, expect) in LANDMARKS {
        let cell = map_grid.nearest_from(world::grid::from_lat_lon(lat, lon), 0);
        let mp = climate.cell_precip_window(cell)?;
        let op = obs_p
            .cells
            .iter()
            .position(|&x| x as usize == cell)
            .map(|k| &obs_p.monthly_mm[k * 12..k * 12 + 12]);
        let summer_share =
            |p: &[f32]| p[5..9].iter().sum::<f32>() as f64 / p.iter().sum::<f32>().max(1e-3) as f64;
        let annual = |p: &[f32]| p.iter().sum::<f32>() as f64;
        let (model, observed, pass) = match expect {
            "dry" => (
                annual(&mp),
                op.map(annual).unwrap_or(f64::NAN),
                annual(&mp) < 250.0,
            ),
            "wet" => (
                annual(&mp),
                op.map(annual).unwrap_or(f64::NAN),
                annual(&mp) > 1500.0,
            ),
            _ => (
                summer_share(&mp),
                op.map(summer_share).unwrap_or(f64::NAN),
                summer_share(&mp) > 0.5 && annual(&mp) >= 500.0,
            ),
        };
        landmarks.push(Landmark {
            name,
            lat,
            lon,
            expect,
            model,
            observed,
            pass,
        });
    }
    let zone_agreement = agree / sw;
    Some(RainScores {
        correlation,
        land_mean_model_mm: mm,
        land_mean_observed_mm: mo,
        zone_agreement,
        zone_share_model: zone_m.iter().map(|v| v / sw).collect(),
        zone_share_observed: zone_o.iter().map(|v| v / sw).collect(),
        pass_landmarks: landmarks.iter().all(|l| l.pass),
        landmarks,
        pass_correlation: correlation >= PASS_RAIN_CORRELATION,
        pass_zones: zone_agreement >= PASS_ZONE_AGREEMENT,
    })
}

/// Agreed Step 4 marks: sea ice extent within 25% of observed.
pub const PASS_SEA_ICE_RELATIVE: f64 = 0.25;

/// Observed 1981-2010 mean sea ice extent, million km² (NSIDC Sea Ice
/// Index v3, monthly climatology): Arctic March and September, Antarctic
/// September and February.
pub const OBSERVED_SEA_ICE: [(&str, usize, bool, f64); 4] = [
    ("Arctic, March", 2, true, 15.5),
    ("Arctic, September", 8, true, 6.5),
    ("Antarctic, September", 8, false, 18.5),
    ("Antarctic, February", 1, false, 3.0),
];

#[derive(Serialize, Clone, Debug)]
pub struct SeaIce {
    pub name: &'static str,
    pub model_million_km2: f64,
    pub observed_million_km2: f64,
    pub pass: bool,
}

/// A contrast between two places that ocean currents create.
#[derive(Serialize, Clone, Debug)]
pub struct Contrast {
    pub name: &'static str,
    /// Warmer place minus cooler place, model and observed, °C.
    pub model_c: f64,
    pub observed_c: f64,
    /// Right sign and at least half the observed difference.
    pub pass: bool,
}

#[derive(Serialize, Clone, Debug)]
pub struct OceanScores {
    pub sea_ice: Vec<SeaIce>,
    pub contrasts: Vec<Contrast>,
}

/// Bergen versus the Labrador coast in January (the North Atlantic's warm
/// water), and Salvador versus Lima over the year (Lima's cold current).
/// Each: name, warmer place (lat, lon), cooler place, month (None = annual).
type ContrastSpec = (&'static str, (f64, f64), (f64, f64), Option<usize>);

const CONTRASTS: [ContrastSpec; 2] = [
    (
        "Bergen warmer than Labrador in January",
        (60.4, 5.3),
        (56.5, -61.7),
        Some(0),
    ),
    (
        "Salvador warmer than Lima over the year",
        (-12.97, -38.5),
        (-12.05, -77.04),
        None,
    ),
];

pub fn ocean_scores(
    climate: &Climate,
    obs: &ObservedTemperature,
    map_grid: &Grid,
) -> Option<OceanScores> {
    let mut sea_ice = Vec::new();
    for (name, month, north, observed) in OBSERVED_SEA_ICE {
        let (n, s) = climate.window_sea_ice_extent_km2(month)?;
        let model = if north { n } else { s } / 1e6;
        sea_ice.push(SeaIce {
            name,
            model_million_km2: model,
            observed_million_km2: observed,
            pass: (model - observed).abs() <= PASS_SEA_ICE_RELATIVE * observed,
        });
    }
    let value = |row: &[f32], month: Option<usize>| match month {
        Some(m) => row[m] as f64,
        None => row.iter().map(|&v| v as f64).sum::<f64>() / 12.0,
    };
    let mut contrasts = Vec::new();
    for (name, warm, cool, month) in CONTRASTS {
        // The observed land cell nearest each place (coastal towns can
        // fall in a mostly-sea cell).
        let nearest = |(lat, lon): (f64, f64)| {
            let p = world::grid::from_lat_lon(lat, lon);
            let k = (0..obs.cells.len())
                .max_by(|&a, &b| {
                    let d = |i: usize| {
                        let c = map_grid.centers[obs.cells[i] as usize];
                        c[0] * p[0] + c[1] * p[1] + c[2] * p[2]
                    };
                    d(a).total_cmp(&d(b))
                })
                .unwrap();
            (
                obs.cells[k] as usize,
                value(&obs.monthly_c[k * 12..k * 12 + 12], month),
            )
        };
        let ((w, ow), (c, oc)) = (nearest(warm), nearest(cool));
        let model_c =
            value(&climate.cell_window(w)?, month) - value(&climate.cell_window(c)?, month);
        let observed_c = ow - oc;
        contrasts.push(Contrast {
            name,
            model_c,
            observed_c,
            pass: model_c >= 0.5 * observed_c,
        });
    }
    Some(OceanScores { sea_ice, contrasts })
}

// ---------------------------------------------------------------------------
// Step 5: eruptions and the El Niño-like cycle.

/// Observed global cooling in the year after Pinatubo (mid-1991): about
/// 0.4 °C at its largest in 1992 (Parker et al. 1996; Hansen et al. 1996);
/// the agreed pass range.
pub const PINATUBO_OBSERVED_C: f64 = 0.4;
pub const PASS_PINATUBO_C: (f64, f64) = (0.3, 0.5);
/// El Niño's rhythm: a cycle of 2-7 years (the agreed range; observed peak
/// 3-5 years) and a typical size within 50% of the observed 0.9 °C.
pub const PASS_ENSO_PERIOD_YEARS: (f64, f64) = (2.0, 7.0);
pub const PASS_ENSO_SIZE_RELATIVE: f64 = 0.5;

#[derive(Serialize, Clone, Debug)]
pub struct VariabilityScores {
    pub pinatubo_cooling_c: Option<f64>,
    pub pinatubo_observed_c: f64,
    pub pinatubo_pass: bool,
    pub enso_size_c: f64,
    pub enso_observed_size_c: f64,
    pub enso_size_pass: bool,
    pub enso_period_years: Option<f64>,
    pub enso_period_pass: bool,
    /// Rain in El Niño months, % of normal, at the far west of the basin
    /// (observed: drier) and on its eastern coast (observed: wetter).
    pub west_rain_percent: Option<f64>,
    pub west_pass: bool,
    pub east_rain_percent: Option<f64>,
    pub east_pass: bool,
}

/// Scores for the comparison window. `pinatubo_cooling_c` comes from a
/// twin run with and without a Pinatubo-sized eruption.
pub fn variability_scores(
    climate: &Climate,
    pinatubo_cooling_c: Option<f64>,
) -> Option<VariabilityScores> {
    let s: Vec<f64> = climate.window_enso.iter().map(|&v| v as f64).collect();
    if s.len() < 120 {
        return None;
    }
    let mean = s.iter().sum::<f64>() / s.len() as f64;
    let size = (s.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / s.len() as f64).sqrt();
    let observed = crate::variability::ENSO_REFERENCE_SIZE_C;
    let period = crate::variability::period_months(&s).map(|m| m / 12.0);
    let rain = climate.window_enso_rain_percent();
    let (lo, hi) = PASS_PINATUBO_C;
    Some(VariabilityScores {
        pinatubo_cooling_c,
        pinatubo_observed_c: PINATUBO_OBSERVED_C,
        pinatubo_pass: pinatubo_cooling_c.is_some_and(|c| (lo..=hi).contains(&c)),
        enso_size_c: size,
        enso_observed_size_c: observed,
        enso_size_pass: (size - observed).abs() <= PASS_ENSO_SIZE_RELATIVE * observed,
        enso_period_years: period,
        enso_period_pass: period
            .is_some_and(|p| (PASS_ENSO_PERIOD_YEARS.0..=PASS_ENSO_PERIOD_YEARS.1).contains(&p)),
        west_rain_percent: rain.map(|r| r[0]),
        west_pass: rain.is_some_and(|r| r[0] < 100.0),
        east_rain_percent: rain.map(|r| r[1]),
        east_pass: rain.is_some_and(|r| r[1] > 100.0),
    })
}
