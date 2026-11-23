//! The present-day check: the model run in today's world, compared with
//! observed 1981-2010 temperatures (CHELSA V2.1, land) and with the
//! observed global mean.

use crate::model::Climate;
use serde::Serialize;
use world::data::ObservedTemperature;
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
