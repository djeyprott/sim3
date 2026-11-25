//! Climate data for the viewer: maps of temperature, rainfall, snow and
//! climate zones, the winds, a place's climate chart, the world temperature
//! history and the present-day check scores.

use axum::body::Bytes;
use axum::extract::{Query, State};
use axum::http::header;
use axum::response::IntoResponse;
use axum::Json;
use serde::Deserialize;

use crate::AppState;

/// Marks "no value" in temperature maps (hundredths of °C as i16).
const NO_VALUE: i16 = i16::MIN;

#[derive(Deserialize)]
pub struct MapQuery {
    /// "temperature" (default), "rain", "snow" or "zone".
    field: Option<String>,
    /// "current", "normal", "observed" or "difference".
    view: String,
    /// 0-11; absent for the annual mean.
    month: Option<usize>,
}

fn encode(values: impl Iterator<Item = Option<f32>>) -> Bytes {
    let mut out = Vec::new();
    for v in values {
        let x = match v {
            Some(t) if t.is_finite() => (t * 100.0).round().clamp(-32767.0, 32767.0) as i16,
            _ => NO_VALUE,
        };
        out.extend_from_slice(&x.to_le_bytes());
    }
    Bytes::from(out)
}

/// Cells compared with observations: at least half land today.
fn compared_land(world: &sim::WorldContext) -> &[u8] {
    &world.data.epochs[world.epoch_for_year(2000)].land_share
}

/// Observed 1981-2010 temperature per map cell for a month (or the year),
/// on land only (cells at least half land today): the agreed comparison.
fn observed(world: &sim::WorldContext, month: Option<usize>) -> Vec<Option<f32>> {
    let mut out = vec![None; world.map_grid.len()];
    let land = compared_land(world);
    if let Some(obs) = &world.data.observed_temperature {
        for (k, &cell) in obs.cells.iter().enumerate() {
            if land[cell as usize] < 128 {
                continue;
            }
            let row = &obs.monthly_c[k * 12..k * 12 + 12];
            out[cell as usize] = Some(match month {
                Some(m) => row[m],
                None => row.iter().sum::<f32>() / 12.0,
            });
        }
    }
    out
}

/// Observed 1981-2010 rainfall (mm) per map cell, for a month or the year
/// total, on land.
fn observed_rain(world: &sim::WorldContext, month: Option<usize>) -> Vec<Option<f32>> {
    let mut out = vec![None; world.map_grid.len()];
    let land = compared_land(world);
    if let Some(obs) = &world.data.observed_precipitation {
        for (k, &cell) in obs.cells.iter().enumerate() {
            if land[cell as usize] < 128 {
                continue;
            }
            let row = &obs.monthly_mm[k * 12..k * 12 + 12];
            out[cell as usize] = Some(match month {
                Some(m) => row[m],
                None => row.iter().sum(),
            });
        }
    }
    out
}

/// Climate zone (index into `climate::koppen::Zone::ALL`) of the
/// observations per map cell, on land.
fn observed_zones(world: &sim::WorldContext, northern: &[bool]) -> Vec<Option<f32>> {
    let mut out = vec![None; world.map_grid.len()];
    let land = compared_land(world);
    let (Some(ot), Some(op)) = (
        &world.data.observed_temperature,
        &world.data.observed_precipitation,
    ) else {
        return out;
    };
    let rain: std::collections::HashMap<u32, usize> =
        op.cells.iter().enumerate().map(|(k, &c)| (c, k)).collect();
    for (k, &cell) in ot.cells.iter().enumerate() {
        let c = cell as usize;
        let Some(&kp) = rain.get(&cell) else { continue };
        if land[c] < 128 {
            continue;
        }
        let t: [f32; 12] = ot.monthly_c[k * 12..k * 12 + 12].try_into().unwrap();
        let p: [f32; 12] = op.monthly_mm[kp * 12..kp * 12 + 12].try_into().unwrap();
        out[c] = Some(climate::koppen::classify(&t, &p, northern[c]).index() as f32);
    }
    out
}

fn northern(world: &sim::WorldContext) -> Vec<bool> {
    world
        .map_grid
        .centers
        .iter()
        .map(|&c| c[2] >= 0.0)
        .collect()
}

/// Per map cell, i16 little-endian; -32768 = no value. Temperature in
/// hundredths of °C; rain and snow in mm (rain differences as 100 ×
/// log10(model / observed)); zones as the zone's index.
pub async fn get_map(
    State(shared): State<AppState>,
    Query(q): Query<MapQuery>,
) -> impl IntoResponse {
    let month = q.month.filter(|&m| m < 12);
    let field = q.field.as_deref().unwrap_or("temperature");
    let bytes = {
        let sim = shared.sim.lock().unwrap();
        let c = &sim.climate;
        let window = sim.scenario.comparison_window().is_some() && c.window_years() > 0;
        let pick = |row: [f32; 12], total: bool| match month {
            Some(m) => row[m],
            None if total => row.iter().sum(),
            None => row.iter().sum::<f32>() / 12.0,
        };
        let world = &shared.world;
        let n = world.map_grid.len();
        match (field, q.view.as_str()) {
            ("rain", "current") => encode_scaled(c.map_precip_current().into_iter().map(Some), 1.0),
            ("rain", "anomaly") => {
                let now = c.map_precip_current();
                let normal = c.map_precip(Some(c.last_month));
                encode_scaled(
                    (0..n).map(|i| Some(100.0 * (now[i].max(1.0) / normal[i].max(1.0)).log10())),
                    1.0,
                )
            }
            ("rain", "observed") => encode_scaled(observed_rain(world, month).into_iter(), 1.0),
            ("rain", "difference") => {
                let obs = observed_rain(world, month);
                encode_scaled(
                    (0..n).map(|i| {
                        let m = c.cell_precip_window(i).map(|r| pick(r, true))?;
                        let o = obs[i]?;
                        Some(100.0 * ((m.max(1.0)) / o.max(1.0)).log10())
                    }),
                    1.0,
                )
            }
            ("rain", _) => encode_scaled(c.map_precip(month).into_iter().map(Some), 1.0),
            ("seaice", _) => {
                encode_scaled(c.map_sea_ice().into_iter().map(|v| Some(v * 100.0)), 1.0)
            }
            ("snow", "current") => encode_scaled(c.map_snow_current().into_iter().map(Some), 1.0),
            ("snow", _) => encode_scaled(c.map_snow(month).into_iter().map(Some), 1.0),
            ("zone", "observed") => {
                encode_scaled(observed_zones(world, &northern(world)).into_iter(), 1.0)
            }
            ("zone", _) => encode_scaled(
                c.map_zones(window, &northern(world))
                    .into_iter()
                    .map(|z| z.map(|z| z.index() as f32)),
                1.0,
            ),
            (_, "current") => encode(c.map_current().into_iter().map(Some)),
            (_, "anomaly") => {
                let now = c.map_current();
                let normal = c.map_normal(Some(c.last_month));
                encode((0..n).map(|i| Some(now[i] - normal[i])))
            }
            (_, "observed") => encode(observed(world, month).into_iter()),
            (_, "difference") => {
                let obs = observed(world, month);
                encode((0..n).map(
                    |i| match (c.cell_window(i).map(|r| pick(r, false)), obs[i]) {
                        (Some(m), Some(o)) => Some(m - o),
                        _ => None,
                    },
                ))
            }
            _ => encode(c.map_normal(month).into_iter().map(Some)),
        }
    };
    ([(header::CONTENT_TYPE, "application/octet-stream")], bytes)
}

/// Like `encode` but values are multiplied by `scale` (not 100).
fn encode_scaled(values: impl Iterator<Item = Option<f32>>, scale: f32) -> Bytes {
    encode(values.map(|v| v.map(|x| x * scale / 100.0)))
}

/// This year's ocean currents: per coarse cell [lat, lon, east m/s, north
/// m/s] (cm/s precision); cells without a current are left out.
pub async fn get_currents(State(shared): State<AppState>) -> Json<serde_json::Value> {
    let sim = shared.sim.lock().unwrap();
    let w: Vec<[f32; 4]> = sim
        .climate
        .ocean
        .currents()
        .into_iter()
        .filter(|c| c.2.is_finite() && c.3.is_finite() && (c.2 != 0.0 || c.3 != 0.0))
        .map(|(a, b, u, v)| {
            [
                a as f32,
                b as f32,
                (u * 1000.0).round() as f32 / 1000.0,
                (v * 1000.0).round() as f32 / 1000.0,
            ]
        })
        .collect();
    Json(serde_json::json!({ "currents": w }))
}

/// This month's winds: per coarse cell [lat, lon, east m/s, north m/s].
pub async fn get_winds(State(shared): State<AppState>) -> Json<serde_json::Value> {
    let sim = shared.sim.lock().unwrap();
    let w: Vec<[f32; 4]> = sim
        .climate
        .atmosphere
        .winds()
        .into_iter()
        .map(|(a, b, u, v)| {
            [
                a as f32,
                b as f32,
                (u * 10.0).round() as f32 / 10.0,
                (v * 10.0).round() as f32 / 10.0,
            ]
        })
        .collect();
    Json(serde_json::json!({ "winds": w }))
}

#[derive(Deserialize)]
pub struct CellQuery {
    cell: usize,
}

pub async fn get_cell(
    State(shared): State<AppState>,
    Query(q): Query<CellQuery>,
) -> Json<serde_json::Value> {
    let sim = shared.sim.lock().unwrap();
    if q.cell >= shared.world.map_grid.len() {
        return Json(serde_json::json!({ "error": "no such cell" }));
    }
    let land = compared_land(&shared.world);
    let compare = sim.scenario.comparison_window().is_some() && land[q.cell] >= 128;
    let find = |cells: &[u32]| cells.iter().position(|&c| c as usize == q.cell);
    let data = &shared.world.data;
    let observed = data
        .observed_temperature
        .as_ref()
        .filter(|_| compare)
        .and_then(|o| find(&o.cells).map(|k| o.monthly_c[k * 12..k * 12 + 12].to_vec()));
    let observed_rain = data
        .observed_precipitation
        .as_ref()
        .filter(|_| compare)
        .and_then(|o| find(&o.cells).map(|k| o.monthly_mm[k * 12..k * 12 + 12].to_vec()));
    let northern = shared.world.map_grid.centers[q.cell][2] >= 0.0;
    let c = &sim.climate;
    let zone = |t: Option<[f32; 12]>, p: Option<[f32; 12]>| {
        Some(climate::koppen::classify(&t?, &p?, northern).code())
    };
    let window = c.cell_window(q.cell);
    let window_rain = c.cell_precip_window(q.cell);
    let years = c.years_simulated();
    let as12 = |v: &Option<Vec<f32>>| v.as_ref().map(|v| <[f32; 12]>::try_from(&v[..]).unwrap());
    Json(serde_json::json!({
        "normal": c.cell_normal(q.cell),
        "window": window,
        "observed": observed,
        "rain_normal": c.cell_precip_normal(q.cell),
        "rain_window": window_rain,
        "rain_observed": observed_rain,
        "snow_normal": c.cell_snow_normal(q.cell),
        "zone": if years > 0 { zone(Some(c.cell_normal(q.cell)), Some(c.cell_precip_normal(q.cell))) } else { None },
        "zone_window": zone(window, window_rain),
        "zone_observed": zone(as12(&observed), as12(&observed_rain)),
        "years": years,
    }))
}

/// Global mean temperature per year: at most ~1,000 points (averaged in
/// groups when the run is longer).
pub async fn get_global(State(shared): State<AppState>) -> Json<serde_json::Value> {
    let sim = shared.sim.lock().unwrap();
    let h = &sim.climate.global_history;
    let group = h.len().div_ceil(1000).max(1);
    let points: Vec<(f64, f64)> = h
        .chunks(group)
        .map(|c| {
            let year = c.iter().map(|p| p.0 as f64).sum::<f64>() / c.len() as f64;
            let t = c.iter().map(|p| p.1 as f64).sum::<f64>() / c.len() as f64;
            (year, t)
        })
        .collect();
    // The El Niño-like cycle in the widest equatorial basin: the last 50
    // years, month by month.
    let e = &sim.climate.enso_history;
    let start = sim.scenario.start_year() as f64;
    let from = e.len().saturating_sub(50 * 12);
    let enso: Vec<(f64, f32)> = (from..e.len())
        .map(|k| (start + k as f64 / 12.0, e[k]))
        .collect();
    let basin = sim
        .climate
        .enso
        .main_basin()
        .map(|b| &sim.climate.enso.basins[b])
        .map(|b| (b.west_lon, b.east_lon, b.width_km));
    Json(serde_json::json!({
        "points": points,
        "start_year": sim.scenario.start_year(),
        "enso": enso,
        "enso_basin": basin,
        "enso_strong_c": climate::variability::ENSO_STRONG_C,
    }))
}

pub async fn get_check(State(shared): State<AppState>) -> Json<serde_json::Value> {
    let sim = shared.sim.lock().unwrap();
    let scores = shared
        .world
        .data
        .observed_temperature
        .as_ref()
        .filter(|_| sim.scenario.comparison_window().is_some())
        .and_then(|obs| {
            climate::check::scores(
                &sim.climate,
                obs,
                &shared.world.map_grid,
                &shared.world.data.epochs[sim.epoch].land_share,
            )
        });
    let rain = (|| {
        sim.scenario.comparison_window()?;
        climate::check::rain_scores(
            &sim.climate,
            shared.world.data.observed_temperature.as_ref()?,
            shared.world.data.observed_precipitation.as_ref()?,
            &shared.world.map_grid,
            &shared.world.data.epochs[sim.epoch].land_share,
        )
    })();
    let ocean = (|| {
        sim.scenario.comparison_window()?;
        climate::check::ocean_scores(
            &sim.climate,
            shared.world.data.observed_temperature.as_ref()?,
            &shared.world.map_grid,
        )
    })();
    let variability = sim
        .scenario
        .comparison_window()
        .and_then(|_| climate::check::variability_scores(&sim.climate, sim.pinatubo_check()));
    Json(serde_json::json!({
        "window": sim.scenario.comparison_window(),
        "scores": scores,
        "rain": rain,
        "ocean": ocean,
        "variability": variability,
        "zones": climate::koppen::Zone::ALL.iter().map(|z| z.code()).collect::<Vec<_>>(),
        "pass_marks": {
            "global_mean_c": climate::check::PASS_GLOBAL_MEAN_C,
            "annual_typical_error_c": climate::check::PASS_ANNUAL_TYPICAL_ERROR_C,
            "seasonal_swing_error_c": climate::check::PASS_SEASONAL_SWING_ERROR_C,
            "rain_correlation": climate::check::PASS_RAIN_CORRELATION,
            "zone_agreement": climate::check::PASS_ZONE_AGREEMENT,
        },
        "observed_global_mean_c": climate::check::OBSERVED_GLOBAL_MEAN_C,
        "observed_global_uncertainty_c": climate::check::OBSERVED_GLOBAL_MEAN_UNCERTAINTY_C,
    }))
}
