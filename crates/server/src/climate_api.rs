//! Climate data for the viewer: temperature maps, a place's climate chart,
//! the world temperature history and the present-day check scores.

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

/// Observed 1981-2010 temperature per map cell for a month (or the year),
/// on land only (cells at least half land today): the agreed comparison.
fn observed(world: &sim::WorldContext, month: Option<usize>) -> Vec<Option<f32>> {
    let mut out = vec![None; world.map_grid.len()];
    let today = world.epoch_for_year(2000);
    let land = &world.data.epochs[today].land_share;
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

/// Per map cell, i16 little-endian, hundredths of °C; -32768 = no value.
pub async fn get_map(
    State(shared): State<AppState>,
    Query(q): Query<MapQuery>,
) -> impl IntoResponse {
    let month = q.month.filter(|&m| m < 12);
    let bytes = {
        let sim = shared.sim.lock().unwrap();
        let c = &sim.climate;
        match q.view.as_str() {
            "current" => encode(c.map_current().into_iter().map(Some)),
            "observed" => encode(observed(&shared.world, month).into_iter()),
            "difference" => {
                let obs = observed(&shared.world, month);
                let n = obs.len();
                let model: Vec<Option<f32>> = (0..n)
                    .map(|i| {
                        c.cell_window(i).map(|row| {
                            month
                                .map(|m| row[m])
                                .unwrap_or(row.iter().sum::<f32>() / 12.0)
                        })
                    })
                    .collect();
                encode((0..n).map(|i| match (model[i], obs[i]) {
                    (Some(m), Some(o)) => Some(m - o),
                    _ => None,
                }))
            }
            _ => encode(c.map_normal(month).into_iter().map(Some)),
        }
    };
    ([(header::CONTENT_TYPE, "application/octet-stream")], bytes)
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
    let land = &shared.world.data.epochs[shared.world.epoch_for_year(2000)].land_share;
    let observed = shared
        .world
        .data
        .observed_temperature
        .as_ref()
        .filter(|_| land[q.cell] >= 128)
        .and_then(|obs| {
            obs.cells
                .iter()
                .position(|&c| c as usize == q.cell)
                .map(|k| obs.monthly_c[k * 12..k * 12 + 12].to_vec())
        });
    Json(serde_json::json!({
        "normal": sim.climate.cell_normal(q.cell),
        "window": sim.climate.cell_window(q.cell),
        "observed": if sim.scenario.comparison_window().is_some() { observed } else { None },
        "years": sim.climate.years_simulated(),
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
    Json(serde_json::json!({ "points": points, "start_year": sim.scenario.start_year() }))
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
    Json(serde_json::json!({
        "window": sim.scenario.comparison_window(),
        "scores": scores,
        "pass_marks": {
            "global_mean_c": climate::check::PASS_GLOBAL_MEAN_C,
            "annual_typical_error_c": climate::check::PASS_ANNUAL_TYPICAL_ERROR_C,
            "seasonal_swing_error_c": climate::check::PASS_SEASONAL_SWING_ERROR_C,
        },
        "observed_global_mean_c": climate::check::OBSERVED_GLOBAL_MEAN_C,
        "observed_global_uncertainty_c": climate::check::OBSERVED_GLOBAL_MEAN_UNCERTAINTY_C,
    }))
}
