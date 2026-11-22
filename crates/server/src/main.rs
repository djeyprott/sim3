//! `sim3`: runs the simulation on a background thread and serves the viewer
//! at http://127.0.0.1:<port> for the browser on this computer.
//!
//! Usage: sim3 [--port N] [--no-browser]

mod map;

use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use axum::extract::{Query, State};
use axum::http::header;
use axum::response::{Html, IntoResponse};
use axum::routing::{get, post};
use axum::{Json, Router};
use serde::{Deserialize, Serialize};
use sim::time::Date;
use sim::Simulation;

const DEFAULT_PORT: u16 = 7340;

/// Playback speeds offered in the viewer, in simulated days per real
/// second: a day, a month, a year, a decade and a century per second.
/// 0 means "as fast as the computer can go".
const SPEEDS: [u32; 6] = [1, 30, 365, 3_652, 36_524, 0];

/// Largest single step the viewer can request in days: 100 years.
const MAX_STEP_DAYS: u64 = 36_525;

struct Control {
    paused: bool,
    days_per_second: u32,
    pending_steps: u64,
}

struct Shared {
    sim: Mutex<Simulation>,
    control: Mutex<Control>,
    map: map::MapData,
}

type AppState = Arc<Shared>;

fn run_simulation(shared: AppState) {
    let mut budget = 0.0_f64;
    let mut last = Instant::now();
    loop {
        let (paused, dps, steps) = {
            let mut c = shared.control.lock().unwrap();
            let s = c.pending_steps;
            c.pending_steps = 0;
            (c.paused, c.days_per_second, s)
        };
        let now = Instant::now();
        let elapsed = now.duration_since(last).as_secs_f64();
        last = now;

        if steps > 0 {
            let mut sim = shared.sim.lock().unwrap();
            for _ in 0..steps {
                sim.step();
            }
        }

        if paused {
            budget = 0.0;
            thread::sleep(Duration::from_millis(20));
            continue;
        }

        if dps == 0 {
            // Flat out: hold the lock in short slices so the viewer can still read.
            let slice_end = Instant::now() + Duration::from_millis(10);
            let mut sim = shared.sim.lock().unwrap();
            while Instant::now() < slice_end {
                sim.step();
            }
            drop(sim);
            thread::yield_now();
        } else {
            budget += elapsed * dps as f64;
            // Never try to catch up more than one second's worth at once.
            budget = budget.min(dps as f64);
            let due = budget.floor() as u32;
            if due > 0 {
                budget -= due as f64;
                let mut sim = shared.sim.lock().unwrap();
                for _ in 0..due {
                    sim.step();
                }
            }
            thread::sleep(Duration::from_millis(5));
        }
    }
}

#[derive(Serialize)]
struct StateResponse {
    tick: u64,
    date: Date,
    years_elapsed: i64,
    date_label: String,
    year_label: String,
    paused: bool,
    days_per_second: u32,
    speeds: &'static [u32],
    event_count: usize,
}

async fn get_state(State(shared): State<AppState>) -> Json<StateResponse> {
    let (paused, dps) = {
        let c = shared.control.lock().unwrap();
        (c.paused, c.days_per_second)
    };
    let sim = shared.sim.lock().unwrap();
    let date = sim.date();
    Json(StateResponse {
        tick: sim.tick(),
        date,
        years_elapsed: date.years_since_start(),
        date_label: date.label(),
        year_label: date.year_label(),
        paused,
        days_per_second: dps,
        speeds: &SPEEDS,
        event_count: sim.log.len(),
    })
}

#[derive(Deserialize)]
#[serde(tag = "action", rename_all = "snake_case")]
enum ControlRequest {
    Pause,
    Resume,
    /// Advance by `count` calendar days, months or years.
    Step {
        unit: StepUnit,
        count: Option<u32>,
    },
    Speed {
        days_per_second: u32,
    },
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
enum StepUnit {
    Day,
    Month,
    Year,
}

async fn post_control(
    State(shared): State<AppState>,
    Json(req): Json<ControlRequest>,
) -> impl IntoResponse {
    let current_tick = shared.sim.lock().unwrap().tick();
    let mut c = shared.control.lock().unwrap();
    match req {
        ControlRequest::Pause => c.paused = true,
        ControlRequest::Resume => c.paused = false,
        ControlRequest::Step { unit, count } => {
            // Clamp the request to at most 100 years, in its own unit, so the
            // step still lands on a matching calendar date.
            let count = count.unwrap_or(1);
            let from_tick = current_tick + c.pending_steps;
            let from = Date::from_tick(from_tick);
            let to = match unit {
                StepUnit::Day => from_tick + (count as u64).min(MAX_STEP_DAYS),
                StepUnit::Month => from.add_months(count.min(1_200)).to_tick(),
                StepUnit::Year => from.add_years(count.min(100)).to_tick(),
            };
            c.pending_steps += to - from_tick;
        }
        ControlRequest::Speed { days_per_second } => {
            if SPEEDS.contains(&days_per_second) {
                c.days_per_second = days_per_second;
            }
        }
    }
    Json(serde_json::json!({ "ok": true }))
}

#[derive(Deserialize)]
struct EventsQuery {
    from: Option<u64>,
    limit: Option<usize>,
}

#[derive(Serialize)]
struct EventView {
    #[serde(flatten)]
    event: history::Event,
    date_label: String,
}

async fn get_events(
    State(shared): State<AppState>,
    Query(q): Query<EventsQuery>,
) -> Json<Vec<EventView>> {
    let sim = shared.sim.lock().unwrap();
    let events = sim
        .log
        .from_id(q.from.unwrap_or(0), q.limit.unwrap_or(200).min(1000))
        .iter()
        .map(|e| EventView {
            date_label: Date::from_tick(e.tick).label(),
            event: e.clone(),
        })
        .collect();
    Json(events)
}

async fn get_map_info(State(shared): State<AppState>) -> Json<serde_json::Value> {
    Json(serde_json::to_value(&shared.map.info).unwrap())
}

fn binary(bytes: axum::body::Bytes) -> impl IntoResponse {
    ([(header::CONTENT_TYPE, "application/octet-stream")], bytes)
}

async fn get_map_raster(State(shared): State<AppState>) -> impl IntoResponse {
    binary(shared.map.raster.clone())
}

#[derive(Deserialize)]
struct EpochQuery {
    epoch: usize,
}

async fn get_map_cells(
    State(shared): State<AppState>,
    Query(q): Query<EpochQuery>,
) -> axum::response::Response {
    match shared.map.cell_blobs.get(q.epoch) {
        Some(blob) => binary(blob.clone()).into_response(),
        None => (axum::http::StatusCode::NOT_FOUND, "no such epoch").into_response(),
    }
}

async fn index() -> Html<&'static str> {
    Html(include_str!("../../../viewer/index.html"))
}

async fn app_js() -> impl IntoResponse {
    (
        [(header::CONTENT_TYPE, "text/javascript; charset=utf-8")],
        include_str!("../../../viewer/app.js"),
    )
}

async fn map_js() -> impl IntoResponse {
    (
        [(header::CONTENT_TYPE, "text/javascript; charset=utf-8")],
        include_str!("../../../viewer/map.js"),
    )
}

async fn style_css() -> impl IntoResponse {
    (
        [(header::CONTENT_TYPE, "text/css; charset=utf-8")],
        include_str!("../../../viewer/style.css"),
    )
}

fn open_browser(url: &str) {
    #[cfg(target_os = "windows")]
    let result = std::process::Command::new("cmd")
        .args(["/C", "start", "", url])
        .spawn();
    #[cfg(target_os = "macos")]
    let result = std::process::Command::new("open").arg(url).spawn();
    #[cfg(all(unix, not(target_os = "macos")))]
    let result = std::process::Command::new("xdg-open").arg(url).spawn();
    if result.is_err() {
        println!("Open {url} in your browser.");
    }
}

#[tokio::main]
async fn main() {
    let mut port = DEFAULT_PORT;
    let mut browser = true;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--port" => {
                port = args
                    .next()
                    .and_then(|p| p.parse().ok())
                    .expect("--port needs a number")
            }
            "--no-browser" => browser = false,
            other => {
                eprintln!("Unknown option: {other}\nUsage: sim3 [--port N] [--no-browser]");
                std::process::exit(2);
            }
        }
    }

    println!("Preparing the world map...");
    let started = Instant::now();
    let map = map::MapData::build();
    println!("Map ready in {:.1} s.", started.elapsed().as_secs_f64());

    let shared: AppState = Arc::new(Shared {
        map,
        sim: Mutex::new(Simulation::new()),
        control: Mutex::new(Control {
            paused: true,
            days_per_second: 365,
            pending_steps: 0,
        }),
    });

    {
        let shared = shared.clone();
        thread::spawn(move || run_simulation(shared));
    }

    let app = Router::new()
        .route("/", get(index))
        .route("/app.js", get(app_js))
        .route("/map.js", get(map_js))
        .route("/style.css", get(style_css))
        .route("/api/state", get(get_state))
        .route("/api/control", post(post_control))
        .route("/api/events", get(get_events))
        .route("/api/map/info", get(get_map_info))
        .route("/api/map/raster", get(get_map_raster))
        .route("/api/map/cells", get(get_map_cells))
        .with_state(shared);

    // Only listen on this computer, never on the network.
    let listener = match tokio::net::TcpListener::bind(("127.0.0.1", port)).await {
        Ok(l) => l,
        Err(e) => {
            eprintln!("Could not use port {port}: {e}\nTry: sim3 --port 7341");
            std::process::exit(1);
        }
    };
    let url = format!("http://127.0.0.1:{port}");
    println!("sim3 is running. Viewer: {url}");
    println!("Close this window or press Ctrl+C to stop.");
    if browser {
        open_browser(&url);
    }

    axum::serve(listener, app)
        .with_graceful_shutdown(async {
            let _ = tokio::signal::ctrl_c().await;
        })
        .await
        .unwrap();
}
