//! `sim3`: runs the simulation on a background thread and serves the viewer
//! at http://127.0.0.1:<port> for the browser on this computer.
//!
//! Usage: sim3 [--port N] [--no-browser]

use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use axum::extract::{Query, State};
use axum::http::header;
use axum::response::{Html, IntoResponse};
use axum::routing::{get, post};
use axum::{Json, Router};
use serde::{Deserialize, Serialize};
use sim::Simulation;

const DEFAULT_PORT: u16 = 7340;

/// Playback speeds offered in the viewer, in simulated months per real
/// second. 0 means "as fast as the computer can go".
const SPEEDS: [u32; 5] = [1, 12, 120, 1200, 0];

struct Control {
    paused: bool,
    months_per_second: u32,
    pending_steps: u32,
}

struct Shared {
    sim: Mutex<Simulation>,
    control: Mutex<Control>,
}

type AppState = Arc<Shared>;

fn run_simulation(shared: AppState) {
    let mut budget = 0.0_f64;
    let mut last = Instant::now();
    loop {
        let (paused, mps, steps) = {
            let mut c = shared.control.lock().unwrap();
            let s = c.pending_steps;
            c.pending_steps = 0;
            (c.paused, c.months_per_second, s)
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

        if mps == 0 {
            // Flat out: hold the lock in short slices so the viewer can still read.
            let slice_end = Instant::now() + Duration::from_millis(10);
            let mut sim = shared.sim.lock().unwrap();
            while Instant::now() < slice_end {
                sim.step();
            }
            drop(sim);
            thread::yield_now();
        } else {
            budget += elapsed * mps as f64;
            // Never try to catch up more than one second's worth at once.
            budget = budget.min(mps as f64);
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
    date: sim::time::Date,
    date_label: String,
    year_label: String,
    paused: bool,
    months_per_second: u32,
    speeds: &'static [u32],
    event_count: usize,
}

async fn get_state(State(shared): State<AppState>) -> Json<StateResponse> {
    let (paused, mps) = {
        let c = shared.control.lock().unwrap();
        (c.paused, c.months_per_second)
    };
    let sim = shared.sim.lock().unwrap();
    let date = sim.date();
    Json(StateResponse {
        tick: sim.tick(),
        date,
        date_label: date.label(),
        year_label: date.year_label(),
        paused,
        months_per_second: mps,
        speeds: &SPEEDS,
        event_count: sim.log.len(),
    })
}

#[derive(Deserialize)]
#[serde(tag = "action", rename_all = "snake_case")]
enum ControlRequest {
    Pause,
    Resume,
    Step { months: Option<u32> },
    Speed { months_per_second: u32 },
}

async fn post_control(
    State(shared): State<AppState>,
    Json(req): Json<ControlRequest>,
) -> impl IntoResponse {
    let mut c = shared.control.lock().unwrap();
    match req {
        ControlRequest::Pause => c.paused = true,
        ControlRequest::Resume => c.paused = false,
        ControlRequest::Step { months } => {
            c.pending_steps = c
                .pending_steps
                .saturating_add(months.unwrap_or(1).min(1200))
        }
        ControlRequest::Speed { months_per_second } => {
            if SPEEDS.contains(&months_per_second) {
                c.months_per_second = months_per_second;
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
            date_label: sim::time::Date::from_tick(e.tick).label(),
            event: e.clone(),
        })
        .collect();
    Json(events)
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

    let shared: AppState = Arc::new(Shared {
        sim: Mutex::new(Simulation::new()),
        control: Mutex::new(Control {
            paused: true,
            months_per_second: 12,
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
        .route("/style.css", get(style_css))
        .route("/api/state", get(get_state))
        .route("/api/control", post(post_control))
        .route("/api/events", get(get_events))
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
