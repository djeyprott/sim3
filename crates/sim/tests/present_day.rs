//! Runs the present-day check (1931-2010) on the real world data and
//! reports the scores and the time the climate took. Needs data/world.bin
//! (format v2 with observed temperatures).
//!
//! cargo test -p sim --release --test present_day -- --nocapture

use std::sync::Arc;
use std::time::Instant;

use sim::{Scenario, Simulation, WorldContext};
use world::data::WorldData;

fn world() -> Option<Arc<WorldContext>> {
    let bytes = std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/world.bin")).ok()?;
    let data = WorldData::parse(&bytes).ok()?;
    data.observed_temperature.as_ref()?;
    Some(Arc::new(WorldContext::new(data)))
}

#[test]
fn present_day_check() {
    let Some(world) = world() else {
        eprintln!("skipped: data/world.bin without observed temperatures");
        return;
    };
    let mut sim = Simulation::new(world.clone(), Scenario::PresentDayCheck);
    let start = Instant::now();
    while sim.date().year <= 2010 {
        sim.step();
    }
    let secs = start.elapsed().as_secs_f64();
    let years = sim.climate.years_simulated();
    println!(
        "{years} years in {secs:.2} s = {:.0} simulated years per second",
        years as f64 / secs
    );
    let obs = world.data.observed_temperature.as_ref().unwrap();
    let s = climate::check::scores(
        &sim.climate,
        obs,
        &world.map_grid,
        &world.data.epochs[sim.epoch].land_share,
    )
    .unwrap();
    println!(
        "global {:.2} (obs {:.1}) | land bias {:.2} | annual error {:.2} | swing error {:.2}",
        s.global_mean_model_c,
        s.global_mean_observed_c,
        s.land_bias_c,
        s.annual_typical_error_c,
        s.seasonal_swing_error_c
    );
    for b in &s.bands {
        println!(
            "{:>4}..{:>4}: annual {:>6.1} vs {:>6.1} | swing {:>5.1} vs {:>5.1}",
            b.lat_from, b.lat_to, b.model_c, b.observed_c, b.model_swing_c, b.observed_swing_c
        );
    }
    // Step 2 expectations: the world average passes; the agreed land pass
    // marks are reported but only required from Step 4.
    assert!(
        s.pass_global_mean,
        "global mean {:.2} °C",
        s.global_mean_model_c
    );
}

#[test]
#[ignore]
fn worst_swing_cells() {
    let world = world().unwrap();
    let mut sim = Simulation::new(world.clone(), Scenario::PresentDayCheck);
    while sim.date().year <= 2010 {
        sim.step();
    }
    let obs = world.data.observed_temperature.as_ref().unwrap();
    let land = &world.data.epochs[sim.epoch].land_share;
    let mut errs = vec![];
    for (k, &cell) in obs.cells.iter().enumerate() {
        let cell = cell as usize;
        if land[cell] < 128 {
            continue;
        }
        let m = sim.climate.cell_window(cell).unwrap();
        let o = &obs.monthly_c[k * 12..k * 12 + 12];
        let sw = |v: Vec<f32>| {
            v.iter().cloned().fold(f32::MIN, f32::max) - v.iter().cloned().fold(f32::MAX, f32::min)
        };
        let (ms, os) = (sw(m.to_vec()), sw(o.to_vec()));
        errs.push(((ms - os).abs(), cell, ms, os, m, o.to_vec()));
    }
    errs.sort_by(|a, b| b.0.total_cmp(&a.0));
    let n = errs.len();
    for p in [0.5, 0.9, 0.99] {
        println!("p{}: {:.1}", p, errs[((1.0 - p) * n as f64) as usize].0);
    }
    for e in errs.iter().take(8) {
        let (lat, lon) = world::grid::to_lat_lon(world.map_grid.centers[e.1]);
        println!(
            "{:.1},{:.1}: swing {:.1} vs {:.1}\n  model {:?}\n  obs   {:?}",
            lat,
            lon,
            e.2,
            e.3,
            e.4.map(|v| v.round()),
            e.5.iter().map(|v| v.round()).collect::<Vec<_>>()
        );
    }
}

/// One-time tuning of the four values in `climate::model::Tuning` against
/// today's world: the observed global mean and the land temperature profile
/// from equator to pole (10° bands, CHELSA 1981-2010). Spin-up 30 years,
/// 20 years averaged. Run with --ignored; results recorded in model.rs.
#[test]
#[ignore]
fn tune() {
    use climate::{Climate, Orbit, Tuning};
    let world = world().unwrap();
    let epoch = world.epoch_for_year(1981);
    let obs = world.data.observed_temperature.as_ref().unwrap();
    let co2 = world.data.co2_1981_2010_ppm.unwrap() as f64;
    let land = &world.data.epochs[epoch].land_share;
    let run = |tuning: Tuning| {
        let mut c = Climate::new(
            &world.climate_grid,
            &world.map_grid,
            world.map_to_climate.clone(),
            &world.data.epochs[epoch],
            &Orbit::present(),
            tuning,
        );
        c.set_window(30, 49);
        let days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
        for y in 0..50 {
            for (m, &d) in days.iter().enumerate() {
                c.step_month(y, m, d, co2);
            }
        }
        let s = climate::check::scores(&c, obs, &world.map_grid, land).unwrap();
        let (mut sq, mut w) = (0.0, 0.0);
        for b in &s.bands {
            let wt = (b.lat_to.to_radians().sin() - b.lat_from.to_radians().sin()).abs();
            sq += wt * (b.model_c - b.observed_c).powi(2);
            w += wt;
        }
        let profile = (sq / w).sqrt();
        let global = s.global_mean_model_c - s.global_mean_observed_c;
        (global.abs() + profile, global, profile, s)
    };
    let mut best: Option<(f64, Tuning)> = None;
    for a0 in [0.30, 0.32, 0.34, 0.36] {
        for a2 in [0.05, 0.10, 0.15, 0.20, 0.25] {
            for ice in [0.62, 0.70, 0.78] {
                for d in [0.4, 0.5, 0.6, 0.7] {
                    let t = Tuning {
                        albedo_a0: a0,
                        albedo_a2: a2,
                        albedo_ice: ice,
                        diffusion: d,
                    };
                    let (score, g, p, _) = run(t);
                    if best.is_none_or(|b| score < b.0) {
                        println!("{t:?}: global {g:+.2}, profile {p:.2}, score {score:.2}");
                        best = Some((score, t));
                    }
                }
            }
        }
    }
    let (_, t) = best.unwrap();
    let (_, g, p, s) = run(t);
    println!(
        "BEST {t:?}: global {g:+.2} profile {p:.2} annual {:.2} swing {:.2} bias {:.2}",
        s.annual_typical_error_c, s.seasonal_swing_error_c, s.land_bias_c
    );
    for b in &s.bands {
        println!(
            "{:>4}..{:>4}: {:>6.1} vs {:>6.1}",
            b.lat_from, b.lat_to, b.model_c, b.observed_c
        );
    }
}

#[test]
#[ignore]
fn profile_speed() {
    use std::time::Instant;
    let world = world().unwrap();
    let mut sim = Simulation::new(world.clone(), Scenario::Holocene);
    for _ in 0..365 * 5 {
        sim.step();
    }
    let t = Instant::now();
    for _ in 0..365 * 20 {
        sim.step();
    }
    println!(
        "whole sim: {:.2} ms per year",
        t.elapsed().as_secs_f64() * 1e3 / 20.0
    );
    let t = Instant::now();
    let mut passes = 0;
    for m in 0..240 {
        sim.climate.step_month(0, m % 12, 30, 265.0);
        passes += sim.climate.solver_passes();
    }
    println!("average solver passes {:.1}", passes as f64 / 240.0);
    println!(
        "climate.step_month: {:.3} ms per month",
        t.elapsed().as_secs_f64() * 1e3 / 240.0
    );
}
