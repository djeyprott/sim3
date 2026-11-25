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
    println!(
        "gyres (Sv): north subtropics {:.0}, south subtropics {:.0}, north subpolar {:.0}",
        sim.climate.ocean.gyre_strength_sv(15.0, 40.0),
        sim.climate.ocean.gyre_strength_sv(-40.0, -15.0),
        sim.climate.ocean.gyre_strength_sv(45.0, 65.0)
    );
    for b in &sim.climate.ocean.basins {
        println!(
            "overturning basin at {:.0},{:.0}: {:.1} Sv, {:.2} PW, salinity diff {:.2}, P-E warm {:.2} cold {:.2} m/yr, T {:.1} / {:.1}",
            b.name_lat_lon.0,
            b.name_lat_lon.1,
            b.flow / 1e6,
            b.heat_w / 1e15,
            b.salinity_diff,
            b.fresh_m_per_year.0,
            b.fresh_m_per_year.1,
            b.temperatures.0,
            b.temperatures.1
        );
    }
    println!("deep water {:.1} °C", sim.climate.ocean.deep_temperature);
    let o = climate::check::ocean_scores(&sim.climate, obs, &world.map_grid).unwrap();
    for i in &o.sea_ice {
        println!(
            "sea ice {}: {:.1} vs {:.1} million km² {}",
            i.name,
            i.model_million_km2,
            i.observed_million_km2,
            if i.pass { "pass" } else { "MISS" }
        );
    }
    for c in &o.contrasts {
        println!(
            "{}: {:.1} vs {:.1} °C {}",
            c.name,
            c.model_c,
            c.observed_c,
            if c.pass { "pass" } else { "MISS" }
        );
    }
    let v = climate::check::variability_scores(&sim.climate, sim.pinatubo_check()).unwrap();
    let pf = |b: bool| if b { "pass" } else { "MISS" };
    println!(
        "Pinatubo cooling next year: {:.2} °C (observed {:.1}) {}",
        v.pinatubo_cooling_c.unwrap_or(f64::NAN),
        v.pinatubo_observed_c,
        pf(v.pinatubo_pass)
    );
    println!(
        "El Niño size {:.2} °C (observed {:.1}) {}; cycle {:.1} years {}",
        v.enso_size_c,
        v.enso_observed_size_c,
        pf(v.enso_size_pass),
        v.enso_period_years.unwrap_or(f64::NAN),
        pf(v.enso_period_pass)
    );
    println!(
        "El Niño months rain: far west {:.0}% {}, east coast {:.0}% {}",
        v.west_rain_percent.unwrap_or(f64::NAN),
        pf(v.west_pass),
        v.east_rain_percent.unwrap_or(f64::NAN),
        pf(v.east_pass)
    );
    for b in &sim.climate.enso.basins {
        println!(
            "El Niño basin {:.0}..{:.0}°E, {:.0} km",
            b.west_lon, b.east_lon, b.width_km
        );
    }
    if let Some(op) = world.data.observed_precipitation.as_ref() {
        let r = climate::check::rain_scores(
            &sim.climate,
            obs,
            op,
            &world.map_grid,
            &world.data.epochs[sim.epoch].land_share,
        )
        .unwrap();
        println!(
            "rain correlation {:.2} | land mean {:.0} vs {:.0} mm | zones agree {:.0}%",
            r.correlation,
            r.land_mean_model_mm,
            r.land_mean_observed_mm,
            100.0 * r.zone_agreement
        );
        for (z, (m, o)) in climate::koppen::Zone::ALL
            .iter()
            .zip(r.zone_share_model.iter().zip(&r.zone_share_observed))
        {
            println!(
                "  {:>2}: {:>4.1}% vs {:>4.1}%",
                z.code(),
                100.0 * m,
                100.0 * o
            );
        }
        for l in &r.landmarks {
            println!(
                "  {:>20} ({}): {:.2} vs {:.2} {}",
                l.name,
                l.expect,
                l.model,
                l.observed,
                if l.pass { "pass" } else { "MISS" }
            );
        }
    }
    // The world average must pass; the agreed land, rain, sea ice and
    // ocean pass marks are reported and become required in the climate
    // polish step (see the plan).
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
    use climate::{Climate, Orbit, Tuning, ATMOSPHERE_TUNED};
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
            ATMOSPHERE_TUNED,
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
        // Score: the world average plus the two typical land errors the
        // agreed marks use (yearly average and summer-winter swing).
        (
            global.abs() + s.annual_typical_error_c + s.seasonal_swing_error_c,
            global,
            profile,
            s,
        )
    };
    // Coordinate search from the current values: try each value one step up
    // and down (both at once, on two threads), keep any improvement, and
    // halve the steps when nothing improves.
    let get = |t: &Tuning, k: usize| [t.albedo_a0, t.albedo_a2, t.albedo_ice, t.diffusion][k];
    let set = |t: &mut Tuning, k: usize, v: f64| match k {
        0 => t.albedo_a0 = v,
        1 => t.albedo_a2 = v,
        2 => t.albedo_ice = v,
        _ => t.diffusion = v,
    };
    let mut steps = [0.02, 0.05, 0.08, 0.1];
    let mut current = climate::TUNED;
    let mut current_score = run(current).0;
    println!("start {current:?}: score {current_score:.2}");
    for _round in 0..3 {
        loop {
            let mut improved = false;
            for (k, &step) in steps.iter().enumerate() {
                let candidates: Vec<Tuning> = [-1.0, 1.0]
                    .iter()
                    .map(|sgn| {
                        let mut t = current;
                        set(&mut t, k, (get(&current, k) + sgn * step).max(0.0));
                        t
                    })
                    .collect();
                let results: Vec<(f64, f64, f64)> = std::thread::scope(|sc| {
                    let hs: Vec<_> = candidates
                        .iter()
                        .map(|&t| {
                            sc.spawn(move || {
                                let r = run(t);
                                (r.0, r.1, r.2)
                            })
                        })
                        .collect();
                    hs.into_iter().map(|h| h.join().unwrap()).collect()
                });
                for (t, r) in candidates.iter().zip(&results) {
                    if r.0 < current_score {
                        current = *t;
                        current_score = r.0;
                        improved = true;
                        println!(
                            "{t:?}: global {:+.2}, profile {:.2}, score {:.2}",
                            r.1, r.2, r.0
                        );
                    }
                }
            }
            if !improved {
                break;
            }
        }
        steps.iter_mut().for_each(|s| *s /= 2.0);
    }
    let t = current;
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
    // Measure after the oceans have settled (about 35 years from the
    // start): that is the speed of nearly all of a long run.
    for _ in 0..365 * 40 {
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
    let before = sim.climate.timing;
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
    let after = sim.climate.timing;
    let [w, h, m] = [0, 1, 2].map(|k| (after[k] - before[k]) * 1e3 / 240.0);
    println!(
        "moisture passes last month {}",
        sim.climate.atmosphere.moisture_sweeps
    );
    println!("  of which winds {w:.2} ms, temperature {h:.2} ms, rain and snow {m:.2} ms");
}
