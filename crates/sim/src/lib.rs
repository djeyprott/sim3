//! The simulation core: time, the scenario, and every simulated system,
//! stepped one day at a time. Systems run at the cadence they need: the
//! climate once a month, later systems daily, seasonally or yearly.

pub mod scenario;
pub mod time;
pub mod world_context;

use std::sync::Arc;

use climate::{Climate, Orbit, TUNED};
use history::{EventKind, EventLog, Tick};
pub use scenario::Scenario;
use time::Date;
pub use world_context::WorldContext;

/// The orbit changes by well under 1% per century, so the sunshine table
/// is recomputed once every hundred years.
const ORBIT_UPDATE_YEARS: i64 = 100;

pub struct Simulation {
    pub scenario: Scenario,
    tick: Tick,
    pub log: EventLog,
    pub world: Arc<WorldContext>,
    pub epoch: usize,
    pub climate: Climate,
    orbit_year: i64,
}

impl Simulation {
    pub fn new(world: Arc<WorldContext>, scenario: Scenario) -> Self {
        let mut log = EventLog::new();
        log.record(0, EventKind::SimulationStarted, vec![]);
        let start = scenario.start_year();
        let epoch = world.epoch_for_year(start);
        let orbit = orbit_for(&world, scenario, start);
        let mut climate = Climate::new(
            &world.climate_grid,
            &world.map_grid,
            world.map_to_climate.clone(),
            &world.data.epochs[epoch],
            &orbit,
            TUNED,
        );
        if let Some((a, b)) = scenario.comparison_window() {
            climate.set_window(a, b);
        }
        let mut sim = Simulation {
            scenario,
            tick: 0,
            log,
            world,
            epoch,
            climate,
            orbit_year: start,
        };
        sim.start_month();
        sim
    }

    /// Whether a run with an end (the present-day check) has reached it.
    pub fn finished(&self) -> bool {
        self.scenario
            .comparison_window()
            .is_some_and(|(_, last)| self.date().year > last)
    }

    pub fn tick(&self) -> Tick {
        self.tick
    }

    pub fn date(&self) -> Date {
        self.date_of(self.tick)
    }

    pub fn date_of(&self, tick: Tick) -> Date {
        Date::days_after_new_year(self.scenario.start_year(), tick)
    }

    /// Advances the world by one day.
    pub fn step(&mut self) {
        self.tick += 1;
        if self.date().day == 1 {
            self.start_month();
        }
    }

    /// Runs the climate for the month that begins today.
    fn start_month(&mut self) {
        let d = self.date();
        if d.month == 1 && d.year - self.orbit_year >= ORBIT_UPDATE_YEARS {
            let orbit = orbit_for(&self.world, self.scenario, d.year);
            self.climate.set_orbit(&orbit);
            self.orbit_year = d.year;
        }
        let co2 = co2_for(&self.world, self.scenario, d.year);
        let days = time::days_in_month(d.year, d.month) as u32;
        self.climate
            .step_month(d.year, d.month as usize - 1, days, co2);
    }
}

/// Years before 1950 (the reference year of paleoclimate records).
fn years_bp(year: i64) -> f64 {
    (1950 - year) as f64
}

fn orbit_for(world: &WorldContext, scenario: Scenario, year: i64) -> Orbit {
    match scenario {
        Scenario::Holocene if !world.data.orbit.is_empty() => {
            Orbit::at(&world.data.orbit, years_bp(year))
        }
        _ => Orbit::present(),
    }
}

/// CO2 (ppm): the ice-core record for the Holocene run, the observed
/// 1981-2010 average for the present-day check.
fn co2_for(world: &WorldContext, scenario: Scenario, year: i64) -> f64 {
    match scenario {
        Scenario::PresentDayCheck => world
            .data
            .co2_1981_2010_ppm
            .map(|v| v as f64)
            .unwrap_or(climate::model::OLR_CO2_PPM),
        Scenario::Holocene => {
            let rec = &world.data.co2_record;
            if rec.is_empty() {
                return climate::model::OLR_CO2_PPM;
            }
            let age = years_bp(year) as f32;
            let i = rec.partition_point(|&(a, _)| a < age);
            if i == 0 {
                return rec[0].1 as f64;
            }
            if i >= rec.len() {
                return rec[rec.len() - 1].1 as f64;
            }
            let ((a0, c0), (a1, c1)) = (rec[i - 1], rec[i]);
            let t = if a1 > a0 { (age - a0) / (a1 - a0) } else { 0.0 };
            (c0 + t * (c1 - c0)) as f64
        }
    }
}
