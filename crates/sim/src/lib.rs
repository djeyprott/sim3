//! The simulation core: time, the scenario, and every simulated system,
//! stepped one day at a time. Systems run at the cadence they need: the
//! climate once a month, later systems daily, seasonally or yearly.

pub mod scenario;
pub mod time;
pub mod world_context;

use std::sync::Arc;

use climate::variability::{EruptionStats, PINATUBO_TGS};
use climate::{Climate, Orbit, ATMOSPHERE_TUNED, TUNED};
use history::{EnsoReport, EventId, EventKind, EventLog, Tick};
pub use scenario::Scenario;
use time::Date;
pub use world_context::WorldContext;

/// The orbit changes by well under 1% per century, so the sunshine table
/// is recomputed once every hundred years.
const ORBIT_UPDATE_YEARS: i64 = 100;

/// Eruptions at least Pinatubo's size go into the Chronicle (smaller ones
/// still dim the sky). Tg S.
const CHRONICLE_ERUPTION_TGS: f64 = PINATUBO_TGS;
/// A year is noticeably cold when the world's average is more than two
/// standard deviations below the previous thirty years' (the usual test
/// for an exceptional year); the cooling is put down to an eruption in the
/// Chronicle within the last three years (a veil lasts about two, the sea
/// takes a little longer to recover).
const COLD_YEAR_SIGMAS: f64 = 2.0;
const COLD_YEAR_BASELINE_YEARS: usize = 30;
const COOLING_CAUSE_YEARS: i64 = 3;

pub struct Simulation {
    pub scenario: Scenario,
    tick: Tick,
    pub log: EventLog,
    pub world: Arc<WorldContext>,
    pub epoch: usize,
    pub climate: Climate,
    orbit_year: i64,
    /// Eruptions in the Chronicle: (event, year, sulfur Tg S, explained a
    /// cold year already).
    recent_eruptions: Vec<(EventId, i64, f64, bool)>,
    /// The Pinatubo test result, worked out once a check run has finished.
    pinatubo: std::sync::OnceLock<f64>,
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
            ATMOSPHERE_TUNED,
        );
        if let Some((a, b)) = scenario.comparison_window() {
            climate.set_window(a, b);
        }
        climate.set_variability(eruption_stats(&world, scenario), run_seed(scenario));
        let mut sim = Simulation {
            scenario,
            tick: 0,
            log,
            world,
            epoch,
            climate,
            orbit_year: start,
            recent_eruptions: Vec::new(),
            pinatubo: std::sync::OnceLock::new(),
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
        let years_before = self.climate.global_history.len();
        self.climate
            .step_month(d.year, d.month as usize - 1, days, co2);
        self.log_climate_events(d.year);
        if self.climate.global_history.len() > years_before {
            self.check_cold_year();
        }
    }

    /// How much a Pinatubo-sized eruption cools the world: the climate is
    /// run on twice for two years from now, with and without one (9 Tg S
    /// at 15°N, Pinatubo's latitude), and the world averages compared over
    /// months 7-18 after it (the year after). Same random history in both.
    pub fn pinatubo_response(&self) -> f64 {
        let mut with = self.climate.clone();
        let mut without = self.climate.clone();
        with.volcanism.inject(15.0, PINATUBO_TGS);
        let mut d = self.date();
        let co2 = co2_for(&self.world, self.scenario, d.year);
        let mut diff = 0.0;
        for k in 0..18 {
            let days = time::days_in_month(d.year, d.month) as u32;
            let m = d.month as usize - 1;
            with.step_month(d.year, m, days, co2);
            without.step_month(d.year, m, days, co2);
            if k >= 6 {
                diff += without.global_mean_now() - with.global_mean_now();
            }
            d = if d.month == 12 {
                Date {
                    year: d.year + 1,
                    month: 1,
                    day: 1,
                }
            } else {
                Date {
                    month: d.month + 1,
                    day: 1,
                    ..d
                }
            };
        }
        diff / 12.0
    }

    /// The Pinatubo test, once a run with an end has finished.
    pub fn pinatubo_check(&self) -> Option<f64> {
        self.finished()
            .then(|| *self.pinatubo.get_or_init(|| self.pinatubo_response()))
    }

    /// Puts this month's eruptions and El Niño years into the Chronicle.
    fn log_climate_events(&mut self, year: i64) {
        for e in std::mem::take(&mut self.climate.eruptions) {
            if e.sulfur_tgs >= CHRONICLE_ERUPTION_TGS {
                let kind = EventKind::VolcanicEruption {
                    lat: e.lat as f32,
                    lon: e.lon as f32,
                    sulfur_tg: e.sulfur_tgs as f32,
                };
                let id = self.log.record(self.tick, kind, vec![]);
                self.recent_eruptions.push((id, year, e.sulfur_tgs, false));
            }
        }
        self.recent_eruptions
            .retain(|&(_, y, _, _)| year - y <= COOLING_CAUSE_YEARS);
        for ev in std::mem::take(&mut self.climate.enso_events) {
            let Some(((dlat, dlon, dpct), (wlat, wlon, wpct))) = self.climate.rain_extremes()
            else {
                continue;
            };
            let b = &self.climate.enso.basins[ev.basin];
            let mid = 0.5 * (b.west_lon + b.east_lon);
            let report = EnsoReport {
                basin_lon: (((mid + 180.0).rem_euclid(360.0)) - 180.0) as f32,
                strength_c: ev.strength_c as f32,
                dry_lat: dlat as f32,
                dry_lon: dlon as f32,
                dry_percent: dpct as f32,
                wet_lat: wlat as f32,
                wet_lon: wlon as f32,
                wet_percent: wpct as f32,
            };
            let kind = if ev.strength_c > 0.0 {
                EventKind::ElNino(report)
            } else {
                EventKind::LaNina(report)
            };
            self.log.record(self.tick, kind, vec![]);
        }
    }

    /// After a year ends: if it was exceptionally cold and a big eruption
    /// came shortly before, the cooling goes into the Chronicle.
    fn check_cold_year(&mut self) {
        let h = &self.climate.global_history;
        let n = h.len();
        if n <= COLD_YEAR_BASELINE_YEARS {
            return;
        }
        let base: Vec<f64> = h[n - 1 - COLD_YEAR_BASELINE_YEARS..n - 1]
            .iter()
            .map(|&(_, t)| t as f64)
            .collect();
        let mean = base.iter().sum::<f64>() / base.len() as f64;
        let sd = (base.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / base.len() as f64).sqrt();
        let cooling = mean - h[n - 1].1 as f64;
        if cooling <= COLD_YEAR_SIGMAS * sd {
            return;
        }
        let cause = self
            .recent_eruptions
            .iter_mut()
            .filter(|e| !e.3)
            .max_by(|a, b| a.2.total_cmp(&b.2));
        if let Some(e) = cause {
            e.3 = true;
            let id = e.0;
            self.log.record(
                self.tick,
                EventKind::WorldCooled {
                    cooling_c: cooling as f32,
                },
                vec![id],
            );
        }
    }
}

/// How often and how big eruptions are: the recorded eruptions for the
/// Holocene run. The present-day check has none, since its observed
/// averages already include the eruptions of 1981-2010 and today's
/// tuning was made without them.
fn eruption_stats(world: &WorldContext, scenario: Scenario) -> Option<EruptionStats> {
    match scenario {
        Scenario::Holocene => world.data.eruptions.as_ref().map(|r| EruptionStats {
            sizes_tgs: r.sizes_tgs.clone(),
            per_month: r.sizes_tgs.len() as f64 / (r.span_years as f64 * 12.0),
        }),
        Scenario::PresentDayCheck => None,
    }
}

/// Each Holocene run gets its own random history; the present-day check
/// is the same every time, so its scores can be compared.
fn run_seed(scenario: Scenario) -> u64 {
    match scenario {
        Scenario::Holocene => std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(1),
        Scenario::PresentDayCheck => 1,
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
