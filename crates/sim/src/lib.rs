//! The simulation core.
//!
//! Phase 0 only advances time and records that the run began. Each later
//! phase adds a system (climate, ecology, populations, ...) that is stepped
//! from [`Simulation::step`].

pub mod time;

use history::{EventKind, EventLog, Tick};
use time::Date;

pub struct Simulation {
    tick: Tick,
    pub log: EventLog,
}

impl Default for Simulation {
    fn default() -> Self {
        Self::new()
    }
}

impl Simulation {
    pub fn new() -> Self {
        let mut log = EventLog::new();
        log.record(0, EventKind::SimulationStarted, vec![]);
        Simulation { tick: 0, log }
    }

    pub fn tick(&self) -> Tick {
        self.tick
    }

    pub fn date(&self) -> Date {
        Date::from_tick(self.tick)
    }

    /// Advances the world by one day.
    ///
    /// Systems run at the cadence they need, not all every day: fast ones
    /// (weather, people's actions) daily, slow ones (births and deaths,
    /// harvests, ice sheets, sea level) on the first day of a season or
    /// year, using the calendar date to decide.
    pub fn step(&mut self) {
        self.tick += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_run_records_its_start() {
        let sim = Simulation::new();
        assert_eq!(sim.log.len(), 1);
        assert_eq!(sim.tick(), 0);
    }

    #[test]
    fn the_first_year_is_365_steps() {
        let mut sim = Simulation::new();
        for _ in 0..365 {
            sim.step();
        }
        assert_eq!(sim.date().year_label(), "9999 BCE");
    }
}
