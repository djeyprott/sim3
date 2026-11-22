//! Simulation time.
//!
//! The simulation advances one month per tick. Calendar years use
//! astronomical numbering internally (year 0 exists, 1 BCE = year 0,
//! 10,000 BCE = year -9999) so arithmetic stays simple; conversion to
//! BCE/CE happens only for display, where there is no year zero.

use history::Tick;
use serde::Serialize;

pub const MONTHS_PER_YEAR: u64 = 12;

/// Astronomical year in which the simulation starts: 10,000 BCE.
pub const START_YEAR: i64 = -9999;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct Date {
    /// Astronomical year (can be zero or negative).
    pub year: i64,
    /// 1..=12
    pub month: u8,
}

impl Date {
    pub fn from_tick(tick: Tick) -> Self {
        let years = (tick / MONTHS_PER_YEAR) as i64;
        let month = (tick % MONTHS_PER_YEAR) as u8 + 1;
        Date {
            year: START_YEAR + years,
            month,
        }
    }

    /// Year as people write it: "10000 BCE", "1 CE".
    pub fn year_label(&self) -> String {
        if self.year <= 0 {
            format!("{} BCE", 1 - self.year)
        } else {
            format!("{} CE", self.year)
        }
    }

    pub fn label(&self) -> String {
        format!("Month {} of {}", self.month, self.year_label())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn starts_in_10000_bce() {
        let d = Date::from_tick(0);
        assert_eq!(d.year_label(), "10000 BCE");
        assert_eq!(d.month, 1);
    }

    #[test]
    fn months_roll_over_into_years() {
        let d = Date::from_tick(13);
        assert_eq!(d.year_label(), "9999 BCE");
        assert_eq!(d.month, 2);
    }

    #[test]
    fn there_is_no_year_zero() {
        // 1 BCE is followed directly by 1 CE.
        let last_bce = Date { year: 0, month: 12 };
        let first_ce = Date { year: 1, month: 1 };
        assert_eq!(last_bce.year_label(), "1 BCE");
        assert_eq!(first_ce.year_label(), "1 CE");
        // 10,000 BCE to 1 CE is 10,000 years.
        let ticks_to_1ce = 10_000 * MONTHS_PER_YEAR;
        assert_eq!(Date::from_tick(ticks_to_1ce), first_ce);
    }
}
