//! Simulation time.
//!
//! The simulation advances one day per tick. Dates are labelled with the
//! proleptic Gregorian calendar (today's calendar extended backwards), which
//! keeps the seasons in place over the whole run: the mean Gregorian year
//! (365.2425 days) differs from the tropical year by about one day in
//! 3,000 years.
//!
//! Years use astronomical numbering internally (year 0 exists, 1 BCE = year 0,
//! 10,000 BCE = year -9999) so arithmetic stays simple; conversion to BCE/CE
//! happens only for display, where there is no year zero.

use history::Tick;
use serde::Serialize;

/// Astronomical year in which the simulation starts: 10,000 BCE.
pub const START_YEAR: i64 = -9999;

/// Mean length of the Gregorian year in days.
pub const DAYS_PER_YEAR: f64 = 365.2425;

const MONTH_NAMES: [&str; 12] = [
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
pub struct Date {
    /// Astronomical year (can be zero or negative).
    pub year: i64,
    /// 1..=12
    pub month: u8,
    /// 1..=31
    pub day: u8,
}

pub fn is_leap_year(year: i64) -> bool {
    year.rem_euclid(4) == 0 && (year.rem_euclid(100) != 0 || year.rem_euclid(400) == 0)
}

pub fn days_in_month(year: i64, month: u8) -> u8 {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if is_leap_year(year) => 29,
        2 => 28,
        _ => panic!("month out of range: {month}"),
    }
}

/// Days since 1 January 1970 (proleptic Gregorian). Algorithm from
/// Howard Hinnant, "chrono-Compatible Low-Level Date Algorithms".
fn days_from_civil(year: i64, month: u8, day: u8) -> i64 {
    let y = if month <= 2 { year - 1 } else { year };
    let era = y.div_euclid(400);
    let yoe = y.rem_euclid(400);
    let m = month as i64;
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + day as i64 - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// Inverse of [`days_from_civil`].
fn civil_from_days(days: i64) -> Date {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = (doy - (153 * mp + 2) / 5 + 1) as u8;
    let month = if mp < 10 { mp + 3 } else { mp - 9 } as u8;
    let year = yoe + era * 400 + if month <= 2 { 1 } else { 0 };
    Date { year, month, day }
}

fn start_days() -> i64 {
    days_from_civil(START_YEAR, 1, 1)
}

impl Date {
    pub fn from_tick(tick: Tick) -> Self {
        civil_from_days(start_days() + tick as i64)
    }

    /// The tick on which this date falls. Dates before the start clamp to 0.
    pub fn to_tick(self) -> Tick {
        (days_from_civil(self.year, self.month, self.day) - start_days()).max(0) as Tick
    }

    /// The same day `n` months later, clamped to the end of shorter months
    /// (31 January + 1 month = 28 or 29 February).
    pub fn add_months(self, n: u32) -> Date {
        let total = self.year * 12 + (self.month as i64 - 1) + n as i64;
        let year = total.div_euclid(12);
        let month = (total.rem_euclid(12) + 1) as u8;
        let day = self.day.min(days_in_month(year, month));
        Date { year, month, day }
    }

    pub fn add_years(self, n: u32) -> Date {
        self.add_months(n * 12)
    }

    /// Year as people write it: "10000 BCE", "1 CE".
    pub fn year_label(&self) -> String {
        if self.year <= 0 {
            format!("{} BCE", 1 - self.year)
        } else {
            format!("{} CE", self.year)
        }
    }

    /// Full date: "14 March 9873 BCE".
    pub fn label(&self) -> String {
        format!(
            "{} {} {}",
            self.day,
            MONTH_NAMES[self.month as usize - 1],
            self.year_label()
        )
    }

    /// Whole calendar years since the simulation started.
    pub fn years_since_start(&self) -> i64 {
        self.year - START_YEAR
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn starts_on_1_january_10000_bce() {
        assert_eq!(Date::from_tick(0).label(), "1 January 10000 BCE");
    }

    #[test]
    fn days_roll_over_into_months_and_years() {
        assert_eq!(Date::from_tick(31).label(), "1 February 10000 BCE");
        // 10000 BCE is astronomical year -9999: not a leap year.
        assert_eq!(Date::from_tick(365).label(), "1 January 9999 BCE");
    }

    #[test]
    fn leap_years_follow_gregorian_rules() {
        assert!(is_leap_year(2000));
        assert!(!is_leap_year(1900));
        assert!(is_leap_year(0)); // 1 BCE
        assert!(is_leap_year(-4)); // 5 BCE
        assert!(!is_leap_year(-100)); // 101 BCE
        assert!(is_leap_year(-9996)); // 9997 BCE
    }

    #[test]
    fn matches_known_dates() {
        assert_eq!(days_from_civil(1970, 1, 1), 0);
        assert_eq!(days_from_civil(2000, 3, 1), 11_017);
        assert_eq!(
            civil_from_days(11_017),
            Date {
                year: 2000,
                month: 3,
                day: 1
            }
        );
    }

    #[test]
    fn round_trips_across_the_whole_run() {
        for tick in (0..5_000_000u64).step_by(9_973) {
            assert_eq!(Date::from_tick(tick).to_tick(), tick);
        }
    }

    #[test]
    fn there_is_no_year_zero() {
        let last_bce = Date {
            year: 0,
            month: 12,
            day: 31,
        };
        let first_ce = Date {
            year: 1,
            month: 1,
            day: 1,
        };
        assert_eq!(last_bce.label(), "31 December 1 BCE");
        assert_eq!(first_ce.label(), "1 January 1 CE");
        assert_eq!(first_ce.to_tick() - last_bce.to_tick(), 1);
    }

    #[test]
    fn mean_year_length_is_gregorian() {
        let span = Date {
            year: 1,
            month: 1,
            day: 1,
        }
        .to_tick();
        let mean = span as f64 / 10_000.0;
        assert!((mean - DAYS_PER_YEAR).abs() < 0.001, "mean year {mean}");
    }

    #[test]
    fn adding_months_clamps_to_month_end() {
        let d = Date {
            year: -9999,
            month: 1,
            day: 31,
        };
        assert_eq!(
            d.add_months(1),
            Date {
                year: -9999,
                month: 2,
                day: 28
            }
        );
        assert_eq!(
            d.add_years(1),
            Date {
                year: -9998,
                month: 1,
                day: 31
            }
        );
    }
}
