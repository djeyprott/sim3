//! The Earth's orbit and the sunlight it delivers.
//!
//! Orbital parameters come from the Berger & Loutre (1991) table in the
//! world data (1,000-year steps), interpolated linearly: over 1,000 years
//! the parameters change smoothly (the perihelion moves about 17°).
//!
//! Daily-mean sunshine at the top of the atmosphere follows Berger (1978),
//! as implemented in climlab (Rose 2018, MIT licence): the calendar day is
//! converted to the Earth's position on its orbit (solar longitude, 0 at the
//! northern spring equinox, fixed at day 80 of the year), then the Sun's
//! declination, distance and day length give the daily average.

use std::f64::consts::PI;
use world::data::OrbitRow;

/// Total solar irradiance at 1 AU, W/m², satellite measurement (Kopp & Lean
/// 2011, Geophys. Res. Lett. 38, L01706).
pub const SOLAR_CONSTANT: f64 = 1361.0;

/// Length of the tropical year in days, which the climate calendar follows.
pub const DAYS_PER_YEAR: f64 = 365.2422;

/// Day of the year (1 = 1 January) on which the northern spring equinox
/// falls in the convention of Berger (1978).
const EQUINOX_DAY: f64 = 80.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orbit {
    pub eccentricity: f64,
    /// Longitude of perihelion measured from the northern spring equinox,
    /// radians (Berger 1978 convention, +180° relative to the table file).
    pub perihelion: f64,
    /// Axial tilt, radians.
    pub obliquity: f64,
}

impl Orbit {
    /// Today's orbit (Berger & Loutre 1991, year 1950).
    pub fn present() -> Orbit {
        Orbit {
            eccentricity: 0.017236,
            perihelion: 281.37_f64.to_radians(),
            obliquity: 23.446_f64.to_radians(),
        }
    }

    /// Orbit `years_bp` years before 1950, from the table.
    pub fn at(table: &[OrbitRow], years_bp: f64) -> Orbit {
        assert!(table.len() >= 2, "orbit table missing from world data");
        let kyr = -years_bp / 1000.0;
        let i = table
            .iter()
            .rposition(|r| r.kyr as f64 <= kyr)
            .unwrap_or(0)
            .min(table.len() - 2);
        let (a, b) = (table[i], table[i + 1]);
        let t = ((kyr - a.kyr as f64) / (b.kyr - a.kyr) as f64).clamp(0.0, 1.0);
        // Perihelion wraps at 360°: interpolate along the shorter way.
        let mut dp = b.perihelion_deg as f64 - a.perihelion_deg as f64;
        if dp > 180.0 {
            dp -= 360.0;
        } else if dp < -180.0 {
            dp += 360.0;
        }
        Orbit {
            eccentricity: a.eccentricity as f64 + t * (b.eccentricity - a.eccentricity) as f64,
            perihelion: (a.perihelion_deg as f64 + t * dp)
                .rem_euclid(360.0)
                .to_radians(),
            obliquity: (a.obliquity_deg as f64 + t * (b.obliquity_deg - a.obliquity_deg) as f64)
                .to_radians(),
        }
    }

    /// Solar longitude (radians, 0 at the northern spring equinox) on a
    /// day of the year (1-based, fractional). Berger (1978), section 3.
    pub fn solar_longitude(&self, day: f64) -> f64 {
        let e = self.eccentricity;
        let p = self.perihelion;
        let delta_lambda = (day - EQUINOX_DAY) * 2.0 * PI / DAYS_PER_YEAR;
        let beta = (1.0 - e * e).sqrt();
        let lambda_m = -2.0
            * ((e / 2.0 + e.powi(3) / 8.0) * (1.0 + beta) * (-p).sin()
                - e * e / 4.0 * (0.5 + beta) * (-2.0 * p).sin()
                + e.powi(3) / 8.0 * (1.0 / 3.0 + beta) * (-3.0 * p).sin())
            + delta_lambda;
        lambda_m
            + (2.0 * e - e.powi(3) / 4.0) * (lambda_m - p).sin()
            + 1.25 * e * e * (2.0 * (lambda_m - p)).sin()
            + 13.0 / 12.0 * e.powi(3) * (3.0 * (lambda_m - p)).sin()
    }

    /// Daily-mean sunshine at the top of the atmosphere (W/m²) at latitude
    /// `lat` (radians) on day `day` of the year.
    pub fn daily_insolation(&self, lat: f64, day: f64, solar_constant: f64) -> f64 {
        self.insolation_at_longitude(lat, self.solar_longitude(day), solar_constant)
    }

    /// Daily-mean sunshine (W/m²) at latitude `lat` when the Earth is at
    /// solar longitude `lambda` (radians).
    pub fn insolation_at_longitude(&self, lat: f64, lambda: f64, solar_constant: f64) -> f64 {
        let delta = (self.obliquity.sin() * lambda.sin()).asin();
        let e = self.eccentricity;
        let irradiance =
            (1.0 + e * (lambda - self.perihelion).cos()).powi(2) / (1.0 - e * e).powi(2);
        // Hour angle of sunset; 0 in polar night, π in polar day.
        let cos_h0 = (-lat.tan() * delta.tan()).clamp(-1.0, 1.0);
        let h0 = cos_h0.acos();
        let coszen = (h0 * lat.sin() * delta.sin() + lat.cos() * delta.cos() * h0.sin()) / PI;
        solar_constant * irradiance * coszen.max(0.0)
    }
}

/// Number of latitude rows in the sunshine table (0.5° each).
pub const LAT_ROWS: usize = 360;

/// Mean daily sunshine (W/m²) for each calendar month and each 0.5°
/// latitude row, for one orbit. Rows run from the south pole northwards.
#[derive(Clone)]
pub struct MonthlySunshine {
    pub table: Vec<[f32; 12]>,
}

impl MonthlySunshine {
    /// `month_days[m]` = (first day of year, number of days) for month m.
    pub fn new(
        orbit: &Orbit,
        solar_constant: f64,
        month_days: &[(f64, f64); 12],
    ) -> MonthlySunshine {
        let table = (0..LAT_ROWS)
            .map(|r| {
                let lat = (-90.0 + (r as f64 + 0.5) * 180.0 / LAT_ROWS as f64).to_radians();
                let mut row = [0.0f32; 12];
                for (m, &(first, len)) in month_days.iter().enumerate() {
                    let n = len.round() as usize;
                    let sum: f64 = (0..n)
                        .map(|d| {
                            orbit.daily_insolation(lat, first + d as f64 + 0.5, solar_constant)
                        })
                        .sum();
                    row[m] = (sum / n as f64) as f32;
                }
                row
            })
            .collect();
        MonthlySunshine { table }
    }

    pub fn at(&self, lat_deg: f64, month: usize) -> f32 {
        self.table[Self::row(lat_deg)][month]
    }

    /// Table row for a latitude in degrees.
    pub fn row(lat_deg: f64) -> usize {
        (((lat_deg + 90.0) / 180.0 * LAT_ROWS as f64) as usize).min(LAT_ROWS - 1)
    }
}

/// First day of the year and length of each month in a 365-day year.
pub fn month_days() -> [(f64, f64); 12] {
    let lens = [
        31.0, 28.25, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0,
    ];
    let mut out = [(0.0, 0.0); 12];
    let mut first = 1.0;
    for (m, &l) in lens.iter().enumerate() {
        out[m] = (first, l);
        first += l;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // The Berger & Loutre (1991) table lists "July" sunshine at 65°N,
    // computed with a solar constant of 1360 W/m² at solar longitude 120°:
    // 426.76 W/m² for 1950 and 467.00 W/m² for 8,000 years earlier.
    #[test]
    fn matches_the_published_65n_july_values() {
        let now = Orbit::present();
        let past = Orbit {
            eccentricity: 0.019674,
            perihelion: (326.82f64 + 180.0).rem_euclid(360.0).to_radians(),
            obliquity: 24.206f64.to_radians(),
        };
        let lat = 65f64.to_radians();
        let july = 120f64.to_radians();
        let (a, b) = (
            now.insolation_at_longitude(lat, july, 1360.0),
            past.insolation_at_longitude(lat, july, 1360.0),
        );
        assert!((a - 426.76).abs() < 0.1, "today {a}");
        assert!((b - 467.00).abs() < 0.1, "8 kyr ago {b}");
    }

    #[test]
    fn exact_geometric_cases() {
        let o = Orbit {
            eccentricity: 0.0,
            perihelion: 0.0,
            obliquity: 23.44f64.to_radians(),
        };
        // Equator at the equinox (day 80): S0 / π.
        let q = o.daily_insolation(0.0, 80.0, SOLAR_CONSTANT);
        assert!((q - SOLAR_CONSTANT / PI).abs() < 1.0, "{q}");
        // North pole at the June solstice (a quarter year later): S0 sin(tilt).
        let q = o.daily_insolation(PI / 2.0 - 1e-9, 80.0 + DAYS_PER_YEAR / 4.0, SOLAR_CONSTANT);
        assert!((q - SOLAR_CONSTANT * o.obliquity.sin()).abs() < 1.0, "{q}");
        // Polar night: 80°N in mid-January gets no sun.
        assert_eq!(
            o.daily_insolation(80f64.to_radians(), 15.0, SOLAR_CONSTANT),
            0.0
        );
    }

    #[test]
    fn annual_global_mean_is_a_quarter_of_the_solar_constant() {
        let o = Orbit::present();
        let mut sum = 0.0;
        let mut w = 0.0;
        for r in 0..180 {
            let lat = (-90.0 + r as f64 + 0.5).to_radians();
            for d in 0..365 {
                sum += o.daily_insolation(lat, d as f64 + 0.5, SOLAR_CONSTANT) * lat.cos();
                w += lat.cos();
            }
        }
        let mean = sum / w;
        // Slightly above S0/4 because of eccentricity: S0 / (4 sqrt(1-e^2)).
        let expected = SOLAR_CONSTANT / 4.0 / (1.0 - 0.017236f64.powi(2)).sqrt();
        assert!((mean - expected).abs() < 0.3, "{mean} vs {expected}");
    }
}
