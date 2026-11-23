//! The kinds of run the simulator can start.

use serde::Serialize;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Scenario {
    /// The normal run: Earth on 1 January 10,000 BCE, with the orbit and
    /// CO2 changing as the record says.
    Holocene,
    /// The present-day check: today's coastlines, ice, orbit and CO2 held
    /// fixed from 1 January 1931; 1931-1980 lets the model settle, and
    /// 1981-2010 is averaged and compared with observations.
    PresentDayCheck,
}

impl Scenario {
    /// Astronomical year in which the run starts.
    pub fn start_year(self) -> i64 {
        match self {
            Scenario::Holocene => crate::time::START_YEAR,
            Scenario::PresentDayCheck => 1931,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Scenario::Holocene => "From 10,000 BCE",
            Scenario::PresentDayCheck => "Present-day check",
        }
    }

    /// Years averaged for comparison with observations, if any.
    pub fn comparison_window(self) -> Option<(i64, i64)> {
        match self {
            Scenario::Holocene => None,
            Scenario::PresentDayCheck => Some((1981, 2010)),
        }
    }
}
