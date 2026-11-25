// Numerical code here loops over cells and indexes several parallel arrays
// by the same cell index; iterator chains would be less clear.
#![allow(clippy::needless_range_loop)]

//! The climate: sunlight and temperature (Step 2), later rain, snow,
//! ocean, variability and ice. Runs once per simulated month on the 40,962
//! cells of the climate grid; see `docs/design.md` and the Phase 1 plan.

pub mod atmosphere;
pub mod check;
pub mod koppen;
pub mod model;
pub mod ocean;
pub mod orbit;
pub mod solver;
pub mod variability;

pub use atmosphere::{AtmosphereTuning, ATMOSPHERE_TUNED};
pub use model::{Climate, Tuning, TUNED};
pub use orbit::Orbit;
