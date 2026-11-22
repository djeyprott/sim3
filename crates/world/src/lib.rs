//! The world: the geodesic grids and the terrain, sea and ice on them.

pub mod data;
pub mod grid;
pub mod terrain;

/// Subdivision level of the map grid: 655,362 cells about 28 km across.
pub const MAP_LEVEL: u32 = 8;
/// Subdivision level of the climate grid: 40,962 cells about 110 km across.
pub const CLIMATE_LEVEL: u32 = 6;
