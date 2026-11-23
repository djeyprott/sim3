//! Everything about the world that doesn't change during a run and is
//! shared between runs: the data file and the grids built from it.

use world::data::WorldData;
use world::grid::{parent_map, Grid};

pub struct WorldContext {
    pub data: WorldData,
    /// 655,362 cells, ~28 km.
    pub map_grid: Grid,
    /// 40,962 cells, ~110 km.
    pub climate_grid: Grid,
    /// Climate cell containing each map cell.
    pub map_to_climate: Vec<u32>,
}

impl WorldContext {
    pub fn new(data: WorldData) -> WorldContext {
        let map_grid = Grid::new(data.grid_level);
        let climate_grid = Grid::new(world::CLIMATE_LEVEL);
        let map_to_climate = parent_map(&map_grid, &climate_grid);
        WorldContext {
            data,
            map_grid,
            climate_grid,
            map_to_climate,
        }
    }

    /// Index of the terrain snapshot for a year: the 10,000 BCE snapshot
    /// for ancient years, today's for recent ones.
    pub fn epoch_for_year(&self, year: i64) -> usize {
        let (mut best, mut dist) = (0, i64::MAX);
        for (i, e) in self.data.epochs.iter().enumerate() {
            let d = (e.year as i64 - year).abs();
            if d < dist {
                best = i;
                dist = d;
            }
        }
        best
    }
}
