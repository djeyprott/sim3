//! Map data for the viewer: everything the browser needs to draw the world,
//! prepared once at startup and served as compact binary arrays.

use axum::body::Bytes;
use serde::Serialize;

use world::grid::{from_lat_lon, Grid};
use world::terrain::{slopes, Terrain};

/// The world data file, built by .github/workflows/data.yml (see build.rs).
pub const WORLD_DATA: &[u8] =
    include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/world.bin"));

/// Size of the lookup image that tells the browser which cell each map pixel
/// belongs to (an equirectangular projection, 0.088° per pixel: about 10 km
/// at the equator, finer than the 28 km cells).
pub const RASTER_W: usize = 4096;
pub const RASTER_H: usize = 2048;

#[derive(Serialize)]
pub struct EpochInfo {
    pub index: usize,
    /// Astronomical year.
    pub year: i32,
    pub label: String,
    /// Share of the Earth's surface above sea level.
    pub land_share_of_earth: f64,
}

#[derive(Serialize)]
pub struct MapInfo {
    pub cells: usize,
    pub raster_width: usize,
    pub raster_height: usize,
    pub epochs: Vec<EpochInfo>,
}

pub struct MapData {
    pub info: MapInfo,
    /// Per epoch, the cell arrays in the layout described at [`cell_blob`].
    pub cell_blobs: Vec<Bytes>,
    /// u32 width, u32 height, then one u32 cell index per pixel, row by row
    /// from the north.
    pub raster: Bytes,
}

impl MapData {
    pub fn build(world: &sim::WorldContext) -> MapData {
        let data = &world.data;
        let grid = &world.map_grid;

        let mut epochs = Vec::new();
        let mut cell_blobs = Vec::new();
        for (index, epoch) in data.epochs.iter().enumerate() {
            let terrain = Terrain::new(grid, epoch);
            epochs.push(EpochInfo {
                index,
                year: epoch.year,
                label: if epoch.year < 0 {
                    format!("1 January {} BCE", 1 - epoch.year)
                } else {
                    "Today".to_string()
                },
                land_share_of_earth: Terrain::land_share_of_earth(grid, epoch),
            });
            cell_blobs.push(Bytes::from(cell_blob(grid, epoch, &terrain)));
        }

        MapData {
            info: MapInfo {
                cells: grid.len(),
                raster_width: RASTER_W,
                raster_height: RASTER_H,
                epochs,
            },
            cell_blobs,
            raster: Bytes::from(raster(grid)),
        }
    }
}

/// Per-cell arrays for one epoch, little-endian, one after another:
/// i16 mean elevation (m), i16 lowest point (m), i16 highest point (m),
/// u16 ice thickness (m), i16 eastward slope and i16 northward slope
/// (metres per km, ×10), u8 land share (0-255), u8 surface kind
/// (0 ocean, 1 land, 2 inland basin).
fn cell_blob(grid: &Grid, epoch: &world::data::Epoch, terrain: &Terrain) -> Vec<u8> {
    let n = grid.len();
    let slope = slopes(grid, &epoch.elev_mean_m);
    let mut out = Vec::with_capacity(n * 14);
    for v in [&epoch.elev_mean_m, &epoch.elev_min_m, &epoch.elev_max_m] {
        v.iter()
            .for_each(|x| out.extend_from_slice(&x.to_le_bytes()));
    }
    epoch
        .ice_thickness_m
        .iter()
        .for_each(|x| out.extend_from_slice(&x.to_le_bytes()));
    let scaled = |x: f64| (x * 1e4).round().clamp(-32768.0, 32767.0) as i16;
    slope
        .iter()
        .for_each(|s| out.extend_from_slice(&scaled(s.0).to_le_bytes()));
    slope
        .iter()
        .for_each(|s| out.extend_from_slice(&scaled(s.1).to_le_bytes()));
    out.extend_from_slice(&epoch.land_share);
    out.extend(terrain.surface.iter().map(|&s| s as u8));
    out
}

fn raster(grid: &Grid) -> Vec<u8> {
    let mut out = Vec::with_capacity(8 + RASTER_W * RASTER_H * 4);
    out.extend_from_slice(&(RASTER_W as u32).to_le_bytes());
    out.extend_from_slice(&(RASTER_H as u32).to_le_bytes());
    let mut guess = 0;
    for y in 0..RASTER_H {
        let lat = 90.0 - (y as f64 + 0.5) * 180.0 / RASTER_H as f64;
        for x in 0..RASTER_W {
            let lon = -180.0 + (x as f64 + 0.5) * 360.0 / RASTER_W as f64;
            guess = grid.nearest_from(from_lat_lon(lat, lon), guess);
            out.extend_from_slice(&(guess as u32).to_le_bytes());
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use world::data::WorldData;
    use world::terrain::Surface;

    fn setup() -> (Grid, WorldData) {
        let data = WorldData::parse(WORLD_DATA).unwrap();
        (Grid::new(data.grid_level), data)
    }

    #[test]
    fn todays_land_is_about_29_percent_of_the_earth() {
        let (grid, data) = setup();
        let today = data.epochs.iter().find(|e| e.year > 0).unwrap();
        let share = Terrain::land_share_of_earth(&grid, today);
        assert!((0.285..0.295).contains(&share), "land share {share}");
    }

    /// Whether two places are joined by land (cells that are mostly land).
    fn connected_by_land(grid: &Grid, terrain: &Terrain, a: (f64, f64), b: (f64, f64)) -> bool {
        let start = grid.nearest_from(world::grid::from_lat_lon(a.0, a.1), 0);
        let goal = grid.nearest_from(world::grid::from_lat_lon(b.0, b.1), 0);
        let mut seen = vec![false; grid.len()];
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(i) = stack.pop() {
            if i == goal {
                return true;
            }
            for &nb in grid.neighbors(i) {
                let nb = nb as usize;
                if !seen[nb] && terrain.surface[nb] == Surface::Land {
                    seen[nb] = true;
                    stack.push(nb);
                }
            }
        }
        false
    }

    #[test]
    fn land_bridges_join_continents_in_10000_bce() {
        let (grid, data) = setup();
        let (then, now) = (&data.epochs[0], &data.epochs[1]);
        assert_eq!(then.year, -9999);
        let (t_then, t_now) = (Terrain::new(&grid, then), Terrain::new(&grid, now));
        for (name, a, b) in [
            ("Britain and France (Doggerland)", (52.5, -1.5), (47.0, 2.5)),
            (
                "Siberia and Alaska (Beringia)",
                (66.0, 170.0),
                (64.0, -155.0),
            ),
            ("Sumatra and Borneo (Sundaland)", (0.0, 102.0), (0.5, 114.0)),
            (
                "New Guinea and Australia (Sahul)",
                (-5.0, 141.0),
                (-20.0, 135.0),
            ),
        ] {
            assert!(
                connected_by_land(&grid, &t_then, a, b),
                "{name}: should be joined in 10,000 BCE"
            );
            assert!(
                !connected_by_land(&grid, &t_now, a, b),
                "{name}: should be separate today"
            );
        }
    }

    #[test]
    fn ice_covers_hudson_bay_in_10000_bce() {
        let (grid, data) = setup();
        let i = grid.nearest_from(world::grid::from_lat_lon(60.0, -85.0), 0);
        assert!(data.epochs[0].ice_thickness_m[i] > 1000);
        assert_eq!(data.epochs[1].ice_thickness_m[i], 0);
    }

    #[test]
    fn caspian_is_an_inland_basin_and_the_black_sea_connects_today() {
        let (grid, data) = setup();
        let today = &data.epochs[1];
        let terrain = Terrain::new(&grid, today);
        let caspian = grid.nearest_from(world::grid::from_lat_lon(42.0, 51.0), 0);
        let black_sea = grid.nearest_from(world::grid::from_lat_lon(43.0, 34.0), 0);
        assert_eq!(terrain.surface[caspian], Surface::InlandBasin);
        assert_eq!(terrain.surface[black_sea], Surface::Ocean);
    }
}
