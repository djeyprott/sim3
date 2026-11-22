//! Terrain on the map grid: which cells are sea, land or ice at a given
//! moment, derived from the world data.

use crate::data::Epoch;
use crate::grid::{Grid, Vec3, EARTH_RADIUS_M};

/// What a cell's surface is, for maps and (later) for every system that
/// cares whether a place is sea or land.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Surface {
    /// Mostly below sea level and connected to the world ocean.
    Ocean = 0,
    /// Mostly above sea level.
    Land = 1,
    /// Mostly below sea level but cut off from the ocean: an inland basin
    /// such as the Caspian or the Dead Sea. Whether it holds a lake or is dry
    /// depends on rainfall and evaporation, which arrive with the climate.
    InlandBasin = 2,
}

pub struct Terrain {
    pub surface: Vec<Surface>,
    /// Cells whose lowest point is below sea level and that connect to the
    /// world ocean through other such cells. Ships can pass these.
    pub sea_connected: Vec<bool>,
}

impl Terrain {
    pub fn new(grid: &Grid, epoch: &Epoch) -> Terrain {
        let n = grid.len();
        assert_eq!(n, epoch.elev_mean_m.len(), "grid and data must match");

        // Flood the ocean outward from the deepest cell, through every cell
        // with at least some ground below sea level. Narrow straits stay
        // open because a cell's lowest point decides, not its average.
        let deepest = (0..n).min_by_key(|&i| epoch.elev_min_m[i]).unwrap();
        let mut sea_connected = vec![false; n];
        let mut stack = vec![deepest];
        sea_connected[deepest] = true;
        while let Some(i) = stack.pop() {
            for &nb in grid.neighbors(i) {
                let nb = nb as usize;
                if !sea_connected[nb] && epoch.elev_min_m[nb] < 0 {
                    sea_connected[nb] = true;
                    stack.push(nb);
                }
            }
        }

        let surface = (0..n)
            .map(|i| {
                let mostly_land = epoch.land_share[i] >= 128;
                if mostly_land {
                    Surface::Land
                } else if sea_connected[i] {
                    Surface::Ocean
                } else {
                    Surface::InlandBasin
                }
            })
            .collect();
        Terrain {
            surface,
            sea_connected,
        }
    }

    /// Share of the Earth's surface that is land (above sea level), weighting
    /// each cell by its area and its own land share.
    pub fn land_share_of_earth(grid: &Grid, epoch: &Epoch) -> f64 {
        let total: f64 = grid.areas_m2.iter().sum();
        let land: f64 = grid
            .areas_m2
            .iter()
            .zip(&epoch.land_share)
            .map(|(a, &s)| a * s as f64 / 255.0)
            .sum();
        land / total
    }
}

/// Unit vectors pointing east and north at `c`.
fn east_north(c: Vec3) -> (Vec3, Vec3) {
    let z = [0.0, 0.0, 1.0];
    let e = [
        z[1] * c[2] - z[2] * c[1],
        z[2] * c[0] - z[0] * c[2],
        z[0] * c[1] - z[1] * c[0],
    ];
    let len = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
    let e = if len < 1e-12 {
        [1.0, 0.0, 0.0]
    } else {
        [e[0] / len, e[1] / len, e[2] / len]
    };
    let nth = [
        c[1] * e[2] - c[2] * e[1],
        c[2] * e[0] - c[0] * e[2],
        c[0] * e[1] - c[1] * e[0],
    ];
    (e, nth)
}

/// Terrain slope of each cell as (eastward, northward) gradient in metres
/// per metre, from a least-squares fit to its neighbours' elevations.
pub fn slopes(grid: &Grid, elevation_m: &[i16]) -> Vec<(f64, f64)> {
    (0..grid.len())
        .map(|i| {
            let c = grid.centers[i];
            let (e, nth) = east_north(c);
            let (mut sxx, mut sxy, mut syy, mut sxz, mut syz) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for &nb in grid.neighbors(i) {
                let p = grid.centers[nb as usize];
                let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
                let x = (d[0] * e[0] + d[1] * e[1] + d[2] * e[2]) * EARTH_RADIUS_M;
                let y = (d[0] * nth[0] + d[1] * nth[1] + d[2] * nth[2]) * EARTH_RADIUS_M;
                let dz = elevation_m[nb as usize] as f64 - elevation_m[i] as f64;
                sxx += x * x;
                sxy += x * y;
                syy += y * y;
                sxz += x * dz;
                syz += y * dz;
            }
            let det = sxx * syy - sxy * sxy;
            if det.abs() < 1e-9 {
                return (0.0, 0.0);
            }
            ((syy * sxz - sxy * syz) / det, (sxx * syz - sxy * sxz) / det)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slope_of_a_tilted_plane() {
        // Elevation rising 1 m per km northward everywhere near the equator.
        let g = Grid::new(5);
        let elev: Vec<i16> = g
            .centers
            .iter()
            .map(|c| (c[2].asin() * EARTH_RADIUS_M / 1000.0).round() as i16)
            .collect();
        let s = slopes(&g, &elev);
        let i = g.nearest_from([1.0, 0.0, 0.0], 0);
        assert!(s[i].0.abs() < 2e-4, "east slope {}", s[i].0);
        assert!((s[i].1 - 1e-3).abs() < 2e-4, "north slope {}", s[i].1);
    }
}
