//! Geodesic grid: hexagonal cells covering the sphere.
//!
//! Built by subdividing an icosahedron `level` times (each triangle split
//! into four, new points pushed out onto the sphere). The cells are the
//! regions around each resulting point: hexagons, except for the 12 original
//! icosahedron corners, which are pentagons.
//!
//! Point (cell) indices are stable across levels: the points of level k are
//! the first points of level k+1, in the same order. A level-k cell index is
//! therefore also a valid level-(k+n) cell index for the same location.
//!
//! Cell count = 10 * 4^level + 2.
//!
//! Cells are not all the same size: at level 8 they range from about 690 to
//! 940 km² (mean 778 km²), larger near the middle of each icosahedron face.
//! Every calculation uses each cell's own area, so this affects only how
//! finely different places are resolved, not any totals.

use std::collections::HashMap;

/// Mean radius of the Earth in metres (IUGG mean radius R1 = 6,371,008.8 m;
/// Moritz 2000, "Geodetic Reference System 1980").
pub const EARTH_RADIUS_M: f64 = 6_371_008.8;

pub type Vec3 = [f64; 3];

pub struct Grid {
    pub level: u32,
    /// Unit vectors to cell centres.
    pub centers: Vec<Vec3>,
    /// Neighbours of cell i are `nbrs[offsets[i]..offsets[i + 1]]`, in
    /// counter-clockwise order seen from outside the sphere.
    pub offsets: Vec<u32>,
    pub nbrs: Vec<u32>,
    /// Corner j of cell i is `corners[corner_ids[offsets[i] + j]]`; corner j
    /// lies between neighbour j and neighbour j + 1.
    pub corner_ids: Vec<u32>,
    /// Unit vectors to cell corners (one per icosphere triangle).
    pub corners: Vec<Vec3>,
    /// Cell areas in square metres.
    pub areas_m2: Vec<f64>,
}

pub fn cell_count(level: u32) -> usize {
    10 * 4usize.pow(level) + 2
}

fn normalize(v: Vec3) -> Vec3 {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Area of the spherical triangle abc on the unit sphere (steradians).
/// Van Oosterom & Strackee (1983).
fn spherical_triangle_area(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let num = dot(a, cross(b, c)).abs();
    let den = 1.0 + dot(a, b) + dot(b, c) + dot(c, a);
    2.0 * num.atan2(den)
}

/// Latitude and longitude in degrees of a unit vector (z = north pole,
/// x = 0° longitude).
pub fn to_lat_lon(v: Vec3) -> (f64, f64) {
    let lat = v[2].clamp(-1.0, 1.0).asin().to_degrees();
    let lon = v[1].atan2(v[0]).to_degrees();
    (lat, lon)
}

pub fn from_lat_lon(lat_deg: f64, lon_deg: f64) -> Vec3 {
    let (lat, lon) = (lat_deg.to_radians(), lon_deg.to_radians());
    [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()]
}

/// Cell corners: the circumcentre of each triangle, i.e. the point equally
/// far from its three cell centres. This makes cells exact Voronoi regions,
/// so "the cell containing a point" is "the nearest centre".
fn circumcenters(points: &[Vec3], faces: &[[u32; 3]]) -> Vec<Vec3> {
    faces
        .iter()
        .map(|f| {
            let (a, b, c) = (
                points[f[0] as usize],
                points[f[1] as usize],
                points[f[2] as usize],
            );
            let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            normalize(cross(ab, ac))
        })
        .collect()
}

fn icosahedron() -> (Vec<Vec3>, Vec<[u32; 3]>) {
    // Vertices of a regular icosahedron, rotated so two corners sit at the
    // poles (pentagon cells then fall at the poles and in open ocean at
    // ±26.57° latitude).
    let mut v = vec![[0.0, 0.0, 1.0]];
    let lat = (0.5_f64).atan(); // 26.565°
    for i in 0..5 {
        let lon = (i as f64) * 72.0_f64.to_radians();
        v.push([lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()]);
    }
    for i in 0..5 {
        let lon = (i as f64 * 72.0 + 36.0).to_radians();
        v.push([lat.cos() * lon.cos(), lat.cos() * lon.sin(), -lat.sin()]);
    }
    v.push([0.0, 0.0, -1.0]);

    let mut f = Vec::new();
    for i in 0..5u32 {
        let a = 1 + i;
        let b = 1 + (i + 1) % 5;
        let c = 6 + i;
        let d = 6 + (i + 1) % 5;
        f.push([0, a, b]); // top cap
        f.push([a, c, b]); // upper band
        f.push([b, c, d]); // lower band
        f.push([c, 11, d]); // bottom cap
    }
    (v, f)
}

impl Grid {
    pub fn new(level: u32) -> Grid {
        let (mut points, mut faces) = icosahedron();
        for _ in 0..level {
            let mut midpoint: HashMap<(u32, u32), u32> =
                HashMap::with_capacity(faces.len() * 3 / 2);
            let mut next = Vec::with_capacity(faces.len() * 4);
            let mut mid = |a: u32, b: u32, points: &mut Vec<Vec3>| -> u32 {
                let key = if a < b { (a, b) } else { (b, a) };
                *midpoint.entry(key).or_insert_with(|| {
                    let (p, q) = (points[a as usize], points[b as usize]);
                    points.push(normalize([p[0] + q[0], p[1] + q[1], p[2] + q[2]]));
                    (points.len() - 1) as u32
                })
            };
            for &[a, b, c] in &faces {
                let ab = mid(a, b, &mut points);
                let bc = mid(b, c, &mut points);
                let ca = mid(c, a, &mut points);
                next.push([a, ab, ca]);
                next.push([b, bc, ab]);
                next.push([c, ca, bc]);
                next.push([ab, bc, ca]);
            }
            faces = next;
        }
        let n = points.len();
        debug_assert_eq!(n, cell_count(level));

        // Make every face counter-clockwise seen from outside.
        for f in faces.iter_mut() {
            let (a, b, c) = (
                points[f[0] as usize],
                points[f[1] as usize],
                points[f[2] as usize],
            );
            let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            if dot(cross(ab, ac), a) < 0.0 {
                f.swap(1, 2);
            }
        }

        let corners = circumcenters(&points, &faces);

        // For each cell: map "neighbour -> (face, next neighbour)" going
        // counter-clockwise, then walk the ring.
        let mut ring: Vec<Vec<(u32, u32, u32)>> = vec![Vec::with_capacity(6); n];
        for (fi, f) in faces.iter().enumerate() {
            for k in 0..3 {
                let (p, q, r) = (f[k], f[(k + 1) % 3], f[(k + 2) % 3]);
                ring[p as usize].push((q, fi as u32, r));
            }
        }
        let mut offsets = Vec::with_capacity(n + 1);
        let mut nbrs = Vec::with_capacity(n * 6);
        let mut corner_ids = Vec::with_capacity(n * 6);
        offsets.push(0u32);
        for r in &ring {
            let mut cur = r[0];
            for _ in 0..r.len() {
                nbrs.push(cur.0);
                corner_ids.push(cur.1);
                let next_nbr = cur.2;
                cur = *r.iter().find(|e| e.0 == next_nbr).expect("closed ring");
            }
            offsets.push(nbrs.len() as u32);
        }

        let r2 = EARTH_RADIUS_M * EARTH_RADIUS_M;
        let areas_m2 = (0..n)
            .map(|i| {
                let (s, e) = (offsets[i] as usize, offsets[i + 1] as usize);
                let c = points[i];
                (s..e)
                    .map(|j| {
                        let j2 = if j + 1 == e { s } else { j + 1 };
                        spherical_triangle_area(
                            c,
                            corners[corner_ids[j] as usize],
                            corners[corner_ids[j2] as usize],
                        )
                    })
                    .sum::<f64>()
                    * r2
            })
            .collect();

        Grid {
            level,
            centers: points,
            offsets,
            nbrs,
            corner_ids,
            corners,
            areas_m2,
        }
    }

    pub fn len(&self) -> usize {
        self.centers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.centers.is_empty()
    }

    pub fn neighbors(&self, i: usize) -> &[u32] {
        &self.nbrs[self.offsets[i] as usize..self.offsets[i + 1] as usize]
    }

    pub fn cell_corners(&self, i: usize) -> impl Iterator<Item = Vec3> + '_ {
        self.corner_ids[self.offsets[i] as usize..self.offsets[i + 1] as usize]
            .iter()
            .map(|&c| self.corners[c as usize])
    }

    /// Index of the cell whose centre is nearest to `p`, found by walking
    /// from `start` towards `p` (cells are Voronoi regions, so the walk
    /// always reaches the right one).
    pub fn nearest_from(&self, p: Vec3, start: usize) -> usize {
        let mut best = start;
        let mut best_d = dot(self.centers[best], p);
        loop {
            let mut improved = false;
            for &nb in self.neighbors(best) {
                let d = dot(self.centers[nb as usize], p);
                if d > best_d {
                    best_d = d;
                    best = nb as usize;
                    improved = true;
                }
            }
            if !improved {
                return best;
            }
        }
    }
}

/// For each cell of `fine`, the nearest cell of `coarse` (a coarser grid
/// built the same way).
pub fn parent_map(fine: &Grid, coarse: &Grid) -> Vec<u32> {
    let mut out = Vec::with_capacity(fine.len());
    let mut guess = 0;
    for &c in &fine.centers {
        guess = coarse.nearest_from(c, guess);
        out.push(guess as u32);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cell_counts_match_formula() {
        for level in 0..5 {
            assert_eq!(Grid::new(level).len(), cell_count(level));
        }
    }

    #[test]
    fn twelve_pentagons_rest_hexagons() {
        let g = Grid::new(4);
        let pent = (0..g.len()).filter(|&i| g.neighbors(i).len() == 5).count();
        let hex = (0..g.len()).filter(|&i| g.neighbors(i).len() == 6).count();
        assert_eq!(pent, 12);
        assert_eq!(pent + hex, g.len());
    }

    #[test]
    fn areas_sum_to_earth_surface() {
        let g = Grid::new(5);
        let total: f64 = g.areas_m2.iter().sum();
        let sphere = 4.0 * std::f64::consts::PI * EARTH_RADIUS_M * EARTH_RADIUS_M;
        assert!((total / sphere - 1.0).abs() < 1e-9, "{total} vs {sphere}");
    }

    #[test]
    fn neighbours_are_mutual() {
        let g = Grid::new(4);
        for i in 0..g.len() {
            for &n in g.neighbors(i) {
                assert!(g.neighbors(n as usize).contains(&(i as u32)));
            }
        }
    }

    #[test]
    fn coarse_centres_lie_in_their_own_fine_cells() {
        // Indices carry over between levels: coarse cell i is near fine cell i.
        let a = Grid::new(3);
        let b = Grid::new(5);
        let parents = parent_map(&b, &a);
        for (i, &p) in parents.iter().take(a.len()).enumerate() {
            assert_eq!(p, i as u32);
        }
    }

    #[test]
    fn nearest_walk_finds_true_nearest() {
        let g = Grid::new(4);
        for k in 0..500 {
            let lat = -89.0 + (k as f64 * 37.3) % 178.0;
            let lon = -180.0 + (k as f64 * 91.7) % 360.0;
            let p = from_lat_lon(lat, lon);
            let brute = (0..g.len())
                .max_by(|&a, &b| dot(g.centers[a], p).total_cmp(&dot(g.centers[b], p)))
                .unwrap();
            assert_eq!(g.nearest_from(p, 0), brute);
        }
    }
}
