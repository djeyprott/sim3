//! The coarse (~440 km) grid on which the air and the open sea's mixed
//! layer are solved: its neighbour geometry, and smoothing weights for
//! handing coarse fields down to the climate cells.
//!
//! Heat carried by the atmosphere and ocean smooths temperature over a
//! length of about R sqrt(D/B) ≈ 3,400 km (R the Earth's radius, D ≈ 0.6
//! W/m²/K the diffusion coefficient on the unit sphere, B ≈ 2.1 W/m²/K the
//! rise in outgoing heat per degree; North et al. 1981), so a ~440 km grid
//! resolves it.

use world::grid::{parent_map, Grid, EARTH_RADIUS_M};

/// Geodesic level on which transport is solved (2,562 cells, ~440 km).
pub const TRANSPORT_LEVEL: u32 = 4;

pub struct Solver {
    offsets: Vec<u32>,
    nbrs: Vec<u32>,
    /// Geometric coupling per neighbour entry: shared edge length /
    /// (centre distance × cell area), 1/m².
    g: Vec<f64>,
    /// For each fine cell: nearby coarse cells and smoothing weights (sum 1),
    /// so fields handed down from the coarse grid vary smoothly instead of
    /// in 440 km blocks.
    blend_idx: Vec<[u32; 7]>,
    blend_w: Vec<[f64; 7]>,
}

fn arc(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    d.clamp(-1.0, 1.0).acos() * EARTH_RADIUS_M
}

impl Solver {
    pub fn new(fine: &Grid) -> Solver {
        let coarse = Grid::new(TRANSPORT_LEVEL);
        let n = coarse.len();
        let parent = parent_map(fine, &coarse);
        // Each coarse cell's area is taken as the sum of its fine cells, so
        // the transport handed down to fine cells sums to exactly zero.
        let mut area = vec![0.0; n];
        for (i, &p) in parent.iter().enumerate() {
            area[p as usize] += fine.areas_m2[i];
        }
        let mut g = Vec::with_capacity(coarse.nbrs.len());
        for i in 0..n {
            let (s, e) = (coarse.offsets[i] as usize, coarse.offsets[i + 1] as usize);
            for j in s..e {
                // The edge shared with neighbour j runs between corners j-1 and j.
                let prev = if j == s { e - 1 } else { j - 1 };
                let edge = arc(
                    coarse.corners[coarse.corner_ids[prev] as usize],
                    coarse.corners[coarse.corner_ids[j] as usize],
                );
                let dist = arc(coarse.centers[i], coarse.centers[coarse.nbrs[j] as usize]);
                g.push(edge / (dist * area[i]));
            }
        }
        // Smoothing weights: the coarse cell and its neighbours, weighted by
        // exp(-(distance / spacing)²), spacing = mean coarse centre distance.
        let spacing = {
            let (mut sum, mut cnt) = (0.0, 0.0);
            for i in 0..n {
                for j in coarse.offsets[i] as usize..coarse.offsets[i + 1] as usize {
                    sum += arc(coarse.centers[i], coarse.centers[coarse.nbrs[j] as usize]);
                    cnt += 1.0;
                }
            }
            sum / cnt
        };
        let mut blend_idx = Vec::with_capacity(fine.len());
        let mut blend_w = Vec::with_capacity(fine.len());
        for i in 0..fine.len() {
            let p = parent[i] as usize;
            let mut idx = [p as u32; 7];
            let mut w = [0.0f64; 7];
            let cands = std::iter::once(p as u32).chain(
                coarse.nbrs[coarse.offsets[p] as usize..coarse.offsets[p + 1] as usize]
                    .iter()
                    .copied(),
            );
            for (k, q) in cands.enumerate() {
                let d = arc(fine.centers[i], coarse.centers[q as usize]) / spacing;
                idx[k] = q;
                w[k] = (-d * d).exp();
            }
            let total: f64 = w.iter().sum();
            w.iter_mut().for_each(|v| *v /= total);
            blend_idx.push(idx);
            blend_w.push(w);
        }
        Solver {
            blend_idx,
            blend_w,
            offsets: coarse.offsets.clone(),
            nbrs: coarse.nbrs.clone(),
            g,
        }
    }

    /// Coarse-grid geometry: neighbour lists (CSR) and coupling g per
    /// neighbour entry (edge / (centre distance × area), 1/m²).
    pub fn coarse(&self) -> (&[u32], &[u32], &[f64]) {
        (&self.offsets, &self.nbrs, &self.g)
    }

    /// Smoothing weights handing a coarse field down to fine cell `i`.
    pub fn blend(&self, i: usize) -> (&[u32; 7], &[f64; 7]) {
        (&self.blend_idx[i], &self.blend_w[i])
    }
}
