//! Implicit monthly heat balance with large-scale heat transport.
//!
//! Each month every climate cell's new temperature T satisfies
//!     c_i T_i = b_i + H_i
//! where c_i collects the cell's own terms (heat capacity / time step plus
//! the increase of outgoing heat with temperature), b_i its forcing, and H_i
//! the heat carried in by the atmosphere and ocean, represented as diffusion:
//!     H = K ∇²T.
//!
//! Diffusion of this strength smooths temperature over a length of about
//! R sqrt(D/B) ≈ 3,400 km (R the Earth's radius, D ≈ 0.6 W/m²/K the
//! diffusion coefficient on the unit sphere, B ≈ 2.1 W/m²/K the rise in
//! outgoing heat per degree; North et al. 1981). Transport is therefore
//! solved on a coarser grid of ~440 km cells, which resolves that length
//! scale, implicitly (stable with month-long steps). Each fine cell then
//! receives its coarse cell's transport per square metre and also exchanges
//! heat with the other fine cells of its coarse cell, at the same coupling
//! strength the coarse cells have with their neighbours. That exchange is
//! solved so it sums to zero within each coarse cell, so energy is conserved
//! exactly. Without it, an island in a mostly-ocean region would get no heat
//! from the sea around it.

use world::grid::{parent_map, Grid, EARTH_RADIUS_M};

/// Geodesic level on which transport is solved (2,562 cells, ~440 km).
pub const TRANSPORT_LEVEL: u32 = 4;

pub struct Solver {
    parent: Vec<u32>,
    /// Fine cells of each coarse cell (CSR).
    child_off: Vec<u32>,
    children: Vec<u32>,
    fine_area: Vec<f64>,
    offsets: Vec<u32>,
    nbrs: Vec<u32>,
    /// Geometric coupling per neighbour entry: shared edge length /
    /// (centre distance × cell area), 1/m².
    g: Vec<f64>,
    area: Vec<f64>,
    c: Vec<f64>,
    b: Vec<f64>,
    t: Vec<f64>,
    /// Last solution for each calendar month: the best starting guess, as
    /// a month resembles the same month a year earlier more than the
    /// previous month.
    t_by_month: Vec<Vec<f64>>,
    /// Mean of Σ_j g_ij over coarse cells: the geometric strength of local
    /// exchange (1/m²).
    local_g: f64,
    /// For each fine cell: nearby coarse cells and smoothing weights (sum 1),
    /// so fields handed down from the coarse grid vary smoothly instead of
    /// in 440 km blocks.
    blend_idx: Vec<[u32; 7]>,
    blend_w: Vec<[f64; 7]>,
    // Reused working storage (fine grid).
    hf: Vec<f64>,
    ts: Vec<f64>,
    /// Passes the last solve needed (for tests and tuning).
    pub last_sweeps: usize,
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
        let local_g = (0..n)
            .map(|i| {
                (coarse.offsets[i] as usize..coarse.offsets[i + 1] as usize)
                    .map(|j| g[j])
                    .sum::<f64>()
            })
            .sum::<f64>()
            / n as f64;
        let mut child_off = vec![0u32; n + 1];
        for &p in &parent {
            child_off[p as usize + 1] += 1;
        }
        for p in 0..n {
            child_off[p + 1] += child_off[p];
        }
        let mut fill = child_off.clone();
        let mut children = vec![0u32; parent.len()];
        for (i, &p) in parent.iter().enumerate() {
            children[fill[p as usize] as usize] = i as u32;
            fill[p as usize] += 1;
        }
        Solver {
            child_off,
            children,
            hf: vec![0.0; fine.len()],
            ts: vec![0.0; fine.len()],
            blend_idx,
            blend_w,
            local_g,
            parent,
            fine_area: fine.areas_m2.clone(),
            offsets: coarse.offsets.clone(),
            nbrs: coarse.nbrs.clone(),
            g,
            area,
            c: vec![0.0; n],
            b: vec![0.0; n],
            t: vec![f64::NAN; n],
            t_by_month: vec![Vec::new(); 12],
            last_sweeps: 0,
        }
    }

    /// Solves the month: `c`, `b` per fine cell; `k` = K in W/K (diffusion
    /// coefficient × R²). Writes new temperatures to `t` and the transport
    /// each cell received (W/m²) to `transport`.
    pub fn solve(
        &mut self,
        month: usize,
        c: &[f64],
        b: &[f64],
        k: f64,
        t: &mut [f64],
        transport: &mut [f64],
    ) {
        if self.t_by_month[month].len() == self.t.len() {
            self.t.copy_from_slice(&self.t_by_month[month]);
        }
        let n = self.area.len();
        self.c.iter_mut().for_each(|v| *v = 0.0);
        self.b.iter_mut().for_each(|v| *v = 0.0);
        for i in 0..c.len() {
            let p = self.parent[i] as usize;
            let a = self.fine_area[i];
            self.c[p] += c[i] * a;
            self.b[p] += b[i] * a;
        }
        for p in 0..n {
            self.c[p] /= self.area[p];
            self.b[p] /= self.area[p];
            if self.t[p].is_nan() {
                self.t[p] = self.b[p] / self.c[p];
            }
        }
        // Over-relaxed Gauss-Seidel from last month's solution until the
        // largest change in a pass is below 0.01 °C.
        const OMEGA: f64 = 1.7;
        let mut sweeps = 0;
        for _ in 0..500 {
            sweeps += 1;
            let mut max_change = 0.0f64;
            for p in 0..n {
                let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
                let (mut sg, mut sgt) = (0.0, 0.0);
                for j in s..e {
                    sg += self.g[j];
                    sgt += self.g[j] * self.t[self.nbrs[j] as usize];
                }
                let gs = (self.b[p] + k * sgt) / (self.c[p] + k * sg);
                let new = self.t[p] + OMEGA * (gs - self.t[p]);
                max_change = max_change.max((new - self.t[p]).abs());
                self.t[p] = new;
            }
            if max_change < 1e-2 {
                break;
            }
        }
        self.last_sweeps = sweeps;
        self.t_by_month[month].clone_from(&self.t);
        let h: Vec<f64> = (0..n)
            .map(|p| {
                let (s, e) = (self.offsets[p] as usize, self.offsets[p + 1] as usize);
                (s..e)
                    .map(|j| k * self.g[j] * (self.t[self.nbrs[j] as usize] - self.t[p]))
                    .sum()
            })
            .collect();
        // Transport per fine cell is its coarse cell's value. (Steps between
        // coarse cells move temperature by only H / kl, about 0.1 °C, because
        // the local exchange below is strong; the reference temperature is
        // what needs smoothing.)
        //
        // Local exchange with a smooth reference temperature T*_i:
        //   T_i = (r_i + kl T*_i) / (c_i + kl),  r_i = b_i + H_i.
        // First T* per coarse cell so the exchange sums to zero there:
        //   T*_p = Σ A r/(c+kl) / Σ A c/(c+kl),
        // then blended to fine cells and shifted per coarse cell (δ_p) so
        // the exchange still sums to zero exactly.
        let mut hf = std::mem::take(&mut self.hf);
        let mut ts = std::mem::take(&mut self.ts);
        let kl = k * self.local_g;
        let (parent, child_off, children, area) = (
            &self.parent,
            &self.child_off,
            &self.children,
            &self.fine_area,
        );
        hf.iter_mut()
            .enumerate()
            .for_each(|(i, v)| *v = h[parent[i] as usize]);
        let (num, den): (Vec<f64>, Vec<f64>) = (0..n)
            .into_iter()
            .map(|p| {
                let (mut nu, mut de) = (0.0, 0.0);
                for &i in &children[child_off[p] as usize..child_off[p + 1] as usize] {
                    let i = i as usize;
                    let r = b[i] + hf[i];
                    nu += area[i] * r / (c[i] + kl);
                    de += area[i] * c[i] / (c[i] + kl);
                }
                (nu, de)
            })
            .unzip();
        let t_star: Vec<f64> = (0..n).map(|p| num[p] / den[p]).collect();
        let (bidx, bw) = (&self.blend_idx, &self.blend_w);
        ts.iter_mut().enumerate().for_each(|(i, v)| {
            let (idx, w) = (&bidx[i], &bw[i]);
            *v = (0..7).map(|k| w[k] * t_star[idx[k] as usize]).sum();
        });
        let delta: Vec<f64> = (0..n)
            .into_iter()
            .map(|p| {
                let mut d = 0.0;
                for &i in &children[child_off[p] as usize..child_off[p + 1] as usize] {
                    let i = i as usize;
                    d += area[i] * (b[i] + hf[i] - c[i] * ts[i]) / (c[i] + kl);
                }
                d / den[p]
            })
            .collect();
        t.iter_mut()
            .zip(transport.iter_mut())
            .enumerate()
            .for_each(|(i, (ti, tr))| {
                let t_ref = ts[i] + delta[parent[i] as usize];
                let r = b[i] + hf[i];
                *ti = (r + kl * t_ref) / (c[i] + kl);
                *tr = hf[i] + kl * (t_ref - *ti);
            });
        self.hf = hf;
        self.ts = ts;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup() -> (Grid, Vec<f64>, Vec<f64>) {
        let g = Grid::new(6);
        let n = g.len();
        // Land-like and ocean-like local terms, sunshine-like forcing.
        let c: Vec<f64> = (0..n)
            .map(|i| {
                if i % 3 == 0 {
                    2.0e8 / 2.6e6 + 2.1
                } else {
                    1.0e7 / 2.6e6 + 2.1
                }
            })
            .collect();
        let b: Vec<f64> = g
            .centers
            .iter()
            .map(|p| 400.0 * (1.0 - p[2] * p[2]) - 200.0)
            .collect();
        (g, c, b)
    }

    #[test]
    fn conserves_energy_and_moves_heat_poleward() {
        let (g, _, b) = setup();
        let c = vec![1.0e8 / 2.6e6 + 2.1; g.len()];
        let mut s = Solver::new(&g);
        let n = g.len();
        let (mut t, mut h) = (vec![0.0; n], vec![0.0; n]);
        s.solve(0, &c, &b, 0.6 * EARTH_RADIUS_M.powi(2), &mut t, &mut h);
        // Transport only moves heat: it sums to zero over the globe.
        let net: f64 =
            (0..n).map(|i| h[i] * g.areas_m2[i]).sum::<f64>() / g.areas_m2.iter().sum::<f64>();
        assert!(net.abs() < 0.01, "net transport {net} W/m²");
        // Tropics lose heat, poles gain it.
        let eq = g.nearest_from([1.0, 0.0, 0.0], 0);
        let pole = g.nearest_from([0.0, 0.0, 1.0], 0);
        assert!(
            h[eq] < -10.0 && h[pole] > 10.0,
            "equator {} pole {}",
            h[eq],
            h[pole]
        );
    }

    #[test]
    fn a_month_is_cheap_when_warm() {
        let (g, c, b) = setup();
        let mut s = Solver::new(&g);
        let n = g.len();
        let (mut t, mut h) = (vec![0.0; n], vec![0.0; n]);
        let k = 0.6 * EARTH_RADIUS_M.powi(2);
        s.solve(0, &c, &b, k, &mut t, &mut h);
        let b2: Vec<f64> = b.iter().map(|v| v + 5.0).collect();
        let start = std::time::Instant::now();
        let reps = 20;
        for r in 0..reps {
            let bb = if r % 2 == 0 { &b2 } else { &b };
            s.solve(r % 12, &c, bb, k, &mut t, &mut h);
        }
        let ms = start.elapsed().as_secs_f64() * 1e3 / reps as f64;
        println!("one month: {ms:.3} ms, {} passes", s.last_sweeps);
        // Budget: a few tenths of a millisecond on a typical PC; allow for
        // slow shared machines (GitHub, Claude's workspace).
        if !cfg!(debug_assertions) {
            assert!(ms < 3.0, "{ms} ms per month");
        }
    }
}
