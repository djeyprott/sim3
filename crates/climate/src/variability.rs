//! Year-to-year variability: volcanic eruptions and the El Niño-like cycle
//! (Phase 1, step 5). Both are random, with the frequencies and sizes of
//! the real record, so each run has its own climate history.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Random numbers

/// A small, fast random number generator (SplitMix64, Steele et al. 2014).
/// Statistical quality is ample for drawing eruptions and weather kicks.
#[derive(Clone)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Standard normal (Box-Muller).
    pub fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }

    /// Poisson-distributed count with the given mean (Knuth's method; the
    /// means used here are well below 1).
    pub fn poisson(&mut self, mean: f64) -> u32 {
        let l = (-mean).exp();
        let mut k = 0;
        let mut p = self.uniform();
        while p > l {
            k += 1;
            p *= self.uniform();
        }
        k
    }
}

// ---------------------------------------------------------------------------
// Volcanic eruptions: the Easy Volcanic Aerosol model (EVA v1.0; Toohey,
// Stevens, Schmidt & Timmreck 2016, Geosci. Model Dev. 9, 4049), the model
// HolVol itself uses to turn ice-core sulfur into veil thickness. Three
// stratospheric boxes (southern extratropics, tropics, northern
// extratropics); sulfur gas turns into sulfate droplets, which mix between
// the boxes and fall out. Seasonal changes in mixing are left out.

/// Box boundaries, degrees latitude (EVA Table 1).
pub const TROPICAL_EDGE_DEG: f64 = 25.0;
/// Effective production of droplets from sulfur gas, months (EVA: 180
/// days).
const TAU_PROD_MONTHS: f64 = 180.0 / 30.44;
/// Loss of droplets from the stratosphere, months (EVA: 330 days).
const TAU_LOSS_MONTHS: f64 = 330.0 / 30.44;
/// Mixing between the tropical and each extratropical box, months (EVA).
const TAU_MIX_MONTHS: f64 = 15.0;
/// One-way transport from the tropics to each extratropical box by the
/// stratosphere's overturning (residual) circulation, months (EVA).
const TAU_RES_MONTHS: f64 = 17.0;
/// Global mean veil thickness (aerosol optical depth at 550 nm) per Tg of
/// sulfur held as droplets, for loads up to Pinatubo's (EVA, fitted to
/// Pinatubo).
const AOD_PER_TGS: f64 = 0.0364;
/// Pinatubo's sulfur injection, Tg S (EVA).
pub const PINATUBO_TGS: f64 = 9.0;
/// Bigger loads make bigger droplets, which dim less per unit mass: above
/// Pinatubo's peak load the veil grows with the 2/3 power of the mass
/// (Crowley & Unterman 2013, as adopted in EVA eq. 13), continuous at that
/// load.
const AOD_LARGE_POWER: f64 = 2.0 / 3.0;
/// Sunlight lost per unit of global veil thickness, W/m² (effective
/// radiative forcing; Marshall et al. 2020, 82 simulations).
const FORCING_PER_AOD: f64 = 20.5;
/// Sunlight absorbed by the Earth on average, W/m² (CERES; the same figure
/// the two-layer model uses).
const ABSORBED_SUNLIGHT: f64 = 239.0;
/// Extratropical veils: Gaussian in latitude centred on 45° with width 14°
/// (EVA).
const EXTRATROPICAL_CENTRE_DEG: f64 = 45.0;
const EXTRATROPICAL_WIDTH_DEG: f64 = 14.0;
/// Tropical veil: Gaussian centred on the equator, its width half the box
/// edge, so the box holds 95% of it (EVA fits this shape to Pinatubo's first
/// months; this matches its spread to about ±25°).
const TROPICAL_WIDTH_DEG: f64 = TROPICAL_EDGE_DEG / 2.0;

/// A recorded set of eruptions, used only for how often and how big.
#[derive(Clone)]
pub struct EruptionStats {
    pub sizes_tgs: Vec<f32>,
    pub per_month: f64,
}

#[derive(Clone)]
pub struct Volcanism {
    stats: Option<EruptionStats>,
    /// Where eruptions can happen: climate cells and cumulative weights.
    sites: Vec<(u32, f64)>,
    site_total: f64,
    /// Sulfur gas and droplets in each box, Tg S.
    gas: [f64; 3],
    droplets: [f64; 3],
    /// Peak droplet load after a Pinatubo-sized eruption, Tg S: the point
    /// above which the veil grows more slowly.
    large_load: f64,
    /// Normalised veil shape of each box on a 1° latitude table.
    shape: [Vec<f64>; 3],
}

/// One eruption, as it happened.
#[derive(Clone, Debug)]
pub struct Eruption {
    pub cell: u32,
    pub lat: f64,
    pub lon: f64,
    pub sulfur_tgs: f64,
}

fn gaussian(x: f64, centre: f64, width: f64) -> f64 {
    (-0.5 * ((x - centre) / width).powi(2)).exp()
}

impl Volcanism {
    /// `sites` pairs climate cells with how likely an eruption is there.
    pub fn new(stats: Option<EruptionStats>, sites: Vec<(u32, f64)>) -> Volcanism {
        let mut cum = Vec::with_capacity(sites.len());
        let mut total = 0.0;
        for (c, w) in sites {
            if w > 0.0 {
                total += w;
                cum.push((c, total));
            }
        }
        // Peak of M(t) = (1 - exp(-t/a)) exp(-t/b): exp(-t/a) = a / (a + b).
        let (a, b) = (TAU_PROD_MONTHS, TAU_LOSS_MONTHS);
        let large_load = PINATUBO_TGS * b / (a + b) * (a / (a + b)).powf(a / b);
        let lats: Vec<f64> = (0..181).map(|k| k as f64 - 90.0).collect();
        let raw = |f: &dyn Fn(f64) -> f64| -> Vec<f64> {
            let v: Vec<f64> = lats.iter().map(|&l| f(l)).collect();
            let wsum: f64 = lats.iter().map(|l| l.to_radians().cos()).sum();
            let mean = v
                .iter()
                .zip(&lats)
                .map(|(x, l)| x * l.to_radians().cos())
                .sum::<f64>()
                / wsum;
            v.into_iter().map(|x| x / mean).collect()
        };
        let shape = [
            raw(&|l| gaussian(l, -EXTRATROPICAL_CENTRE_DEG, EXTRATROPICAL_WIDTH_DEG)),
            raw(&|l| gaussian(l, 0.0, TROPICAL_WIDTH_DEG)),
            raw(&|l| gaussian(l, EXTRATROPICAL_CENTRE_DEG, EXTRATROPICAL_WIDTH_DEG)),
        ];
        Volcanism {
            stats,
            sites: cum,
            site_total: total,
            gas: [0.0; 3],
            droplets: [0.0; 3],
            large_load,
            shape,
        }
    }

    fn box_of(lat: f64) -> usize {
        if lat < -TROPICAL_EDGE_DEG {
            0
        } else if lat > TROPICAL_EDGE_DEG {
            2
        } else {
            1
        }
    }

    /// Adds an eruption's sulfur to the stratosphere above it.
    pub fn inject(&mut self, lat: f64, sulfur_tgs: f64) {
        self.gas[Self::box_of(lat)] += sulfur_tgs;
    }

    /// Draws this month's eruptions (none if the run has no record).
    pub fn draw(&mut self, rng: &mut Rng, lat_lon: impl Fn(u32) -> (f64, f64)) -> Vec<Eruption> {
        let Some(stats) = &self.stats else {
            return Vec::new();
        };
        if self.sites.is_empty() || stats.sizes_tgs.is_empty() {
            return Vec::new();
        }
        let k = rng.poisson(stats.per_month);
        let mut out = Vec::new();
        for _ in 0..k {
            let size = stats.sizes_tgs[(rng.uniform() * stats.sizes_tgs.len() as f64) as usize];
            let x = rng.uniform() * self.site_total;
            let i = self
                .sites
                .partition_point(|&(_, c)| c < x)
                .min(self.sites.len() - 1);
            let cell = self.sites[i].0;
            let (lat, lon) = lat_lon(cell);
            out.push(Eruption {
                cell,
                lat,
                lon,
                sulfur_tgs: size as f64,
            });
        }
        for e in &out {
            self.inject(e.lat, e.sulfur_tgs);
        }
        out
    }

    /// Advances the veil by one month (EVA's monthly box model).
    pub fn step(&mut self) {
        let prod = 1.0 - (-1.0 / TAU_PROD_MONTHS).exp();
        let keep = (-1.0 / TAU_LOSS_MONTHS).exp();
        let mix = 1.0 / TAU_MIX_MONTHS;
        let res = 1.0 / TAU_RES_MONTHS;
        let new: [f64; 3] = std::array::from_fn(|b| self.gas[b] * prod);
        for b in 0..3 {
            self.gas[b] = (self.gas[b] - new[b]) * keep;
        }
        let d = self.droplets;
        let t = d[1];
        self.droplets[1] =
            (t + new[1]) * keep - (t - d[0]) * mix - (t - d[2]) * mix - 2.0 * t * res;
        for b in [0, 2] {
            self.droplets[b] = (d[b] + new[b]) * keep + (t - d[b]) * mix + t * res;
        }
        for v in self.droplets.iter_mut().chain(self.gas.iter_mut()) {
            *v = v.max(0.0);
        }
    }

    /// Global mean veil thickness for a droplet load (Tg S).
    fn global_aod(&self, load: f64) -> f64 {
        if load <= self.large_load {
            AOD_PER_TGS * load
        } else {
            AOD_PER_TGS * self.large_load * (load / self.large_load).powf(AOD_LARGE_POWER)
        }
    }

    pub fn global_aod_now(&self) -> f64 {
        self.global_aod(self.droplets.iter().sum())
    }

    /// Veil thickness at a latitude.
    pub fn aod(&self, lat: f64) -> f64 {
        let load: f64 = self.droplets.iter().sum();
        if load <= 1e-9 {
            return 0.0;
        }
        let scale = self.global_aod(load) / load;
        let k = ((lat + 90.0).round() as usize).min(180);
        scale
            * (0..3)
                .map(|b| self.droplets[b] * self.shape[b][k])
                .sum::<f64>()
    }

    /// Share of sunlight the veil takes away at a latitude.
    pub fn dimming(&self, lat: f64) -> f64 {
        (FORCING_PER_AOD * self.aod(lat) / ABSORBED_SUNLIGHT).min(0.9)
    }
}

// ---------------------------------------------------------------------------
// The El Niño-like cycle: the recharge oscillator (Jin 1997) in its simplest
// form, fitted to observations by Burgers, Jin & van Oldenborgh (2005):
//
//   dT/dt = -2γ T + ω0 h + noise,   dh/dt = -ω0 T
//
// T is the sea surface temperature anomaly of the eastern equatorial basin
// and h the (scaled) heat content of the equatorial band. Observed Pacific:
// damping 1/γ = 24 months, period 2π/ω = 37 months with ω² = ω0² - γ².
// The timescales come from how long equatorial ocean waves take to cross the
// basin (Jin 1997), so they scale with the basin's width.

/// Observed Pacific fit (Burgers et al. 2005), months.
const ENSO_DAMPING_MONTHS: f64 = 24.0;
const ENSO_PERIOD_MONTHS: f64 = 37.0;
/// Width of the Pacific along the equator the fit refers to (Burgers et
/// al.'s heat content index spans 130°E-80°W, 150° of longitude), km.
const ENSO_REFERENCE_WIDTH_KM: f64 = 150.0 * 111.2;
/// Typical size (standard deviation) of the observed eastern Pacific
/// anomaly, °C (Niño3 index, 1981-2010; NOAA Climate Prediction Center).
pub const ENSO_REFERENCE_SIZE_C: f64 = 0.9;
/// Narrowest basin that can hold the oscillator: two equatorial
/// deformation radii (sqrt(c/β) with c = 2.8 m/s for the first baroclinic
/// mode, ≈ 350 km each), as the cycle relies on equatorial ocean waves.
const ENSO_MIN_WIDTH_KM: f64 = 2.0 * 350.0;
/// The anomaly spreads north and south of the equator with this width
/// (observed El Niño sea surface patterns reach about ±10° latitude).
pub const ENSO_LAT_WIDTH_DEG: f64 = 10.0;
/// NOAA's thresholds on the 3-month mean anomaly: an event ends below 0.5
/// °C; above 1.5 °C it is classed strong.
const ENSO_EVENT_END_C: f64 = 0.5;
pub const ENSO_STRONG_C: f64 = 1.5;

#[derive(Clone)]
pub struct EnsoBasin {
    /// Western and eastern edge along the equator (degrees east; the
    /// eastern may exceed 180 when the basin crosses the date line).
    pub west_lon: f64,
    pub east_lon: f64,
    pub width_km: f64,
    gamma: f64,
    omega0: f64,
    noise: f64,
    pub t: f64,
    h: f64,
    last3: [f64; 3],
    /// Strongest 3-month mean of the current event, and whether it has
    /// been reported.
    peak: f64,
    reported: bool,
}

/// A strong El Niño or La Niña at its peak.
#[derive(Clone, Debug)]
pub struct EnsoEvent {
    pub basin: usize,
    /// 3-month mean anomaly, °C (positive: El Niño).
    pub strength_c: f64,
}

#[derive(Clone, Default)]
pub struct Enso {
    pub basins: Vec<EnsoBasin>,
}

impl Enso {
    /// Finds the basins from which stretches of the equator are sea.
    /// `wet[c]` is the equator at longitude -179.5 + c (1° steps).
    pub fn new(wet: &[bool]) -> Enso {
        let n = wet.len();
        let mut basins = Vec::new();
        let mut segs = Vec::new();
        if wet.iter().all(|&w| w) {
            segs.push((0, n));
        } else if let Some(start) = (0..n).find(|&c| !wet[c]) {
            let mut k = 0;
            let mut c = start;
            while k < n {
                if wet[c] {
                    let s0 = c;
                    let mut len = 0;
                    while wet[c] && len < n {
                        c = (c + 1) % n;
                        len += 1;
                        k += 1;
                    }
                    segs.push((s0, len));
                } else {
                    c = (c + 1) % n;
                    k += 1;
                }
            }
        }
        for (s0, len) in segs {
            let width_km = len as f64 * 111.2;
            if width_km < ENSO_MIN_WIDTH_KM {
                continue;
            }
            let scale = width_km / ENSO_REFERENCE_WIDTH_KM;
            let gamma = 1.0 / (ENSO_DAMPING_MONTHS * scale);
            let omega = 2.0 * PI / (ENSO_PERIOD_MONTHS * scale);
            let omega0 = (omega * omega + gamma * gamma).sqrt();
            // The coupling between winds and sea grows with the stretch of
            // warm water the winds blow over, so the typical size scales
            // with width too. With noise only in T, the steady variance is
            // noise² / (4γ) (a short calculation from the equations above).
            let size = ENSO_REFERENCE_SIZE_C * scale.min(1.5);
            basins.push(EnsoBasin {
                west_lon: -180.0 + s0 as f64,
                east_lon: -180.0 + (s0 + len) as f64,
                width_km,
                gamma,
                omega0,
                noise: size * 2.0 * gamma.sqrt(),
                t: 0.0,
                h: 0.0,
                last3: [0.0; 3],
                peak: 0.0,
                reported: false,
            });
        }
        Enso { basins }
    }

    /// The widest basin (the one the viewer's El Niño line follows).
    pub fn main_basin(&self) -> Option<usize> {
        (0..self.basins.len())
            .max_by(|&a, &b| self.basins[a].width_km.total_cmp(&self.basins[b].width_km))
    }

    /// Advances one month; returns strong events that just peaked.
    pub fn step(&mut self, rng: &mut Rng) -> Vec<EnsoEvent> {
        let mut out = Vec::new();
        for (k, b) in self.basins.iter_mut().enumerate() {
            // Semi-implicit step of the damped oscillator (stable for any
            // month length), plus a random kick.
            let t = (b.t + b.omega0 * b.h + b.noise * rng.normal()) / (1.0 + 2.0 * b.gamma);
            b.h -= b.omega0 * t;
            b.t = t;
            b.last3 = [b.last3[1], b.last3[2], t];
            let m3 = b.last3.iter().sum::<f64>() / 3.0;
            let same_sign = m3 * b.peak > 0.0;
            if m3.abs() < ENSO_EVENT_END_C {
                b.peak = 0.0;
                b.reported = false;
            } else if !same_sign {
                // A new event (or a flip straight from one to the other).
                b.peak = m3;
                b.reported = false;
            } else if m3.abs() > b.peak.abs() {
                b.peak = m3;
            } else if !b.reported && b.peak.abs() >= ENSO_STRONG_C {
                // The 3-month mean has turned: the event has peaked.
                b.reported = true;
                out.push(EnsoEvent {
                    basin: k,
                    strength_c: b.peak,
                });
            }
        }
        out
    }

    /// Sea surface anomaly at a place, °C. Along the equator the eastern
    /// end carries the full anomaly and the western end the opposite at half
    /// strength (the observed El Niño pattern: the warm pool shifts east).
    pub fn anomaly(&self, lat: f64, lon: f64) -> f64 {
        let mut a = 0.0;
        for b in &self.basins {
            let mut x = lon;
            while x < b.west_lon {
                x += 360.0;
            }
            if x > b.east_lon {
                continue;
            }
            let f = (x - b.west_lon) / (b.east_lon - b.west_lon);
            a += b.t * (1.5 * f - 0.5) * gaussian(lat, 0.0, ENSO_LAT_WIDTH_DEG);
        }
        a
    }
}

/// The cycle's length, months: twice the lag at which a monthly series is
/// least like itself (for a damped oscillator the autocorrelation goes as
/// e^(-γτ) cos(ωτ), lowest at half a period). Lags from 3 months to 8
/// years are searched.
pub fn period_months(series: &[f64]) -> Option<f64> {
    let n = series.len();
    if n < 120 {
        return None;
    }
    let mean = series.iter().sum::<f64>() / n as f64;
    let var: f64 = series.iter().map(|x| (x - mean).powi(2)).sum();
    if var <= 0.0 {
        return None;
    }
    let ac = |lag: usize| {
        (0..n - lag)
            .map(|i| (series[i] - mean) * (series[i + lag] - mean))
            .sum::<f64>()
            / var
    };
    let (lag, v) = (3..=48.min(n / 2))
        .map(|l| (l, ac(l)))
        .min_by(|a, b| a.1.total_cmp(&b.1))?;
    (v < 0.0).then_some(2.0 * lag as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinatubo_veil_peaks_and_fades() {
        let mut v = Volcanism::new(None, vec![]);
        v.inject(15.0, PINATUBO_TGS);
        let mut peak: f64 = 0.0;
        let mut at = 0;
        for m in 0..48 {
            v.step();
            let a = v.global_aod_now();
            if a > peak {
                peak = a;
                at = m;
            }
        }
        // Observed: global mean peak about 0.1-0.15, some 6-12 months after.
        assert!((0.08..0.2).contains(&peak), "peak {peak}");
        assert!((3..14).contains(&at), "peak month {at}");
        assert!(v.global_aod_now() < 0.02);
    }

    #[test]
    fn oscillator_has_the_observed_rhythm() {
        let mut wet = vec![false; 360];
        for w in wet.iter_mut().take(150) {
            *w = true;
        }
        let mut e = Enso::new(&wet);
        assert_eq!(e.basins.len(), 1);
        let mut rng = Rng::new(7);
        let n = 12 * 3000;
        let series: Vec<f64> = (0..n)
            .map(|_| {
                e.step(&mut rng);
                e.basins[0].t
            })
            .collect();
        let sd = (series.iter().map(|t| t * t).sum::<f64>() / n as f64).sqrt();
        let period_years = period_months(&series).unwrap() / 12.0;
        assert!((0.7..1.1).contains(&sd), "sd {sd}");
        assert!((2.0..7.0).contains(&period_years), "period {period_years}");
    }
}
