//! The world data file: terrain, sea and ice for every map cell at two
//! moments, 1 January 10,000 BCE (the start of a run) and today. Built by
//! `tools/data/build_world.py`; see that file for sources and methods.
//!
//! Layout (little-endian): b"SIM3WRLD", u32 format version, u32 grid level,
//! u32 cell count, u32 epoch count, then per epoch: i32 astronomical year,
//! i16[n] mean elevation (m), i16[n] lowest point (m), i16[n] highest point
//! (m), u8[n] share of the cell above sea level (0-255 = 0-1), u16[n] mean
//! ice thickness (m). Elevations are relative to sea level at that time.
//!
//! Format version 2 adds tagged sections after the epochs (4-byte tag, u32
//! byte length, payload):
//! - `OBST`: observed 1981-2010 monthly mean temperature (CHELSA V2.1) for
//!   mostly-land cells: u32 count, u32[count] cell indices, i16[count * 12]
//!   in hundredths of °C, 12 months per cell.
//! - `CO2 `: ice-core CO2 (Bereiter et al. 2015): u32 count, then f32 pairs
//!   (years before 1950, ppm), oldest last.
//! - `CO2P`: f32, mean Mauna Loa CO2 1981-2010 (ppm).
//! - `ORBT`: Berger & Loutre (1991) orbit: u32 count, then f32 quads (kyr
//!   relative to 1950, eccentricity, longitude of perihelion from the
//!   northern spring equinox in degrees, obliquity in degrees).

pub struct Epoch {
    /// Astronomical year (-9999 = 10,000 BCE).
    pub year: i32,
    pub elev_mean_m: Vec<i16>,
    pub elev_min_m: Vec<i16>,
    pub elev_max_m: Vec<i16>,
    /// Share of the cell above sea level, 0-255 for 0-1.
    pub land_share: Vec<u8>,
    pub ice_thickness_m: Vec<u16>,
}

pub struct ObservedTemperature {
    pub cells: Vec<u32>,
    /// °C, 12 values per cell in `cells` order.
    pub monthly_c: Vec<f32>,
}

#[derive(Clone, Copy, Debug)]
pub struct OrbitRow {
    pub kyr: f32,
    pub eccentricity: f32,
    pub perihelion_deg: f32,
    pub obliquity_deg: f32,
}

pub struct WorldData {
    pub grid_level: u32,
    pub cells: usize,
    pub epochs: Vec<Epoch>,
    pub observed_temperature: Option<ObservedTemperature>,
    /// (years before 1950, ppm), sorted by age.
    pub co2_record: Vec<(f32, f32)>,
    pub co2_1981_2010_ppm: Option<f32>,
    pub orbit: Vec<OrbitRow>,
}

#[derive(Debug)]
pub struct FormatError(pub String);

impl std::fmt::Display for FormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "world data: {}", self.0)
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], FormatError> {
        let end = self.pos.checked_add(n).filter(|&e| e <= self.bytes.len());
        let end = end.ok_or_else(|| FormatError("file is truncated".into()))?;
        let s = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(s)
    }
    fn u32(&mut self) -> Result<u32, FormatError> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn i32(&mut self) -> Result<i32, FormatError> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn i16s(&mut self, n: usize) -> Result<Vec<i16>, FormatError> {
        Ok(self
            .take(n * 2)?
            .as_chunks::<2>()
            .0
            .iter()
            .map(|c| i16::from_le_bytes(*c))
            .collect())
    }
    fn f32s(&mut self, n: usize) -> Result<Vec<f32>, FormatError> {
        Ok(self
            .take(n * 4)?
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| f32::from_le_bytes(*c))
            .collect())
    }
    fn u32s(&mut self, n: usize) -> Result<Vec<u32>, FormatError> {
        Ok(self
            .take(n * 4)?
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| u32::from_le_bytes(*c))
            .collect())
    }
    fn u16s(&mut self, n: usize) -> Result<Vec<u16>, FormatError> {
        Ok(self
            .take(n * 2)?
            .as_chunks::<2>()
            .0
            .iter()
            .map(|c| u16::from_le_bytes(*c))
            .collect())
    }
}

impl WorldData {
    pub fn parse(bytes: &[u8]) -> Result<WorldData, FormatError> {
        let mut r = Reader { bytes, pos: 0 };
        if r.take(8)? != b"SIM3WRLD" {
            return Err(FormatError("not a sim3 world file".into()));
        }
        let version = r.u32()?;
        if version != 1 && version != 2 {
            return Err(FormatError(format!("unsupported format version {version}")));
        }
        let grid_level = r.u32()?;
        let cells = r.u32()? as usize;
        if cells != crate::grid::cell_count(grid_level) {
            return Err(FormatError("cell count does not match grid level".into()));
        }
        let epoch_count = r.u32()?;
        let mut epochs = Vec::new();
        for _ in 0..epoch_count {
            epochs.push(Epoch {
                year: r.i32()?,
                elev_mean_m: r.i16s(cells)?,
                elev_min_m: r.i16s(cells)?,
                elev_max_m: r.i16s(cells)?,
                land_share: r.take(cells)?.to_vec(),
                ice_thickness_m: r.u16s(cells)?,
            });
        }
        let mut out = WorldData {
            grid_level,
            cells,
            epochs,
            observed_temperature: None,
            co2_record: Vec::new(),
            co2_1981_2010_ppm: None,
            orbit: Vec::new(),
        };
        while r.pos < bytes.len() {
            let tag: [u8; 4] = r.take(4)?.try_into().unwrap();
            let len = r.u32()? as usize;
            let mut s = Reader {
                bytes: r.take(len)?,
                pos: 0,
            };
            match &tag {
                b"OBST" => {
                    let n = s.u32()? as usize;
                    let cells = s.u32s(n)?;
                    let monthly_c = s
                        .i16s(n * 12)?
                        .into_iter()
                        .map(|v| v as f32 / 100.0)
                        .collect();
                    out.observed_temperature = Some(ObservedTemperature { cells, monthly_c });
                }
                b"CO2 " => {
                    let n = s.u32()? as usize;
                    let v = s.f32s(n * 2)?;
                    out.co2_record = v.as_chunks::<2>().0.iter().map(|p| (p[0], p[1])).collect();
                }
                b"CO2P" => out.co2_1981_2010_ppm = Some(s.f32s(1)?[0]),
                b"ORBT" => {
                    let n = s.u32()? as usize;
                    let v = s.f32s(n * 4)?;
                    out.orbit = v
                        .as_chunks::<4>()
                        .0
                        .iter()
                        .map(|q| OrbitRow {
                            kyr: q[0],
                            eccentricity: q[1],
                            perihelion_deg: q[2],
                            obliquity_deg: q[3],
                        })
                        .collect();
                }
                _ => {} // unknown sections are skipped, so older programs can read newer files
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_the_synthetic_test_file_sections() {
        // Built by tools/data/build_world.py from synthetic inputs, if present.
        let Ok(bytes) = std::fs::read("/tmp/dtest/world.bin") else {
            return;
        };
        let d = WorldData::parse(&bytes).unwrap();
        assert!(d.observed_temperature.is_some());
        assert_eq!(d.orbit.len(), 16);
        assert!(d.co2_record.len() > 50);
        assert!(d.co2_1981_2010_ppm.unwrap() > 300.0);
    }
}
