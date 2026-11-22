//! The world data file: terrain, sea and ice for every map cell at two
//! moments, 1 January 10,000 BCE (the start of a run) and today. Built by
//! `tools/data/build_world.py`; see that file for sources and methods.
//!
//! Layout (little-endian): b"SIM3WRLD", u32 format version, u32 grid level,
//! u32 cell count, u32 epoch count, then per epoch: i32 astronomical year,
//! i16[n] mean elevation (m), i16[n] lowest point (m), i16[n] highest point
//! (m), u8[n] share of the cell above sea level (0-255 = 0-1), u16[n] mean
//! ice thickness (m). Elevations are relative to sea level at that time.

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

pub struct WorldData {
    pub grid_level: u32,
    pub cells: usize,
    pub epochs: Vec<Epoch>,
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
            .chunks_exact(2)
            .map(|c| i16::from_le_bytes([c[0], c[1]]))
            .collect())
    }
    fn u16s(&mut self, n: usize) -> Result<Vec<u16>, FormatError> {
        Ok(self
            .take(n * 2)?
            .chunks_exact(2)
            .map(|c| u16::from_le_bytes([c[0], c[1]]))
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
        if version != 1 {
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
        Ok(WorldData {
            grid_level,
            cells,
            epochs,
        })
    }
}
