//! Writes the map grid's cell centres for the data pipeline.
//!
//! Usage: gridgen <level> <out-file>
//! Output: little-endian f64 triples (latitude °, longitude °, area m²), one
//! per cell, in cell-index order. Also prints a summary of cell areas.

use std::io::Write;
use world::grid::{to_lat_lon, Grid};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: gridgen <level> <out-file>");
        std::process::exit(2);
    }
    let level: u32 = args[1].parse().expect("level must be a number");
    let grid = Grid::new(level);

    let mut out = std::io::BufWriter::new(std::fs::File::create(&args[2]).expect("create output"));
    for (&c, &area) in grid.centers.iter().zip(&grid.areas_m2) {
        let (lat, lon) = to_lat_lon(c);
        for v in [lat, lon, area] {
            out.write_all(&v.to_le_bytes()).unwrap();
        }
    }

    let (min, max) = grid
        .areas_m2
        .iter()
        .fold((f64::MAX, 0.0_f64), |(a, b), &x| (a.min(x), b.max(x)));
    let mean = grid.areas_m2.iter().sum::<f64>() / grid.len() as f64;
    println!(
        "level {level}: {} cells, area km² min {:.1} mean {:.1} max {:.1}",
        grid.len(),
        min / 1e6,
        mean / 1e6,
        max / 1e6
    );
}
