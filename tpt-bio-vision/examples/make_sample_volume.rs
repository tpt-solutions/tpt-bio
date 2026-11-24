//! Writes a synthetic three-channel 3D volume as an OME-TIFF file, for the viewer example.
//!
//! ```text
//! cargo run -p tpt-bio-vision --example make_sample_volume -- sample-volume.ome.tif
//! ```

// Voxel code uses x, y, z and c for coordinates and channel, as the library does.
#![allow(clippy::many_single_char_names)]
//!
//! The volume is 96 x 96 x 32 voxels with three channels. Channel 0 holds bright nuclei,
//! channel 1 holds membrane shells around the same cells, and channel 2 holds a faint
//! diffuse marker that brightens inside the cells. Values are 16-bit. The noise is
//! deterministic, so the same command always writes the same file.

use std::error::Error;
use std::fs;

use tpt_bio_vision::image::Image;
use tpt_bio_vision::io::ome::encode_ome_tiff;

const WIDTH: usize = 96;
const HEIGHT: usize = 96;
const DEPTH: usize = 32;
const CHANNELS: usize = 3;

/// Cell centres `(x, y, z)` and radius, in voxels.
const CELLS: [(f64, f64, f64, f64); 6] = [
    (24.0, 24.0, 12.0, 11.0),
    (66.0, 26.0, 16.0, 10.0),
    (30.0, 66.0, 20.0, 12.0),
    (70.0, 70.0, 10.0, 9.0),
    (48.0, 48.0, 22.0, 13.0),
    (14.0, 84.0, 16.0, 7.0),
];

/// A small deterministic hash giving noise in `0..=255`.
fn noise(x: usize, y: usize, z: usize, c: usize) -> u16 {
    let mut h = (x as u64).wrapping_mul(0x9E37_79B9)
        ^ (y as u64).wrapping_mul(0x85EB_CA6B)
        ^ (z as u64).wrapping_mul(0xC2B2_AE35)
        ^ (c as u64).wrapping_mul(0x27D4_EB2F);
    h ^= h >> 15;
    h = h.wrapping_mul(0x2C1B_3C6D);
    h ^= h >> 12;
    (h & 0xFF) as u16
}

/// Value of one voxel in one channel.
fn voxel(c: usize, x: usize, y: usize, z: usize) -> u16 {
    let (fx, fy, fz) = (x as f64, y as f64, z as f64);
    let base = noise(x, y, z, c) / 4;
    let mut value = match c {
        0 => 200,
        1 => 100,
        _ => 400 + 8 * x as u16,
    };
    for &(cx, cy, cz, r) in &CELLS {
        let d = ((fx - cx).powi(2) + (fy - cy).powi(2) + (fz - cz).powi(2)).sqrt();
        match c {
            0 if d <= 0.5 * r => value = 3000,
            1 if d <= r && d >= 0.8 * r => value = 2000,
            2 if d <= r => value = value.max(1500),
            _ => {}
        }
    }
    value + base
}

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args_os()
        .nth(1)
        .map_or_else(|| "sample-volume.ome.tif".into(), std::path::PathBuf::from);

    // OME-TIFF stores planes in XYZCT order: Z varies fastest, then C. The volume is
    // therefore laid out channel by channel, then slice by slice.
    let mut data = Vec::with_capacity(WIDTH * HEIGHT * DEPTH * CHANNELS);
    for c in 0..CHANNELS {
        for z in 0..DEPTH {
            for y in 0..HEIGHT {
                for x in 0..WIDTH {
                    data.push(voxel(c, x, y, z));
                }
            }
        }
    }
    let planes = Image::new_3d(WIDTH, HEIGHT, DEPTH * CHANNELS, data)?;
    let bytes = encode_ome_tiff(&planes, (DEPTH, CHANNELS, 1), Some("synthetic cells"), false)?;
    fs::write(&path, bytes)?;
    println!("wrote {} ({WIDTH} x {HEIGHT} x {DEPTH} voxels, {CHANNELS} channels)", path.display());
    Ok(())
}
