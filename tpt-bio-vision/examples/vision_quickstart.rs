//! Vision quickstart: load a PNG, smooth it, threshold it, label the objects, filter
//! them, measure each one and write the results. Run with:
//!
//! ```text
//! cargo run -p tpt-bio-vision --example vision_quickstart [output-dir]
//! ```
//!
//! The input image is generated here, so no data files are needed. The same steps are
//! available from the command line as `tptbio vision segment` and `tptbio vision features`.

use std::error::Error;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use tpt_bio_vision::features::{intensity_features, shape_features_2d};
use tpt_bio_vision::image::Image;
use tpt_bio_vision::io::Planes;
use tpt_bio_vision::io::png::read_png;
use tpt_bio_vision::io::write::{save_png, save_tiff};
use tpt_bio_vision::label::{Connectivity2, label_2d};
use tpt_bio_vision::objects::{ObjectFilter, filter_objects};
use tpt_bio_vision::source::FileSource;
use tpt_bio_vision::spatial::gaussian_smooth;
use tpt_bio_vision::threshold::{apply_threshold, otsu_threshold};

type Failure = Box<dyn Error>;

/// A 96 x 64 grey image: dim background with noise, and bright discs of several sizes.
fn make_scene() -> Image<u8> {
    let (w, h) = (96usize, 64usize);
    let discs = [(14.0, 16.0, 6.0), (40.0, 20.0, 4.0), (70.0, 30.0, 7.0), (22.0, 48.0, 3.0)];
    let mut state: u32 = 12_345;
    let mut data = Vec::with_capacity(w * h);
    for y in 0..h {
        for x in 0..w {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            let noise = ((state >> 16) % 9) as u8;
            let inside = discs.iter().any(|&(cx, cy, r)| {
                let (dx, dy) = (x as f64 - cx, y as f64 - cy);
                dx * dx + dy * dy <= r * r
            });
            data.push(if inside { 190 + noise } else { 10 + noise });
        }
    }
    Image::new_2d(w, h, data).expect("data length matches the shape")
}

fn main() -> Result<(), Failure> {
    let dir: PathBuf = std::env::args_os()
        .nth(1)
        .map_or_else(|| std::env::temp_dir().join("tptbio-vision-quickstart"), PathBuf::from);
    fs::create_dir_all(&dir)?;

    // 1. Generate a PNG and load it back. Channel 0 is the only channel.
    let png = dir.join("scene.png");
    save_png(&png, &make_scene())?;
    let channel = match read_png(&FileSource::open(&png)?)? {
        Planes::U8(mut channels) => channels.swap_remove(0),
        _ => return Err("expected an 8-bit PNG".into()),
    };

    // 2. Smooth with a Gaussian (sigma 1) and convert to f64 for the analysis functions.
    let intensity: Image<f64> = channel.map(|&v| f64::from(v));
    let smooth = gaussian_smooth(&intensity, 1.0)?;

    // 3. Otsu threshold, then the foreground mask.
    let threshold = otsu_threshold(&smooth, 256)?.ok_or("image is empty or constant")?;
    let mask = apply_threshold(&smooth, threshold);
    println!("otsu threshold: {threshold:.1}");

    // 4. Label 8-connected foreground components, then drop small and border objects.
    let labels = label_2d(&mask, Connectivity2::Eight)?;
    println!("components before filtering: {}", labels.count);
    let kept = filter_objects(
        &labels.labels,
        &ObjectFilter { min_area: 10, max_area: None, remove_border: true },
    );
    println!("objects kept: {}", kept.count);

    // 5. Measure shape for every object, and intensity from the raw image.
    let shapes = shape_features_2d(&kept.labels)?;
    let stats = intensity_features(&kept.labels, &intensity)?;

    // 6. Write the label image as 16-bit TIFF, and one CSV row per object.
    let label_image: Image<u16> = kept.labels.map(|&l| u16::try_from(l).unwrap_or(u16::MAX));
    save_tiff(&dir.join("labels.tif"), &label_image, None, false)?;

    let mut csv = String::from("label,area,perimeter,eccentricity,solidity,mean,std,integrated\n");
    for (s, i) in shapes.iter().zip(&stats) {
        writeln!(
            csv,
            "{},{},{:.3},{:.4},{:.4},{:.2},{:.2},{:.1}",
            s.label, s.area, s.perimeter, s.eccentricity, s.solidity, i.mean, i.std, i.integrated
        )?;
    }
    fs::write(dir.join("objects.csv"), csv)?;
    println!("wrote {}", dir.join("objects.csv").display());
    Ok(())
}
