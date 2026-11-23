//! Helpers shared by the integration tests: scratch directories and generated inputs.

#![allow(dead_code)]

use std::fs;
use std::path::{Path, PathBuf};

use tpt_bio_vision::image::Image;
use tpt_bio_vision::io::write::save_png;

/// A fresh, empty scratch directory unique to this test process and `name`.
pub(crate) fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("tpt-pipeline-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("create scratch directory");
    dir
}

/// A 64 by 64 black image with filled squares. Each entry is `(x, y, side, value)`.
pub(crate) fn squares_image(
    width: usize,
    height: usize,
    squares: &[(usize, usize, usize, u8)],
) -> Image<u8> {
    let mut data = vec![0u8; width * height];
    for &(x0, y0, side, value) in squares {
        for y in y0..y0 + side {
            for x in x0..x0 + side {
                data[y * width + x] = value;
            }
        }
    }
    Image::new_2d(width, height, data).expect("shape matches data")
}

/// Writes `image` as a grey PNG at `path`.
pub(crate) fn write_png(path: &Path, image: &Image<u8>) {
    save_png(path, image).expect("write PNG");
}

/// Quotes a path as a TOML literal string. Literal strings take backslashes verbatim.
pub(crate) fn toml_path(path: &Path) -> String {
    let text = path.display().to_string();
    assert!(!text.contains('\''), "path may not contain a single quote: {text}");
    format!("'{text}'")
}
