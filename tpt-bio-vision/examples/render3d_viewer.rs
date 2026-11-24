//! CPU reference render of a multi-channel OME-TIFF volume.
//!
//! ```text
//! cargo run --release -p tpt-bio-vision --example make_sample_volume -- sample-volume.ome.tif
//! cargo run --release -p tpt-bio-vision --example render3d_viewer -- sample-volume.ome.tif viewer-out
//! ```
//!
//! Writes three binary PPM images (`P6`, readable by most image viewers and converters):
//!
//! - `volume.ppm`: the channels composited by ray casting.
//! - `labels.ppm`: the same render with the segmented nuclei overlaid in colour.
//! - `mesh.ppm`: the nuclei surface as a flat-shaded triangle mesh.
//!
//! Nothing is opened in a window. The renderer runs on the CPU and the output is a
//! still image, so this is a reference for the library's `render3d` module, not a viewer.

use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use tpt_bio_vision::image::Image;
use tpt_bio_vision::io::ome::OmeTiff;
use tpt_bio_vision::label::{Connectivity3, label_3d};
use tpt_bio_vision::render3d::Rgba8;
use tpt_bio_vision::render3d::camera::{Camera, Projection};
use tpt_bio_vision::render3d::isosurface::isosurface;
use tpt_bio_vision::render3d::overlay::{LabelOverlay, render_labelled};
use tpt_bio_vision::render3d::rasterise::rasterise_mesh;
use tpt_bio_vision::render3d::raycast::{RenderSettings, render_volume};
use tpt_bio_vision::render3d::transfer::{ChannelTransfer, TransferFunction};
use tpt_bio_vision::render3d::volume::Volume;
use tpt_bio_vision::source::FileSource;

type Failure = Box<dyn Error>;

/// Reads every z slice of every channel of the first series into one 3D image per channel.
fn read_channels(path: &Path) -> Result<Vec<Image<u16>>, Failure> {
    let ome = OmeTiff::open(FileSource::open(path)?)?;
    let pixels = &ome.images().first().ok_or("no image series in file")?.pixels;
    let (w, h, depth) = (pixels.size_x, pixels.size_y, pixels.size_z);
    let mut channels = Vec::with_capacity(pixels.size_c);
    for c in 0..pixels.size_c {
        let mut data = Vec::with_capacity(w * h * depth);
        for z in 0..depth {
            data.extend_from_slice(ome.read_plane::<u16>(0, z, c, 0)?.data());
        }
        channels.push(Image::new_3d(w, h, depth, data)?);
    }
    Ok(channels)
}

/// Writes an RGBA image as binary PPM, dropping the alpha channel.
fn write_ppm(path: &Path, img: &Image<Rgba8>) -> Result<(), Failure> {
    let mut out = BufWriter::new(File::create(path)?);
    write!(out, "P6\n{} {}\n255\n", img.width(), img.height())?;
    for px in img.data() {
        out.write_all(&px[..3])?;
    }
    out.flush()?;
    Ok(())
}

fn main() -> Result<(), Failure> {
    let mut args = std::env::args_os().skip(1);
    let input =
        PathBuf::from(args.next().ok_or("usage: render3d_viewer <volume.ome.tif> [out-dir]")?);
    let out_dir = args.next().map_or_else(|| PathBuf::from("viewer-out"), PathBuf::from);
    fs::create_dir_all(&out_dir)?;

    // 1. Load the channels and build a volume with anisotropic spacing (z is 2 units).
    let channels = read_channels(&input)?;
    let volume = Volume::from_u16(&channels, [1.0, 1.0, 2.0])?;
    println!("volume: {:?} voxels, {} channels", volume.shape(), volume.channel_count());

    // 2. One colour and window per channel. The window maps intensity to opacity.
    let colours = [[0.25, 0.45, 1.0], [0.3, 1.0, 0.35], [1.0, 0.35, 0.25]];
    let mut transfers = Vec::new();
    for (i, colour) in colours.iter().enumerate().take(volume.channel_count()) {
        let (lo, hi) = volume.intensity_range(i).ok_or("channel has no samples")?;
        transfers.push(ChannelTransfer::new(*colour, (lo, hi)));
    }
    let tf = TransferFunction { channels: transfers, clip: None };

    // 3. Frame the whole volume from an oblique angle.
    let extent = volume.extent();
    let centre = [extent[0] / 2.0, extent[1] / 2.0, extent[2] / 2.0];
    let camera = Camera {
        target: centre,
        yaw: 0.6,
        pitch: 0.4,
        distance: 2.6 * extent.iter().copied().fold(0.0, f64::max),
        projection: Projection::Perspective { fov_y: 0.8 },
        viewport: (320, 240),
    };
    let settings = RenderSettings::new(1.0);

    // 4. Plain volume render.
    let image = render_volume(&volume, &camera, &tf, &settings)?;
    write_ppm(&out_dir.join("volume.ppm"), &image)?;

    // 5. Segment the nuclei channel by thresholding and 3D labelling, then overlay the labels.
    let nuclei = volume.channel(0).ok_or("volume has no channel 0")?;
    let mask = nuclei.map(|&v| v >= 1500.0);
    let labels = label_3d(&mask, Connectivity3::Six);
    println!("nuclei objects: {}", labels.count);
    let overlay = LabelOverlay { opacity: 0.6 };
    let labelled = render_labelled(&volume, &labels.labels, &camera, &tf, &settings, &overlay)?;
    write_ppm(&out_dir.join("labels.ppm"), &labelled)?;

    // 6. Extract the nuclei surface at the threshold and rasterise it.
    let field: Image<f32> = nuclei.clone();
    let mesh = isosurface(&field, volume.spacing(), 1500.0)?;
    println!("nuclei mesh: {} vertices, {} triangles", mesh.vertices.len(), mesh.triangles.len());
    let surface = rasterise_mesh(&mesh, &camera, [200, 205, 220], [0, 0, 0])?;
    write_ppm(&out_dir.join("mesh.ppm"), &surface)?;

    println!("wrote volume.ppm, labels.ppm and mesh.ppm to {}", out_dir.display());
    Ok(())
}
