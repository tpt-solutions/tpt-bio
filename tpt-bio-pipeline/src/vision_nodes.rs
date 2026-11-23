//! Built-in vision nodes: image loading, smoothing, Otsu masks, labelling, seeds,
//! watershed, object filtering and per-object measurement.
//!
//! Every node works on 2D images. Images are held as `f64`, so loaded 8- and 16-bit
//! samples are converted once at the load step.

use std::path::Path;

use tpt_bio_vision::features::{intensity_features, shape_features_2d};
use tpt_bio_vision::image::Image;
use tpt_bio_vision::io::Planes;
use tpt_bio_vision::io::{png, tiff};
use tpt_bio_vision::label::{Connectivity2, Labeling, label_2d};
use tpt_bio_vision::objects::{ObjectFilter, filter_objects};
use tpt_bio_vision::seeds::{log_seeds, seeds_to_markers};
use tpt_bio_vision::source::FileSource;
use tpt_bio_vision::{spatial, threshold, watershed};

use crate::definition::{Node, bool_param, float_param, int_param, str_param};
use crate::execute::{Frame, Value};

/// Column names of the frame written by `measure`, in order.
const MEASURE_COLUMNS: [&str; 13] = [
    "label",
    "area",
    "perimeter",
    "eccentricity",
    "solidity",
    "mean",
    "std",
    "min",
    "max",
    "q25",
    "median",
    "q75",
    "integrated",
];

/// Evaluates `node` if it is one of this module's kinds. Returns `None` otherwise.
pub(crate) fn evaluate(node: &Node, inputs: &[&Value]) -> Option<Result<Value, String>> {
    let result = match node.kind.as_str() {
        "load_image" => load_image(node),
        "gaussian" => image_input(inputs, 0).and_then(|img| gaussian(node, img)),
        "subtract_background" => {
            image_input(inputs, 0).and_then(|img| subtract_background(node, img))
        }
        "otsu" => image_input(inputs, 0).and_then(otsu),
        "label" => mask_input(inputs, 0).and_then(|m| label(node, m)),
        "filter_objects" => objects_input(inputs, 0).and_then(|o| filter(node, o)),
        "log_seeds" => image_input(inputs, 0).and_then(|img| seeds(node, img)),
        "watershed" => watershed_node(inputs),
        "measure" => measure(inputs),
        _ => return None,
    };
    Some(result)
}

fn image_err(e: impl std::fmt::Display) -> String {
    e.to_string()
}

fn load_image(node: &Node) -> Result<Value, String> {
    let path = match str_param(&node.params, "path")? {
        Some(p) if !p.is_empty() => p,
        _ => return Err("missing required parameter 'path'".to_string()),
    };
    let channel = usize::try_from(int_param(&node.params, "channel", 0)?.unwrap_or(0))
        .map_err(|_| "parameter 'channel' must be non-negative".to_string())?;

    let source =
        FileSource::open(Path::new(path)).map_err(|e| format!("cannot open '{path}': {e}"))?;
    let extension = Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_ascii_lowercase)
        .unwrap_or_default();
    let planes = match extension.as_str() {
        "png" => png::read_png(&source),
        "tif" | "tiff" => tiff::read_tiff(&source),
        _ => return Err(format!("'{path}' is not a PNG or TIFF file")),
    }
    .map_err(|e| format!("cannot decode '{path}': {e}"))?;

    let img = channel_of(&planes, channel)?;
    if img.depth() != 1 {
        return Err(format!(
            "'{path}' is a volume (depth {}); only 2D images are supported",
            img.depth()
        ));
    }
    Ok(Value::Image(img))
}

fn channel_of(planes: &Planes, channel: usize) -> Result<Image<f64>, String> {
    let out_of_range = || {
        format!(
            "channel {channel} is out of range; the file has {} channel(s)",
            planes.channel_count()
        )
    };
    match planes {
        Planes::U8(v) => v.get(channel).map(|i| i.map(|&x| f64::from(x))),
        Planes::U16(v) => v.get(channel).map(|i| i.map(|&x| f64::from(x))),
        Planes::F32(v) => v.get(channel).map(|i| i.map(|&x| f64::from(x))),
    }
    .ok_or_else(out_of_range)
}

fn gaussian(node: &Node, img: &Image<f64>) -> Result<Value, String> {
    let sigma = float_param(&node.params, "sigma", 0.0)?
        .ok_or_else(|| "missing required parameter 'sigma'".to_string())?;
    spatial::gaussian_smooth(img, sigma).map(Value::Image).map_err(image_err)
}

fn subtract_background(node: &Node, img: &Image<f64>) -> Result<Value, String> {
    let radius = usize::try_from(
        int_param(&node.params, "radius", 0)?
            .ok_or_else(|| "missing required parameter 'radius'".to_string())?,
    )
    .map_err(|_| "parameter 'radius' must be non-negative".to_string())?;
    spatial::subtract_background(img, radius).map(Value::Image).map_err(image_err)
}

fn otsu(img: &Image<f64>) -> Result<Value, String> {
    // Otsu is undefined for an empty or constant image; the mask is then all background.
    let mask = if let Some(t) = threshold::otsu_threshold(img, 256).map_err(image_err)? {
        threshold::apply_threshold(img, t)
    } else {
        let (w, h, d) = img.shape();
        Image::new_3d(w, h, d, vec![false; w * h * d]).map_err(image_err)?
    };
    Ok(Value::Mask(mask))
}

fn label(node: &Node, mask: &Image<bool>) -> Result<Value, String> {
    let conn = match int_param(&node.params, "connectivity", 0)?.unwrap_or(8) {
        4 => Connectivity2::Four,
        _ => Connectivity2::Eight,
    };
    label_2d(mask, conn).map(Value::Objects).map_err(image_err)
}

fn filter(node: &Node, objects: &Labeling) -> Result<Value, String> {
    let min_area = int_param(&node.params, "min_area", 0)?.unwrap_or(0);
    let max_area = int_param(&node.params, "max_area", 0)?;
    let filter = ObjectFilter {
        min_area: usize::try_from(min_area)
            .map_err(|_| "min_area must be non-negative".to_string())?,
        max_area: max_area
            .map(|m| usize::try_from(m).map_err(|_| "max_area must be non-negative".to_string()))
            .transpose()?,
        remove_border: bool_param(&node.params, "remove_border")?.unwrap_or(false),
    };
    Ok(Value::Objects(filter_objects(&objects.labels, &filter)))
}

fn seeds(node: &Node, img: &Image<f64>) -> Result<Value, String> {
    let sigma = float_param(&node.params, "sigma", 0.0)?.unwrap_or(2.0);
    if sigma <= 0.0 {
        return Err("parameter 'sigma' must be greater than 0".to_string());
    }
    let min_distance = usize::try_from(int_param(&node.params, "min_distance", 0)?.unwrap_or(3))
        .map_err(|_| "min_distance must be non-negative".to_string())?;
    let min_response = float_param(&node.params, "min_response", f64::NEG_INFINITY)?.unwrap_or(0.0);
    let found = log_seeds(img, sigma, min_distance, min_response).map_err(image_err)?;
    let (w, h, d) = img.shape();
    let labels = seeds_to_markers(w, h, d, &found).map_err(image_err)?;
    let count =
        u32::try_from(found.len()).map_err(|_| "too many seeds for u32 labels".to_string())?;
    Ok(Value::Objects(Labeling { labels, count }))
}

fn watershed_node(inputs: &[&Value]) -> Result<Value, String> {
    let img = image_input(inputs, 0)?;
    let markers = objects_input(inputs, 1)?;
    let mask = mask_input(inputs, 2)?;
    let labels = watershed::watershed(img, &markers.labels, Some(mask)).map_err(image_err)?;
    Ok(Value::Objects(Labeling { labels, count: markers.count }))
}

fn measure(inputs: &[&Value]) -> Result<Value, String> {
    let objects = objects_input(inputs, 0)?;
    let img = image_input(inputs, 1)?;
    let shape = shape_features_2d(&objects.labels).map_err(image_err)?;
    let intensity = intensity_features(&objects.labels, img).map_err(image_err)?;
    if shape.len() != intensity.len() {
        return Err("shape and intensity measurements disagree on the object count".to_string());
    }
    let rows = shape
        .iter()
        .zip(&intensity)
        .map(|(s, i)| {
            vec![
                f64::from(s.label),
                s.area as f64,
                s.perimeter,
                s.eccentricity,
                s.solidity,
                i.mean,
                i.std,
                i.min,
                i.max,
                i.q25,
                i.median,
                i.q75,
                i.integrated,
            ]
        })
        .collect();
    Ok(Value::Frame(Frame {
        columns: MEASURE_COLUMNS.iter().map(ToString::to_string).collect(),
        rows,
    }))
}

fn image_input<'a>(inputs: &[&'a Value], slot: usize) -> Result<&'a Image<f64>, String> {
    match inputs.get(slot) {
        Some(Value::Image(img)) => Ok(img),
        Some(other) => Err(format!("expected an image input, got {}", other.type_tag())),
        None => Err("missing image input".to_string()),
    }
}

fn mask_input<'a>(inputs: &[&'a Value], slot: usize) -> Result<&'a Image<bool>, String> {
    match inputs.get(slot) {
        Some(Value::Mask(m)) => Ok(m),
        Some(other) => Err(format!("expected a mask input, got {}", other.type_tag())),
        None => Err("missing mask input".to_string()),
    }
}

fn objects_input<'a>(inputs: &[&'a Value], slot: usize) -> Result<&'a Labeling, String> {
    match inputs.get(slot) {
        Some(Value::Objects(o)) => Ok(o),
        Some(other) => Err(format!("expected an objects input, got {}", other.type_tag())),
        None => Err("missing objects input".to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn channel_out_of_range_is_reported() {
        let planes = Planes::U8(vec![Image::new_2d(1, 1, vec![9u8]).unwrap()]);
        let img = channel_of(&planes, 0).unwrap();
        assert_eq!(img.data(), &[9.0]);
        let err = channel_of(&planes, 3).unwrap_err();
        assert!(err.contains("channel 3"), "{err}");
    }

    #[test]
    fn otsu_of_constant_image_is_all_background() {
        let img = Image::new_2d(2, 2, vec![5.0; 4]).unwrap();
        let Value::Mask(m) = otsu(&img).unwrap() else { panic!("expected mask") };
        assert!(m.data().iter().all(|&b| !b));
    }
}
