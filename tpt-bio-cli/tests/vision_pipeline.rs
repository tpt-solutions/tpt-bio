//! End-to-end tests for `vision segment`, `vision features`, `pipeline run` and `view`,
//! run against the built `tptbio` binary. Inputs are generated in a scratch directory
//! that is removed afterwards.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use tpt_bio_vision::image::Image;
use tpt_bio_vision::io::Planes;
use tpt_bio_vision::io::tiff::read_tiff;
use tpt_bio_vision::io::write::{save_png, save_tiff};
use tpt_bio_vision::source::FileSource;

const EXIT_OK: i32 = 0;
const EXIT_RUNTIME: i32 = 1;
const EXIT_USAGE: i32 = 2;

/// Background intensity and object intensity of the generated scenes.
const BACKGROUND: u8 = 10;
const OBJECT: u8 = 200;

/// A unique scratch directory under the system temp dir, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Self {
        let nanos = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
        let dir = std::env::temp_dir()
            .join(format!("tptbio-vision-{tag}-{}-{nanos}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    /// Writes a 40 x 20 PNG: background 10, with filled rectangles `(x0, x1, y0, y1)`
    /// (half-open) at 200.
    fn scene_png(&self, name: &str, rects: &[(usize, usize, usize, usize)]) -> PathBuf {
        let path = self.path(name);
        save_png(&path, &scene(rects)).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn scene(rects: &[(usize, usize, usize, usize)]) -> Image<u8> {
    let (w, h) = (40, 20);
    let mut data = vec![BACKGROUND; w * h];
    for &(x0, x1, y0, y1) in rects {
        for y in y0..y1 {
            for x in x0..x1 {
                data[y * w + x] = OBJECT;
            }
        }
    }
    Image::new_2d(w, h, data).unwrap()
}

/// Two 6 x 6 squares, 4-connectivity apart.
const TWO_SQUARES: &[(usize, usize, usize, usize)] = &[(4, 10, 7, 13), (24, 30, 7, 13)];

fn tptbio(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_tptbio")).args(args).output().unwrap()
}

fn code(out: &Output) -> i32 {
    out.status.code().unwrap_or(-1)
}

fn stdout(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn s(path: &Path) -> &str {
    path.to_str().unwrap()
}

/// Runs `vision segment` and asserts success.
fn segment(input: &Path, output: &Path, extra: &[&str]) -> String {
    let mut args = vec!["vision", "segment", s(input), "-o", s(output)];
    args.extend_from_slice(extra);
    let out = tptbio(&args);
    assert_eq!(code(&out), EXIT_OK, "stderr: {}", stderr(&out));
    stdout(&out)
}

/// Largest label and number of labelled pixels in a 16-bit label TIFF.
fn label_summary(path: &Path) -> (u16, usize) {
    let source = FileSource::open(path).unwrap();
    let Planes::U16(channels) = read_tiff(&source).unwrap() else {
        panic!("label image is not 16-bit");
    };
    let data = channels[0].data();
    (data.iter().copied().max().unwrap_or(0), data.iter().filter(|&&v| v != 0).count())
}

#[test]
fn segment_counts_two_squares() {
    let dir = Scratch::new("segment-two");
    let input = dir.scene_png("in.png", TWO_SQUARES);
    let labels = dir.path("labels.tif");
    let line = segment(&input, &labels, &[]);
    assert!(line.contains("2 objects found, 2 kept"), "stdout: {line}");
    assert_eq!(label_summary(&labels), (2, 72));
}

#[test]
fn segment_min_size_drops_small_objects() {
    let dir = Scratch::new("segment-min");
    // A single-pixel dot at (17, 3) is object 3 and falls below --min-size 4.
    let mut rects = TWO_SQUARES.to_vec();
    rects.push((17, 18, 3, 4));
    let input = dir.scene_png("in.png", &rects);
    let labels = dir.path("labels.tif");
    let line = segment(&input, &labels, &["--min-size", "4"]);
    assert!(line.contains("3 objects found, 2 kept"), "stdout: {line}");
    assert_eq!(label_summary(&labels), (2, 72));
}

#[test]
fn segment_remove_border_drops_edge_objects() {
    let dir = Scratch::new("segment-border");
    // The square at the top-left corner touches the image border.
    let mut rects = TWO_SQUARES.to_vec();
    rects.push((0, 6, 0, 6));
    let input = dir.scene_png("in.png", &rects);
    let labels = dir.path("labels.tif");
    let line = segment(&input, &labels, &["--remove-border"]);
    assert!(line.contains("3 objects found, 2 kept"), "stdout: {line}");
}

#[test]
fn segment_eight_connectivity_joins_diagonal_neighbours() {
    let dir = Scratch::new("segment-conn");
    // The two squares meet only at a corner.
    let input = dir.scene_png("in.png", &[(4, 8, 7, 11), (8, 12, 11, 15)]);
    let labels = dir.path("labels.tif");
    let four = segment(&input, &labels, &[]);
    assert!(four.contains("2 objects found, 2 kept"), "stdout: {four}");
    let eight = segment(&input, &labels, &["--connectivity", "8"]);
    assert!(eight.contains("1 objects found, 1 kept"), "stdout: {eight}");
}

#[test]
fn features_write_shape_and_intensity_columns() {
    let dir = Scratch::new("features-csv");
    let input = dir.scene_png("in.png", TWO_SQUARES);
    let labels = dir.path("labels.tif");
    segment(&input, &labels, &[]);
    let out_csv = dir.path("features.csv");
    let out = tptbio(&[
        "vision",
        "features",
        "--labels",
        s(&labels),
        "--intensity",
        s(&input),
        "-o",
        s(&out_csv),
    ]);
    assert_eq!(code(&out), EXIT_OK, "stderr: {}", stderr(&out));
    let text = fs::read_to_string(&out_csv).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(
        lines[0],
        "label,area,perimeter,eccentricity,solidity,intensity_mean,intensity_std,\
         intensity_min,intensity_max,intensity_median,intensity_integrated"
    );
    assert_eq!(lines.len(), 3, "one row per object: {text}");
    // A 6 x 6 square: area 36, boundary 24, mean 200, integrated 36 * 200.
    assert!(lines[1].starts_with("1,36,24,"), "row: {}", lines[1]);
    assert!(lines[1].contains(",200,0,200,200,200,7200"), "row: {}", lines[1]);
}

#[test]
fn features_without_intensity_write_shape_columns_only() {
    let dir = Scratch::new("features-shape");
    let input = dir.scene_png("in.png", TWO_SQUARES);
    let labels = dir.path("labels.tif");
    segment(&input, &labels, &[]);
    let out_csv = dir.path("shape.csv");
    let out = tptbio(&["vision", "features", "--labels", s(&labels), "-o", s(&out_csv)]);
    assert_eq!(code(&out), EXIT_OK, "stderr: {}", stderr(&out));
    let text = fs::read_to_string(&out_csv).unwrap();
    assert_eq!(text.lines().next(), Some("label,area,perimeter,eccentricity,solidity"));
    assert_eq!(text.lines().count(), 3);
}

#[test]
fn features_follow_output_format_jsonl() {
    let dir = Scratch::new("features-jsonl");
    let input = dir.scene_png("in.png", TWO_SQUARES);
    let labels = dir.path("labels.tif");
    segment(&input, &labels, &[]);
    let out_file = dir.path("shape.jsonl");
    let out = tptbio(&[
        "--output-format",
        "jsonl",
        "vision",
        "features",
        "--labels",
        s(&labels),
        "-o",
        s(&out_file),
    ]);
    assert_eq!(code(&out), EXIT_OK, "stderr: {}", stderr(&out));
    let text = fs::read_to_string(&out_file).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].contains("\"area\":36"), "line: {}", lines[0]);
}

#[test]
fn view_summarizes_png() {
    let dir = Scratch::new("view-png");
    let input = dir.scene_png("in.png", TWO_SQUARES);
    let out = tptbio(&["view", s(&input)]);
    assert_eq!(code(&out), EXIT_OK, "stderr: {}", stderr(&out));
    let text = stdout(&out);
    assert!(text.contains("size: 40 x 20 x 1"), "stdout: {text}");
    assert!(text.contains("pixel type: u8"), "stdout: {text}");
    assert!(text.contains("channels: 1"), "stdout: {text}");
    // 72 object pixels at 200 and 728 background pixels at 10: mean 27.1.
    assert!(text.contains("channel 0: min 10, max 200, mean 27.100"), "stdout: {text}");
}

#[test]
fn view_summarizes_tiff() {
    let dir = Scratch::new("view-tiff");
    let path = dir.path("in.tif");
    let img: Image<u16> = scene(TWO_SQUARES).map(|&v| u16::from(v) * 10);
    save_tiff(&path, &img, None, false).unwrap();
    let out = tptbio(&["view", s(&path)]);
    assert_eq!(code(&out), EXIT_OK, "stderr: {}", stderr(&out));
    let text = stdout(&out);
    assert!(text.contains("pixel type: u16"), "stdout: {text}");
    assert!(text.contains("channel 0: min 100, max 2000"), "stdout: {text}");
}

#[test]
fn pipeline_run_prints_tables_and_label_summary() {
    let dir = Scratch::new("pipeline-ok");
    let def = dir.path("pipeline.toml");
    fs::write(
        &def,
        r#"
[[nodes]]
name = "data"
type = "table_source"
params = { values = [[0.0, 0.0], [0.0, 1.0], [10.0, 10.0], [10.0, 11.0]] }

[[nodes]]
name = "scaled"
type = "normalize"
params = { method = "zscore" }
inputs = ["data"]

[[nodes]]
name = "reduced"
type = "pca"
params = { components = 1 }
inputs = ["scaled"]

[[nodes]]
name = "clusters"
type = "kmeans"
params = { k = 2, seed = 1 }
inputs = ["reduced"]
"#,
    )
    .unwrap();
    let out = tptbio(&["pipeline", "run", s(&def)]);
    assert_eq!(code(&out), EXIT_OK, "stderr: {}", stderr(&out));
    let text = stdout(&out);
    assert!(text.contains("c1,c2"), "stdout: {text}");
    assert!(text.contains("labels: 4 rows, 2 clusters"), "stdout: {text}");
}

#[test]
fn pipeline_definition_error_names_the_node() {
    let dir = Scratch::new("pipeline-bad");
    let def = dir.path("bad.toml");
    fs::write(&def, "[[nodes]]\nname = \"oops\"\ntype = \"no_such_kind\"\n").unwrap();
    let out = tptbio(&["pipeline", "run", s(&def)]);
    assert_eq!(code(&out), EXIT_USAGE);
    assert!(stderr(&out).contains("oops"), "stderr: {}", stderr(&out));
}

#[test]
fn exit_codes_for_usage_and_runtime_errors() {
    let dir = Scratch::new("exit-codes");
    let input = dir.scene_png("in.png", TWO_SQUARES);
    let labels = dir.path("labels.tif");
    segment(&input, &labels, &[]);
    let missing = dir.path("missing.png");
    let jpg = dir.path("photo.jpg");
    fs::write(&jpg, b"not an image").unwrap();
    let out_csv = dir.path("out.csv");
    let txt = dir.path("notes.txt");
    fs::write(&txt, b"text").unwrap();

    let cases: Vec<(&str, Vec<&str>, i32)> = vec![
        // Usage errors: exit 2.
        ("segment without -o", vec!["vision", "segment", s(&input)], EXIT_USAGE),
        (
            "segment with unsupported extension",
            vec!["vision", "segment", s(&jpg), "-o", s(&labels)],
            EXIT_USAGE,
        ),
        (
            "features without --labels",
            vec!["vision", "features", s(&input), "-o", s(&out_csv)],
            EXIT_USAGE,
        ),
        (
            "features with both input and --intensity",
            vec![
                "vision",
                "features",
                s(&input),
                "--intensity",
                s(&input),
                "--labels",
                s(&labels),
                "-o",
                s(&out_csv),
            ],
            EXIT_USAGE,
        ),
        ("vision without subcommand", vec!["vision"], EXIT_USAGE),
        ("view on an unsupported file", vec!["view", s(&txt)], EXIT_USAGE),
        ("pipeline run without a file", vec!["pipeline", "run"], EXIT_USAGE),
        // Runtime errors: exit 1.
        (
            "segment with a missing channel",
            vec!["vision", "segment", s(&input), "-o", s(&labels), "--channel", "5"],
            EXIT_RUNTIME,
        ),
        (
            "segment with a missing input file",
            vec!["vision", "segment", s(&missing), "-o", s(&labels)],
            EXIT_RUNTIME,
        ),
        (
            "features with a missing label image",
            vec!["vision", "features", "--labels", s(&missing), "-o", s(&out_csv)],
            EXIT_RUNTIME,
        ),
        ("view on a missing file", vec!["view", s(&missing)], EXIT_RUNTIME),
        ("pipeline run on a missing file", vec!["pipeline", "run", s(&missing)], EXIT_RUNTIME),
    ];
    for (name, args, expected) in cases {
        let out = tptbio(&args);
        assert_eq!(code(&out), expected, "{name}: stderr: {}", stderr(&out));
    }
}
