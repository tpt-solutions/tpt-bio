//! Runs the spec section 7 example definition on a generated plate image.

mod common;

use std::fs;
use std::path::Path;

use tpt_bio_pipeline::{Pipeline, Value};

use common::{scratch, squares_image, write_png};

/// Four small dim cells and two large bright cells, well separated.
const CELLS: [(usize, usize, usize, u8); 6] = [
    (10, 10, 6, 120),
    (40, 10, 6, 120),
    (70, 10, 6, 120),
    (10, 60, 6, 120),
    (40, 55, 10, 240),
    (70, 55, 10, 240),
];

#[test]
fn spec7_example_runs_on_generated_plate() {
    let dir = scratch("spec7");
    write_png(&dir.join("plate.png"), &squares_image(96, 96, &CELLS));

    let example = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/spec7_screen.toml");
    let text = fs::read_to_string(&example).expect("read example definition");
    // The example writes its inputs and outputs under 'screen/'. Point them at the scratch directory.
    let text = text.replace("'screen/", &format!("'{}/", dir.display()));

    let pipeline = Pipeline::from_toml(&text).expect("example definition is valid");
    let out = pipeline.run().expect("example pipeline runs");
    assert_eq!(out.len(), 14);

    // Step 4: one feature row per cell, every cell of the plate found.
    let features = fs::read_to_string(dir.join("features.csv")).expect("feature CSV written");
    let rows: Vec<&str> = features.lines().skip(1).collect();
    assert_eq!(rows.len(), CELLS.len(), "{features}");

    // Step 5: the phenotype output has one label per cell, in two groups.
    let Value::Labels(labels) = &out[11] else { panic!("expected labels, got {:?}", out[11]) };
    assert_eq!(labels.len(), CELLS.len());
    let mut distinct = labels.clone();
    distinct.sort_unstable();
    distinct.dedup();
    assert_eq!(distinct.len(), 2, "{labels:?}");

    // Step 6: the label export mirrors the labels.
    let phenotypes = fs::read_to_string(dir.join("phenotypes.csv")).expect("label CSV written");
    assert_eq!(phenotypes.lines().count(), CELLS.len() + 1);
}
