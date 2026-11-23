//! Integration tests for the built-in sequence and vision nodes.

mod common;

use std::fmt::Write as _;
use std::fs;

use tpt_bio_pipeline::{Pipeline, Value};
use tpt_bio_seq::fastq::FastqBuf;

use common::{scratch, squares_image, toml_path, write_png};

/// Six bright squares on black: five interior 6x6 objects, one 6x6 object touching
/// the left border, and a single-pixel dot.
fn plate_squares() -> Vec<(usize, usize, usize, u8)> {
    vec![
        (8, 8, 6, 200),
        (24, 8, 6, 200),
        (40, 8, 6, 200),
        (8, 40, 6, 200),
        (40, 40, 6, 200),
        (0, 28, 6, 200),
        (30, 30, 1, 200),
    ]
}

fn table_value(v: &Value) -> i64 {
    let Value::Table(m) = v else { panic!("expected a table, got {v:?}") };
    assert_eq!((m.rows(), m.cols()), (1, 1));
    m.get(0, 0).round() as i64
}

#[test]
fn threshold_label_filter_counts_objects_in_generated_png() {
    let dir = scratch("vision-count");
    let png = dir.join("plate.png");
    write_png(&png, &squares_image(64, 64, &plate_squares()));

    let text = format!(
        r#"
[[nodes]]
name = "raw"
type = "load_image"
params = {{ path = {} }}

[[nodes]]
name = "mask"
type = "otsu"
inputs = ["raw"]

[[nodes]]
name = "labels"
type = "label"
params = {{ connectivity = 8 }}
inputs = ["mask"]

[[nodes]]
name = "kept"
type = "filter_objects"
params = {{ min_area = 4, remove_border = true }}
inputs = ["labels"]

[[nodes]]
name = "count"
type = "count_objects"
inputs = ["kept"]
"#,
        toml_path(&png)
    );
    let out = Pipeline::from_toml(&text).unwrap().run().unwrap();
    assert_eq!(table_value(out.last().unwrap()), 5);

    // Without the filter, all seven components are counted.
    let unfiltered = text.replace("params = { min_area = 4, remove_border = true }", "params = {}");
    let out = Pipeline::from_toml(&unfiltered).unwrap().run().unwrap();
    assert_eq!(table_value(out.last().unwrap()), 7);
}

#[test]
fn fastq_pipeline_reads_and_counts_generated_file() {
    let dir = scratch("seq-fastq");
    let fastq = dir.join("reads.fastq");
    let mut text = String::new();
    for i in 0..25 {
        write!(text, "@read{i}\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    }
    fs::write(&fastq, text).unwrap();

    let definition = format!(
        r#"
[[nodes]]
name = "reads"
type = "read_fastq"
params = {{ path = {} }}

[[nodes]]
name = "count"
type = "count_reads"
inputs = ["reads"]
"#,
        toml_path(&fastq)
    );
    let out = Pipeline::from_toml(&definition).unwrap().run().unwrap();
    assert_eq!(table_value(&out[1]), 25);
    let Value::Fastq(reads) = &out[0] else { panic!("expected reads") };
    assert_eq!(reads.len(), 25);
    assert_eq!(
        reads[0],
        FastqBuf { name: b"read0".to_vec(), seq: b"ACGTACGT".to_vec(), qual: b"IIIIIIII".to_vec() }
    );
}

#[test]
fn fastq_cannot_feed_a_records_sink() {
    let definition = r#"
[[nodes]]
name = "reads"
type = "read_fastq"
params = { path = "reads.fastq" }

[[nodes]]
name = "out"
type = "write_records_csv"
params = { path = "out.csv" }
inputs = ["reads"]
"#;
    let err = Pipeline::from_toml(definition).unwrap_err();
    assert!(err.to_string().contains("'out'"), "{err}");
}

#[test]
fn bam_written_in_test_is_read_counted_and_exported() {
    use tpt_bio_seq::{BamWriter, Record, Reference};

    let dir = scratch("seq-bam");
    let header =
        tpt_bio_seq::Header::from_references(&[Reference { name: "chr1".to_string(), len: 1000 }])
            .unwrap();
    let mut writer = BamWriter::new(Vec::new(), header).unwrap();
    for i in 0..3 {
        let rec = Record {
            name: format!("r{i}").into_bytes(),
            ref_id: 0,
            pos: 10 * i64::from(i),
            mapq: 60,
            seq: b"ACGT".to_vec(),
            ..Record::default()
        };
        writer.write_record(&rec).unwrap();
    }
    let unmapped = Record { name: b"u".to_vec(), ref_id: -1, pos: -1, ..Record::default() };
    writer.write_record(&unmapped).unwrap();
    let bytes = writer.finish().unwrap();
    let bam = dir.join("reads.bam");
    fs::write(&bam, bytes).unwrap();
    let csv = dir.join("records.csv");

    let definition = format!(
        r#"
[[nodes]]
name = "records"
type = "read_bam"
params = {{ path = {} }}

[[nodes]]
name = "count"
type = "count_records"
inputs = ["records"]

[[nodes]]
name = "export"
type = "write_records_csv"
params = {{ path = {} }}
inputs = ["records"]
"#,
        toml_path(&bam),
        toml_path(&csv)
    );
    let out = Pipeline::from_toml(&definition).unwrap().run().unwrap();
    assert_eq!(table_value(&out[1]), 4);
    assert_eq!(out[2], Value::File(csv.display().to_string()));
    let text = fs::read_to_string(&csv).unwrap();
    assert!(text.starts_with("qname,"), "{text}");
    assert!(text.contains("\nr2,"), "{text}");
}

#[test]
fn measure_feeds_frame_csv_with_one_row_per_object() {
    let dir = scratch("vision-measure");
    let png = dir.join("plate.png");
    write_png(&png, &squares_image(64, 64, &plate_squares()));
    let csv = dir.join("features.csv");

    let definition = format!(
        r#"
[[nodes]]
name = "raw"
type = "load_image"
params = {{ path = {} }}

[[nodes]]
name = "mask"
type = "otsu"
inputs = ["raw"]

[[nodes]]
name = "cells"
type = "label"
inputs = ["mask"]

[[nodes]]
name = "features"
type = "measure"
inputs = ["cells", "raw"]

[[nodes]]
name = "export"
type = "write_frame_csv"
params = {{ path = {} }}
inputs = ["features"]
"#,
        toml_path(&png),
        toml_path(&csv)
    );
    let out = Pipeline::from_toml(&definition).unwrap().run().unwrap();
    assert!(matches!(out[4], Value::File(_)));
    let text = fs::read_to_string(&csv).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        "label,area,perimeter,eccentricity,solidity,mean,std,min,max,q25,median,q75,integrated"
    );
    assert_eq!(lines.count(), 7);
}
