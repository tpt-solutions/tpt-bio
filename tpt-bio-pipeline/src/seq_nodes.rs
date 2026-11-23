//! Built-in sequence nodes: FASTQ and BAM readers, record counts and CSV writers.
//!
//! Readers return whole files in memory, which matches the in-memory executor.
//! Sinks write their file and pass its path on as a [`Value::File`].

use std::fmt::Write as _;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use tpt_bio_core::stats::Matrix;
use tpt_bio_seq::fastq::{self, FastqBuf, FastqReader};
use tpt_bio_seq::{BamReader, records_to_csv};

use crate::definition::{Node, str_param};
use crate::execute::{Frame, Records, Value};

/// Evaluates `node` if it is one of this module's kinds. Returns `None` otherwise.
pub(crate) fn evaluate(node: &Node, inputs: &[&Value]) -> Option<Result<Value, String>> {
    let result = match node.kind.as_str() {
        "read_fastq" => path(node).and_then(|p| read_fastq(&p)),
        "read_bam" => path(node).and_then(|p| read_bam(&p)),
        "write_records_csv" => {
            records_input(inputs).and_then(|r| path(node).and_then(|p| write_records_csv(r, &p)))
        }
        "write_frame_csv" => {
            frame_input(inputs).and_then(|f| path(node).and_then(|p| write_frame_csv(f, &p)))
        }
        "write_labels_csv" => {
            labels_input(inputs).and_then(|l| path(node).and_then(|p| write_labels_csv(l, &p)))
        }
        "count_reads" => fastq_input(inputs).and_then(|r| count(r.len())),
        "count_records" => records_input(inputs).and_then(|r| count(r.records.len())),
        "count_objects" => objects_input(inputs).and_then(|o| count(o.count as usize)),
        _ => return None,
    };
    Some(result)
}

fn path(node: &Node) -> Result<String, String> {
    match str_param(&node.params, "path")? {
        Some(p) if !p.is_empty() => Ok(p.to_string()),
        _ => Err("missing required parameter 'path'".to_string()),
    }
}

fn read_fastq(path: &str) -> Result<Value, String> {
    let read = |path: &str| -> Result<Vec<FastqBuf>, tpt_bio_seq::SeqError> {
        let src = fastq::open_reader(Path::new(path))?;
        let mut reader = FastqReader::new(src);
        let mut out = Vec::new();
        loop {
            let mut buf = FastqBuf::default();
            if !reader.read_record(&mut buf)? {
                break;
            }
            out.push(buf);
        }
        Ok(out)
    };
    read(path).map(Value::Fastq).map_err(|e| format!("cannot read '{path}': {e}"))
}

fn read_bam(path: &str) -> Result<Value, String> {
    let read = |path: &str| -> Result<Records, tpt_bio_seq::SeqError> {
        let mut reader = BamReader::new(BufReader::new(File::open(path)?))?;
        let header = reader.header().clone();
        let mut records = Vec::new();
        while let Some(rec) = reader.next_record()? {
            records.push(rec);
        }
        Ok(Records { header, records })
    };
    read(path).map(Value::Records).map_err(|e| format!("cannot read '{path}': {e}"))
}

fn write_records_csv(records: &Records, path: &str) -> Result<Value, String> {
    let text = records_to_csv(&records.records, &records.header).map_err(|e| e.to_string())?;
    write_text(path, &text)
}

fn write_frame_csv(frame: &Frame, path: &str) -> Result<Value, String> {
    let mut text = frame.columns.join(",");
    text.push('\n');
    for row in &frame.rows {
        let cells: Vec<String> = row.iter().map(ToString::to_string).collect();
        text.push_str(&cells.join(","));
        text.push('\n');
    }
    write_text(path, &text)
}

fn write_labels_csv(labels: &[usize], path: &str) -> Result<Value, String> {
    let mut text = String::from("row,label\n");
    for (row, label) in labels.iter().enumerate() {
        writeln!(text, "{row},{label}").expect("writing to a String cannot fail");
    }
    write_text(path, &text)
}

fn write_text(path: &str, text: &str) -> Result<Value, String> {
    std::fs::write(path, text).map_err(|e| format!("cannot write '{path}': {e}"))?;
    Ok(Value::File(path.to_string()))
}

fn count(n: usize) -> Result<Value, String> {
    Matrix::new(1, 1, vec![n as f64]).map(Value::Table).map_err(|e| e.to_string())
}

fn fastq_input<'a>(inputs: &[&'a Value]) -> Result<&'a [FastqBuf], String> {
    match inputs.first() {
        Some(Value::Fastq(r)) => Ok(r),
        Some(other) => Err(format!("expected fastq input, got {}", other.type_tag())),
        None => Err("missing fastq input".to_string()),
    }
}

fn records_input<'a>(inputs: &[&'a Value]) -> Result<&'a Records, String> {
    match inputs.first() {
        Some(Value::Records(r)) => Ok(r),
        Some(other) => Err(format!("expected records input, got {}", other.type_tag())),
        None => Err("missing records input".to_string()),
    }
}

fn frame_input<'a>(inputs: &[&'a Value]) -> Result<&'a Frame, String> {
    match inputs.first() {
        Some(Value::Frame(f)) => Ok(f),
        Some(other) => Err(format!("expected frame input, got {}", other.type_tag())),
        None => Err("missing frame input".to_string()),
    }
}

fn labels_input<'a>(inputs: &[&'a Value]) -> Result<&'a [usize], String> {
    match inputs.first() {
        Some(Value::Labels(l)) => Ok(l),
        Some(other) => Err(format!("expected labels input, got {}", other.type_tag())),
        None => Err("missing labels input".to_string()),
    }
}

fn objects_input<'a>(inputs: &[&'a Value]) -> Result<&'a tpt_bio_vision::label::Labeling, String> {
    match inputs.first() {
        Some(Value::Objects(o)) => Ok(o),
        Some(other) => Err(format!("expected objects input, got {}", other.type_tag())),
        None => Err("missing objects input".to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_is_a_one_by_one_table() {
        let Value::Table(m) = count(7).unwrap() else { panic!("expected table") };
        assert_eq!((m.rows(), m.cols()), (1, 1));
        assert!((m.get(0, 0) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn unrelated_kinds_are_not_claimed() {
        let node = Node {
            name: "x".into(),
            kind: "pca".into(),
            params: toml::Table::new(),
            inputs: Vec::new(),
        };
        assert!(evaluate(&node, &[]).is_none());
    }
}
