//! Writers for row-oriented results in the layouts selected by `--output-format`.
//!
//! Rows are text cells under a header. JSON Lines writes a cell as a number when it
//! is a plain decimal literal, and as a string otherwise.

use std::io::{self, Write};

use serde_json::{Map, Value};

use crate::config::OutputFormat;

/// Writes a header and rows in the chosen format.
///
/// # Errors
/// Returns an I/O error from `out`.
pub(crate) fn write_rows(
    out: &mut dyn Write,
    format: OutputFormat,
    header: &[String],
    rows: &[Vec<String>],
) -> io::Result<()> {
    match format {
        OutputFormat::Csv => write_csv(out, header, rows),
        OutputFormat::Jsonl => write_jsonl(out, header, rows),
        OutputFormat::Table => write_aligned(out, header, rows),
    }
}

fn csv_line(out: &mut dyn Write, cells: &[String]) -> io::Result<()> {
    for (i, cell) in cells.iter().enumerate() {
        if i > 0 {
            out.write_all(b",")?;
        }
        if cell.contains([',', '"', '\n', '\r']) {
            write!(out, "\"{}\"", cell.replace('"', "\"\""))?;
        } else {
            out.write_all(cell.as_bytes())?;
        }
    }
    out.write_all(b"\n")
}

fn write_csv(out: &mut dyn Write, header: &[String], rows: &[Vec<String>]) -> io::Result<()> {
    csv_line(out, header)?;
    for row in rows {
        csv_line(out, row)?;
    }
    Ok(())
}

/// A cell as a JSON value: a number for plain decimal text, otherwise a string.
fn json_cell(cell: &str) -> Value {
    let plain = !cell.is_empty()
        && cell
            .bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'+' | b'-' | b'.' | b'e' | b'E'));
    if plain {
        if let Ok(i) = cell.parse::<i64>() {
            return Value::from(i);
        }
        if let Some(n) = cell.parse::<f64>().ok().and_then(serde_json::Number::from_f64) {
            return Value::Number(n);
        }
    }
    Value::String(cell.to_owned())
}

fn write_jsonl(out: &mut dyn Write, header: &[String], rows: &[Vec<String>]) -> io::Result<()> {
    for row in rows {
        let object: Map<String, Value> =
            header.iter().zip(row).map(|(h, c)| (h.clone(), json_cell(c))).collect();
        serde_json::to_writer(&mut *out, &Value::Object(object))?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

fn write_aligned(out: &mut dyn Write, header: &[String], rows: &[Vec<String>]) -> io::Result<()> {
    let mut widths: Vec<usize> = header.iter().map(|h| h.chars().count()).collect();
    for row in rows {
        for (w, cell) in widths.iter_mut().zip(row) {
            *w = (*w).max(cell.chars().count());
        }
    }
    let line = |out: &mut dyn Write, cells: &[String]| -> io::Result<()> {
        let mut text = String::new();
        for (i, cell) in cells.iter().enumerate() {
            if i > 0 {
                text.push_str("  ");
            }
            text.push_str(cell);
            let width = widths.get(i).copied().unwrap_or(0);
            for _ in cell.chars().count()..width {
                text.push(' ');
            }
        }
        writeln!(out, "{}", text.trim_end())
    };
    line(out, header)?;
    for row in rows {
        line(out, row)?;
    }
    Ok(())
}
