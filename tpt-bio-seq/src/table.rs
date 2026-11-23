//! Tabular export of alignment records through `tpt-strata`.
//!
//! `tpt-strata` is the in-process columnar engine that replaces Arrow for this
//! project. Its core crate has no dependencies. Parquet files would need
//! `tpt-strata-parquet`, which depends on Arrow and so is not used here.

use std::sync::Arc;

use tpt_strata::{ArrayRef, Batch, DataType, Field, I64Array, Scalar, Schema, StrArray, Table};

use crate::error::{Result, SeqError};
use crate::header::Header;
use crate::record::Record;

/// Column names of the record table, in order.
pub const COLUMNS: [&str; 11] =
    ["qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual"];

fn text(bytes: &[u8]) -> Scalar {
    if bytes.is_empty() {
        Scalar::Null
    } else {
        Scalar::Str(String::from_utf8_lossy(bytes).into_owned())
    }
}

/// Builds a columnar table with one row per record.
///
/// Positions are 1-based, as in SAM. Missing values (empty names, unset positions,
/// absent sequences) are NULL. Quality is written as Phred+33 text.
///
/// # Errors
/// Returns [`SeqError::Unsupported`] if the table cannot be built, which should not
/// happen for records that came from this crate.
pub fn records_to_table(records: &[Record], header: &Header) -> Result<Table> {
    let n = records.len();
    let mut qname = Vec::with_capacity(n);
    let mut flag = Vec::with_capacity(n);
    let mut rname = Vec::with_capacity(n);
    let mut pos = Vec::with_capacity(n);
    let mut mapq = Vec::with_capacity(n);
    let mut cigar = Vec::with_capacity(n);
    let mut rnext = Vec::with_capacity(n);
    let mut pnext = Vec::with_capacity(n);
    let mut tlen = Vec::with_capacity(n);
    let mut seq = Vec::with_capacity(n);
    let mut qual = Vec::with_capacity(n);

    for rec in records {
        qname.push(text(&rec.name));
        flag.push(Scalar::I64(i64::from(rec.flags.0)));
        rname
            .push(header.ref_name(rec.ref_id).map_or(Scalar::Null, |s| Scalar::Str(s.to_string())));
        pos.push(if rec.pos < 0 { Scalar::Null } else { Scalar::I64(rec.pos + 1) });
        mapq.push(Scalar::I64(i64::from(rec.mapq)));
        cigar.push(if rec.cigar.is_empty() {
            Scalar::Null
        } else {
            Scalar::Str(crate::record::format_cigar(&rec.cigar))
        });
        rnext.push(if rec.next_ref_id == rec.ref_id && rec.ref_id >= 0 {
            Scalar::Str("=".to_string())
        } else {
            header.ref_name(rec.next_ref_id).map_or(Scalar::Null, |s| Scalar::Str(s.to_string()))
        });
        pnext.push(if rec.next_pos < 0 { Scalar::Null } else { Scalar::I64(rec.next_pos + 1) });
        tlen.push(Scalar::I64(rec.tlen));
        seq.push(text(&rec.seq));
        qual.push(if rec.qual.is_empty() {
            Scalar::Null
        } else {
            Scalar::Str(rec.qual.iter().map(|&q| char::from(q + 33)).collect())
        });
    }

    // The schema is declared, not inferred: a column that is entirely NULL (for
    // example RNEXT when every record has `*`) has no value to infer a type from.
    let utf8 = |values: Vec<Scalar>| {
        ArrayRef::Utf8(StrArray::from_nulls(
            values
                .into_iter()
                .map(|v| match v {
                    Scalar::Str(s) => Some(Arc::from(s.as_str())),
                    _ => None,
                })
                .collect(),
        ))
    };
    let int64 = |values: Vec<Scalar>| {
        ArrayRef::Int64(I64Array::from_nulls(
            values
                .into_iter()
                .map(|v| match v {
                    Scalar::I64(x) => Some(x),
                    _ => None,
                })
                .collect(),
        ))
    };
    let columns: [(Vec<Scalar>, DataType); 11] = [
        (qname, DataType::Utf8),
        (flag, DataType::Int64),
        (rname, DataType::Utf8),
        (pos, DataType::Int64),
        (mapq, DataType::Int64),
        (cigar, DataType::Utf8),
        (rnext, DataType::Utf8),
        (pnext, DataType::Int64),
        (tlen, DataType::Int64),
        (seq, DataType::Utf8),
        (qual, DataType::Utf8),
    ];
    let mut fields = Vec::with_capacity(columns.len());
    let mut arrays = Vec::with_capacity(columns.len());
    for ((values, data_type), name) in columns.into_iter().zip(COLUMNS) {
        fields.push(Field::new(name, data_type, true));
        arrays.push(match data_type {
            DataType::Int64 => int64(values),
            _ => utf8(values),
        });
    }
    let schema = Arc::new(Schema::new(fields));
    let batch = Batch::try_new(schema.clone(), arrays)
        .map_err(|e| SeqError::Unsupported(format!("record table: {e}")))?;
    Ok(Table::from_batches(schema, vec![batch]))
}

/// Renders records as CSV with a header row, via `tpt-strata`.
///
/// # Errors
/// Returns an error if the table cannot be built. See [`records_to_table`].
pub fn records_to_csv(records: &[Record], header: &Header) -> Result<String> {
    Ok(tpt_strata::csv::to_csv(&records_to_table(records, header)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::SamReader;

    const SAM: &str = "@SQ\tSN:chr1\tLN:1000\n\
r1\t99\tchr1\t101\t60\t4M\t=\t201\t150\tACGT\tIIII\n\
r2\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t*\n";

    #[test]
    fn table_has_one_row_per_record_and_expected_columns() {
        let mut reader = SamReader::new(SAM.as_bytes()).unwrap();
        let header = reader.header().clone();
        let mut records = Vec::new();
        while let Some(r) = reader.next_record().unwrap() {
            records.push(r);
        }
        let table = records_to_table(&records, &header).unwrap();
        assert_eq!(table.row_count(), 2);
        assert_eq!(table.column_names(), COLUMNS.to_vec());
    }

    #[test]
    fn csv_export_uses_one_based_positions_and_nulls() {
        let mut reader = SamReader::new(SAM.as_bytes()).unwrap();
        let header = reader.header().clone();
        let mut records = Vec::new();
        while let Some(r) = reader.next_record().unwrap() {
            records.push(r);
        }
        let csv = records_to_csv(&records, &header).unwrap();
        let mut lines = csv.lines();
        assert_eq!(lines.next().unwrap(), COLUMNS.join(","));
        assert_eq!(lines.next().unwrap(), "r1,99,chr1,101,60,4M,=,201,150,ACGT,IIII");
        // The unmapped record has empty (NULL) reference and position fields.
        assert_eq!(lines.next().unwrap(), "r2,4,,,0,,,,0,AC,");
    }
}

/// Writes `s` as a JSON string literal, escaping quotes, backslashes and control bytes.
fn json_string(out: &mut String, s: &[u8]) {
    use std::fmt::Write as _;
    out.push('"');
    for c in String::from_utf8_lossy(s).chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn json_text_or_null(out: &mut String, s: Option<&[u8]>) {
    match s {
        Some(bytes) if !bytes.is_empty() => json_string(out, bytes),
        _ => out.push_str("null"),
    }
}

/// Renders records as JSON Lines: one object per record, with the same fields as
/// [`records_to_table`]. Missing values are `null`. The output can be read back
/// with `tpt_strata::json::parse_json`.
#[must_use]
pub fn records_to_jsonl(records: &[Record], header: &Header) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    for rec in records {
        let rname = header.ref_name(rec.ref_id).map(str::as_bytes);
        let rnext = if rec.next_ref_id == rec.ref_id && rec.ref_id >= 0 {
            Some(&b"="[..])
        } else {
            header.ref_name(rec.next_ref_id).map(str::as_bytes)
        };
        let qual: Vec<u8> = rec.qual.iter().map(|&q| q + 33).collect();
        let cigar = (!rec.cigar.is_empty()).then(|| crate::record::format_cigar(&rec.cigar));

        out.push('{');
        out.push_str("\"qname\":");
        json_text_or_null(&mut out, Some(&rec.name));
        let _ = write!(out, ",\"flag\":{}", rec.flags.0);
        out.push_str(",\"rname\":");
        json_text_or_null(&mut out, rname);
        out.push_str(",\"pos\":");
        json_position(&mut out, rec.pos);
        let _ = write!(out, ",\"mapq\":{}", rec.mapq);
        out.push_str(",\"cigar\":");
        json_text_or_null(&mut out, cigar.as_deref().map(str::as_bytes));
        out.push_str(",\"rnext\":");
        json_text_or_null(&mut out, rnext);
        out.push_str(",\"pnext\":");
        json_position(&mut out, rec.next_pos);
        let _ = write!(out, ",\"tlen\":{}", rec.tlen);
        out.push_str(",\"seq\":");
        json_text_or_null(&mut out, Some(&rec.seq));
        out.push_str(",\"qual\":");
        json_text_or_null(&mut out, Some(&qual));
        out.push_str("}\n");
    }
    out
}

/// Writes a 0-based position as a 1-based JSON number, or `null` when unset.
fn json_position(out: &mut String, zero_based: i64) {
    use std::fmt::Write as _;
    if zero_based < 0 {
        out.push_str("null");
    } else {
        let _ = write!(out, "{}", zero_based + 1);
    }
}

#[cfg(test)]
mod jsonl_tests {
    use super::*;
    use crate::sam::SamReader;

    #[test]
    fn jsonl_round_trips_through_tpt_strata_json_reader() {
        let sam = "@SQ\tSN:chr1\tLN:1000\n\
r\"1\t99\tchr1\t101\t60\t4M\t=\t201\t150\tACGT\tIIII\n\
r2\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t*\n";
        let mut reader = SamReader::new(sam.as_bytes()).unwrap();
        let header = reader.header().clone();
        let mut records = Vec::new();
        while let Some(r) = reader.next_record().unwrap() {
            records.push(r);
        }
        let text = records_to_jsonl(&records, &header);
        assert!(text.contains("\"qname\":\"r\\\"1\""), "{text}");
        let table = tpt_strata::json::parse_json(&text).unwrap();
        assert_eq!(table.row_count(), 2);
        assert!(table.column_names().contains(&"pos"));
    }
}
