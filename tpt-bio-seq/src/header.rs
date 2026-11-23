//! SAM/BAM header: reference dictionary plus the raw header text.

use std::collections::HashMap;
use std::fmt::Write as _;

use crate::error::{Result, SeqError};

/// A reference sequence (`@SQ` line).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    /// Sequence name (`SN`).
    pub name: String,
    /// Sequence length (`LN`).
    pub len: u64,
}

/// A tagged header line (`@HD`, `@RG` or `@PG`) with its `TAG:value` fields in order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HeaderRecord {
    /// The two-letter record type without the `@`, such as `"RG"`.
    pub kind: String,
    /// The `TAG:value` fields, in file order. Tags are two characters.
    pub fields: Vec<(String, String)>,
}

impl HeaderRecord {
    /// The value of the first field with this two-character tag.
    #[must_use]
    pub fn get(&self, tag: &str) -> Option<&str> {
        self.fields.iter().find(|(t, _)| t == tag).map(|(_, v)| v.as_str())
    }
}

/// A parsed header. `text` is kept verbatim so SAM and BAM output round-trip.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Header {
    text: String,
    references: Vec<Reference>,
    index: HashMap<String, i32>,
    records: Vec<HeaderRecord>,
    comments: Vec<String>,
}

/// Parses the `@HD`, `@RG`, `@PG` and `@CO` lines. Reference lines are ignored.
fn parse_meta(text: &str) -> Result<(Vec<HeaderRecord>, Vec<String>)> {
    let mut records = Vec::new();
    let mut comments = Vec::new();
    for line in text.lines() {
        if let Some(comment) = line.strip_prefix("@CO	") {
            comments.push(comment.to_string());
            continue;
        }
        let Some(kind) = ["HD", "RG", "PG"].into_iter().find(|k| {
            line.strip_prefix('@').is_some_and(|rest| rest.starts_with(k))
        }) else {
            continue;
        };
        let mut fields = Vec::new();
        for field in line.split('\t').skip(1) {
            let Some((tag, value)) = field.split_once(':') else {
                return Err(SeqError::malformed(
                    "SAM header",
                    format!("bad header field in line: {line}"),
                ));
            };
            fields.push((tag.to_string(), value.to_string()));
        }
        records.push(HeaderRecord { kind: kind.to_string(), fields });
    }
    Ok((records, comments))
}

impl Header {
    /// Parses header text (lines starting with `@`).
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] for an `@SQ` line without a valid `SN` or `LN`,
    /// or for a duplicate reference name.
    pub fn parse(text: &str) -> Result<Self> {
        let mut references = Vec::new();
        let mut index = HashMap::new();
        let (records, comments) = parse_meta(text)?;
        for line in text.lines() {
            let Some(rest) = line.strip_prefix("@SQ") else { continue };
            let mut name = None;
            let mut len = None;
            for field in rest.split('\t').skip(1) {
                if let Some(v) = field.strip_prefix("SN:") {
                    name = Some(v.to_string());
                } else if let Some(v) = field.strip_prefix("LN:") {
                    len = v.parse::<u64>().ok();
                }
            }
            let (Some(name), Some(len)) = (name, len) else {
                return Err(SeqError::malformed("SAM header", format!("bad @SQ line: {line}")));
            };
            if index.insert(name.clone(), references.len() as i32).is_some() {
                return Err(SeqError::malformed(
                    "SAM header",
                    format!("duplicate reference {name}"),
                ));
            }
            references.push(Reference { name, len });
        }
        let text = text.to_string();
        Ok(Self { text, references, index, records, comments })
    }

    /// Combines header text with an authoritative reference dictionary, as read
    /// from a BAM file. The text is kept as-is even if its `@SQ` lines disagree.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] on duplicate reference names.
    pub fn from_text_and_references(text: &str, references: Vec<Reference>) -> Result<Self> {
        let mut index = HashMap::with_capacity(references.len());
        for (i, r) in references.iter().enumerate() {
            if index.insert(r.name.clone(), i as i32).is_some() {
                return Err(SeqError::malformed(
                    "BAM header",
                    format!("duplicate reference {}", r.name),
                ));
            }
        }
        // Metadata is best-effort here: a BAM header is kept verbatim even when a
        // metadata line is malformed, so reading never fails on it.
        let (records, comments) = parse_meta(text).unwrap_or_default();
        Ok(Self { text: text.to_string(), references, index, records, comments })
    }

    /// Builds a header from references only, with no other lines.
    ///
    /// # Errors
    /// Returns [`SeqError::Malformed`] on duplicate reference names.
    pub fn from_references(references: &[Reference]) -> Result<Self> {
        let mut text = String::new();
        for r in references {
            // Writing to a String cannot fail.
            let _ = writeln!(text, "@SQ\tSN:{}\tLN:{}", r.name, r.len);
        }
        Self::parse(&text)
    }

    /// The raw header text, as stored.
    #[must_use]
    pub fn text(&self) -> &str {
        &self.text
    }

    /// The reference dictionary, in order.
    #[must_use]
    pub fn references(&self) -> &[Reference] {
        &self.references
    }

    /// The `@HD` line, if present. It is kept as a [`HeaderRecord`].
    #[must_use]
    pub fn hd(&self) -> Option<&HeaderRecord> {
        self.records.iter().find(|r| r.kind == "HD")
    }

    /// All `@RG` read-group lines, in order.
    pub fn read_groups(&self) -> impl Iterator<Item = &HeaderRecord> {
        self.records.iter().filter(|r| r.kind == "RG")
    }

    /// All `@PG` program lines, in order.
    pub fn programs(&self) -> impl Iterator<Item = &HeaderRecord> {
        self.records.iter().filter(|r| r.kind == "PG")
    }

    /// The text of all `@CO` comment lines, in order, without the `@CO	` prefix.
    #[must_use]
    pub fn comments(&self) -> &[String] {
        &self.comments
    }

    /// The index of a reference by name.
    #[must_use]
    pub fn ref_id(&self, name: &str) -> Option<i32> {
        self.index.get(name).copied()
    }

    /// The name of a reference by index.
    #[must_use]
    pub fn ref_name(&self, id: i32) -> Option<&str> {
        usize::try_from(id).ok().and_then(|i| self.references.get(i)).map(|r| r.name.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_references_and_keeps_text() {
        let text = "@HD	VN:1.6	SO:coordinate\n@SQ	SN:chr1	LN:248956422\n@SQ	SN:chr2	LN:242193529\n@RG	ID:rg1\n";
        let h = Header::parse(text).unwrap();
        assert_eq!(h.text(), text);
        assert_eq!(h.references().len(), 2);
        assert_eq!(h.ref_id("chr2"), Some(1));
        assert_eq!(h.ref_name(0), Some("chr1"));
        assert_eq!(h.ref_name(-1), None);
    }

    #[test]
    fn parses_header_metadata() {
        let text = "@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:10
@RG	ID:rg1	SM:s1
@RG	ID:rg2	SM:s2
@PG	ID:bwa	PN:bwa
@CO	free text	kept
";
        let h = Header::parse(text).unwrap();
        assert_eq!(h.hd().and_then(|r| r.get("SO")), Some("coordinate"));
        let rgs: Vec<_> = h.read_groups().map(|r| r.get("SM")).collect();
        assert_eq!(rgs, [Some("s1"), Some("s2")]);
        assert_eq!(h.programs().next().and_then(|r| r.get("PN")), Some("bwa"));
        assert_eq!(h.comments(), ["free text	kept"]);
        assert_eq!(h.references().len(), 1);
    }

    #[test]
    fn rejects_bad_and_duplicate_sq() {
        assert!(Header::parse("@SQ	SN:chr1\n").is_err());
        assert!(Header::parse("@SQ	SN:chr1	LN:5\n@SQ	SN:chr1	LN:5\n").is_err());
    }
}
