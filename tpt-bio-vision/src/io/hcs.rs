//! OME-NGFF high-content screening (HCS) metadata: plates, wells, fields and channels.
//!
//! Follows the OME-NGFF plate and well sections (`plate-md`, `well-md`) of the specification
//! in two layouts:
//!
//! - NGFF 0.4 `.zattrs`, with `plate`, `well`, `multiscales` and `omero` at the top level.
//! - NGFF 0.5 `zarr.json`, with the same keys under `attributes.ome`.
//!
//! A plate is a group whose `plate` object lists rows, columns, wells and acquisitions. Each
//! well is a group whose `well` object lists its fields of view (`images`). Each field is an
//! image group, whose `multiscales` axes and `omero` channels describe its channel layout.
//!
//! The parsers take the JSON as a string, so the caller chooses how to read it from a
//! [`ZarrStore`](super::zarr::ZarrStore) or a file. The `validate` methods check the MUST
//! rules of the specification, which the parsers themselves do not enforce.

use std::collections::HashSet;

use serde::Deserialize;
use serde_json::Value;

use super::{FormatError, FormatResult};

/// A parsed well name: a row label and a column number, such as `"B03"` giving
/// `row = "B"` and `column = 3`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WellName {
    /// The letters before the digits. Case is kept as written.
    pub row: String,
    /// The number after the letters. Leading zeros are dropped.
    pub column: u32,
}

/// Parses a well name of the form letters then digits, such as `"A1"` or `"AA12"`.
///
/// # Errors
/// Returns [`FormatError::Invalid`] if the name is not letters followed by digits.
pub fn parse_well_name(name: &str) -> FormatResult<WellName> {
    let split = name.find(|c: char| !c.is_ascii_alphabetic()).unwrap_or(name.len());
    let (row, digits) = name.split_at(split);
    let invalid = || FormatError::Invalid(format!("well name {name:?} is not letters then digits"));
    if row.is_empty() || digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(invalid());
    }
    let column: u32 = digits.parse().map_err(|_| invalid())?;
    Ok(WellName { row: row.to_owned(), column })
}

/// One `acquisition` of a plate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Acquisition {
    /// Acquisition identifier, referred to by wells.
    pub id: u32,
    /// Optional display name.
    pub name: Option<String>,
    /// Largest number of fields acquired in any well of this acquisition.
    pub maximum_field_count: Option<u32>,
    /// Optional description.
    pub description: Option<String>,
    /// Optional start time, as an integer epoch timestamp.
    pub start_time: Option<i64>,
    /// Optional end time, as an integer epoch timestamp.
    pub end_time: Option<i64>,
}

/// One well of a plate, with its position in the plate's row and column lists.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlateWell {
    /// Path of the well group, such as `"B/03"`.
    pub path: String,
    /// Index into [`PlateMeta::rows`].
    pub row_index: usize,
    /// Index into [`PlateMeta::columns`].
    pub column_index: usize,
}

/// Typed OME-NGFF plate metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlateMeta {
    /// Optional plate name.
    pub name: Option<String>,
    /// NGFF version of the plate metadata.
    pub version: Option<String>,
    /// Row labels in plate order, such as `["A", "B", "C"]`.
    pub rows: Vec<String>,
    /// Column labels in plate order, such as `["1", "2", "3"]`.
    pub columns: Vec<String>,
    /// The wells present on the plate.
    pub wells: Vec<PlateWell>,
    /// Maximum number of fields per well, if the plate declares it.
    pub field_count: Option<u32>,
    /// The acquisitions, if any.
    pub acquisitions: Vec<Acquisition>,
}

/// Typed OME-NGFF well metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WellMeta {
    /// NGFF version of the well metadata.
    pub version: Option<String>,
    /// The fields of view of the well.
    pub images: Vec<WellImage>,
}

/// One field of view of a well.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WellImage {
    /// Path of the field's image group, relative to the well.
    pub path: String,
    /// Acquisition this field belongs to.
    pub acquisition: Option<u32>,
}

/// The channel layout of one field image, from its `multiscales` axes and `omero` channels.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImageChannels {
    /// Axis names of the first multiscale, in array order, such as `["t", "c", "z", "y", "x"]`.
    pub axes: Vec<String>,
    /// Index into [`axes`](Self::axes) of the channel axis, if the image declares one.
    pub channel_axis: Option<usize>,
    /// Channel labels in channel order, from `omero.channels`. Empty when there is no `omero`
    /// block. A channel without a label is `None`.
    pub labels: Vec<Option<String>>,
}

#[derive(Deserialize)]
struct Named {
    name: String,
}

#[derive(Deserialize)]
struct RawWell {
    path: String,
    #[serde(rename = "rowIndex")]
    row_index: usize,
    #[serde(rename = "columnIndex")]
    column_index: usize,
}

#[derive(Deserialize)]
struct RawAcquisition {
    id: u32,
    name: Option<String>,
    #[serde(rename = "maximumfieldcount")]
    maximum_field_count: Option<u32>,
    description: Option<String>,
    starttime: Option<i64>,
    endtime: Option<i64>,
}

#[derive(Deserialize)]
struct RawPlate {
    name: Option<String>,
    version: Option<String>,
    #[serde(default)]
    rows: Vec<Named>,
    #[serde(default)]
    columns: Vec<Named>,
    #[serde(default)]
    wells: Vec<RawWell>,
    field_count: Option<u32>,
    #[serde(default)]
    acquisitions: Vec<RawAcquisition>,
}

#[derive(Deserialize)]
struct RawWellImage {
    path: String,
    acquisition: Option<u32>,
}

#[derive(Deserialize)]
struct RawWellBody {
    version: Option<String>,
    #[serde(default)]
    images: Vec<RawWellImage>,
}

/// Parses the JSON text of a group's metadata file into a value.
fn parse_root(json: &str, what: &str) -> FormatResult<Value> {
    serde_json::from_str(json).map_err(|e| FormatError::Decode(format!("{what} metadata: {e}")))
}

/// The NGFF version stored at the `ome` level. NGFF 0.5 puts it there, rather than in the
/// object itself, so this is the fallback for a well with no `version` of its own.
fn ome_version(root: &Value) -> Option<String> {
    root.get("version")
        .or_else(|| root.get("attributes").and_then(|a| a.get("ome")).and_then(|o| o.get("version")))
        .and_then(Value::as_str)
        .map(str::to_owned)
}

/// Finds an NGFF object such as `plate` or `well`. NGFF 0.4 `.zattrs` has it at the top
/// level. NGFF 0.5 `zarr.json` has it under `attributes.ome`.
fn ome_object<'a>(root: &'a Value, key: &str) -> Option<&'a Value> {
    root.get(key)
        .or_else(|| root.get("attributes").and_then(|a| a.get("ome")).and_then(|o| o.get(key)))
        .or_else(|| root.get("ome").and_then(|o| o.get(key)))
}

/// Finds the `(row_index, column_index)` of a well name in the given label lists.
///
/// Row labels match case-insensitively. Column labels match by number, so the
/// well `"B03"` matches a column labelled `"3"` or `"03"`.
///
/// # Errors
/// Returns [`FormatError::Invalid`] if the name is malformed or its row or column
/// is not in the lists.
pub fn well_index(rows: &[String], columns: &[String], name: &str) -> FormatResult<(usize, usize)> {
    let parsed = parse_well_name(name)?;
    let row = rows
        .iter()
        .position(|r| r.eq_ignore_ascii_case(&parsed.row))
        .ok_or_else(|| FormatError::Invalid(format!("row {:?} not in plate", parsed.row)))?;
    let column = columns
        .iter()
        .position(|c| c.trim().parse::<u32>().ok() == Some(parsed.column))
        .ok_or_else(|| FormatError::Invalid(format!("column {} not in plate", parsed.column)))?;
    Ok((row, column))
}

/// True for NGFF names, which must be non-empty and contain only ASCII letters and digits.
fn is_ngff_name(name: &str) -> bool {
    !name.is_empty() && name.bytes().all(|b| b.is_ascii_alphanumeric())
}

fn check_labels(kind: &str, names: &[String]) -> FormatResult<()> {
    let mut seen = HashSet::new();
    for n in names {
        if !is_ngff_name(n) {
            return Err(FormatError::Invalid(format!("{kind} name {n:?} is not alphanumeric")));
        }
        if !seen.insert(n.as_str()) {
            return Err(FormatError::Invalid(format!("{kind} name {n:?} is duplicated")));
        }
    }
    Ok(())
}

impl PlateMeta {
    /// Parses the metadata of a plate group, in NGFF 0.4 `.zattrs` or 0.5 `zarr.json` form.
    ///
    /// # Errors
    /// Returns a decode error if the JSON has no `plate` object or if its fields are
    /// malformed. Call [`PlateMeta::validate`] to check the specification's rules.
    pub fn from_zattrs(json: &str) -> FormatResult<Self> {
        let root = parse_root(json, "plate")?;
        let obj = ome_object(&root, "plate")
            .ok_or_else(|| FormatError::Decode("plate metadata has no plate object".into()))?;
        let raw: RawPlate = serde_json::from_value(obj.clone())
            .map_err(|e| FormatError::Decode(format!("plate object: {e}")))?;
        Ok(Self {
            name: raw.name,
            version: raw.version,
            rows: raw.rows.into_iter().map(|n| n.name).collect(),
            columns: raw.columns.into_iter().map(|n| n.name).collect(),
            wells: raw
                .wells
                .into_iter()
                .map(|w| PlateWell {
                    path: w.path,
                    row_index: w.row_index,
                    column_index: w.column_index,
                })
                .collect(),
            field_count: raw.field_count,
            acquisitions: raw
                .acquisitions
                .into_iter()
                .map(|a| Acquisition {
                    id: a.id,
                    name: a.name,
                    maximum_field_count: a.maximum_field_count,
                    description: a.description,
                    start_time: a.starttime,
                    end_time: a.endtime,
                })
                .collect(),
        })
    }

    /// Checks the MUST rules of the plate specification.
    ///
    /// Row and column names must be alphanumeric and unique. The version must be present.
    /// Every well's row and column indices must be in range, and its path must be
    /// `row/column` for those indices. No two wells may share a position. Acquisition ids
    /// must be unique, and counts must be positive.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] naming the first rule that fails.
    pub fn validate(&self) -> FormatResult<()> {
        if self.version.is_none() {
            return Err(FormatError::Invalid("plate must declare a version".into()));
        }
        check_labels("row", &self.rows)?;
        check_labels("column", &self.columns)?;
        if self.field_count == Some(0) {
            return Err(FormatError::Invalid("field_count must be positive".into()));
        }
        let mut positions = HashSet::new();
        for w in &self.wells {
            let (Some(row), Some(col)) = (self.rows.get(w.row_index), self.columns.get(w.column_index))
            else {
                return Err(FormatError::Invalid(format!(
                    "well {:?} has an index outside the plate",
                    w.path
                )));
            };
            let expected = format!("{row}/{col}");
            if w.path != expected {
                return Err(FormatError::Invalid(format!(
                    "well path {:?} must be {expected:?} for its indices",
                    w.path
                )));
            }
            if !positions.insert((w.row_index, w.column_index)) {
                return Err(FormatError::Invalid(format!("well {:?} is listed twice", w.path)));
            }
        }
        let mut ids = HashSet::new();
        for a in &self.acquisitions {
            if !ids.insert(a.id) {
                return Err(FormatError::Invalid(format!("acquisition id {} is duplicated", a.id)));
            }
            if a.maximum_field_count == Some(0) {
                return Err(FormatError::Invalid(format!(
                    "acquisition {} maximumfieldcount must be positive",
                    a.id
                )));
            }
        }
        Ok(())
    }

    /// Maps a well name such as `"B03"` to `(row_index, column_index)`.
    ///
    /// # Errors
    /// Same as [`well_index`].
    pub fn well_index(&self, name: &str) -> FormatResult<(usize, usize)> {
        well_index(&self.rows, &self.columns, name)
    }

    /// The plate entry for a well name, if the plate lists that well.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] if the name does not map to a row and column
    /// of the plate. Returns `Ok(None)` if the position is valid but the plate does
    /// not list a well there.
    pub fn find_well(&self, name: &str) -> FormatResult<Option<&PlateWell>> {
        let (row, column) = self.well_index(name)?;
        Ok(self.wells.iter().find(|w| w.row_index == row && w.column_index == column))
    }
}

impl WellMeta {
    /// Parses the metadata of a well group, in NGFF 0.4 `.zattrs` or 0.5 `zarr.json` form.
    ///
    /// # Errors
    /// Returns a decode error if the JSON has no `well` object.
    pub fn from_zattrs(json: &str) -> FormatResult<Self> {
        let root = parse_root(json, "well")?;
        let obj = ome_object(&root, "well")
            .ok_or_else(|| FormatError::Decode("well metadata has no well object".into()))?;
        let body: RawWellBody = serde_json::from_value(obj.clone())
            .map_err(|e| FormatError::Decode(format!("well object: {e}")))?;
        Ok(Self {
            version: body.version.or_else(|| ome_version(&root)),
            images: body
                .images
                .into_iter()
                .map(|i| WellImage { path: i.path, acquisition: i.acquisition })
                .collect(),
        })
    }

    /// Checks the MUST rules of the well specification, against `plate` when it is given.
    ///
    /// Field paths must be alphanumeric and unique. When the plate has more than one
    /// acquisition, every field must name one, and it must be an acquisition of the plate.
    /// The number of fields must not exceed the plate's `field_count`, and each acquisition's
    /// fields must not exceed its `maximumfieldcount`.
    ///
    /// # Errors
    /// Returns [`FormatError::Invalid`] naming the first rule that fails.
    pub fn validate(&self, plate: Option<&PlateMeta>) -> FormatResult<()> {
        let mut paths = HashSet::new();
        for img in &self.images {
            if !is_ngff_name(&img.path) {
                return Err(FormatError::Invalid(format!("field path {:?} is not alphanumeric", img.path)));
            }
            if !paths.insert(img.path.as_str()) {
                return Err(FormatError::Invalid(format!("field path {:?} is duplicated", img.path)));
            }
        }
        let Some(plate) = plate else { return Ok(()) };
        if let Some(max) = plate.field_count
            && self.images.len() > max as usize
        {
            return Err(FormatError::Invalid(format!(
                "well has {} fields, plate field_count is {max}",
                self.images.len()
            )));
        }
        let multi = plate.acquisitions.len() > 1;
        let mut per_acquisition = std::collections::HashMap::<u32, usize>::new();
        for img in &self.images {
            match img.acquisition {
                Some(id) => {
                    if !plate.acquisitions.is_empty() && !plate.acquisitions.iter().any(|a| a.id == id) {
                        return Err(FormatError::Invalid(format!(
                            "field {:?} names acquisition {id}, which the plate does not define",
                            img.path
                        )));
                    }
                    *per_acquisition.entry(id).or_default() += 1;
                }
                None if multi => {
                    return Err(FormatError::Invalid(format!(
                        "field {:?} must name an acquisition, the plate has several",
                        img.path
                    )));
                }
                None => {}
            }
        }
        for a in &plate.acquisitions {
            if let (Some(max), Some(&n)) = (a.maximum_field_count, per_acquisition.get(&a.id))
                && n > max as usize
            {
                return Err(FormatError::Invalid(format!(
                    "acquisition {} has {n} fields, maximumfieldcount is {max}",
                    a.id
                )));
            }
        }
        Ok(())
    }
}

impl ImageChannels {
    /// Reads the channel layout of a field image from its `multiscales` and `omero` metadata.
    ///
    /// # Errors
    /// Returns a decode error if the JSON has no `multiscales` entry, or if its axes are
    /// malformed.
    pub fn from_zattrs(json: &str) -> FormatResult<Self> {
        let root = parse_root(json, "image")?;
        let first = ome_object(&root, "multiscales")
            .and_then(Value::as_array)
            .and_then(|a| a.first())
            .ok_or_else(|| FormatError::Decode("image has no multiscales entry".into()))?;
        let axes_json = first
            .get("axes")
            .and_then(Value::as_array)
            .ok_or_else(|| FormatError::Decode("multiscales has no axes list".into()))?;
        let mut axes = Vec::with_capacity(axes_json.len());
        let mut channel_axis = None;
        for (i, axis) in axes_json.iter().enumerate() {
            // NGFF 0.3 lists axes as plain names. 0.4 and later use objects.
            let (name, kind) = match axis {
                Value::String(s) => (s.clone(), None),
                Value::Object(o) => {
                    let name = o
                        .get("name")
                        .and_then(Value::as_str)
                        .ok_or_else(|| FormatError::Decode("axis object has no name".into()))?;
                    (name.to_owned(), o.get("type").and_then(Value::as_str))
                }
                _ => return Err(FormatError::Decode("axis is neither a name nor an object".into())),
            };
            if channel_axis.is_none() && (kind == Some("channel") || (kind.is_none() && name == "c")) {
                channel_axis = Some(i);
            }
            axes.push(name);
        }
        let labels = ome_object(&root, "omero")
            .and_then(|o| o.get("channels"))
            .and_then(Value::as_array)
            .map(|chs| {
                chs.iter().map(|c| c.get("label").and_then(Value::as_str).map(str::to_owned)).collect()
            })
            .unwrap_or_default();
        Ok(Self { axes, channel_axis, labels })
    }
}

/// Checks that every field of a well (or plate) has the same channel layout, and returns it.
///
/// Two fields agree when their axes, channel axis and channel labels are all equal. Images in
/// one well are usually acquired with the same channels, so a mismatch is worth reporting.
///
/// # Errors
/// Returns [`FormatError::Invalid`] if `fields` is empty or if any field differs from the first.
pub fn check_same_channels(fields: &[ImageChannels]) -> FormatResult<&ImageChannels> {
    let first = fields
        .first()
        .ok_or_else(|| FormatError::Invalid("no fields to compare channels for".into()))?;
    for (i, f) in fields.iter().enumerate().skip(1) {
        if f != first {
            return Err(FormatError::Invalid(format!(
                "field {i} has channel layout {:?}, field 0 has {:?}",
                f.labels, first.labels
            )));
        }
    }
    Ok(first)
}

#[cfg(test)]
mod tests {
    use super::*;

    const PLATE: &str = r#"{
        "plate": {
            "name": "screen-1",
            "version": "0.4",
            "field_count": 2,
            "rows": [{"name": "A"}, {"name": "B"}, {"name": "C"}],
            "columns": [{"name": "1"}, {"name": "2"}, {"name": "3"}, {"name": "4"}],
            "wells": [
                {"path": "A/1", "rowIndex": 0, "columnIndex": 0},
                {"path": "B/3", "rowIndex": 1, "columnIndex": 2},
                {"path": "C/4", "rowIndex": 2, "columnIndex": 3}
            ],
            "acquisitions": [
                {"id": 0, "name": "first pass", "maximumfieldcount": 2},
                {"id": 1, "description": "second pass"}
            ]
        }
    }"#;

    /// The 0.5 specification's example plate: two acquisitions, six wells, four fields.
    const SPEC_05_PLATE: &str = r#"{
        "zarr_format": 3, "node_type": "group",
        "attributes": { "ome": { "version": "0.5", "plate": {
            "acquisitions": [
                { "id": 1, "maximumfieldcount": 2, "name": "Meas_01(2012-07-31_10-41-12)", "starttime": 1343731272000 },
                { "id": 2, "maximumfieldcount": 2, "name": "Meas_02(201207-31_11-56-41)", "starttime": 1343735801000 }
            ],
            "columns": [ { "name": "1" }, { "name": "2" }, { "name": "3" } ],
            "field_count": 4, "name": "test",
            "rows": [ { "name": "A" }, { "name": "B" } ],
            "version": "0.5",
            "wells": [
                { "path": "A/1", "rowIndex": 0, "columnIndex": 0 },
                { "path": "A/2", "rowIndex": 0, "columnIndex": 1 },
                { "path": "A/3", "rowIndex": 0, "columnIndex": 2 },
                { "path": "B/1", "rowIndex": 1, "columnIndex": 0 },
                { "path": "B/2", "rowIndex": 1, "columnIndex": 1 },
                { "path": "B/3", "rowIndex": 1, "columnIndex": 2 }
            ]
        } } }
    }"#;

    /// The 0.5 specification's example well: four fields over two acquisitions.
    const SPEC_05_WELL: &str = r#"{
        "zarr_format": 3, "node_type": "group",
        "attributes": { "ome": { "version": "0.5", "well": { "images": [
            { "acquisition": 1, "path": "0" },
            { "acquisition": 1, "path": "1" },
            { "acquisition": 2, "path": "2" },
            { "acquisition": 2, "path": "3" }
        ] } } }
    }"#;

    #[test]
    fn parses_plate_metadata() {
        let plate = PlateMeta::from_zattrs(PLATE).unwrap();
        assert_eq!(plate.name.as_deref(), Some("screen-1"));
        assert_eq!(plate.version.as_deref(), Some("0.4"));
        assert_eq!(plate.rows, ["A", "B", "C"]);
        assert_eq!(plate.columns.len(), 4);
        assert_eq!(plate.field_count, Some(2));
        assert_eq!(plate.wells.len(), 3);
        assert_eq!(plate.wells[1].path, "B/3");
        assert_eq!(plate.acquisitions[0].maximum_field_count, Some(2));
        assert_eq!(plate.acquisitions[1].name, None);
        assert_eq!(plate.acquisitions[1].description.as_deref(), Some("second pass"));
        plate.validate().unwrap();
    }

    #[test]
    fn parses_ngff_05_plate_layout() {
        let plate = PlateMeta::from_zattrs(SPEC_05_PLATE).unwrap();
        assert_eq!(plate.name.as_deref(), Some("test"));
        assert_eq!(plate.version.as_deref(), Some("0.5"));
        assert_eq!(plate.wells.len(), 6);
        assert_eq!(plate.acquisitions[1].start_time, Some(1_343_735_801_000));
        assert_eq!(plate.acquisitions[0].end_time, None);
        plate.validate().unwrap();
        assert_eq!(plate.find_well("B03").unwrap().unwrap().path, "B/3");
    }

    #[test]
    fn validate_rejects_spec_violations() {
        let mut plate = PlateMeta::from_zattrs(SPEC_05_PLATE).unwrap();
        plate.wells[0].path = "A/2".into();
        assert!(plate.validate().is_err(), "path must match indices");

        let mut plate = PlateMeta::from_zattrs(SPEC_05_PLATE).unwrap();
        plate.wells[5].row_index = 7;
        assert!(plate.validate().is_err(), "index out of range");

        let mut plate = PlateMeta::from_zattrs(SPEC_05_PLATE).unwrap();
        let dup = plate.wells[0].clone();
        plate.wells.push(dup);
        assert!(plate.validate().is_err(), "duplicate well");

        let mut plate = PlateMeta::from_zattrs(SPEC_05_PLATE).unwrap();
        plate.rows[1].push('-');
        assert!(plate.validate().is_err(), "non-alphanumeric row");

        let mut plate = PlateMeta::from_zattrs(SPEC_05_PLATE).unwrap();
        plate.acquisitions[1].id = 1;
        assert!(plate.validate().is_err(), "duplicate acquisition id");

        let mut plate = PlateMeta::from_zattrs(SPEC_05_PLATE).unwrap();
        plate.version = None;
        assert!(plate.validate().is_err(), "version is required");
    }

    #[test]
    fn maps_well_names_to_indices() {
        let plate = PlateMeta::from_zattrs(PLATE).unwrap();
        assert_eq!(plate.well_index("A01").unwrap(), (0, 0));
        assert_eq!(plate.well_index("B03").unwrap(), (1, 2));
        assert_eq!(plate.well_index("c4").unwrap(), (2, 3));
        assert!(plate.well_index("D01").is_err());
        assert!(plate.well_index("A9").is_err());
        assert!(plate.well_index("03").is_err());
    }

    #[test]
    fn finds_listed_wells() {
        let plate = PlateMeta::from_zattrs(PLATE).unwrap();
        assert_eq!(plate.find_well("B03").unwrap().unwrap().path, "B/3");
        assert!(plate.find_well("A02").unwrap().is_none());
    }

    #[test]
    fn parses_well_names() {
        assert_eq!(parse_well_name("AA12").unwrap(), WellName { row: "AA".into(), column: 12 });
        assert_eq!(parse_well_name("B03").unwrap(), WellName { row: "B".into(), column: 3 });
        assert!(parse_well_name("B").is_err());
        assert!(parse_well_name("1B").is_err());
        assert!(parse_well_name("B3x").is_err());
    }

    #[test]
    fn parses_well_metadata() {
        let json = r#"{"well": {"version": "0.4", "images": [
            {"path": "0", "acquisition": 0}, {"path": "1", "acquisition": 0}]}}"#;
        let well = WellMeta::from_zattrs(json).unwrap();
        assert_eq!(well.images.len(), 2);
        assert_eq!(well.images[1], WellImage { path: "1".into(), acquisition: Some(0) });
        assert!(WellMeta::from_zattrs(r#"{"plate": {}}"#).is_err());
    }

    #[test]
    fn parses_ngff_05_well_and_validates_against_plate() {
        let plate = PlateMeta::from_zattrs(SPEC_05_PLATE).unwrap();
        let well = WellMeta::from_zattrs(SPEC_05_WELL).unwrap();
        assert_eq!(well.images.len(), 4);
        well.validate(Some(&plate)).unwrap();
        well.validate(None).unwrap();
    }

    #[test]
    fn validate_well_rejects_spec_violations() {
        let plate = PlateMeta::from_zattrs(SPEC_05_PLATE).unwrap();
        let mut well = WellMeta::from_zattrs(SPEC_05_WELL).unwrap();
        well.images[0].acquisition = Some(9);
        assert!(well.validate(Some(&plate)).is_err(), "unknown acquisition");

        let mut well = WellMeta::from_zattrs(SPEC_05_WELL).unwrap();
        well.images[0].acquisition = None;
        assert!(well.validate(Some(&plate)).is_err(), "acquisition required with several");

        let mut well = WellMeta::from_zattrs(SPEC_05_WELL).unwrap();
        well.images.push(WellImage { path: "4".into(), acquisition: Some(2) });
        assert!(well.validate(Some(&plate)).is_err(), "more fields than field_count");

        let mut well = WellMeta::from_zattrs(SPEC_05_WELL).unwrap();
        well.images[1].path = "0".into();
        assert!(well.validate(None).is_err(), "duplicate field path");

        let mut well = WellMeta::from_zattrs(SPEC_05_WELL).unwrap();
        well.images[0].acquisition = Some(2);
        well.images[1].acquisition = Some(2);
        well.images[2].acquisition = Some(2);
        assert!(well.validate(Some(&plate)).is_err(), "acquisition over maximumfieldcount");
    }

    #[test]
    fn rejects_plate_without_plate_object() {
        assert!(PlateMeta::from_zattrs(r#"{"multiscales": []}"#).is_err());
    }

    const FIELD_04: &str = r#"{
        "multiscales": [{ "axes": [
            {"name": "c", "type": "channel"}, {"name": "y", "type": "space"}, {"name": "x", "type": "space"}
        ] }],
        "omero": { "channels": [ {"label": "DAPI"}, {"label": "GFP"} ] }
    }"#;

    const FIELD_05: &str = r#"{
        "zarr_format": 3, "node_type": "group",
        "attributes": { "ome": { "version": "0.5",
            "multiscales": [{ "axes": [
                {"name": "t", "type": "time"}, {"name": "c", "type": "channel"},
                {"name": "y", "type": "space"}, {"name": "x", "type": "space"}
            ] }],
            "omero": { "channels": [ {"label": "DAPI"}, {"label": "GFP"} ] }
        } }
    }"#;

    #[test]
    fn reads_channel_layout_in_both_layouts() {
        let v04 = ImageChannels::from_zattrs(FIELD_04).unwrap();
        assert_eq!(v04.axes, ["c", "y", "x"]);
        assert_eq!(v04.channel_axis, Some(0));
        assert_eq!(v04.labels, [Some("DAPI".into()), Some("GFP".into())]);

        let v05 = ImageChannels::from_zattrs(FIELD_05).unwrap();
        assert_eq!(v05.axes, ["t", "c", "y", "x"]);
        assert_eq!(v05.channel_axis, Some(1));
        assert_eq!(v05.labels.len(), 2);
    }

    #[test]
    fn channel_axis_is_absent_without_a_channel_axis() {
        let json = r#"{"multiscales": [{"axes": ["y", "x"]}]}"#;
        let img = ImageChannels::from_zattrs(json).unwrap();
        assert_eq!(img.axes, ["y", "x"]);
        assert_eq!(img.channel_axis, None);
        assert!(img.labels.is_empty());
    }

    #[test]
    fn check_same_channels_detects_mismatch() {
        let a = ImageChannels::from_zattrs(FIELD_04).unwrap();
        let b = ImageChannels::from_zattrs(FIELD_05).unwrap();
        assert_eq!(check_same_channels(&[a.clone(), a.clone()]).unwrap(), &a);
        assert!(check_same_channels(&[a.clone(), b]).is_err());
        assert!(check_same_channels(&[]).is_err());
    }

    // Reference fixtures written by ome-zarr 0.21.0 (see tests/fixtures/hcs/README.md).
    const REF_04_PLATE: &str = include_str!("../../tests/fixtures/hcs/0.4/.zattrs");
    const REF_04_WELL_A1: &str = include_str!("../../tests/fixtures/hcs/0.4/A/1/.zattrs");
    const REF_04_WELL_B2: &str = include_str!("../../tests/fixtures/hcs/0.4/B/2/.zattrs");
    const REF_04_FIELD_A1_0: &str = include_str!("../../tests/fixtures/hcs/0.4/A/1/0/.zattrs");
    const REF_04_FIELD_A1_1: &str = include_str!("../../tests/fixtures/hcs/0.4/A/1/1/.zattrs");
    const REF_04_FIELD_B2_0: &str = include_str!("../../tests/fixtures/hcs/0.4/B/2/0/.zattrs");
    const REF_05_PLATE: &str = include_str!("../../tests/fixtures/hcs/0.5/zarr.json");
    const REF_05_WELL_A1: &str = include_str!("../../tests/fixtures/hcs/0.5/A/1/zarr.json");
    const REF_05_WELL_B2: &str = include_str!("../../tests/fixtures/hcs/0.5/B/2/zarr.json");
    const REF_05_FIELD_A1_0: &str = include_str!("../../tests/fixtures/hcs/0.5/A/1/0/zarr.json");
    const REF_05_FIELD_A1_1: &str = include_str!("../../tests/fixtures/hcs/0.5/A/1/1/zarr.json");

    fn check_reference_plate(plate_json: &str, wells: &[(&str, &str)], fields: &[&str]) {
        let plate = PlateMeta::from_zattrs(plate_json).unwrap();
        plate.validate().unwrap();
        assert_eq!(plate.name.as_deref(), Some("screen"));
        assert_eq!(plate.rows, ["A", "B"]);
        assert_eq!(plate.columns, ["1", "2", "3"]);
        assert_eq!(plate.field_count, Some(2));
        assert_eq!(plate.acquisitions.len(), 2);
        assert_eq!(plate.find_well("B02").unwrap().unwrap().path, "B/2");
        for (name, json) in wells {
            let well = WellMeta::from_zattrs(json).unwrap();
            well.validate(Some(&plate)).unwrap_or_else(|e| panic!("well {name}: {e}"));
            assert_eq!(well.images.len(), 2, "well {name}");
            assert_eq!(well.version.as_deref(), Some(if plate.version.as_deref() == Some("0.4") { "0.4" } else { "0.5" }));
        }
        let layouts: Vec<ImageChannels> =
            fields.iter().map(|j| ImageChannels::from_zattrs(j).unwrap()).collect();
        let shared = check_same_channels(&layouts).unwrap();
        assert_eq!(shared.axes, ["c", "y", "x"]);
        assert_eq!(shared.channel_axis, Some(0));
        assert_eq!(shared.labels, [Some("DAPI".into()), Some("GFP".into()), Some("RFP".into())]);
    }

    #[test]
    fn reference_plate_ngff_04_parses_and_validates() {
        check_reference_plate(
            REF_04_PLATE,
            &[("A/1", REF_04_WELL_A1), ("B/2", REF_04_WELL_B2)],
            &[REF_04_FIELD_A1_0, REF_04_FIELD_A1_1, REF_04_FIELD_B2_0],
        );
    }

    #[test]
    fn reference_plate_ngff_05_parses_and_validates() {
        check_reference_plate(
            REF_05_PLATE,
            &[("A/1", REF_05_WELL_A1), ("B/2", REF_05_WELL_B2)],
            &[REF_05_FIELD_A1_0, REF_05_FIELD_A1_1],
        );
    }

    #[test]
    fn reference_well_04_maps_fields_to_acquisitions() {
        let well = WellMeta::from_zattrs(REF_04_WELL_B2).unwrap();
        assert!(well.images.iter().all(|i| i.acquisition == Some(1)));
    }
}
