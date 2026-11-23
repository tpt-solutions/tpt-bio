//! OME-NGFF high-content screening (HCS) metadata: plates and wells.
//!
//! A plate is a Zarr group whose `.zattrs` holds a `plate` object with `rows`,
//! `columns`, `wells`, `field_count` and `acquisitions`. Each well is a group
//! whose `.zattrs` holds a `well` object with its `images` (fields). This module
//! parses both into typed structures, and maps well names such as `"B03"` to
//! row and column indices.
//!
//! The attribute JSON is passed in as a string, so the caller chooses how to read
//! it from a [`ZarrStore`](super::zarr::ZarrStore) or a file.

use serde::Deserialize;

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
    /// Number of fields per well, if the plate declares it.
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
struct PlateAttrs {
    plate: RawPlate,
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

#[derive(Deserialize)]
struct WellAttrs {
    well: RawWellBody,
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
        .ok_or_else(|| {
        FormatError::Invalid(format!("column {} not in plate", parsed.column))
    })?;
    Ok((row, column))
}

impl PlateMeta {
    /// Parses the `.zattrs` JSON of a plate group.
    ///
    /// # Errors
    /// Returns a decode error if the JSON does not have the NGFF `plate` object,
    /// or if the wells or acquisitions are malformed.
    pub fn from_zattrs(json: &str) -> FormatResult<Self> {
        let attrs: PlateAttrs = serde_json::from_str(json)
            .map_err(|e| FormatError::Decode(format!("plate .zattrs: {e}")))?;
        let raw = attrs.plate;
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
                })
                .collect(),
        })
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
    /// Parses the `.zattrs` JSON of a well group.
    ///
    /// # Errors
    /// Returns a decode error if the JSON does not have the NGFF `well` object.
    pub fn from_zattrs(json: &str) -> FormatResult<Self> {
        let attrs: WellAttrs = serde_json::from_str(json)
            .map_err(|e| FormatError::Decode(format!("well .zattrs: {e}")))?;
        Ok(Self {
            version: attrs.well.version,
            images: attrs
                .well
                .images
                .into_iter()
                .map(|i| WellImage { path: i.path, acquisition: i.acquisition })
                .collect(),
        })
    }
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
    fn rejects_plate_without_plate_object() {
        assert!(PlateMeta::from_zattrs(r#"{"multiscales": []}"#).is_err());
    }
}
