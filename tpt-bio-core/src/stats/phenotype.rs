//! Phenotype grouping and per-group aggregation, such as per-well summaries.
//!
//! A phenotype or well is just a label per measurement. [`aggregate`] reduces the
//! measurements sharing a label to one number, and [`group_indices`] lists the
//! row indices behind each label.

use alloc::collections::BTreeMap;
use alloc::vec::Vec;

use super::StatsError;
use super::normalize::median;

/// The reduction applied to each group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Aggregate {
    /// Arithmetic mean.
    Mean,
    /// Median. Averages the two middle values for even counts.
    Median,
}

/// Maps each label to the indices of the rows that carry it, in ascending order.
#[must_use]
pub fn group_indices<K: Ord + Clone>(labels: &[K]) -> BTreeMap<K, Vec<usize>> {
    let mut groups: BTreeMap<K, Vec<usize>> = BTreeMap::new();
    for (i, label) in labels.iter().enumerate() {
        groups.entry(label.clone()).or_default().push(i);
    }
    groups
}

/// Reduces `values` to one number per distinct key in `keys`.
///
/// For per-well aggregation, pass the well identifier as the key. For phenotype
/// summaries, pass the phenotype label.
///
/// # Errors
/// Returns [`StatsError::Shape`] when the two slices differ in length, and
/// [`StatsError::Empty`] when there are no measurements.
pub fn aggregate<K: Ord + Clone>(
    keys: &[K],
    values: &[f64],
    how: Aggregate,
) -> Result<BTreeMap<K, f64>, StatsError> {
    if keys.len() != values.len() {
        return Err(StatsError::Shape("one key is required per value"));
    }
    if keys.is_empty() {
        return Err(StatsError::Empty);
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err(StatsError::NonFinite);
    }
    let out = group_indices(keys)
        .into_iter()
        .map(|(key, idx)| {
            let vals: Vec<f64> = idx.iter().map(|&i| values[i]).collect();
            let v = match how {
                Aggregate::Mean => vals.iter().sum::<f64>() / vals.len() as f64,
                Aggregate::Median => median(&vals),
            };
            (key, v)
        })
        .collect();
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn per_well_mean_and_median() {
        let wells = ["A1", "A2", "A1", "A1", "A2"];
        let vals = [1.0, 10.0, 3.0, 8.0, 20.0];
        let mean = aggregate(&wells, &vals, Aggregate::Mean).unwrap();
        assert!((mean["A1"] - 4.0).abs() < 1e-12);
        assert!((mean["A2"] - 15.0).abs() < 1e-12);
        let med = aggregate(&wells, &vals, Aggregate::Median).unwrap();
        assert!((med["A1"] - 3.0).abs() < 1e-12);
        assert!((med["A2"] - 15.0).abs() < 1e-12);
    }

    #[test]
    fn grouping_lists_indices_in_order() {
        let phen = ["mitotic", "normal", "mitotic"];
        let g = group_indices(&phen);
        assert_eq!(g["mitotic"], vec![0, 2]);
        assert_eq!(g["normal"], vec![1]);
    }

    #[test]
    fn rejects_mismatched_lengths_and_empty_input() {
        assert!(matches!(aggregate(&[1u8, 2], &[1.0], Aggregate::Mean), Err(StatsError::Shape(_))));
        let empty: [u8; 0] = [];
        assert_eq!(aggregate(&empty, &[], Aggregate::Mean), Err(StatsError::Empty));
    }
}
