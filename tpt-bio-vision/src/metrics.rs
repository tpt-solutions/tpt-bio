//! Object-level validation of a predicted label image against ground truth.
//!
//! Objects are matched one-to-one by intersection-over-union (`IoU`). Candidate pairs
//! with `IoU >= threshold` are taken greedily in descending `IoU` order (ties by
//! predicted then true id), so each object matches at most once. A predicted object
//! with no match is a false positive; an unmatched ground-truth object is a false
//! negative.
//!
//! Precision = TP / (TP + FP), recall = TP / (TP + FN), and `F1` = 2TP / (2TP + FP + FN).
//! Each ratio with a zero denominator is 0, except that `F1` is 1 when both images have
//! no objects, since they agree exactly.

use std::collections::HashMap;

use crate::image::{Image, ImageError};
use crate::label::relabel_sequential;

/// One matched pair of objects.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObjectMatch {
    /// Consecutive id of the predicted object.
    pub predicted: u32,
    /// Consecutive id of the ground-truth object.
    pub truth: u32,
    /// Intersection over union of the pair.
    pub iou: f64,
}

/// Result of matching a prediction against ground truth.
#[derive(Debug, Clone, PartialEq)]
pub struct MatchStats {
    /// Matched pairs, in the order they were accepted (descending `IoU`).
    pub matches: Vec<ObjectMatch>,
    /// Number of matched pairs.
    pub true_positives: usize,
    /// Predicted objects with no match.
    pub false_positives: usize,
    /// Ground-truth objects with no match.
    pub false_negatives: usize,
    /// `TP / (TP + FP)`.
    pub precision: f64,
    /// `TP / (TP + FN)`.
    pub recall: f64,
    /// `2TP / (2TP + FP + FN)`.
    pub f1: f64,
}

/// Matches the objects of `predicted` to those of `truth` by `IoU`, and computes precision,
/// recall and `F1`. Label values are treated as object identities; `0` is background.
///
/// `iou_threshold` should be in `(0, 1]`. A value at or below 0 accepts any overlap.
///
/// # Errors
/// Returns [`ImageError::ShapeMismatch`] if the two label images differ in shape.
pub fn match_objects(
    predicted: &Image<u32>,
    truth: &Image<u32>,
    iou_threshold: f64,
) -> Result<MatchStats, ImageError> {
    predicted.check_same_shape(truth)?;
    let p = relabel_sequential(predicted);
    let t = relabel_sequential(truth);
    let mut pred_area = vec![0usize; p.count as usize];
    let mut true_area = vec![0usize; t.count as usize];
    let mut overlap: HashMap<(u32, u32), usize> = HashMap::new();
    for (&pl, &tl) in p.labels.data().iter().zip(t.labels.data()) {
        if pl != 0 {
            pred_area[pl as usize - 1] += 1;
        }
        if tl != 0 {
            true_area[tl as usize - 1] += 1;
        }
        if pl != 0 && tl != 0 {
            *overlap.entry((pl, tl)).or_insert(0) += 1;
        }
    }

    let mut candidates: Vec<ObjectMatch> = overlap
        .iter()
        .filter_map(|(&(pl, tl), &inter)| {
            let union = pred_area[pl as usize - 1] + true_area[tl as usize - 1] - inter;
            let iou = inter as f64 / union as f64;
            (iou >= iou_threshold).then_some(ObjectMatch { predicted: pl, truth: tl, iou })
        })
        .collect();
    candidates.sort_by(|a, b| {
        b.iou.total_cmp(&a.iou).then(a.predicted.cmp(&b.predicted)).then(a.truth.cmp(&b.truth))
    });

    let mut pred_used = vec![false; p.count as usize];
    let mut true_used = vec![false; t.count as usize];
    let mut matches = Vec::new();
    for c in candidates {
        let (pi, ti) = (c.predicted as usize - 1, c.truth as usize - 1);
        if !pred_used[pi] && !true_used[ti] {
            pred_used[pi] = true;
            true_used[ti] = true;
            matches.push(c);
        }
    }

    let tp = matches.len();
    let fp = p.count as usize - tp;
    let fn_ = t.count as usize - tp;
    let precision = ratio(tp, tp + fp);
    let recall = ratio(tp, tp + fn_);
    let f1 = if p.count == 0 && t.count == 0 { 1.0 } else { ratio(2 * tp, 2 * tp + fp + fn_) };
    Ok(MatchStats {
        matches,
        true_positives: tp,
        false_positives: fp,
        false_negatives: fn_,
        precision,
        recall,
        f1,
    })
}

/// `IoU` of two specific objects, identified by their label values in each image.
/// Returns `None` if either object is absent.
///
/// # Errors
/// Returns [`ImageError::ShapeMismatch`] if the two label images differ in shape.
pub fn object_iou(
    a: &Image<u32>,
    a_label: u32,
    b: &Image<u32>,
    b_label: u32,
) -> Result<Option<f64>, ImageError> {
    a.check_same_shape(b)?;
    let (mut inter, mut area_a, mut area_b) = (0usize, 0usize, 0usize);
    for (&x, &y) in a.data().iter().zip(b.data()) {
        let in_a = x == a_label;
        let in_b = y == b_label;
        area_a += usize::from(in_a);
        area_b += usize::from(in_b);
        inter += usize::from(in_a && in_b);
    }
    let union = area_a + area_b - inter;
    Ok((union > 0 && area_a > 0 && area_b > 0).then(|| inter as f64 / union as f64))
}

fn ratio(num: usize, den: usize) -> f64 {
    if den == 0 { 0.0 } else { num as f64 / den as f64 }
}

#[allow(clippy::float_cmp)]
#[cfg(test)]
mod tests {
    use super::*;

    fn labels(w: usize, data: Vec<u32>) -> Image<u32> {
        Image::new_2d(w, data.len() / w, data).unwrap()
    }

    #[test]
    fn perfect_match_has_f1_one() {
        let a = labels(4, vec![1, 1, 0, 0, 0, 0, 2, 2]);
        let s = match_objects(&a, &a, 0.5).unwrap();
        assert_eq!(s.true_positives, 2);
        assert_eq!(s.f1, 1.0);
        assert!(s.matches.iter().all(|m| (m.iou - 1.0).abs() < 1e-12));
    }

    #[test]
    fn hand_computed_iou_and_counts() {
        // Truth: object A = 4 pixels (x 0..2, y 0..2). Object B = 2 pixels at (4,0),(5,0).
        // Prediction: object X = 3 pixels, covering 3 of A's pixels -> IoU = 3 / (3+4-3) = 0.75.
        // Object Y = 1 pixel at (4,0), which is half of B -> IoU = 1 / (1+2-1) = 0.5.
        // Object Z = 1 pixel at (7,0), outside everything -> false positive.
        let w = 8;
        let mut truth = vec![0u32; 16];
        for &(x, y) in &[(0, 0), (1, 0), (0, 1), (1, 1), (4, 0), (5, 0)] {
            truth[y * w + x] = if x < 4 { 1 } else { 2 };
        }
        let mut pred = vec![0u32; 16];
        for &(x, y) in &[(0, 0), (1, 0), (0, 1)] {
            pred[y * w + x] = 1;
        }
        pred[4] = 2;
        pred[7] = 3;
        let t = Image::new_2d(w, 2, truth).unwrap();
        let p = Image::new_2d(w, 2, pred).unwrap();

        let s = match_objects(&p, &t, 0.4).unwrap();
        assert_eq!(s.true_positives, 2);
        assert_eq!(s.false_positives, 1);
        assert_eq!(s.false_negatives, 0);
        // Greedy: X (0.75) takes A first; Y (0.5) takes B.
        assert_eq!(s.matches[0].predicted, 1);
        assert!((s.matches[0].iou - 0.75).abs() < 1e-12);
        assert!((s.matches[1].iou - 0.5).abs() < 1e-12);
        // F1 = 2*2 / (2*2 + 1 + 0) = 0.8.
        assert!((s.f1 - 0.8).abs() < 1e-12);
        assert!((s.precision - 2.0 / 3.0).abs() < 1e-12);
        assert_eq!(s.recall, 1.0);
    }

    #[test]
    fn threshold_rejects_weak_overlap() {
        let t = labels(4, vec![1, 1, 1, 1]);
        let p = labels(4, vec![1, 0, 0, 0]);
        // IoU = 1/4 = 0.25 < 0.5.
        let s = match_objects(&p, &t, 0.5).unwrap();
        assert_eq!(s.true_positives, 0);
        assert_eq!(s.f1, 0.0);
        assert_eq!((s.false_positives, s.false_negatives), (1, 1));
    }

    #[test]
    fn empty_images_agree() {
        let e = labels(2, vec![0, 0]);
        assert_eq!(match_objects(&e, &e, 0.5).unwrap().f1, 1.0);
    }

    #[test]
    fn object_iou_for_specific_labels() {
        let a = labels(4, vec![1, 1, 0, 0]);
        let b = labels(4, vec![0, 1, 1, 0]);
        // Intersection 1, union 3.
        assert!((object_iou(&a, 1, &b, 1).unwrap().unwrap() - 1.0 / 3.0).abs() < 1e-12);
        assert_eq!(object_iou(&a, 9, &b, 1).unwrap(), None);
    }

    #[test]
    fn rejects_shape_mismatch() {
        let a = labels(2, vec![0, 0]);
        let b = labels(1, vec![0, 0]);
        assert!(match_objects(&a, &b, 0.5).is_err());
    }
}
