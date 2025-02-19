pub mod intervaltree;
use crate::intervaltree::IntervalTree;
use std::collections::{BTreeSet, HashMap};


/// Represents the coverage level.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoverageLevel {
    Low,
    Medium,
    High,
}

/// Represents the type of data (treatment or control).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    Treatment,
    Control,
}

/// Represents the strand of the genomic region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

/// A struct to hold count metadata for a genomic region
///
/// # Fields
/// 
/// - `treatment_counts`, `control_counts`: vector of replicate counts in same order
///    as input
/// - `treatment_coverage`, `control_coverage`: Coverage level categories.
#[derive(Debug, Clone)]
pub struct CountableRegionMetadata {
    pub treatment_counts: Vec<u32>,
    pub treatment_coverage: Option<CoverageLevel>,
    pub control_counts: Vec<u32>,
    pub control_coverage: Option<CoverageLevel>,
}

/// Represents an interval in the interval tree. GenomicIntervals are 0 indexed, half open,
/// e.g. [0, 10) is 0-9.
/// 
/// # Fields
/// 
/// - `chr`: The chromosome name.
/// - `start`: The start position of the interval.
/// - `end`: The end position of the interval.
/// - `strand`: The strand of the genomic region, represented as a Strand enum.
/// - `data`: An optional field to hold additional data associated with the interval.
/// 
/// # Example
/// 
/// ```rust
/// 
/// use tag_counter::{CountableRegionMetadata, GenomicInterval, Strand};
/// 
/// let interval1 = GenomicInterval {
///   chr: "chr1".to_string(),
///   start: 100,
///   end: 200,
///   strand: Some(Strand::Plus),
///   data: Some(CountableRegionMetadata {
///     treatment_counts: vec![0, 3],
///     treatment_coverage: None,
///     control_counts: vec![0, 1],
///     control_coverage: None,
///   }),
/// };
/// 
/// let interval2: GenomicInterval<CountableRegionMetadata> = GenomicInterval{
///   chr: "chr13".to_string(),
///   start: 60000,
///   end: 777777,
///   strand: None,
///   data: None
/// };
/// ```
#[derive(Debug, Clone)]
pub struct GenomicInterval<T> {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub strand: Option<Strand>,
    pub data: Option<T>,
}

/// Represents a collection of genomic regions organized by chromosome.
/// 
/// ## Fields
/// 
/// - `trees`: A map from chromosome names to interval trees containing regions.
/// - `region_map`: A map from chromosome names to vectors of regions. Eg,
///    'chr1': [(0, 300), (5000, 6000), (5001, 6001), ...] Note: this map
///    is for the regions explicitly provided by the user.
#[derive(Debug)]
pub struct CountableRegionTree {
    /// Stores preprocessed, non-overlapping intervals for efficient querying.
    pub trees: HashMap<String, IntervalTree<u32, GenomicInterval<CountableRegionMetadata>>>,
    pub num_treatment_replicates: u8,
    pub num_control_replicates: u8,
}

impl CountableRegionTree {
    /// Creates a new `CountableRegionTree`.
    pub fn new(
        num_treatment_replicates: u8, 
        num_control_replicates: u8
    ) -> Self {
        Self {
            trees: HashMap::new(),
            num_treatment_replicates,
            num_control_replicates,
        }
    }

    /// Constructs the interval tree from the region map using `BTreeSet`.
    ///
    /// This method ensures that the interval tree contains **non-overlapping** intervals.
    /// It processes intervals from `region_map`, extracts their breakpoints, and builds
    /// the `IntervalTree` for each chromosome.
    ///
    /// # Example
    ///
    /// ```rust
    /// use std::collections::{HashMap, HashSet, BTreeSet};
    /// use tag_counter::{CountableRegionTree, GenomicInterval, CountableRegionMetadata};
    ///
    /// let mut region_map: HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>> = HashMap::new();
    /// 
    /// // Create a new CountableRegionTree
    /// let mut tree = CountableRegionTree::new(3, 2);
    ///
    /// // Manually add regions to the region map
    /// region_map.insert(
    ///     "chr1".to_string(),
    ///     vec![
    ///         GenomicInterval { start: 10, end: 50, chr: "chr1".to_string(), strand: None, data: None },
    ///         GenomicInterval { start: 30, end: 70, chr: "chr1".to_string(), strand: None, data: None },
    ///         GenomicInterval { start: 80, end: 90, chr: "chr1".to_string(), strand: None, data: None },
    ///     ],
    /// );
    ///
    /// // Construct the interval tree from region_map
    /// tree.construct(&region_map);
    ///
    /// // Retrieve the constructed interval tree for "chr1"
    /// let interval_tree = tree.trees.get("chr1").expect("Interval tree should exist");
    ///
    /// // Verify that the tree contains the expected non-overlapping intervals
    /// let expected_intervals = vec![(10, 30), (30, 50), (50, 70), (70, 80), (80, 90)];
    /// let actual_intervals: Vec<_> = interval_tree
    ///     .find(0..u32::MAX)  // Fetch all stored intervals
    ///     .map(|entry| (entry.interval().start, entry.interval().end)) // Ensure direct field access
    ///     .collect();
    ///
    /// assert_eq!(
    ///     actual_intervals.iter().collect::<HashSet<_>>(),
    ///     expected_intervals.iter().collect::<HashSet<_>>()
    /// );
    /// ```
    pub fn construct(&mut self, region_map: &HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>>) {
        for (chr, intervals) in region_map {
            let mut btreeset = BTreeSet::new();

            // Insert all intervals into the BTreeSet to get non-overlapping regions
            for interval in intervals {
                btreeset.insert(interval.start);
                btreeset.insert(interval.end);
            }

            // Convert BTreeSet into an IntervalTree
            let tree = self.trees.entry(chr.clone()).or_insert_with(IntervalTree::new);

            let mut iter = btreeset.iter().peekable();
            while let Some(&seg_start) = iter.next() {
                if let Some(&seg_end) = iter.peek() {
                    tree.insert(
                        seg_start..*seg_end, 
                        GenomicInterval {
                            chr: chr.clone(),
                            start: seg_start,
                            end: *seg_end,
                            strand: None,
                            data: Some(CountableRegionMetadata {
                                treatment_counts: vec![0; self.num_treatment_replicates as usize],
                                treatment_coverage: None,
                                control_counts: vec![0; self.num_control_replicates as usize],
                                control_coverage: None,
                            }),
                        }
                    );
                }
            }
        }
    }

    /// Adds count data to overlapping genomic regions in the interval tree.
    ///
    /// This function updates the count of either treatment or control data for a specific replicate
    /// in all regions overlapping the given range.
    ///
    /// # Example
    ///
    /// ```rust
    /// use std::collections::HashMap;
    /// use tag_counter::{CountableRegionTree, GenomicInterval, CountableRegionMetadata, DataType};
    ///
    /// 
    /// let mut region_map: HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>> = HashMap::new();
    /// // Initialize the region tree with 3 treatment replicates and 2 control replicates.
    /// let mut region_tree = CountableRegionTree::new(3, 2);
    ///
    /// // Manually add a genomic region to the region map.
    /// let interval = GenomicInterval {
    ///     chr: "chr1".to_string(),
    ///     start: 100,
    ///     end: 200,
    ///     strand: None,
    ///     data: Some(CountableRegionMetadata {
    ///         treatment_counts: vec![0; 3],
    ///         treatment_coverage: None,
    ///         control_counts: vec![0; 2],
    ///         control_coverage: None,
    ///     }),
    /// };
    ///
    /// // Insert the interval into the region map.
    /// region_map.insert("chr1".to_string(), vec![interval]);
    ///
    /// // Construct the interval tree.
    /// region_tree.construct(&region_map);
    ///
    /// // Add a treatment count of 5 to the first treatment replicate (index 0).
    /// region_tree.add_counts("chr1", 120, 180, 5, DataType::Treatment, 0).unwrap();
    ///
    /// // Verify that the count was updated correctly.
    /// let overlaps = region_tree.trees.get("chr1").unwrap()
    ///     .find(120..180)
    ///     .collect::<Vec<_>>();
    ///
    /// assert_eq!(overlaps.len(), 1);
    /// assert_eq!(overlaps[0].data().data.as_ref().unwrap().treatment_counts[0], 5);
    ///
    /// // Add a control count of 3 to the first control replicate (index 0).
    /// region_tree.add_counts("chr1", 150, 190, 3, DataType::Control, 0).unwrap();
    ///
    /// // Verify that the control count was updated correctly.
    /// let overlaps = region_tree.trees.get("chr1").unwrap()
    ///     .find(150..190)
    ///     .collect::<Vec<_>>();
    ///
    /// assert_eq!(overlaps.len(), 1);
    /// assert_eq!(overlaps[0].data().data.as_ref().unwrap().control_counts[0], 3);
    /// ```
    pub fn add_counts(&mut self, chrom: &str, start: u32, end: u32, count: u32, data_type: DataType, replicate_index: u8) -> Result<(), String> {

        // Check if the chromosome exists in the tree
        let Some(tree) = self.trees.get_mut(chrom) else {
            return Err(format!("Chromosome '{}' not found in the interval tree.", chrom));
        };

        // Raise error if data_type is not Treatment or Control
        if data_type != DataType::Treatment && data_type != DataType::Control {
            return Err(format!("Invalid data type: {:?}", data_type));
        }

        // Raise error if replicate_index is out of bounds
        // (data_type is unsigned, so a negative can't be passed. Same re: count)
        if replicate_index >= self.num_treatment_replicates && data_type == DataType::Treatment {
            return Err(format!("Replicate index {} out of bounds for treatment counts.", replicate_index));
        } else if replicate_index >= self.num_control_replicates && data_type == DataType::Control {
            return Err(format!("Replicate index {} out of bounds for control counts.", replicate_index));
        }

        // Find overlapping regions
        let overlapping_regions: Vec<_> = tree.find(start..end).collect();

        // If the length of the overlapping region is greater than 1, raise an error
        if overlapping_regions.len() > 1 {
            return Err(format!("Multiple overlapping regions found for {}: {:?}", chrom, overlapping_regions));
        }

        // Update counts for each overlapping region
        for mut entry in tree.find_mut(start..end) {
            if let Some(region) = entry.data().data.as_mut() {
                match data_type {
                    DataType::Treatment => {
                        region.treatment_counts[replicate_index as usize] += count;
                    }
                    DataType::Control => {
                        region.control_counts[replicate_index as usize] += count;
                    }
                }
            }
        }

    Ok(())
    }


    /// Returns overlapping regions for a query position.
    ///
    /// # Arguments
    ///
    /// * `chrom` - The name of the chromosome.
    /// * `query_start` - The start position of the query.
    /// * `query_end` - The end position of the query.
    ///
    /// # Returns
    ///
    /// - If the chromosome is **not found**, returns `Err(String)`.
    /// - If no overlaps are found, returns `Ok(None)`.
    /// - If overlaps exist, returns `Ok(Some(Vec<GenomicInterval<CountableRegionMetadata>>))`.
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut region_tree = CountableRegionTree::new(3, 2);
    ///
    /// // Add overlapping intervals
    /// region_tree.add_region("chr1", 100, 200, None, true).unwrap();
    /// region_tree.add_region("chr1", 150, 250, None, true).unwrap();
    ///
    /// // Query for overlaps within the full range
    /// let overlaps = region_tree.find_overlaps("chr1", 100, 250).unwrap().unwrap();
    ///
    /// // Assert the number of intervals returned
    /// assert_eq!(overlaps.len(), 3);
    ///
    /// // Assert the specific intervals returned
    /// assert_eq!(overlaps[0].start, 100);
    /// assert_eq!(overlaps[0].end, 150);
    ///
    /// assert_eq!(overlaps[1].start, 150);
    /// assert_eq!(overlaps[1].end, 200);
    ///
    /// assert_eq!(overlaps[2].start, 200);
    /// assert_eq!(overlaps[2].end, 250);
    /// ```
    pub fn find_overlaps<'a>(&'a self, chrom: &str, query_start: u32, query_end: u32) -> Result<Option<Vec<&'a GenomicInterval<CountableRegionMetadata>>>, String> {
        let Some(tree) = self.trees.get(chrom) else {
            return Err(format!("Chromosome '{}' not found in the interval tree.", chrom));
        };

        let results: Vec<_> = tree.find(query_start..query_end).map(|entry| entry.data()).collect();
        if results.is_empty() {
            Ok(None)
        } else {
            Ok(Some(results))
        }
    }


}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_add_region_splits_intervals_correctly() {
//         let mut region_tree = CountableRegionTree::new(3, 2);

//         // Add overlapping regions
//         region_tree.add_region("chr1", 100, 200, Some(Strand::Plus), true).unwrap();
//         region_tree.add_region("chr1", 150, 250, Some(Strand::Plus), true).unwrap();

//         // Query for overlaps
//         let overlaps = region_tree.find_overlaps("chr1", 100, 250).unwrap().unwrap();

//         // Expected intervals: [100, 150], [150, 200], [200, 250]
//         let expected_intervals = vec![(100, 150), (150, 200), (200, 250)];
//         let actual_intervals: Vec<_> = overlaps.iter().map(|region| (region.start, region.end)).collect();

//         assert_eq!(actual_intervals, expected_intervals);
//     }

//     #[test]
//     fn test_add_region_with_no_overlap() {
//         let mut region_tree = CountableRegionTree::new(3, 2);

//         // Add two non-overlapping regions
//         region_tree.add_region("chr1", 100, 200, Some(Strand::Minus), true).unwrap();
//         region_tree.add_region("chr1", 300, 400, Some(Strand::Minus), true).unwrap();

//         // Check that the regions exist
//         let overlaps_1 = region_tree.find_overlaps("chr1", 100, 200).unwrap().unwrap();
//         let overlaps_2 = region_tree.find_overlaps("chr1", 300, 400).unwrap().unwrap();

//         assert_eq!(overlaps_1.len(), 1);
//         assert_eq!(overlaps_1[0].start, 100);
//         assert_eq!(overlaps_1[0].end, 200);

//         assert_eq!(overlaps_2.len(), 1);
//         assert_eq!(overlaps_2[0].start, 300);
//         assert_eq!(overlaps_2[0].end, 400);
//     }

//     #[test]
//     fn test_add_counts_updates_correct_replicate() {
//         let mut region_tree = CountableRegionTree::new(3, 2);

//         // Add an interval
//         region_tree.add_region("chr1", 100, 200, Some(Strand::Plus), true).unwrap();

//         // Add count to the first treatment replicate
//         region_tree.add_counts("chr1", 100, 200, 5, DataType::Treatment, 0).unwrap();

//         // Verify the count was updated
//         let overlaps = region_tree.find_overlaps("chr1", 100, 200).unwrap().unwrap();
//         assert_eq!(overlaps.len(), 1);
//         assert_eq!(overlaps[0].data.as_ref().unwrap().treatment_counts[0], 5);
//     }

//     #[test]
//     fn test_add_counts_errors_on_invalid_replicate_index() {
//         let mut region_tree = CountableRegionTree::new(3, 2);

//         // Add an interval
//         region_tree.add_region("chr1", 100, 200, Some(Strand::Plus), true).unwrap();

//         // Try adding counts to an out-of-bounds replicate index
//         let result = region_tree.add_counts("chr1", 100, 200, 5, DataType::Treatment, 3);
//         assert!(result.is_err());
//         assert_eq!(result.unwrap_err(), "Replicate index 3 out of bounds for treatment counts.");
//     }

//     #[test]
//     fn test_find_overlaps_returns_none_when_no_overlap() {
//         let mut region_tree = CountableRegionTree::new(3, 2);

//         // Add an interval
//         region_tree.add_region("chr1", 100, 200, Some(Strand::Plus), true).unwrap();

//         // Query a non-overlapping range
//         let result = region_tree.find_overlaps("chr1", 300, 400);

//         // Ensure it returns Ok(None)
//         assert!(result.is_ok());
//         assert!(result.unwrap().is_none());
//     }

//     #[test]
//     fn test_find_overlaps_returns_error_for_missing_chromosome() {
//         let region_tree = CountableRegionTree::new(3, 2);

//         // Query a chromosome that hasn't been added
//         let result = region_tree.find_overlaps("chrX", 100, 200);

//         // Ensure it returns an error
//         assert!(result.is_err());
//         assert_eq!(result.unwrap_err(), "Chromosome 'chrX' not found in the interval tree.");
//     }
// }

