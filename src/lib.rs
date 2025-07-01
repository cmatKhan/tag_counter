pub mod intervaltree;
pub mod parsers;
pub mod quantification;

use crate::intervaltree::IntervalTree;
use std::collections::{BTreeSet, HashMap};

/// Represents the coverage level.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoverageLevel {
    Low,
    Medium,
    High,
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
/// - `counts`: vector of replicate counts in same order
///   as input
/// - `coverage`: Coverage level categories.
#[derive(Debug, Clone)]
pub struct CountableRegionMetadata {
    pub counts: Vec<u32>,
    pub coverage: Option<CoverageLevel>,
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
///     counts: vec![0, 3],
///     coverage: None,
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
///   'chr1': [(0, 300), (5000, 6000), (5001, 6001), ...] Note: this map
///   is for the regions explicitly provided by the user.
#[derive(Debug)]
pub struct CountableRegionTree {
    /// Stores preprocessed, non-overlapping intervals for efficient querying.
    pub trees: HashMap<String, IntervalTree<u32, GenomicInterval<CountableRegionMetadata>>>,
    pub n_replicates: u8,
}

impl CountableRegionTree {
    /// Creates a new `CountableRegionTree`.
    pub fn new(n_replicates: u8) -> Self {
        Self {
            trees: HashMap::new(),
            n_replicates,
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
    /// let mut tree = CountableRegionTree::new(3);
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
    /// let expected_intervals = vec![(10, 30), (30, 50), (50, 70), (80, 90)];
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
    pub fn construct(
        &mut self,
        region_map: &HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>>,
    ) {
        for (chr, intervals) in region_map {
            let mut btreeset = BTreeSet::new();

            // Insert all intervals into the BTreeSet to get non-overlapping regions
            for interval in intervals {
                btreeset.insert(interval.start);
                btreeset.insert(interval.end);
            }

            // Convert BTreeSet into an IntervalTree
            let tree = self.trees.entry(chr.clone()).or_default();

            let mut iter = btreeset.iter().peekable();
            while let Some(&seg_start) = iter.next() {
                if let Some(&seg_end) = iter.peek() {
                    let is_in_regions = intervals
                        .iter()
                        .any(|orig| seg_start >= orig.start && *seg_end <= orig.end);
                    if is_in_regions {
                        tree.insert(
                            seg_start..*seg_end,
                            GenomicInterval {
                                chr: chr.clone(),
                                start: seg_start,
                                end: *seg_end,
                                strand: None,
                                data: Some(CountableRegionMetadata {
                                    counts: vec![0; self.n_replicates as usize],
                                    coverage: None,
                                }),
                            },
                        );
                    }
                }
            }
        }
    }

    /// Adds count data to a single, non-overlapping genomic interval in the
    /// interval tree
    ///
    /// This function updates the count of a single genomic interval. If the interval
    /// over which the count applies overlaps multiple regions, an error is raised.
    ///
    /// # Example
    ///
    /// ```rust
    /// use std::collections::HashMap;
    /// use tag_counter::{CountableRegionTree, GenomicInterval, CountableRegionMetadata};
    ///
    ///
    /// let mut region_map: HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>> = HashMap::new();
    /// // Initialize the region tree with 3 treatment replicates and 2 control replicates.
    /// let mut region_tree = CountableRegionTree::new(3);
    ///
    /// // Manually add a genomic region to the region map.
    /// let interval = GenomicInterval {
    ///     chr: "chr1".to_string(),
    ///     start: 100,
    ///     end: 200,
    ///     strand: None,
    ///     data: Some(CountableRegionMetadata {
    ///         counts: vec![0; 3],
    ///         coverage: None,
    ///     }),
    /// };
    ///
    /// // Insert the interval into the region map.
    /// region_map.insert("chr1".to_string(), vec![interval]);
    ///
    /// // Construct the interval tree.
    /// region_tree.construct(&region_map);
    ///
    /// // Add a count of 5 to the first treatment replicate (index 0).
    /// region_tree.add_counts("chr1", 120, 180, 5, 0).unwrap();
    ///
    /// // Add a count of 10 to [99, 100). This should not add to the interval [100, 200).
    /// region_tree.add_counts("chr1", 99, 100, 10, 0).unwrap();
    ///
    /// // Verify that the count was updated correctly.
    /// let overlaps = region_tree.trees.get("chr1").unwrap()
    ///     .find(120..180)
    ///     .collect::<Vec<_>>();
    ///
    /// assert_eq!(overlaps.len(), 1);
    /// assert_eq!(overlaps[0].data().data.as_ref().unwrap().counts[0], 5);
    ///
    /// // Add a control count of 3 to the first control replicate (index 0).
    /// region_tree.add_counts("chr1", 199, 200, 3, 0).unwrap();
    ///
    /// // Add count of 10 to index 200 -- this should not add to the interval [100, 200).
    /// region_tree.add_counts("chr1", 200, 201, 10, 0).unwrap();
    ///
    /// // Verify that the control count was updated correctly.
    /// let overlaps = region_tree.trees.get("chr1").unwrap()
    ///     .find(150..190)
    ///     .collect::<Vec<_>>();
    ///
    /// assert_eq!(overlaps.len(), 1);
    /// assert_eq!(overlaps[0].data().data.as_ref().unwrap().counts[0], 8);
    /// ```
    pub fn add_counts(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        count: u32,
        replicate_index: u8,
    ) -> Result<u32, String> {
        let Some(tree) = self.trees.get_mut(chrom) else {
            return Err(format!(
                "CountableRegionTree.add_counts() error: Chromosome '{}' not found in the interval tree.",
                chrom
            ));
        };

        if replicate_index >= self.n_replicates {
            return Err(format!(
                "Replicate index {} out of bounds.",
                replicate_index
            ));
        }

        let mut overlapping: Vec<_> = tree.find_mut(start..end).collect();

        if overlapping.len() > 1 {
            return Err(format!(
                "Multiple overlapping regions found for {}:[{}, {}): {:?}",
                chrom, start, end, overlapping
            ));
        }

        if let Some(mut entry) = overlapping.pop() {
            if let Some(region) = entry.data().data.as_mut() {
                region.counts[replicate_index as usize] += count;
                // count successfully added
                return Ok(count);
            }
        }
        // no overlapping region found, so no count added. However, no error
        Ok(0)
    }

    /// Returns overlapping regions for a query position.
    ///
    /// This function retrieves genomic intervals from the interval tree that overlap the given query range.
    ///
    /// # Example
    ///
    /// ```rust
    /// use std::collections::{HashMap, HashSet};
    /// use tag_counter::{CountableRegionTree, GenomicInterval, CountableRegionMetadata};
    ///
    /// let mut region_map: HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>> = HashMap::new();
    ///
    /// // Create a new CountableRegionTree
    /// let mut region_tree = CountableRegionTree::new(3);
    ///
    /// // Manually add genomic regions to the region map
    /// region_map.insert(
    ///     "chr1".to_string(),
    ///     vec![
    ///         GenomicInterval {
    ///             chr: "chr1".to_string(),
    ///             start: 100,
    ///             end: 200,
    ///             strand: None,
    ///             data: Some(CountableRegionMetadata {
    ///                 counts: vec![0; 3],
    ///                 coverage: None,
    ///             }),
    ///         },
    ///         GenomicInterval {
    ///             chr: "chr1".to_string(),
    ///             start: 150,
    ///             end: 250,
    ///             strand: None,
    ///             data: Some(CountableRegionMetadata {
    ///                 counts: vec![0; 3],
    ///                 coverage: None,
    ///             }),
    ///         },
    ///     ],
    /// );
    ///
    /// // Construct the interval tree from the region map
    /// region_tree.construct(&region_map);
    ///
    /// // Query for overlaps within the full range
    /// let overlaps = region_tree.find_overlaps("chr1", 100, 200).unwrap().unwrap();
    ///
    /// // Verify that two intervals are found. Note that this should be 3 intervals
    /// // since intervals stored in the tree are non-overlapping.
    /// assert_eq!(overlaps.len(), 2);
    ///
    /// // Verify the exact regions returned
    /// let expected_intervals: HashSet<_> = vec![
    ///     (100, 150),
    ///     (150, 200),
    /// ].into_iter().collect();
    ///
    /// let actual_intervals: HashSet<_> = overlaps
    ///     .iter()
    ///     .map(|region| (region.start, region.end))
    ///     .collect();
    /// assert_eq!(actual_intervals, expected_intervals);
    /// ```
    pub fn find_overlaps<'a>(
        &'a self,
        chrom: &str,
        query_start: u32,
        query_end: u32,
    ) -> Result<Option<Vec<&'a GenomicInterval<CountableRegionMetadata>>>, String> {
        let Some(tree) = self.trees.get(chrom) else {
            return Err(format!(
                "CountableRegionTree.find_overlaps error: Chromosome '{}' not found in the interval tree.",
                chrom
            ));
        };

        let results: Vec<_> = tree
            .find(query_start..query_end)
            .map(|entry| entry.data())
            .collect();
        if results.is_empty() {
            Ok(None)
        } else {
            Ok(Some(results))
        }
    }

    /// Get a total of counts from each replicate over a specified region
    pub fn get_counts(&self, chrom: &str, start: u32, end: u32) -> Result<Vec<u32>, String> {
        let mut counts = vec![0; self.n_replicates as usize];

        match self.find_overlaps(chrom, start, end)? {
            Some(overlapping_regions) => {
                for region in overlapping_regions {
                    if let Some(data) = region.data.as_ref() {
                        for i in 0..self.n_replicates {
                            counts[i as usize] += data.counts[i as usize];
                        }
                    }
                }
                Ok(counts)
            }
            None => Ok(counts),
        }
    }

    /// Returns the total counts across all chromosomes and replicates.
    ///
    /// This is the sum of the values returned by `get_chromosome_totals()`.
    ///
    /// # Example
    ///
    /// ```rust
    /// use std::collections::HashMap;
    /// use tag_counter::{CountableRegionTree, GenomicInterval, CountableRegionMetadata};
    ///
    /// let mut region_map: HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>> = HashMap::new();
    ///
    /// // Create regions for chr1 and chr2
    /// region_map.insert(
    ///     "chr1".to_string(),
    ///     vec![GenomicInterval {
    ///         chr: "chr1".to_string(),
    ///         start: 100,
    ///         end: 200,
    ///         strand: None,
    ///         data: Some(CountableRegionMetadata {
    ///             counts: vec![0; 3],
    ///             coverage: None,
    ///         }),
    ///     }]
    /// );
    ///
    /// region_map.insert(
    ///     "chr2".to_string(),
    ///     vec![GenomicInterval {
    ///         chr: "chr2".to_string(),
    ///         start: 0,
    ///         end: 10,
    ///         strand: None,
    ///         data: Some(CountableRegionMetadata {
    ///             counts: vec![0; 3],
    ///             coverage: None,
    ///         }),
    ///     }]
    /// );
    ///
    /// let mut region_tree = CountableRegionTree::new(3);
    /// region_tree.construct(&region_map);
    ///
    /// region_tree.add_counts("chr1", 100, 120, 5, 0).unwrap();
    /// region_tree.add_counts("chr1", 199, 200, 3, 2).unwrap();
    /// region_tree.add_counts("chr2", 3, 5, 10, 0).unwrap();
    ///
    /// let total = region_tree.get_total_counts();
    /// assert_eq!(total, vec![15, 0, 3]);
    /// ```
    pub fn get_chromosome_totals(&self) -> HashMap<String, Vec<u32>> {
        let mut result = HashMap::new();

        for chr in self.trees.keys() {
            match self.get_counts(chr, 0, u32::MAX) {
                Ok(counts) => {
                    result.insert(chr.clone(), counts);
                }
                Err(e) => {
                    eprintln!("Warning: failed to get counts for {}: {}", chr, e);
                }
            }
        }

        result
    }

    /// Returns the total counts across all chromosomes and replicates.
    /// This is the sum of the values in `get_chromosome_totals`.
    pub fn get_total_counts(&self) -> Vec<u32> {
        let mut total = vec![0u32; self.n_replicates as usize];

        for chr_counts in self.get_chromosome_totals().values() {
            for (i, &count) in chr_counts.iter().enumerate() {
                total[i] += count;
            }
        }

        total
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_countable_region_tree() {
        // Step 1: Create a region_map with 3 chromosomes, each having 3 regions
        let mut region_map: HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>> =
            HashMap::new();

        let chromosomes = vec!["chr1", "chr2", "chr3"];
        for &chr in &chromosomes {
            region_map.insert(
                chr.to_string(),
                vec![
                    GenomicInterval {
                        chr: chr.to_string(),
                        start: 100,
                        end: 200,
                        strand: None,
                        data: Some(CountableRegionMetadata {
                            counts: vec![0; 3],
                            coverage: None,
                        }),
                    },
                    GenomicInterval {
                        chr: chr.to_string(),
                        start: 150,
                        end: 250,
                        strand: None,
                        data: Some(CountableRegionMetadata {
                            counts: vec![0; 3],
                            coverage: None,
                        }),
                    },
                    GenomicInterval {
                        chr: chr.to_string(),
                        start: 300,
                        end: 400,
                        strand: None,
                        data: Some(CountableRegionMetadata {
                            counts: vec![0; 3],
                            coverage: None,
                        }),
                    },
                ],
            );
        }

        // Step 2: Instantiate two CountableRegionTrees (control & treatment)
        let mut control_tree = CountableRegionTree::new(3);
        let mut treatment_tree = CountableRegionTree::new(3);

        // Step 3: Construct trees from the region_map
        control_tree.construct(&region_map);
        treatment_tree.construct(&region_map);

        // Step 4: Add single-base-pair counts to control and treatment trees
        let test_counts = vec![
            (100, 101, 5), // First base of first region
            (150, 151, 3), // First base of second overlapping region
            (200, 201, 7), // Just outside the first region, should not be counted
            (250, 251, 2), // First base of the end of second region
            (300, 301, 6), // First base of the third region
        ];

        for (start, end, count) in &test_counts {
            control_tree
                .add_counts("chr1", *start, *end, *count, 0)
                .unwrap();
            treatment_tree
                .add_counts("chr1", *start, *end, *count + 2, 0) // Treatment gets 2 more counts
                .unwrap();
        }

        // Step 5: Extract counts from the region_map and compare
        let expected_counts = vec![
            (100, 200, (8, 12)),  // Expected control count in region [100, 200)
            (150, 250, (10, 14)), // Expected control count in region [150, 250)
            (300, 400, (6, 8)),   // Expected control count in region [300, 400)
        ];

        for (start, end, expected) in expected_counts {
            let control_counts = control_tree.get_counts("chr1", start, end).unwrap();
            let treatment_counts = treatment_tree.get_counts("chr1", start, end).unwrap();

            // Ensure control counts match expected
            assert_eq!(control_counts[0], expected.0);

            // Ensure treatment is control +2 for each count
            assert_eq!(treatment_counts[0], expected.1);
        }

        // Step 6: Assert total counts per chromosome
        let control_chr_totals = control_tree.get_chromosome_totals();
        let treatment_chr_totals = treatment_tree.get_chromosome_totals();

        // note that this is 21/31 b/c the 2 counts from [250, 251) are not in
        // the region_map
        assert_eq!(control_chr_totals["chr1"], vec![21, 0, 0]);
        assert_eq!(treatment_chr_totals["chr1"], vec![29, 0, 0]);

        // Other chromosomes had no counts added
        for chr in ["chr2", "chr3"] {
            assert_eq!(control_chr_totals[chr], vec![0, 0, 0]);
            assert_eq!(treatment_chr_totals[chr], vec![0, 0, 0]);
        }

        // Step 7: Assert overall totals
        assert_eq!(control_tree.get_total_counts(), vec![21, 0, 0]);
        assert_eq!(treatment_tree.get_total_counts(), vec![29, 0, 0]);
    }

    #[test]
    fn test_single_region_single_count() {
        let mut tree = CountableRegionTree::new(1);
        let mut region_map = HashMap::new();

        region_map.insert(
            "chr1".to_string(),
            vec![GenomicInterval {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                strand: None,
                data: Some(CountableRegionMetadata {
                    counts: vec![0; 1],
                    coverage: None,
                }),
            }],
        );

        tree.construct(&region_map);
        tree.add_counts("chr1", 150, 160, 10, 0).unwrap();

        assert_eq!(tree.get_total_counts(), vec![10]);
    }

    #[test]
    fn test_no_matching_region() {
        let mut tree = CountableRegionTree::new(1);
        let mut region_map = HashMap::new();

        region_map.insert(
            "chr1".to_string(),
            vec![GenomicInterval {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                strand: None,
                data: Some(CountableRegionMetadata {
                    counts: vec![0; 1],
                    coverage: None,
                }),
            }],
        );

        tree.construct(&region_map);

        let result = tree.add_counts("chr1", 250, 300, 5, 0);
        assert_eq!(result.unwrap(), 0);
        assert_eq!(tree.get_total_counts(), vec![0]);
    }

    #[test]
    fn test_multiple_overlapping_regions_error() {
        let mut tree = CountableRegionTree::new(1);
        let mut region_map = HashMap::new();

        region_map.insert(
            "chr1".to_string(),
            vec![
                GenomicInterval {
                    chr: "chr1".to_string(),
                    start: 100,
                    end: 200,
                    strand: None,
                    data: Some(CountableRegionMetadata {
                        counts: vec![0; 1],
                        coverage: None,
                    }),
                },
                GenomicInterval {
                    chr: "chr1".to_string(),
                    start: 150,
                    end: 250,
                    strand: None,
                    data: Some(CountableRegionMetadata {
                        counts: vec![0; 1],
                        coverage: None,
                    }),
                },
            ],
        );

        tree.construct(&region_map);
        let result = tree.add_counts("chr1", 140, 180, 5, 0);
        assert!(result.is_err());
    }

    #[test]
    fn test_multiple_replicates() {
        let mut tree = CountableRegionTree::new(2);
        let mut region_map = HashMap::new();

        region_map.insert(
            "chr1".to_string(),
            vec![GenomicInterval {
                chr: "chr1".to_string(),
                start: 0,
                end: 100,
                strand: None,
                data: Some(CountableRegionMetadata {
                    counts: vec![0; 2],
                    coverage: None,
                }),
            }],
        );

        tree.construct(&region_map);
        tree.add_counts("chr1", 10, 20, 5, 0).unwrap();
        tree.add_counts("chr1", 30, 40, 3, 1).unwrap();

        assert_eq!(tree.get_total_counts(), vec![5, 3]);
    }

    #[test]
    fn test_interval_edges() {
        let mut tree = CountableRegionTree::new(1);
        let mut region_map = HashMap::new();

        region_map.insert(
            "chr1".to_string(),
            vec![GenomicInterval {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                strand: None,
                data: Some(CountableRegionMetadata {
                    counts: vec![0; 1],
                    coverage: None,
                }),
            }],
        );

        tree.construct(&region_map);

        assert_eq!(tree.add_counts("chr1", 100, 100, 1, 0).unwrap(), 0);
        assert_eq!(tree.add_counts("chr1", 99, 100, 1, 0).unwrap(), 0);
        assert_eq!(tree.add_counts("chr1", 100, 101, 1, 0).unwrap(), 1);

        assert_eq!(tree.get_total_counts(), vec![1]);
    }
}
