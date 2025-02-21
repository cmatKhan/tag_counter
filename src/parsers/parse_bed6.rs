use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

use crate::{CountableRegionMetadata, GenomicInterval};

/// Reads a BED file and parses it into a `region_map`. Note that this currently
/// ignores strand
///
/// The BED file is expected to be tab-separated with at least 3 columns:
/// `chrom`, `start`, `end`. Additional columns are ignored.
///
/// # Arguments
///
/// * `bed_path` - Path to the BED file.
///
/// # Returns
///
/// * `Ok(region_map)` - A hashmap mapping chromosome names to vectors of `GenomicInterval`.
/// * `Err(String)` - If the file cannot be read or if there is a formatting issue.
pub fn parse_bed6<P: AsRef<Path>>(
    bed_path: P,
) -> Result<HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>>, String> {
    let file = File::open(&bed_path).map_err(|e| format!("Failed to open file: {}", e))?;
    let reader = io::BufReader::new(file);
    let mut region_map: HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>> =
        HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading line {}: {}", line_num + 1, e))?;
        let fields: Vec<&str> = line.split_whitespace().collect();

        if fields.len() < 3 {
            return Err(format!(
                "Invalid BED format at line {}: {:?}",
                line_num + 1,
                line
            ));
        }

        let chr = fields[0].to_string();
        let start: u32 = fields[1]
            .parse()
            .map_err(|_| format!("Invalid start at line {}", line_num + 1))?;
        let end: u32 = fields[2]
            .parse()
            .map_err(|_| format!("Invalid end at line {}", line_num + 1))?;

        // Ensure BED is 0-indexed, half-open.
        if start >= end {
            return Err(format!(
                "Invalid interval at line {}: start >= end",
                line_num + 1
            ));
        }

        let interval = GenomicInterval {
            chr: chr.clone(),
            start,
            end,
            strand: None,
            data: Some(CountableRegionMetadata {
                counts: vec![0; 3], // Placeholder: Modify based on replicates
                coverage: None,
            }),
        };

        region_map.entry(chr).or_default().push(interval);
    }

    Ok(region_map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_parse_bed6() {
        let bed_path = Path::new("tests/data/synthetic_genome_promoters.bed");

        // Parse the BED file
        let region_map = parse_bed6(bed_path).expect("Failed to parse BED file");

        // Check that at least 3 chromosomes exist in the parsed data
        assert!(region_map.len() >= 3, "Expected at least 3 chromosomes");

        // Expected chromosome names
        let expected_chromosomes = vec!["chr1", "chr2", "chr3"];
        for chr in &expected_chromosomes {
            assert!(
                region_map.contains_key(*chr),
                "Missing expected chromosome: {}",
                chr
            );
        }

        // Check that each chromosome has at least 3 regions
        for chr in expected_chromosomes {
            let intervals = region_map.get(chr).unwrap();
            assert!(
                intervals.len() >= 2 && intervals.len() <= 3,
                "Chromosome {} should have at least 3 regions, found {}",
                chr,
                intervals.len()
            );

            // Ensure correct parsing of coordinates
            for interval in intervals {
                assert!(
                    interval.start < interval.end,
                    "Invalid interval detected on {}: {}-{}",
                    interval.chr,
                    interval.start,
                    interval.end
                );
            }
        }

        println!("test_parse_bed6 passed!");
    }
}
