use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

/// Parses a TSV file containing (chr, start, end, count) into a HashMap.
///
/// # Arguments
/// - `file_path`: Path to the input TSV file.
///
/// # Returns
/// - `Ok(HashMap<(String, u32, u32), u32>)`: A mapping from (chr, start, end) → count.
/// - `Err(String)`: Error message if parsing fails.
pub fn region_count_parser<P: AsRef<Path>>(
    file_path: P,
) -> Result<HashMap<(String, u32, u32), u32>, String> {
    let file = File::open(&file_path).map_err(|e| format!("Failed to open file: {}", e))?;
    let reader = io::BufReader::new(file);
    let mut region_counts: HashMap<(String, u32, u32), u32> = HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading line {}: {}", line_num + 1, e))?;
        let fields: Vec<&str> = line.split_whitespace().collect();

        if fields.len() != 4 {
            return Err(format!(
                "Invalid format at line {}: expected 4 columns, found {}",
                line_num + 1,
                fields.len()
            ));
        }

        let chr = fields[0].to_string();
        let start: u32 = fields[1]
            .parse()
            .map_err(|_| format!("Invalid start at line {}", line_num + 1))?;
        let end: u32 = fields[2]
            .parse()
            .map_err(|_| format!("Invalid end at line {}", line_num + 1))?;
        let count: u32 = fields[3]
            .parse()
            .map_err(|_| format!("Invalid count at line {}", line_num + 1))?;

        // Ensure valid genomic intervals
        if start >= end {
            return Err(format!(
                "Invalid interval at line {}: start ({}) >= end ({})",
                line_num + 1,
                start,
                end
            ));
        }

        region_counts.insert((chr, start, end), count);
    }

    Ok(region_counts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_region_count_parser_valid() {
        let file_path = Path::new("tests/data/treatment_sorted_r1_5p_cov_region_counts.tsv");

        let result = region_count_parser(file_path);
        assert!(result.is_ok(), "Failed to parse valid region counts file");

        let region_counts = result.unwrap();

        // Ensure correct number of regions parsed
        // NOTE: there are 8 regions in the test file, but one is a duplicate.
        assert_eq!(
            region_counts.len(),
            7,
            "Expected 7 regions, found {}",
            region_counts.len()
        );

        // Check specific values
        assert_eq!(
            region_counts.get(&("chr2".to_string(), 1336, 1836)),
            Some(&614)
        );
        assert_eq!(
            region_counts.get(&("chr2".to_string(), 2296, 2796)),
            Some(&632)
        );
        assert_eq!(
            region_counts.get(&("chr2".to_string(), 696, 1196)),
            Some(&617)
        );
        assert_eq!(
            region_counts.get(&("chr3".to_string(), 2949, 3449)),
            Some(&604)
        );
        assert_eq!(
            region_counts.get(&("chr3".to_string(), 3904, 4404)),
            Some(&486)
        );
        assert_eq!(
            region_counts.get(&("chr1".to_string(), 463, 963)),
            Some(&629)
        );
        assert_eq!(
            region_counts.get(&("chr1".to_string(), 1264, 1764)),
            Some(&670)
        );

        println!("test_region_count_parser_valid passed!");
    }

    #[test]
    fn test_region_count_parser_invalid_format() {
        let file_path = Path::new("tests/data/invalid_region_counts.tsv");

        let result = region_count_parser(file_path);
        assert!(result.is_err(), "Expected failure for invalid file format");

        println!("test_region_count_parser_invalid_format passed!");
    }

    #[test]
    fn test_region_count_parser_invalid_numbers() {
        let file_path = Path::new("tests/data/invalid_numbers_region_counts.tsv");

        let result = region_count_parser(file_path);
        assert!(
            result.is_err(),
            "Expected failure due to non-numeric values"
        );

        println!("test_region_count_parser_invalid_numbers passed!");
    }
}
