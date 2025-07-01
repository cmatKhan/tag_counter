use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

/// A parser for Bedtools genomecov -5 -dz output.
///
/// This struct implements an **iterator** over a BEDTools genomecov file.
/// Each iteration produces a tuple: `(chrom, start, start+1, count)`.
pub struct TagCountParser<R: BufRead> {
    reader: R,
}

impl TagCountParser<io::BufReader<File>> {
    /// Creates a new `TagCountParser` from a file path.
    ///
    /// # Arguments
    /// * `path` - Path to the input genomecov file.
    ///
    /// # Returns
    /// * `Ok(TagCountParser)` if the file is opened successfully.
    /// * `Err(String)` if there is an issue opening the file.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let file = File::open(&path).map_err(|e| {
            let path = path.as_ref();
            format!("Failed to open file {}: {}", path.display(), e)
        })?;
        let reader = io::BufReader::new(file);
        Ok(TagCountParser { reader }) // Ensure type is inferred properly
    }
}

impl<R: BufRead> Iterator for TagCountParser<R> {
    type Item = Result<(String, u32, u32, u32), String>;

    /// Parses each line lazily and returns `(chr, start, start+1, count)`.
    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();

        // Read the next line
        match self.reader.read_line(&mut line) {
            Ok(0) => None, // EOF
            Ok(_) => {
                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.len() != 3 {
                    return Some(Err(format!("Invalid line format: {}", line.trim())));
                }

                let chr = fields[0].to_string();
                let start = match fields[1].parse::<u32>() {
                    Ok(val) => val,
                    Err(_) => return Some(Err(format!("Invalid start value: {}", fields[1]))),
                };
                let count = match fields[2].parse::<u32>() {
                    Ok(val) => val,
                    Err(_) => return Some(Err(format!("Invalid count value: {}", fields[2]))),
                };

                Some(Ok((chr, start, start + 1, count))) // Return as expected tuple
            }
            Err(e) => Some(Err(format!("Error reading file: {}", e))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_parse_tag_count() {
        let path = Path::new("tests/data/control_sorted_r1_5p_cov.txt");

        // Initialize parser
        let parser = TagCountParser::from_path(path).expect("Failed to open test file");

        let expected_first_five = vec![
            ("chr1".to_string(), 249, 250, 1),
            ("chr1".to_string(), 253, 254, 1),
            ("chr1".to_string(), 255, 256, 1),
            ("chr1".to_string(), 256, 257, 1),
            ("chr1".to_string(), 257, 258, 1),
        ];

        let expected_last_five = vec![
            ("chr3".to_string(), 4315, 4316, 2),
            ("chr3".to_string(), 4316, 4317, 1),
            ("chr3".to_string(), 4317, 4318, 1),
            ("chr3".to_string(), 4318, 4319, 1),
            ("chr3".to_string(), 4320, 4321, 2),
        ];

        // Read the first five lines
        let parsed_results = parser.take(5).collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(parsed_results, expected_first_five);

        // Reopen the parser to read the entire file again to get the last 5 lines
        let parser = TagCountParser::from_path(path).expect("Failed to reopen test file");

        // Read all lines and collect the last 5
        let last_five = parser
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
            .into_iter()
            .rev()
            .take(5)
            .collect::<Vec<_>>();

        // Reverse again to maintain the original order
        let last_five: Vec<_> = last_five.into_iter().rev().collect();

        assert_eq!(last_five, expected_last_five);
    }
}
