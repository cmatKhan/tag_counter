use clap::Parser;
use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use tag_counter::parsers::parse_bed6::parse_bed6;
use tag_counter::parsers::region_count_parser::region_count_parser;
use tag_counter::parsers::tag_count_parser::TagCountParser;
use tag_counter::quantification::region_stats::RegionStats;
use tag_counter::{CountableRegionMetadata, CountableRegionTree, GenomicInterval};

/// Command-line arguments parser
#[derive(Parser, Debug)]
#[command(about = "Process tag count data and compute region counts")]
struct Args {
    /// Path to the BED file of regions of interest
    #[arg(short, long)]
    bed_file: PathBuf,

    /// Path to the JSON file specifying tag count file paths
    #[arg(short, long)]
    json_file: PathBuf,

    /// Output directory
    #[arg(short, long, default_value = "tag_count_output")]
    output_dir: PathBuf,

    /// Path to background region count file
    #[arg(long)]
    background_counts: Option<PathBuf>,

    /// Total number of background tags
    #[arg(long)]
    background_tag_total: Option<u32>,
}

/// Load JSON file and parse it into a HashMap<String, Vec<String>>
fn load_json<P: AsRef<Path>>(json_path: P) -> Result<HashMap<String, Vec<String>>, Box<dyn Error>> {
    let file = File::open(json_path)?;
    let data: HashMap<String, Vec<String>> = serde_json::from_reader(file)?;
    Ok(data)
}

/// Generate output file name by replacing the suffix and appending "_region_counts.tsv"
fn generate_output_filename(input_path: &str) -> String {
    let path = Path::new(input_path);
    let filename = path.file_stem().unwrap().to_string_lossy();
    format!("{}_region_counts.tsv", filename.replace("_r1_5p_cov", ""))
}

/// Write the count data for each replicate and optionally the combined counts if applicable.
fn write_counts(
    output_dir: &Path,
    replicate_paths: &[String],
    region_map: &HashMap<String, Vec<GenomicInterval<CountableRegionMetadata>>>,
    tree: &CountableRegionTree,
    total_counts_map: &HashMap<String, u64>,
    background_counts: Option<&HashMap<(String, u32, u32), u32>>,
    background_tag_total: Option<u32>,
) -> Result<(), io::Error> {
    const PSEUDOCOUNT: f64 = 0.1;

    // Create output file handles for each replicate
    let mut output_files: Vec<(String, File)> = replicate_paths
        .iter()
        .map(|path| {
            let output_filename = generate_output_filename(path);
            File::create(output_dir.join(&output_filename)).map(|file| (path.to_string(), file))
            // Store filename along with File
        })
        .collect::<Result<Vec<_>, _>>()?;

    // If applicable, add a file handle for the combined output
    let mut combined_file = if replicate_paths.len() > 1 {
        Some(File::create(output_dir.join("combined_region_counts.tsv"))?)
    } else {
        None
    };

    // iterate over the region_map (ie promoter regions)
    for (chr, regions) in region_map {
        for region in regions {
            let counts = match tree.get_counts(chr, region.start, region.end) {
                Ok(c) => c,
                Err(e) => {
                    eprintln!("Warning: {}", e);
                    continue;
                }
            };

            let background_count = background_counts
                .and_then(|bc| bc.get(&(chr.clone(), region.start, region.end)))
                .copied()
                .unwrap_or(0);

            // NOTE: output_files has the same length and order as the array of counts
            // from get_counts() above, so we can safely zip them together into
            // (filename, file, count) for each replicate
            for ((filename, file), count) in output_files.iter_mut().zip(counts.iter()) {
                // if background_counts and background_tag_total are both some, then
                // instantiate a RegionStats struct with the background_tag_total and
                // treatment_tag_total with a pseudocount of 0.1
                // if background_counts and background_tag_total are both none, then
                // do not instantiate a RegionStats struct
                let region_stats =
                    background_tag_total
                        .zip(background_counts)
                        .map(|(bg_total, _)| {
                            RegionStats::new(
                                bg_total,
                                total_counts_map.get(filename).copied().unwrap_or(0) as u32,
                                PSEUDOCOUNT,
                            )
                        });

                // if region_stats is some, then compute the enrichment and poisson p-value
                // and write them to the file
                if let Some(stats) = region_stats {
                    let enrichment_score = stats.enrichment(background_count, *count);
                    let poisson_p = stats.poisson_pval(background_count, *count);

                    // if enrichment_score is an error, then set it to 0.0
                    // if poisson_p is an error, then set it to 1.0
                    let enrichment_score = match enrichment_score {
                        Ok(score) => score,
                        Err(_) => f64::NAN, // Now using NaN instead of 0.0
                    };
                    let poisson_p = match poisson_p {
                        Ok(p) => p,
                        Err(_) => f64::NAN, // Now using NaN instead of 1.0
                    };

                    writeln!(
                        file,
                        "{}\t{}\t{}\t{}\t{}\t{:.15}\t{:.15}",
                        chr,
                        region.start,
                        region.end,
                        count,
                        background_count,
                        enrichment_score,
                        poisson_p
                    )?;
                } else {
                    writeln!(file, "{}\t{}\t{}\t{}", chr, region.start, region.end, count)?;
                }
            }

            // Write combined counts if applicable
            if let Some(file) = combined_file.as_mut() {
                let combined_count: u32 = counts.iter().sum();

                let region_stats =
                    background_tag_total
                        .zip(background_counts)
                        .map(|(bg_total, _)| {
                            RegionStats::new(
                                bg_total,
                                total_counts_map
                                    .values()
                                    .copied()
                                    .sum::<u64>()
                                    .try_into()
                                    .expect(
                                        "Error: Total treatment tags exceed maximum value for u32.",
                                    ),
                                PSEUDOCOUNT,
                            )
                        });

                if let Some(stats) = region_stats {
                    let enrichment_score = stats.enrichment(background_count, combined_count);
                    let poisson_p = stats.poisson_pval(background_count, combined_count);

                    // if enrichment_score is an error, then set it to 0.0
                    // if poisson_p is an error, then set it to 1.0
                    let enrichment_score = match enrichment_score {
                        Ok(score) => score,
                        Err(_) => f64::NAN, // Now using NaN instead of 0.0
                    };
                    let poisson_p = match poisson_p {
                        Ok(p) => p,
                        Err(_) => f64::NAN, // Now using NaN instead of 1.0
                    };
                    writeln!(
                        file,
                        "{}\t{}\t{}\t{}\t{}\t{:.15}\t{:.15}",
                        chr,
                        region.start,
                        region.end,
                        combined_count,
                        background_count,
                        enrichment_score,
                        poisson_p
                    )?;
                } else {
                    writeln!(
                        file,
                        "{}\t{}\t{}\t{}",
                        chr, region.start, region.end, combined_count
                    )?;
                }
            }
        }
    }

    Ok(())
}

/// Write total counts per replicate and combined to a file inside the TF-specific subdirectory
fn write_total_counts(
    output_dir: &Path,
    totals: &HashMap<String, u64>,
    region_total: &HashMap<String, u32>,
) -> io::Result<()> {
    let output_path = output_dir.join("total_tag_counts.tsv");
    let mut file = File::create(output_path)?;

    for (filename, total) in totals {
        let region_total = region_total.get(filename).copied().unwrap_or(0);
        let output_filename = generate_output_filename(filename);
        writeln!(file, "{}\t{}\t{}", output_filename, total, region_total)?;
    }

    // if the number in totals is greater than 1, then add another line to the file
    // which is the sum of all the totals
    if totals.len() > 1 {
        let sum_total: u64 = totals.values().sum();
        let sum_region_total: u32 = region_total.values().sum();
        writeln!(file, "combined\t{}\t{}", sum_total, sum_region_total)?;
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse command-line arguments
    let args = Args::parse();

    if args.background_counts.is_some() != args.background_tag_total.is_some() {
        eprintln!(
            "Error: Both --background_counts and --background_tag_total must be provided together."
        );
        std::process::exit(1);
    }

    let background_counts = if let Some(background_file) = &args.background_counts {
        Some(region_count_parser(background_file)?)
    } else {
        None
    };

    let background_tag_total = args.background_tag_total;

    // Parse BED file to create the region_map
    let region_map = parse_bed6(&args.bed_file)?;

    // Load the JSON file containing TF -> replicate count file mappings
    let tf_replicates = load_json(&args.json_file)?;

    // Create output directory if not exists
    fs::create_dir_all(&args.output_dir)?;

    for (tf, replicate_paths) in tf_replicates.iter() {
        let tf_output_dir = args.output_dir.join(tf);
        fs::create_dir_all(&tf_output_dir)?;
        let mut total_counts_map: HashMap<String, u64> = HashMap::new();
        let mut region_counts_map: HashMap<String, u32> = HashMap::new();

        // Create a new CountableRegionTree for this TF
        let mut tree = CountableRegionTree::new(replicate_paths.len() as u8);
        tree.construct(&region_map);

        for (chr, tree) in &tree.trees {
            let mut intervals: Vec<_> = tree
                .find(0..u32::MAX)
                .map(|entry| entry.interval().clone())
                .collect();

            intervals.sort_by_key(|iv| iv.start);

            for window in intervals.windows(2) {
                let a = &window[0];
                let b = &window[1];
                assert!(
                    a.end <= b.start,
                    "Overlapping canonical intervals in tree for {chr}: {:?} and {:?}",
                    a,
                    b
                );
            }
        }

        // Iterate over each replicate and add counts
        // for (rep_index, path) in replicate_paths.iter().enumerate() {
        //     let parser = TagCountParser::from_path(path)?;
        //     total_counts_map.insert(path.clone(), 0);

        //     for record in parser {
        //         let (chr, start, end, count) = record?;
        //         match tree.add_counts(&chr, start, end, count, rep_index as u8) {
        //             Ok(_) => *total_counts_map.entry(path.to_string()).or_insert(0) += count as u64,
        //             Err(e) => {
        //                 eprintln!("Error: {}. This count is entirely discarded. Count cannot overlap multiple countable intervals.", e);
        //                 continue;
        //             }
        //         }
        //     }
        // }
        for (rep_index, path) in replicate_paths.iter().enumerate() {
            let parser = TagCountParser::from_path(path)?;
            total_counts_map.insert(path.clone(), 0);
            region_counts_map.insert(path.clone(), 0);

            for record in parser {
                let (chr, start, end, count) = record?;

                // Always increment the total count

                match tree.add_counts(&chr, start, end, count, rep_index as u8) {
                    Ok(added) if added > 0 => {
                        *total_counts_map.entry(path.clone()).or_insert(0) += count as u64;
                        *region_counts_map.entry(path.clone()).or_insert(0) += added;
                    }
                    Ok(_) => {
                        *total_counts_map.entry(path.clone()).or_insert(0) += count as u64;
                    }
                    Err(e) => {
                        eprintln!("Error: {}. This count is entirely discarded. Count cannot overlap multiple countable intervals.", e);
                    }
                }
            }
        }

        // // zip replicate_paths together with tree.get_total_counts() to region_totals_by_replicate
        // let region_totals_by_replicate: HashMap<String, u32> = replicate_paths
        //     .iter()
        //     .cloned()
        //     .zip(tree.get_total_counts())
        //     .collect();

        // Write replicate and combined counts
        write_counts(
            &tf_output_dir,
            replicate_paths,
            &region_map,
            &tree,
            &total_counts_map,
            background_counts.as_ref(),
            background_tag_total,
        )?;

        // Write total tag counts for all replicates
        write_total_counts(&tf_output_dir, &total_counts_map, &region_counts_map)?;
    }

    println!(
        "Processing complete. Output stored in {:?}",
        args.output_dir
    );
    Ok(())
}
