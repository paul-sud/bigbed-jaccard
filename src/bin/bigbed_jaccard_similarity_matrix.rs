use bigbed_jaccard::bed::{get_offset_data, IntervalQuery};
use bigbed_jaccard::bottom_k::{compute_k_minhashes, jaccard, BoundedPriorityQueue};
use bigbed_jaccard::request::download_to_tempfile;

use bigtools::bbiread::BBIRead;
use bigtools::bigbedread::BigBedRead;

use csv::{Writer, WriterBuilder};

use itertools::Itertools;

use serde::Serialize;

use std::error::Error;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;
use std::string::String;
use std::time::Instant;

use structopt::StructOpt;

use rayon::prelude::*;

const QUEUE_SIZE: usize = 100;

#[derive(Debug, Serialize)]
struct JaccardResult<'a> {
    id1: &'a str,
    id2: &'a str,
    jaccard: f64,
}

#[derive(StructOpt)]
struct Cli {
    #[structopt(parse(from_os_str))]
    infile: std::path::PathBuf,
    #[structopt(parse(from_os_str))]
    outfile: std::path::PathBuf,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::from_args();
    let mut file = File::open(&args.infile)?;
    let mut file_contents = String::new();
    file.read_to_string(&mut file_contents)?;
    let bigbed_paths = file_contents.lines().collect::<Vec<_>>();
    let minhashes = bigbed_paths
        .par_iter()
        .map(|&item| genome_wide_minhash(item))
        .collect::<Vec<_>>();

    let outfile = File::create(&args.outfile)?;
    let mut writer = get_writer(&outfile);
    // Note no replacement, so we don't compute jaccard of diagonal which is always 1
    for pair in minhashes.iter().combinations(2) {
        let start = Instant::now();
        let jaccard = jaccard(&pair[0].1, &pair[1].1);
        println!("Jaccard of {} and {}: {}", pair[0].0, pair[1].0, &jaccard,);
        println!("Time elapsed in jaccard() is: {:?}", start.elapsed());
        writer.serialize(JaccardResult {
            id1: Path::new(pair[0].0).file_stem().unwrap().to_str().unwrap(),
            id2: Path::new(pair[1].0).file_stem().unwrap().to_str().unwrap(),
            jaccard,
        })?;
    }
    Ok(())
}

fn genome_wide_minhash(bigbed_path: &str) -> (&str, BoundedPriorityQueue) {
    let download = download_to_tempfile(bigbed_path).unwrap();
    let mut reader =
        BigBedRead::from_file_and_attach(download.path().to_str().unwrap().to_string()).unwrap();
    let chroms = reader.get_chroms();
    let queries = chroms
        .iter()
        .map(|chrom| IntervalQuery::new(chrom.name.to_string(), 0, chrom.length))
        .collect::<Vec<_>>();
    let mut start = Instant::now();
    let data = get_offset_data(&mut reader, &queries).unwrap();
    let mut duration = start.elapsed();
    println!("Time elapsed in getting_data() is: {:?}", duration);
    start = Instant::now();
    let mut minhashes = compute_k_minhashes(&data, QUEUE_SIZE);
    duration = start.elapsed();
    println!("Time elapsed in compute_k_minhashes() is: {:?}", duration);
    minhashes.shrink_to_queue_size();
    (bigbed_path, minhashes)
}

fn get_writer<W: Write>(wtr: W) -> Writer<W> {
    WriterBuilder::new().delimiter(b'\t').from_writer(wtr)
}
