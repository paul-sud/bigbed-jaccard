use bigbed_jaccard::bed::{get_offset_data, IntervalQuery};
use bigbed_jaccard::lsh::Lsh;
use bigbed_jaccard::oph::OnePermutationHasher;
use bigbed_jaccard::request::download_to_tempfile;

use bigtools::bbiread::BBIRead;
use bigtools::bigbedread::BigBedRead;

use rayon::prelude::*;

use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::string::String;
use std::time::Instant;

use structopt::StructOpt;

const OPH_NUM_BINS: usize = 100;

#[derive(StructOpt)]
struct Cli {
    #[structopt(parse(from_os_str))]
    infile: std::path::PathBuf,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::from_args();
    let mut file = File::open(&args.infile)?;
    let mut file_contents = String::new();
    file.read_to_string(&mut file_contents)?;
    let bigbed_paths = file_contents.lines().collect::<Vec<_>>();
    let sketches = bigbed_paths
        .par_iter()
        .map(|&item| genome_wide_sketch(item))
        .collect::<Vec<_>>();

    let mut lsh: Lsh<u32, &str> = Lsh::new(10_usize, 10_usize);
    let mut start = Instant::now();
    for (path, sketch) in sketches.iter() {
        lsh.insert(sketch, path);
    }
    let mut duration = start.elapsed();
    println!("Time inserting into LSH is: {:?}", duration);
    start = Instant::now();
    let query_results = lsh.query(&sketches[0].1).unwrap();
    duration = start.elapsed();
    println!("Time querying LSH is: {:?}", duration);
    print!("Query results are {:?}", query_results);
    Ok(())
}

fn genome_wide_sketch(bigbed_path: &str) -> (&str, Vec<u32>) {
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
    let mut oph = OnePermutationHasher::new(OPH_NUM_BINS);
    start = Instant::now();
    let minhashes = oph.dense_sketch(&data).unwrap();
    duration = start.elapsed();
    println!("Time elapsed in oph.dense_sketch() is: {:?}", duration);
    (bigbed_path, minhashes)
}
