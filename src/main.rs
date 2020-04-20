use bigbed_jaccard::request::download_to_tempfile;
use bigbed_jaccard::{compute_k_minhashes, get_offset_data, jaccard, BoundedPriorityQueue, IntervalQuery};

use bigtools::bbiread::BBIRead;
use bigtools::bigbedread::BigBedRead;

use itertools::Itertools;

use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::string::String;
use std::time::Instant;

use rayon::prelude::*;

const QUEUE_SIZE: usize = 100;

fn main() -> Result<(), Box<dyn Error>> {
    let mut file = File::open("files.txt")?;
    let mut file_contents = String::new();
    file.read_to_string(&mut file_contents)?;
    let bigbed_paths = file_contents.lines().collect::<Vec<_>>();
    let minhashes = bigbed_paths
        .par_iter()
        .map(|&item| {
            genome_wide_minhash(item)
        })
        .collect::<Vec<_>>();

    for pair in minhashes.iter().combinations(2) {
        let start = Instant::now();
        println!(
            "Jaccard of {} and {}: {}",
            pair[0].0,
            pair[1].0,
            jaccard(&pair[0].1, &pair[1].1)
        );
        println!("Time elapsed in jaccard() is: {:?}", start.elapsed());
    }
    Ok(())
}

fn genome_wide_minhash(bigbed_path: &str) -> (&str, BoundedPriorityQueue) {
    let download = download_to_tempfile(bigbed_path).unwrap();
    let mut reader =
        BigBedRead::from_file_and_attach(download.path().to_str().unwrap().to_string())
            .unwrap();
    let chroms = reader.get_chroms();
    let queries = chroms.iter().map(|chrom| IntervalQuery::new(chrom.name.to_string(), 0, chrom.length)).collect::<Vec<_>>();
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
