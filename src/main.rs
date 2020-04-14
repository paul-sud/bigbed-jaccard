use bigbed_jaccard::{
    compute_k_minhashes, get_chrom_end, get_data, jaccard, BigBedReader, IntervalQuery,
};

use bigtools::bigbedread::BigBedRead;
use bigtools::remote_file::RemoteFile;

use itertools::Itertools;

use std::time::Instant;

use rayon::prelude::*;

const QUEUE_SIZE: usize = 100;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let domain = "https://encode-public.s3.amazonaws.com";
    let bigbed_paths = [
        "/2020/01/17/7d2573b1-86f4-4592-a68a-ac3d5d0372d6/ENCFF592UJG.bigBed",
        "/2016/11/13/87d48fc0-c813-4d07-a55b-b60ccde13b25/ENCFF187BSA.bigBed",
        // "/2019/09/07/52581224-261e-4812-a345-0939fdaa84e6/ENCFF512NUN.bigBed"
    ];
    let minhashes = bigbed_paths
        .par_iter()
        .map(|bigbed_path| {
            // This should really be a library function
            let bigbed = RemoteFile::new(&format!("{}{}", domain, bigbed_path));
            let mut reader = BigBedRead::from(bigbed).unwrap();
            let chrom_end = get_chrom_end(&mut reader, "chr1");
            let mut start = Instant::now();
            let data = get_data(
                &mut reader,
                &IntervalQuery::new("chr1".to_string(), 0, chrom_end),
            );
            let mut duration = start.elapsed();
            println!("Time elapsed in getting_data() is: {:?}", duration);
            start = Instant::now();
            let mut minhashes = compute_k_minhashes(&data, QUEUE_SIZE);
            duration = start.elapsed();
            println!("Time elapsed in compute_k_minhashes() is: {:?}", duration);
            minhashes.shrink_to_queue_size();
            (*bigbed_path, minhashes)
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

    // let bigbed_path = "/2020/01/17/7d2573b1-86f4-4592-a68a-ac3d5d0372d6/ENCFF592UJG.bigBed";
    // let mut reader = BigBedReader::new(format!("{}{}", domain, bigbed_path));
    // let mut reader = BigBedReader::new("ENCFF592UJG.bigBed".to_string());
    // let chrom_end = reader.get_chrom_end("chr1")?;
    // let data = reader.get_data(
    //     &IntervalQuery::new("chr1".to_string(), 1000000, 1100000),
    // );
    // dbg!(data);
    Ok(())
}
