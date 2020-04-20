# bigbed-jaccard

[![MIT License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Overview

`bigbed-jaccard` is a tool designed to quickly estimate pairwise Jaccard similarities between large numbers of bigBed files by using minhashing. If you only want to compare a couple of bed files, it will probably be faster to download them and compute the Jaccard distances directly using [`bedtools jaccard`](https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html).

See [Details](#details) for more information on the algorithm.

## Installation

Clone this repo, `cd` into it, and run `cargo build --release`. The binary will then be available at `target/release/bigbed-jaccard`

This will require that you have `cargo` installed, I would recommend using [rustup](rustup.rs) to install the Rust toolchain. Build has been tested on `rustc` version `1.42.0`

## Usage

Here, `infile` is simply a file containing a list of file URLs, one per line, and `outfile` is the name of the file to write the comparisons to.

```bash
$ target/release/bigbed-jaccard [infile] [outfile]
```

There is some utility code in this repository to facilitate data aquisition from the [ENCODE portal](https://www.encodeproject.org). Namely, [this script](scripts/get_urls_from_report.py) will allow you to flatten ENCODE file report TSVs containing JSON objects in the `cloud_metadata` field, requires `pandas>=1.0`. You can then convert the output of that script into a file list by extracting the `url` field and deleting the header.

## Benchmarks

Running pairwise comparisons on the 2,677 bigBed files in [this file](files.txt) took ~13 minutes on a Mac with 2.6 GHz 6-Core Intel i7 processor and an SSD, resulting in 3,581,827 output comparisons. Memory usage hovered around 1.5 GB.

## Details

The minhash algorithm variant used here is bottom-_k_ minhashing, also known as one-permutation hashing (OPH), using a single 32-bit hash function described in [1].

## Future Work

I'd like to enable the ability to specify regions of interest for which to compute the jaccard. The code should determine based on the number of query intervals and the query range whether to download to a tempfile (currently hardwired behavior) or to query the bigBed remotely. More benchmarking is needed to see how those decisions should be made.

I'd also like to make the underlying code available as a Rust library and also generate Python (possibly NumPy) bindings for the Rust functions.

Once this is accomplished, it should be more straightforward to interoperate with Jupyter notebooks for further analysis like clustering, multi-dimensional scaling (MDS), dimensionality reduction, etc.

## References

1. Mikkel Thorup. 2013. Bottom-k and priority sampling, set similarity and subset sums with minimal independence. In Proceedings of the forty-fifth annual ACM symposium on Theory of Computing (STOC ’13). Association for Computing Machinery, New York, NY, USA, 371–380. DOI:https://doi.org/10.1145/2488608.2488655
