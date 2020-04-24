# bigbed-jaccard

[![MIT License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Overview

`bigbed-jaccard` is a tool designed to quickly estimate pairwise Jaccard similarities between large numbers of bigBed files by using minhashing. If you only want to compare a couple of bed files, it will probably be faster to download them and compute the Jaccard distances directly using [`bedtools jaccard`](https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html).

See [Details](#details) for more information on the algorithm.

## Installation

Clone this repo, `cd` into it, and run `cargo build --release`. The binary will then be available at `target/release/bigbed-jaccard`

This will require that you have `cargo` installed, I would recommend using [rustup](https://rustup.rs) to install the Rust toolchain. Build has been tested on `rustc` version `1.42.0`

## Usage

Here, `infile` is simply a file containing a list of file URLs, one per line, and `outfile` is the name of the file to write the comparisons to.

```bash
$ target/release/bigbed-jaccard [infile] [outfile]
```

There is some utility code in this repository to facilitate data aquisition from the [ENCODE portal](https://www.encodeproject.org). Namely, [this script](scripts/get_urls_from_report.py) will allow you to flatten ENCODE file report TSVs containing JSON objects in the `cloud_metadata` field, requires `pandas>=1.0`. You can then convert the output of that script into a file list by extracting the `url` field and deleting the header.

## Benchmarks

Running pairwise comparisons on the 2,677 bigBed files in [this file](files.txt) took ~13 minutes on a Mac with 2.6 GHz 6-Core Intel i7 processor and an SSD, resulting in 3,581,827 output comparisons. Memory usage hovered around 1.5 GB.

## Details

The minhash algorithm variant used here is bottom-_k_ minhashing, using a single 32-bit hash function described in [1].

## Future Work

I'd like to make the underlying code available as a Rust library and also generate Python (possibly NumPy) bindings for the Rust functions.

## References

1. Mikkel Thorup. 2013. Bottom-k and priority sampling, set similarity and subset sums with minimal independence. In Proceedings of the forty-fifth annual ACM symposium on Theory of Computing (STOC ’13). Association for Computing Machinery, New York, NY, USA, 371–380.

2. Anshumali Shrivastava and Ping Li. 2014. Improved densification of one permutation hashing. In Proceedings of the Thirtieth Conference on Uncertainty in Artificial Intelligence (UAI’14). AUAI Press, Arlington, Virginia, USA, 732–741.
