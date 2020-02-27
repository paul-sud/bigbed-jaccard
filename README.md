# bigbed-jaccard

[![MIT License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Overview

`bigbed-jaccard` is a tool designed to quickly estimate pairwise Jaccard similarities between large numbers of bigBed files by using minhashing. If you only want to compare a couple of bed files, it will probably be faster to download them and compute the Jaccard distances directly using [`bedtools jaccard`](https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html).

See [Details](#details) for more information on the algorithm.

## Installation

Clone this repo, `cd` into it, and run `cargo build --release`. The binary will then be available at `target/release/bigbed-jaccard`

## Usage

TODO: add CLI + docs

## Details

The minhash algorithm variant used here is bottom-k minhashing, using a single 32-bit hash function described in [1].

## References

1. Mikkel Thorup. 2013. Bottom-k and priority sampling, set similarity and subset sums with minimal independence. In Proceedings of the forty-fifth annual ACM symposium on Theory of Computing (STOC ’13). Association for Computing Machinery, New York, NY, USA, 371–380. DOI:https://doi.org/10.1145/2488608.2488655
