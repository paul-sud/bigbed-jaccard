"""
Resources:
Great explanantion of minhashing with worked examples:
http://infolab.stanford.edu/~ullman/mmds/ch3.pdf

High level overview of using many hash functions, Python code example:
https://mccormickml.com/2015/06/12/minhash-tutorial-with-python-code/

This is the approach using multiple hash functions, and is slow due to large number of
elements in set (number of genomic positions, in the hundreds of millions). Instead, we
probably want to use a single hash function (so only one hash computed per element)
https://en.wikipedia.org/wiki/MinHash#Variant_with_a_single_hash_function
"""

import pyBigWig
from multiprocessing import Pool
from typing import List, Tuple
import random

# Smallest prime larger than chr1 length (248,956,422 on GRCh38)
CHR_1_PRIME = 248_956_429
CHR_1_SIZE = 248_956_422
NUM_HASHES = 100


random.seed(30)


def main():
    domain = "https://encode-public.s3.amazonaws.com"
    urls = [
        f"{domain}/2020/01/17/7d2573b1-86f4-4592-a68a-ac3d5d0372d6/ENCFF592UJG.bigBed",
        f"{domain}/2016/11/13/87d48fc0-c813-4d07-a55b-b60ccde13b25/ENCFF187BSA.bigBed",
    ]
    coeffs = [get_coeffs(NUM_HASHES, CHR_1_SIZE + 1)] * 2
    with Pool(processes=2) as pool:
        sig_matrices = pool.starmap(get_all, zip(urls, coeffs), chunksize=1)
    print(jaccard(*sig_matrices))


def get_all(url: str, coeffs: List[Tuple[int, int]]) -> List[int]:
    data = get_data(url)
    sig_matrix = compute_sig_matrix(data, coeffs)
    return sig_matrix


def jaccard(sig_a: List[int], sig_b: List[int]) -> float:
    """
    We approximate the jaccard by computing the fraction of rows in the signature matrix
    that are the same.
    """
    same = 0
    total = len(sig_a)
    for a, b in zip(sig_a, sig_b):
        if a == b:
            same += 1
    return same / total


def compute_sig_matrix(
    data: List[Tuple[int, int]],
    coeffs: List[Tuple[int, int]]
) -> List[int]:
    """
    Compute the signature matrix for a bigbed (minhash). For each hash function defined
    by coeff pair, check if the hash of the current entry (genomic position) is smaller
    than the current entry in the signature vector, and update the signature if it is.

    In the minhash algorithm, we skip rows of the characteristic matrix that have a
    value of 0, that is to say the element corresponding to that row is not present.
    With a bigBed, we only recieve the coordinates of elements that are present in that
    set, so we can just iterate through them.

    For each element we need to decide what to do with it. The start positions are in
    sorted order, so the current start will always be >= the previous. The current end
    could be any value greater than current start coordinate though, and potentially
    less than a previous end position. If the current interval is completely contained
    in a previous interval, then we can skip computing the hashes entirely since it
    doesn't contain any new elements. If its end position is greater than any previous
    end position we have seen then we only need to compute the hashes for the new
    elements (there are two further sub cases in here, where the new start does/does not
    lie in a previous interval).

    We can keep track of the seen positions with a single virtual interval that we keep
    extending if consecutive intervals overlap. If we encounter a new interval that
    doesn't overlap, then we can set the virtual interval equal to the new one and
    forget the old one.
    """

    # Initialize signature vector to value greater than the max of the hash functions.
    sig = [CHR_1_PRIME + 1 for i in range(len(coeffs))]

    # Need to initialize the virtual interval, and compute the hashes for it since it is
    # assumed that hashes have already been calculated for the virtual interval
    virtual_start, virtual_end = data[0]
    for i in range(virtual_start, virtual_end):
        for j, (coeff_a, coeff_b) in enumerate(coeffs):
            current_hash = (coeff_a * i + coeff_b) % CHR_1_PRIME
            if current_hash < sig[j]:
                sig[j] = current_hash

    # Now we can iterate over the rest of the data since the virtual interval has had
    # its hashes computed
    for start, end in data[1:]:
        # Case 1: Interval is completely contained in the virutal interval
        if end <= virtual_end:
            continue
        # Case 2: The interval lies completely outside the virtual interval
        elif start >= virtual_end:
            start_at, end_at = start, end
            virtual_start, virtual_end = start, end
        # Case 3: The start of the interval lies within the virtual interval, but the
        # end lies outside. Need to extend the virtual interval here.
        else:
            start_at, end_at = virtual_end, end
            virtual_end = end

        # Compute the hashes of the elements
        for i in range(start_at, end_at):
            for j, (coeff_a, coeffb) in enumerate(coeffs):
                current_hash = (coeff_a * i + coeff_b) % CHR_1_PRIME
                if current_hash < sig[j]:
                    sig[j] = current_hash
    return sig


def get_coeffs(k: int, max_id: int) -> List[Tuple[int, int]]:
    """
    Get (a, b) coefficients for hash function f(x) = (ax + b) % N
    """
    coeffs_a = pick_random_coeffs(k, max_id)
    coeffs_b = pick_random_coeffs(k, max_id)
    return list(zip(coeffs_a, coeffs_b))


def pick_random_coeffs(k: int, max_id: int) -> List[int]:
    # Create a list of 'k' random values.
    coeffs: List[int] = []
    while k > 0:
        # Get a random shingle ID.
        rand_coeff = random.randint(0, max_id)

        # Ensure that each random number is unique.
        while rand_coeff in coeffs:
            rand_coeff = random.randint(0, max_id)

        # Add the random number to the list.
        coeffs.append(rand_coeff)
        k = k - 1

    return coeffs


def get_data(url: str, chrom: str = "chr1") -> List[Tuple[int, int]]:
    """
    Returns list of [start, stop) feature coordinates. These are always in sorted order
    by the start coordinate, half open, 1-indexed.
    """
    bb = pyBigWig.open(url)
    chroms = bb.chroms()
    return bb.entries(chrom, 1, chroms[chrom], withString=False)


if __name__ == '__main__':
    main()
