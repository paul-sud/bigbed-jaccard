"""
Unlike the other script, here we use a single hash function (so only one hash computed
per element), store the k elements with the lowest hashes (priority queue) for each set,
then subsquently compute the jaccard from these signatures (see following for details)
https://en.wikipedia.org/wiki/MinHash#Variant_with_a_single_hash_function. This ends up
being much faster that computing k hashes for every element of the set.

The __future__ annotations buisiness is to avoid complaints when annotating a method
returning a new instance of a class (at least for python 3.7)
"""

from __future__ import annotations
import heapq
from multiprocessing import Pool
import os
import random
from typing import List, Tuple
from timeit import default_timer as timer

import attr
import pyBigWig


# Smallest prime larger than chr1 length (248,956,422 on GRCh38)
CHR_1_PRIME = 248_956_429
CHR_1_SIZE = 248_956_422
# Number of elements to maintain in priority queue
QUEUE_SIZE = 100

# Generate coeffs for hash function
random.seed(30)
A = random.randint(0, CHR_1_SIZE + 1)
B = random.randint(0, CHR_1_SIZE + 1)


def main():
    domain = "https://encode-public.s3.amazonaws.com"
    urls = [
        f"{domain}/2020/01/17/7d2573b1-86f4-4592-a68a-ac3d5d0372d6/ENCFF592UJG.bigBed",
        f"{domain}/2016/11/13/87d48fc0-c813-4d07-a55b-b60ccde13b25/ENCFF187BSA.bigBed",
    ]
    with Pool(processes=2) as pool:
        sig_matrices = pool.map(get_all, urls, chunksize=1)
    print(jaccard(*sig_matrices))


@attr.s(auto_attribs=True)
class PriorityQueue:
    """
    Bounded min queue implementation. We bound it so that the memory footprint remains
    small (queue size < 1000 for our application). However, checking the length is slow.

    Note that we assume two elements will not have the same priority, i.e. no hash
    collisions between the minhashes.
    """
    queue_size: int
    heap: List[Tuple[int, int]] = []
    queue_at_capacity: bool = False

    def __len__(self):
        return len(self.heap)

    def shrink_to_queue_size(self):
        """
        Take the queue_size number of greatest items from the heap. We only need this
        when merging two heaps together.
        """
        self.heap = sorted(self.heap)[-(self.queue_size + 1):-1]
        self.queue_at_capacity = True

    def push(self, elem: Tuple[int, int]):
        """
        We add the new element to the queue, and if we exceed the queue_size we need to
        pop an element, where the element is a tuple of (hash(item), item). Note that
        the first item in the order is popped (minimum hash). Comparing the length to
        the queue size each time was slow, so instead we use a flag. Peformance now
        bounded by heapq.heappushpop
        """
        if self.queue_at_capacity:
            heapq.heappushpop(self.heap, elem)
        else:
            heapq.heappush(self.heap, elem)
            if len(self) == self.queue_size:
                self.queue_at_capacity = True

    def merge(self, other: PriorityQueue) -> PriorityQueue:
        """
        Merge two priority queues and return the merged queue, preserving the original.
        If two values are the same, they are deduped (set intersection).
        """
        new_heap = sorted(list(set(heapq.merge(self.heap, other.heap))))
        merged_queue = PriorityQueue(self.queue_size, new_heap)
        merged_queue.shrink_to_queue_size()
        return merged_queue

    def as_set(self):
        """
        Return a set containing the elements of the heap without the priorities
        """
        return {elem for _, elem in self.heap}


def get_all(url: str) -> PriorityQueue:
    start = timer()
    data = get_data(url)
    end = timer()
    print(f"Time getting data in process {os.getpid()}: {end - start}")
    start = timer()
    minhashes = compute_k_minhashes(data, QUEUE_SIZE)
    end = timer()
    print(f"Time computing in process {os.getpid()}: {end - start}")
    return minhashes


def jaccard(sig_a: PriorityQueue, sig_b: PriorityQueue) -> float:
    """
    We approximate the jaccard by finding the elements of the union with the k lowest
    hashes, then taking the union of the elements of that set, the elements of sig_a,
    and the elements of sig_b.
    """
    merged = sig_a.merge(sig_b)
    shared = merged.as_set().intersection(sig_a.as_set()).intersection(sig_b.as_set())
    return len(shared) / sig_a.queue_size


def compute_k_minhashes(data: List[Tuple[int, int]], k: int) -> PriorityQueue:
    """
    Compute the k smallest minhashes for elements of a bigbed. We store the elements
    with the lowest values in a priority queue, where the hash is used to denote the
    priorty. When the queue grows beyond the size of k we evict high hashes from the
    queue. Note that unlike when we use multiple hash functions, we store the element
    itself in addition to the hash so that we can select the bottom k hashes from the
    union of the sets later.

    In the minhash algorithm, we skip rows of the characteristic matrix that have a
    value of 0, that is to say the element corresponding to that row is not present.
    With a bigBed, we only recieve the coordinates of elements that are present in that
    set, so we can just iterate through them.

    For each element we need to decide what to do with it. The start positions are in
    sorted order, so the current start will always be >= the previous. The current end
    could be any value greater than current start coordinate though, and potentially
    less than a previous end position. If the current interval is completely contained
    in a previous interval, then we can skip computing the hash entirely since it
    doesn't contain any new elements. If its end position is greater than any previous
    end position we have seen then we only need to compute the hash for the new
    elements (there are two further sub cases in here, where the new start does/does not
    lie in a previous interval).

    We can keep track of the seen positions with a single virtual interval that we keep
    extending if consecutive intervals overlap. If we encounter a new interval that
    doesn't overlap, then we can set the virtual interval equal to the new one and
    forget the old one.

    Returns: A sorted list of (hash, element) pairs
    """
    minhashes = PriorityQueue(queue_size=k)

    # Need to initialize the virtual interval, and compute the hashes for it since it is
    # assumed that hashes have already been calculated for the virtual interval
    virtual_start, virtual_end = data[0]
    for i in range(virtual_start, virtual_end):
        current_hash = -1 * hasher(i)
        minhashes.push((current_hash, i))

    # Now we can iterate over the rest of the data since the virtual interval has had
    # its hashes computed
    for start, end in data[1:]:
        # Case 1: Interval is completely contained in the virtual interval
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
            current_hash = -1 * hasher(i)
            minhashes.push((current_hash, i))

    return minhashes


def hasher(x: int) -> int:
    """
    This is a 2-independent hash, required for error bounds for the bottom-k approach to
    hold. See the following for details. This method is not ideal because we are
    required to know a prime before hand, eventually will want to use right shift
    instead of moduloing the primes.

    `if the elements are 32-bit keys, we pick two random 64-bit numbers a and b. The
    hash of key x is computed with the C-code (a * x + b) >> 32, where ∗ is 64-bit
    multiplication which as usual discards overflow, and >> is a right shift.`

    https://en.wikipedia.org/wiki/K-independent_hashing#Polynomials_with_random_coefficients
    https://arxiv.org/pdf/1303.5479v2.pdf
    """
    return (A * x + B) % CHR_1_PRIME


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
