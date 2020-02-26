use bigtools::bbiread::BBIRead;
use bigtools::bigbedread::BigBedRead;
use bigtools::bigwig::BedEntry;
use bigtools::seekableread::{Reopen, SeekableRead};

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};
use std::iter::FromIterator;

const A: u64 = 11_927_359_292_700_924_260;
const B: u64 = 6_512_515_406_574_399_413;
const CAPACITY: usize = 1_000_000;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct HashWithValue {
    hash: u32,
    value: u32,
}

// We do not enforce the queue size on push (which requires a push and a pop), since
// both push and pop incur a cost O(log n). Instead we use the queue size to determine
// how many elements to read off. We want this to be a min queue so we use Reverse
pub struct BoundedPriorityQueue {
    queue_size: usize,
    heap: BinaryHeap<Reverse<HashWithValue>>,
}

impl BoundedPriorityQueue {
    fn new(queue_size: usize, capacity: usize) -> Self {
        BoundedPriorityQueue {
            queue_size,
            heap: BinaryHeap::with_capacity(capacity),
        }
    }

    fn push(&mut self, elem: HashWithValue) {
        self.heap.push(Reverse(elem));
    }

    pub fn shrink_to_queue_size(&mut self) {
        let mut new_heap: BinaryHeap<Reverse<HashWithValue>> =
            BinaryHeap::with_capacity(self.queue_size);
        while new_heap.len() < self.queue_size {
            match self.heap.pop() {
                Some(x) => new_heap.push(x),
                None => break,
            }
        }
        self.heap = new_heap;
    }

    fn merge(&self, other: &Self) -> Self {
        let mut merged_heaps = self
            .heap
            .iter()
            .chain(other.heap.iter())
            .copied()
            .collect::<Vec<_>>();
        merged_heaps.sort();
        merged_heaps.dedup();
        let new_heap = BinaryHeap::from_iter(merged_heaps);
        let mut new_queue = BoundedPriorityQueue {
            queue_size: self.queue_size,
            heap: new_heap,
        };
        new_queue.shrink_to_queue_size();
        new_queue
    }
}

#[allow(clippy::implicit_hasher)]
impl From<&BoundedPriorityQueue> for HashSet<Reverse<HashWithValue>> {
    fn from(queue: &BoundedPriorityQueue) -> Self {
        HashSet::from_iter(queue.heap.clone())
    }
}

pub fn compute_k_minhashes(data: &[BedEntry], k: usize) -> BoundedPriorityQueue {
    let mut minhashes = BoundedPriorityQueue::new(k, CAPACITY);
    let BedEntry {
        start: virtual_start,
        end: mut virtual_end,
        ..
    } = data[0];
    for i in virtual_start..virtual_end {
        let current_hash = hash(i);
        minhashes.push(HashWithValue {
            hash: current_hash,
            value: i,
        })
    }
    let (mut start_at, mut end_at);
    for datum in data.iter().skip(1) {
        let BedEntry { start, end, .. } = *datum;
        if end <= virtual_end {
            continue;
        } else if start >= virtual_end {
            start_at = start;
            end_at = end;
            virtual_end = end;
        } else {
            start_at = virtual_end;
            end_at = end;
            virtual_end = end;
        }
        for i in start_at..end_at {
            let current_hash = hash(i);
            minhashes.push(HashWithValue {
                hash: current_hash,
                value: i,
            })
        }
    }
    minhashes
}

pub fn jaccard(sig_a: &BoundedPriorityQueue, sig_b: &BoundedPriorityQueue) -> f64 {
    let merged: HashSet<Reverse<HashWithValue>> = (&sig_a.merge(sig_b)).into();
    let sig_a_set: HashSet<Reverse<HashWithValue>> = sig_a.into();
    let sig_b_set: HashSet<Reverse<HashWithValue>> = sig_b.into();
    let a_intersect_b = HashSet::from_iter(sig_a_set.intersection(&sig_b_set).copied());
    let shared = merged.intersection(&a_intersect_b);
    shared.count() as f64 / sig_a.queue_size as f64
}

// from https://arxiv.org/pdf/1303.5479v2.pdf
fn hash(value: u32) -> u32 {
    (((A.overflowing_mul(value as u64).0).overflowing_add(B).0) >> 32) as u32
}

// TODO: PR fn against bigtools
pub fn get_data<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    reader: &mut BigBedRead<R, S>,
    chrom: &str,
    start: u32,
    end: u32,
) -> Vec<BedEntry> {
    let remote_intervals: Vec<_> = reader
        .get_interval(chrom, start, end)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    remote_intervals
}

// TODO: PR fn against bigtools
pub fn get_chrom_end<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    reader: &mut BigBedRead<R, S>,
    chrom: &str,
) -> u32 {
    let chroms = reader.get_chroms();
    let chrom = chroms.iter().find(|v| v.name == chrom).unwrap();
    chrom.length
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_queue_push() {
        let mut queue = BoundedPriorityQueue::new(2, 3);
        queue.push(HashWithValue { hash: 1, value: 2 });
        queue.push(HashWithValue { hash: 2, value: 4 });
        queue.push(HashWithValue { hash: 4, value: 8 });
        assert_eq!(
            queue.heap.peek(),
            Some(&Reverse(HashWithValue { hash: 1, value: 2 }))
        );
    }

    #[test]
    fn test_queue_merge() {
        let mut queue = BoundedPriorityQueue::new(2, 3);
        queue.push(HashWithValue { hash: 1, value: 2 });
        queue.push(HashWithValue { hash: 5, value: 4 });
        let mut other_queue = BoundedPriorityQueue::new(2, 3);
        other_queue.push(HashWithValue { hash: 2, value: 2 });
        other_queue.push(HashWithValue { hash: 3, value: 2 });
        let merged_queue = queue.merge(&other_queue);
        assert_eq!(
            merged_queue.heap.peek(),
            Some(&Reverse(HashWithValue { hash: 1, value: 2 }))
        );
    }

    #[test]
    fn test_queue_shrink_to_queue_size() {
        let mut queue = BoundedPriorityQueue::new(2, 2);
        queue.push(HashWithValue { hash: 1, value: 2 });
        queue.push(HashWithValue { hash: 5, value: 4 });
        queue.push(HashWithValue { hash: 6, value: 3 });
        queue.shrink_to_queue_size();
        assert_eq!(queue.heap.len(), 2);
    }

    #[test]
    fn test_compute_k_minhashes() {
        let data = [BedEntry {
            start: 1,
            end: 4,
            rest: String::from(""),
        }];
        let queue = compute_k_minhashes(&data, 2);
        assert_eq!(
            queue.heap.peek(),
            Some(&Reverse(HashWithValue {
                hash: 1_257_542_785,
                value: 3
            }))
        );
        assert_eq!(queue.heap.len(), 3);
    }

    // #[test]
    // fn test_jaccard() {
    //     let mut queue = BoundedPriorityQueue::new(2, 3);
    //     queue.push(HashWithValue { hash: 1, value: 2 });
    //     queue.push(HashWithValue { hash: 2, value: 4 });
    //     queue.push(HashWithValue { hash: 4, value: 8 });
    //     let mut other_queue = BoundedPriorityQueue::new(2, 3);
    //     other_queue.push(HashWithValue { hash: 3, value: 8 });
    //     other_queue.push(HashWithValue { hash: 6, value: 9 });
    //     other_queue.push(HashWithValue { hash: 7, value: 10 });
    //     let result = jaccard(&queue, &other_queue);
    //     assert_eq!(result, 0_f64);
    // }

    #[test]
    fn test_hash() {
        let result = hash(100);
        assert_eq!(result, 48_913_034);
    }
}
