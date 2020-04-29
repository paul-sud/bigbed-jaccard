use crate::hash::hash;

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};
use std::iter::FromIterator;

// Min queue capacity
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

pub fn compute_k_minhashes(data: &[u32], k: usize) -> BoundedPriorityQueue {
    let mut minhashes = BoundedPriorityQueue::new(k, CAPACITY);
    for &i in data.iter() {
        let current_hash = hash(i);
        minhashes.push(HashWithValue {
            hash: current_hash,
            value: i,
        })
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

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::approx_eq;

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
        let data = [1, 2, 3];
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

    #[test]
    fn test_jaccard_no_overlap() {
        let mut queue = BoundedPriorityQueue::new(2, 3);
        queue.push(HashWithValue { hash: 1, value: 2 });
        queue.push(HashWithValue { hash: 2, value: 4 });
        queue.push(HashWithValue { hash: 4, value: 8 });
        let mut other_queue = BoundedPriorityQueue::new(2, 3);
        other_queue.push(HashWithValue { hash: 3, value: 8 });
        other_queue.push(HashWithValue { hash: 6, value: 9 });
        other_queue.push(HashWithValue { hash: 7, value: 10 });
        let result = jaccard(&queue, &other_queue);
        assert!(approx_eq!(f64, result, 0.));
    }

    #[test]
    fn test_jaccard() {
        let mut queue = BoundedPriorityQueue::new(2, 3);
        queue.push(HashWithValue { hash: 1, value: 2 });
        queue.push(HashWithValue { hash: 2, value: 4 });
        queue.push(HashWithValue { hash: 4, value: 8 });
        let mut other_queue = BoundedPriorityQueue::new(2, 3);
        other_queue.push(HashWithValue { hash: 1, value: 2 });
        other_queue.push(HashWithValue { hash: 6, value: 9 });
        other_queue.push(HashWithValue { hash: 7, value: 10 });
        let result = jaccard(&queue, &other_queue);
        assert!(approx_eq!(f64, result, 0.5));
    }
}
