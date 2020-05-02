use crate::oph::OnePermutationHasher;

use rayon::prelude::*;

use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::hash::Hash;
use std::sync::{Arc, Mutex};

#[derive(Debug)]
pub struct IncorrectlyChunkedInputSketchError {
    message: String,
}

impl fmt::Display for IncorrectlyChunkedInputSketchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl Error for IncorrectlyChunkedInputSketchError {}

pub struct Lsh<K, V> {
    /// K should be a hashable type, and V should be something small like a String, as
    /// opposed to the entire sketch vector. The first index is the bin number, the
    /// second is the actual hashmap. This implementation is not optimal because the
    /// entire structure is locked globally, if we are modifying data in independent
    /// bins then those operations can occur without synchronization. It is only if we
    /// try to concurrently update or concurrently read and update data in the same bin
    /// that we end up with the potential for a data race. However, we would still need
    /// the index of the independent Arc<Mutex<HashMap>> maps to be thread safe. A
    /// better way to manage them would be using a channel.
    hashers: Arc<Mutex<HashMap<usize, HashMap<Vec<K>, Vec<V>>>>>,
    /// Corresponds to L in the literature, L independent hash tables over groups of K
    /// values.
    n_hashes: usize,
    /// Corresponds to K in the literature, number of values (hashes) per L (bin). NOT
    /// the same K used as the generic type parameter for this struct.
    n_values: usize,
    /// Table mapping from sketch ID to actual data, so we only need to store IDs in the
    /// hashers
    sketches: HashMap<V, Vec<K>>,
}

impl<K: Clone + Eq + Hash + Send + Sync, V: Clone + Copy + Hash + Ord + Send + Sync> Lsh<K, V> {
    /// # Arguments
    ///
    /// * `n_hashes` - L in the literature. The number of hash tables to construct
    /// * `n_values` - K in the literature. The length of the vector for each table to hash.
    pub fn new(n_hashes: usize, n_values: usize) -> Self {
        let mut hashers = HashMap::new();
        for i in 0..n_hashes {
            hashers.insert(i, HashMap::new());
        }
        Self {
            hashers: Arc::new(Mutex::new(hashers)),
            n_hashes,
            n_values,
            sketches: HashMap::new(),
        }
    }

    fn chunk_sketch(
        &self,
        sketch: &[K],
    ) -> Result<Vec<Vec<K>>, IncorrectlyChunkedInputSketchError> {
        let expected = self.n_hashes * self.n_values;
        let sketch_len = sketch.len();
        if sketch_len != expected {
            return Err(IncorrectlyChunkedInputSketchError {
                message: format!(
                    "sketch does not have right number of elements, expected {} but found {}",
                    expected, sketch_len
                ),
            });
        };
        let mut result = Vec::with_capacity(self.n_hashes);
        for chunk in sketch.chunks(self.n_values) {
            result.push(chunk.to_vec());
        }
        Ok(result)
    }

    /// Insert data for bins into the appropriate map in parallel. Will panic if the
    /// right bin is not found. If the data was already indexed, then will return
    /// without inserting.
    pub fn insert(&mut self, sketch: &[K], sketch_id: V) {
        if self.sketches.contains_key(&sketch_id) {
            return;
        }
        self.sketches.insert(sketch_id, sketch.to_vec());
        let chunked_sketch = self.chunk_sketch(sketch).unwrap();
        chunked_sketch
            .par_iter()
            .enumerate()
            .map(|(i, vals)| {
                let hashers = self.hashers.clone();
                let mut map = hashers.lock().unwrap();
                match map.get_mut(&i) {
                    Some(map) => match map.get_mut(vals) {
                        Some(val) => val.push(sketch_id),
                        None => {
                            map.insert(vals.to_vec(), vec![sketch_id]);
                        }
                    },
                    None => {
                        return Err(IncorrectlyChunkedInputSketchError {
                            message: "Data had more bins than expected".to_string(),
                        })
                    }
                }
                Ok(())
            })
            .collect::<Vec<_>>();
    }

    /// Given a query sketch vector, return the list of IDs that returned hits in any of
    /// the hash tables. If no results are found, then will return None. To get ranked
    /// results, see `ranked_query`.
    pub fn query(&self, sketch: &[K]) -> Option<Vec<V>> {
        let chunked_sketch = self.chunk_sketch(sketch).unwrap();
        let hashers = self.hashers.clone();
        let map = hashers.lock().unwrap();
        let mut results = chunked_sketch
            .par_iter()
            .enumerate()
            .map(|(i, chunk)| {
                // Index i should always be in the map. However there are no such
                // guarantees for the chunk, so we can't unwrap.
                match map.get(&i).unwrap().get(chunk) {
                    Some(v) => v.to_owned(),
                    None => vec![],
                }
            })
            .flatten()
            .collect::<Vec<_>>();

        results.sort();
        results.dedup();
        if results.is_empty() {
            return None;
        }
        Some(results)
    }

    /// Given a vector of queries, retrieve results in parallel
    pub fn query_bulk(&self, sketches: &[&[K]]) -> Vec<Option<Vec<V>>> {
        sketches
            .par_iter()
            .map(|sketch| self.query(sketch))
            .collect::<Vec<_>>()
    }

    /// Given a set of putative nearest neighbors, compute the actual distance between
    /// them and return the ranked results as (result, similarity) tuples.
    fn get_top_hits(&self, results: Option<Vec<V>>, query: &[K]) -> Option<Vec<(V, f64)>> {
        results.as_ref()?;
        let res = results.unwrap();
        let ranking = res
            .par_iter()
            .map(|result| {
                let jac = OnePermutationHasher::jaccard(self.sketches.get(result).unwrap(), &query)
                    .unwrap();
                (jac, result)
            })
            .collect::<Vec<_>>();
        let mut sorted = ranking.as_slice().to_owned();
        sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
        Some(sorted.iter().map(|(jac, &result)| (result, *jac)).collect())
    }

    pub fn ranked_query(&self, query: &[K]) -> Option<Vec<(V, f64)>> {
        let results = self.query(query);
        self.get_top_hits(results, query)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::approx_eq;

    #[test]
    fn test_lsh_new() {
        let lsh: Lsh<u32, u32> = Lsh::new(3_usize, 2_usize);
        let result = lsh.hashers.lock().unwrap().keys().len();
        assert_eq!(result, 3)
    }

    #[test]
    fn test_chunk_sketch() {
        let lsh: Lsh<u32, u32> = Lsh::new(3_usize, 2_usize);
        assert_eq!(
            lsh.chunk_sketch(&[1, 2, 3, 4, 5, 6]).unwrap(),
            vec![vec![1, 2], vec![3, 4], vec![5, 6]]
        );
    }

    #[test]
    #[should_panic]
    fn test_chunk_sketch_wrong_size_should_panic() {
        let lsh: Lsh<u32, u32> = Lsh::new(3_usize, 2_usize);
        lsh.chunk_sketch(&[1, 2, 3, 4, 5]).unwrap();
    }

    #[test]
    fn test_lsh_query_some() {
        let mut lsh: Lsh<u32, u32> = Lsh::new(3_usize, 2_usize);
        lsh.insert(&[1, 2, 3, 4, 5, 6], 12);
        let result = lsh.query(&[1, 2, 17, 18, 5, 6]).unwrap();
        assert_eq!(result, vec![12])
    }

    #[test]
    fn test_lsh_query_no_hits_returns_none() {
        let mut lsh: Lsh<u32, u32> = Lsh::new(2_usize, 2_usize);
        lsh.insert(&[1, 2, 3, 4], 12);
        let result = lsh.query(&[17, 18, 5, 6]);
        assert!(result.is_none())
    }

    #[test]
    fn test_lsh_query_bulk() {
        let mut lsh: Lsh<u32, u32> = Lsh::new(2_usize, 2_usize);
        lsh.insert(&[1, 2, 3, 4], 12);
        lsh.insert(&[5, 6, 7, 8], 13);
        let result = lsh.query_bulk(&[&[17, 18, 3, 4], &[5, 6, 10, 11], &[9, 10, 11, 12]]);
        assert_eq!(result, vec![Some(vec![12]), Some(vec![13]), None])
    }

    #[test]
    fn test_lsh_insert() {
        let mut lsh: Lsh<u32, &str> = Lsh::new(3_usize, 2_usize);
        lsh.insert(&[1, 2, 3, 4, 5, 6], "foo");
        let hashers = lsh.hashers.lock().unwrap();
        assert_eq!(
            hashers.get(&0_usize).unwrap().get(&vec![1, 2]).unwrap(),
            &vec!["foo"]
        );
        assert_eq!(
            hashers.get(&1_usize).unwrap().get(&vec![3, 4]).unwrap(),
            &vec!["foo"]
        );
        assert_eq!(
            hashers.get(&2_usize).unwrap().get(&vec![5, 6]).unwrap(),
            &vec!["foo"]
        );
        assert_eq!(lsh.sketches.get("foo").unwrap(), &vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn test_lsh_get_top_hits() {
        let mut lsh: Lsh<u32, u32> = Lsh::new(2_usize, 2_usize);
        lsh.insert(&[1, 2, 3, 4], 12);
        lsh.insert(&[1, 2, 8, 9], 13);
        let query = vec![1, 2, 3, 5];
        let results = Some(vec![12, 13]);
        let hits = lsh.get_top_hits(results, &query).unwrap();
        let (hit_0_id, hit_0_jac) = hits[0];
        assert_eq!(hit_0_id, 12);
        assert!(approx_eq!(f64, hit_0_jac, 0.75));
        let (hit_1_id, hit_1_jac) = hits[1];
        assert_eq!(hit_1_id, 13);
        assert!(approx_eq!(f64, hit_1_jac, 0.5));
    }
}
