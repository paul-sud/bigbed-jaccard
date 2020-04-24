use bigtools::bigwig::BedEntry;

use crate::hash::hash;

use rand::prelude::*;

use std::error::Error;
use std::fmt;

#[derive(Debug)]
pub struct EmptySketchError;

impl fmt::Display for EmptySketchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "All bins of sketch were empty, could not densify.")
    }
}

impl Error for EmptySketchError {}

// Seed for generating OPH indicators for densification
// Generated with this code:
// use rand::prelude::*;
// fn main() {
//     let mut rng = rand::thread_rng();
//     let result: u64 = rng.gen();
//     dbg!(result);
// }
const OPH_SEED: u64 = 12_963_566_592_248_599_713;

pub struct OnePermutationHasher {
    num_bins: usize,
    offset: u32,
    rng: StdRng,
    indicators: Option<Vec<bool>>,
}

impl OnePermutationHasher {
    pub fn new(num_bins: usize) -> Self {
        let rng: StdRng = SeedableRng::seed_from_u64(OPH_SEED);
        let offset = std::u32::MAX / (num_bins as u32) + 1;
        OnePermutationHasher {
            num_bins,
            offset,
            rng,
            indicators: None,
        }
    }

    pub fn get_indicators(&mut self) -> &Vec<bool> {
        match &self.indicators {
            Some(_) => {}
            None => {
                let mut indicators = vec![false; self.num_bins];
                for elem in indicators.iter_mut() {
                    *elem = self.rng.gen_bool(0.5);
                }
                self.indicators = Some(indicators);
            }
        }
        self.indicators.as_ref().unwrap()
    }

    /// Given a sketch containing potentially empty values (None), densify it according
    /// to the improved densification scheme of Shrivastava and Li 2014
    /// http://www.auai.org/uai2014/proceedings/individuals/225.pdf
    /// Basically, to fill empty bins, we find the next non-empty bin in the direction
    /// depending on the value for the indicator for that bin. Then we fill the bin with
    /// the value of the neighboring non-empty bin plus a constant times the number of
    /// bins traveled to reach the non-empty bin (always at least 1).
    fn densify(&mut self, sketch: &[Option<u32>]) -> Result<Vec<u32>, EmptySketchError> {
        let mut dense = vec![0; self.num_bins];
        let offset = self.offset;
        let indicators = self.get_indicators();
        for (i, value) in sketch.iter().enumerate() {
            dense[i] = match value {
                Some(val) => *val,
                None => {
                    let search_result = if indicators[i] {
                        // Look for closest non-empty bin to the right
                        sketch[i + 1..sketch.len()]
                            .iter()
                            .chain(sketch[0..i].iter())
                            .enumerate()
                            .find(|(_, x)| !x.is_none())
                    } else {
                        // Look to the left, note the reversal of the iterator order
                        sketch[i + 1..sketch.len()]
                            .iter()
                            .chain(sketch[0..i].iter())
                            .rev()
                            .enumerate()
                            .find(|(_, x)| !x.is_none())
                    };
                    match search_result {
                        // We already checked above that value is Some.
                        Some((index, Some(value))) => value + (index as u32 + 1) * offset,
                        _ => return Err(EmptySketchError),
                    }
                }
            };
        }
        Ok(dense)
    }

    /// Return an un-densified sketch of the data using the one-permutation hashing (OPH).
    /// TODO: this should be decoupled from BedEntry specific things to just accept a
    /// bag of integers. The overlapping BedEntry logic should be a separate method.
    /// When generators hit stable then it will be easy to convert, for now converting
    /// to an iterator would require storing the iterator state on a struct manually.
    fn sketch(&self, data: &[BedEntry]) -> Vec<Option<u32>> {
        let mut output = vec![None; self.num_bins];
        let BedEntry {
            start: virtual_start,
            end: mut virtual_end,
            ..
        } = data[0];
        for i in virtual_start..virtual_end {
            let current_hash = hash(i);
            let bin = (current_hash % self.num_bins as u32) as usize;
            match output[bin] {
                Some(val) => {
                    if current_hash < val {
                        output[bin] = Some(current_hash);
                    }
                }
                None => output[bin] = Some(current_hash),
            }
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
                let bin = (current_hash % self.num_bins as u32) as usize;
                match output[bin] {
                    Some(val) => {
                        if current_hash < val {
                            output[bin] = Some(current_hash);
                        }
                    }
                    None => output[bin] = Some(current_hash),
                }
            }
        }
        output
    }

    pub fn dense_sketch(&mut self, data: &[BedEntry]) -> Result<Vec<u32>, EmptySketchError> {
        let sketch = self.sketch(data);
        self.densify(&sketch)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test using the constant RNG seed.
    #[test]
    fn test_one_permutation_hasher_new() {
        let mut oph = OnePermutationHasher::new(3);
        assert_eq!(oph.get_indicators(), &vec![true, false, true]);
    }

    #[test]
    fn test_one_permutation_hasher_densify() {
        let mut oph = OnePermutationHasher::new(4);
        let data = [None, None, None, Some(2)];
        let dense = oph.densify(&data).unwrap();
        assert_eq!(dense, vec![2 + 3 * oph.offset, 2 + 2 * oph.offset, 2 + oph.offset, 2]);
    }
}
