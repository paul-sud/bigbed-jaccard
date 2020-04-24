use bigtools::bbiread::{BBIRead, ChromAndSize};
use bigtools::bigbedread::BigBedRead;
use bigtools::bigwig::BedEntry;
use bigtools::seekableread::{Reopen, SeekableRead};

use std::collections::HashMap;
use std::error::Error;
use std::fmt;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct ChromData {
    length: u32,
    offset: u32,
}

impl ChromData {
    fn new(length: u32, offset: u32) -> Self {
        ChromData { length, offset }
    }
}

#[derive(Eq, Ord, PartialOrd)]
pub struct ChromInfo {
    name: std::string::String,
    data: ChromData,
}

impl PartialEq for ChromInfo {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

impl From<&ChromAndSize> for ChromInfo {
    fn from(chrom: &ChromAndSize) -> Self {
        ChromInfo {
            name: chrom.name.clone(),
            data: ChromData::new(chrom.length, 0),
        }
    }
}

pub struct IntervalQuery {
    chrom: std::string::String,
    start: u32,
    end: u32,
}

impl IntervalQuery {
    pub fn new(chrom: std::string::String, start: u32, end: u32) -> Self {
        IntervalQuery { chrom, start, end }
    }
}

#[derive(Debug)]
pub struct ChromNotFoundError {
    chrom: std::string::String,
}

impl fmt::Display for ChromNotFoundError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Could not find data for chrom {}. Ensure that the files contain the same set of chromosomes you are querying.", self.chrom)
    }
}

impl Error for ChromNotFoundError {}

pub fn get_data<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    reader: &mut BigBedRead<R, S>,
    query: &IntervalQuery,
) -> Vec<BedEntry> {
    let remote_intervals: Vec<_> = reader
        .get_interval(&query.chrom, query.start, query.end)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    remote_intervals
}

pub fn get_chrom_end<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    reader: &mut BigBedRead<R, S>,
    chrom: &str,
) -> u32 {
    let chroms = reader.get_chroms();
    let chrom = chroms.iter().find(|v| v.name == chrom).unwrap();
    chrom.length
}

pub fn get_chrom_data(chroms: &[ChromAndSize]) -> HashMap<std::string::String, ChromData> {
    let mut chrom_and_data = HashMap::new();
    let mut sorted_chroms = chroms.iter().map(ChromInfo::from).collect::<Vec<_>>();
    sorted_chroms.sort();
    let mut offset = 0;
    for chrom in sorted_chroms.iter_mut() {
        chrom.data.offset = offset;
        chrom_and_data.insert(chrom.name.clone(), chrom.data.clone());
        offset += chrom.data.length;
    }
    chrom_and_data
}

fn offset_data(data: &mut [BedEntry], offset: u32) {
    for entry in data.iter_mut() {
        entry.start += offset;
        entry.end += offset;
    }
}

/// Given a set of query intervals, fetch results from all of the queries and
/// aggregate them into a single vector with offsets added to the positions such
/// that the positions will all be unique.
fn get_offset_intervals<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    reader: &mut BigBedRead<R, S>,
    queries: &[IntervalQuery],
) -> Result<Vec<BedEntry>, ChromNotFoundError> {
    let chrom_data = reader.get_chroms();
    let chrom_data = get_chrom_data(&chrom_data);
    let mut all_data = Vec::new();
    for query in queries.iter() {
        let current_chrom_data = match chrom_data.get(&(*query.chrom).to_string()) {
            Some(current_chrom_data) => current_chrom_data,
            None => {
                return Err(ChromNotFoundError {
                    chrom: (*query.chrom).to_string(),
                })
            }
        };
        let mut data = get_data(reader, &query);
        offset_data(&mut data, current_chrom_data.offset);
        all_data.extend(data.into_iter());
    }
    Ok(all_data)
}

pub fn get_offset_data<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    reader: &mut BigBedRead<R, S>,
    queries: &[IntervalQuery],
) -> Result<Vec<u32>, ChromNotFoundError> {
    match get_offset_intervals(reader, queries) {
        Ok(intervals) => Ok(intervals_to_vec(&intervals)),
        Err(i) => Err(i),
    }
}

fn intervals_to_vec(data: &[BedEntry]) -> Vec<u32> {
    let mut result = Vec::with_capacity(data.len());
    let BedEntry {
        start: virtual_start,
        end: mut virtual_end,
        ..
    } = data[0];
    for i in virtual_start..virtual_end {
        result.push(i);
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
            result.push(i);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_chrom_data() {
        let chroms = [
            ChromAndSize {
                name: "chr1".to_string(),
                length: 10,
            },
            ChromAndSize {
                name: "chr2".to_string(),
                length: 20,
            },
        ];
        let result = get_chrom_data(&chroms);
        assert_eq!(
            result[&"chr1".to_string()],
            ChromData {
                length: 10,
                offset: 0
            }
        );
        assert_eq!(
            result[&"chr2".to_string()],
            ChromData {
                length: 20,
                offset: 10
            }
        );
    }

    #[test]
    fn test_offset_data() {
        let mut data = [
            BedEntry {
                start: 1,
                end: 2,
                rest: "".to_string(),
            },
            BedEntry {
                start: 4,
                end: 10,
                rest: "".to_string(),
            },
        ];
        offset_data(&mut data, 3);
        assert_eq!(data[0].start, 4);
        assert_eq!(data[0].end, 5);
        assert_eq!(data[1].start, 7);
        assert_eq!(data[1].end, 13);
    }
}
