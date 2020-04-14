use bigtools::bbiread::{BBIRead, ChromAndSize};
use bigtools::bigbedread::BigBedRead;
use bigtools::bigwig::BedEntry;
use bigtools::remote_file::RemoteFile;

use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::string::String;

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
    name: String,
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
    chrom: String,
    start: u32,
    end: u32,
}

impl IntervalQuery {
    pub fn new(chrom: String, start: u32, end: u32) -> Self {
        IntervalQuery { chrom, start, end }
    }

    fn is_valid(&self) -> bool {
        self.end > self.start
    }

    /// Valid is defined as the end coordinate being strictly greater than the start. If
    /// this is not the case then we assume the input coordinates were just swapped, and
    /// we swap them back.
    pub fn make_valid(&mut self) {
        if !self.is_valid() {
            let old_start = self.start;
            self.start = self.end;
            self.end = old_start;
        }
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

enum FileOrUrlBigBedRead<T, U> {
    UrlRead(T),
    FileRead(U),
}

pub struct BigBedReader<T, U> {
    reader: FileOrUrlBigBedRead<T, U>,
    pub uri: String,
    chrom_data: Option<HashMap<String, ChromData>>,
}

type FileOrUrlBigBedReader = BigBedReader<
    BigBedRead<RemoteFile, RemoteFile>,
    BigBedRead<bigtools::seekableread::ReopenableFile, std::fs::File>,
>;

/// URI can refer to a local file or a remote file accessible over HTTP, will return a
/// reader of the appropriate type, albeit wrapped in an enum. Remote file URIs must
/// begin with http, otherwise it will be assumed the file is local.
impl FileOrUrlBigBedReader {
    pub fn new(uri: String) -> FileOrUrlBigBedReader {
        if uri.starts_with("http") {
            let bigbed = RemoteFile::new(&uri);
            BigBedReader {
                reader: FileOrUrlBigBedRead::UrlRead(BigBedRead::from(bigbed).unwrap()),
                uri,
                chrom_data: None,
            }
        } else {
            BigBedReader {
                reader: FileOrUrlBigBedRead::FileRead(
                    BigBedRead::from_file_and_attach(uri.clone()).unwrap(),
                ),
                uri,
                chrom_data: None,
            }
        }
    }

    fn get_chroms(&mut self) -> Vec<ChromAndSize> {
        match &self.reader {
            FileOrUrlBigBedRead::FileRead(rdr) => rdr.get_chroms(),
            FileOrUrlBigBedRead::UrlRead(rdr) => rdr.get_chroms(),
        }
    }

    pub fn get_data(&mut self, query: &IntervalQuery) -> Vec<BedEntry> {
        // Can only make valid if query is declared as mut: not the cleanest API.
        // Probably want to make a mutable copy internally.
        // query.make_valid();
        let results: Vec<_> = match &mut self.reader {
            FileOrUrlBigBedRead::UrlRead(rdr) => rdr
                .get_interval(&query.chrom, query.start, query.end)
                .unwrap()
                .collect::<Result<_, _>>()
                .unwrap(),
            FileOrUrlBigBedRead::FileRead(rdr) => rdr
                .get_interval(&query.chrom, query.start, query.end)
                .unwrap()
                .collect::<Result<_, _>>()
                .unwrap(),
        };
        results
    }

    /// If the chrom data has not been fetched yet, we pull it from the reader.
    /// Then in both cases we can extract the HashMap from the Option and return a
    /// reference to it.
    fn get_chrom_data(&mut self) -> &HashMap<String, ChromData> {
        if self.chrom_data.is_none() {
            let chroms = self.get_chroms();
            let mut chrom_and_data = HashMap::new();
            let mut sorted_chroms = chroms.iter().map(ChromInfo::from).collect::<Vec<_>>();
            sorted_chroms.sort();
            let mut offset = 0;
            for chrom in sorted_chroms.iter_mut() {
                chrom.data.offset = offset;
                chrom_and_data.insert(chrom.name.clone(), chrom.data.clone());
                offset += chrom.data.length;
            }
            self.chrom_data = Some(chrom_and_data);
        }
        self.chrom_data.as_ref().unwrap()
    }

    pub fn get_chrom_end(&mut self, chrom: &str) -> Result<u32, ChromNotFoundError> {
        let chrom_data = self.get_chrom_data();
        match chrom_data.get(chrom) {
            Some(data) => Ok(data.length),
            None => Err(ChromNotFoundError {
                chrom: chrom.to_string(),
            }),
        }
    }

    fn offset_data(data: &mut [BedEntry], offset: u32) {
        for entry in data.iter_mut() {
            entry.start += offset;
            entry.end += offset;
        }
    }

    // Given a set of query intervals, fetch results from all of the queries and
    // aggregate them into a single vector with offsets added to the positions such
    // that the positions will all be unique.
    //
    // This code doesn't currently work, although it should once Polonius borrow
    // checker hits stable: https://github.com/rust-lang/rust/issues/68117
    // fn get_offset_data(&mut self, queries: &[IntervalQuery]) -> Result<Vec<BedEntry>, ChromNotFoundError> {
    //     let chrom_data = self.get_chrom_data();
    //     let mut all_data = Vec::new();
    //     for query in queries.iter() {
    //         let current_chrom_data = match chrom_data.get(&(*query.chrom).to_string()) {
    //             Some(current_chrom_data) => current_chrom_data,
    //             None => {
    //                 return Err(ChromNotFoundError {
    //                     chrom: (*query.chrom).to_string(),
    //                 })
    //             }
    //         };
    //         let mut data = self.get_data(&query);  // cannot borrow `*self` as mutable more than once at a time
    //         Self::offset_data(&mut data, current_chrom_data.offset);
    //         all_data.extend(data.into_iter());
    //     }
    //     Ok(all_data)
    // }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_query_is_valid() {
        let query = IntervalQuery::new("chr1".to_string(), 10, 1);
        assert!(!query.is_valid())
    }

    #[test]
    fn test_interval_query_make_valid_already_valid() {
        let mut query = IntervalQuery::new("chr1".to_string(), 1, 10);
        query.make_valid();
        assert!(query.is_valid());
        assert_eq!(query.start, 1);
        assert_eq!(query.end, 10);
    }

    #[test]
    fn test_interval_query_make_valid_from_invalid() {
        let mut query = IntervalQuery::new("chr1".to_string(), 10, 1);
        query.make_valid();
        assert!(query.is_valid());
        assert_eq!(query.start, 1);
        assert_eq!(query.end, 10);
    }

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
