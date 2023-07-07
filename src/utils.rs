use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufReader, Read};
use std::ops::Range;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::alphabets::dna::complement;
use bio::io::fasta;
use bio::io::fastq::{self, Record};
use flate2::read::MultiGzDecoder;

use crate::blockalign::{block_align_read, ReadBlockAlign};
use crate::blockinfo::BlockInfo;

// fastq reader, file as arg, decide based on extension
pub fn get_reader(path: &str) -> Box<dyn Read + Send> {
    if path.ends_with(".gz") {
        let f = File::open(path).unwrap();
        // Box::new(bufread::MultiGzDecoder::new(BufReader::new(f)))
        Box::new(MultiGzDecoder::new(BufReader::new(f)))
    } else {
        let f = File::open(path).unwrap();
        Box::new(BufReader::new(f))
    }
}

pub fn read_fasta(fa_file: &str) -> Result<HashMap<String, Vec<u8>>> {
    let index_reader = fasta::Reader::from_file(fa_file).unwrap();
    let mut index_hash: HashMap<String, Vec<u8>> = HashMap::new();
    for record_r in index_reader.records() {
        let record = record_r.unwrap();
        let sample_name = record.id().to_string();
        index_hash.insert(sample_name.to_owned(), record.seq().to_vec());
    }
    Ok(index_hash)
}

pub fn write_fasta(seq_map: &HashMap<String, Vec<u8>>, out_file: &str) -> Result<()> {
    let handle = io::BufWriter::new(fs::File::create(out_file).unwrap());
    let mut writer = fasta::Writer::new(handle);
    seq_map.iter().for_each(|(read_id, seq)| {
        let record = fasta::Record::with_attrs(read_id, None, seq);
        let write_result = writer.write_record(&record);
    });

    Ok(())
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut rc_seq: Vec<u8> = Vec::new();
    rc_seq.extend(seq.iter().rev().cloned().map(complement));
    rc_seq
}

#[test]
fn test_revcomp() {
    let seq = b"AATTCCGG";
    let rc_seq = revcomp(seq);
    // dbg!(std::str::from_utf8(&rc_seq).unwrap());
    assert_eq!(&rc_seq, b"CCGGAATT")

    // dbg!(rc_seq);
}

/// union multiple Range
pub fn union(ranges: Vec<std::ops::Range<i32>>) -> Vec<std::ops::Range<i32>> {
    let mut sorted_ranges = ranges.clone();
    sorted_ranges.sort_by(|a, b| a.start.cmp(&b.start));

    let mut result = Vec::new();
    let mut current = sorted_ranges[0].clone();

    for range in sorted_ranges {
        if range.start <= current.end {
            current.end = current.end.max(range.end);
        } else {
            result.push(current);
            current = range;
        }
    }

    result.push(current);
    result
}

/// range difference
fn difference(
    range1: std::ops::Range<i32>,
    ranges: Vec<std::ops::Range<i32>>,
) -> Vec<std::ops::Range<i32>> {
    let mut sorted_ranges = ranges.clone();
    sorted_ranges.sort_by(|a, b| a.start.cmp(&b.start));

    let mut result = Vec::new();
    let mut current = range1.clone();

    for range in sorted_ranges {
        if range.start <= current.end {
            if current.start < range.start {
                result.push(current.start..range.start);
            }
            if current.end > range.end {
                current.start = range.end;
            } else {
                return result;
            }
        }
    }

    result.push(current.start..current.end);
    result
}

pub fn dna_to_spans(dna: &str, ranges: &[(Range<usize>, String)], other_class: &str) -> String {
    let mut sorted_ranges = ranges.to_vec();
    sorted_ranges.sort_by(|a, b| a.0.start.cmp(&b.0.start));

    let mut result = String::new();
    let mut current_end = 0;

    for (range, class) in sorted_ranges {
        if range.start > current_end {
            result.push_str(&format!(
                "<span class='{}'>{}</span>",
                other_class,
                &dna[current_end..range.start]
            ));
        }
        result.push_str(&format!(
            "<span class='{}'>{}</span>",
            class,
            &dna[range.start..range.end]
        ));
        current_end = range.end;
    }

    if current_end < dna.len() {
        result.push_str(&format!(
            "<span class='{}'>{}</span>",
            other_class,
            &dna[current_end..]
        ));
    }

    result
}

#[test]
fn test_dna_to_spans() {
    let dna = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let ranges = [
        (2..6, "A".to_string()),
        (8..12, "B".to_string()),
        (14..18, "C".to_string()),
    ];

    let output = dna_to_spans(dna, &ranges, "other");
    assert_eq!(
        &output,
        "<span class='other'>AC</span><span class='A'>GTAC</span>\
<span class='other'>GT</span><span class='B'>ACGT</span><span class='other'>AC</span>\
<span class='C'>GTAC</span><span class='other'>GTACGTACGTACGT</span>"
    )
}
