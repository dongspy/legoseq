use std::collections::HashMap;
use std::fs::{File, self};
use std::io::{BufReader, Read, self};
use std::path::{PathBuf, Path};

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

pub fn write_fasta(seq_map: &HashMap<String, Vec<u8>>, out_file: &str) -> Result<()>{
    let handle = io::BufWriter::new(
        fs::File::create(out_file).unwrap()
    );
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