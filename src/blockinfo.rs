#![allow(unused)]
use once_cell::sync::Lazy;
use serde::Deserialize;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::hash::Hash;
use std::sync::{Arc, Mutex};

use crate::aligner::BAligner;
use crate::utils::read_fasta;

pub static BLOCKFLAGS: Lazy<Arc<Mutex<HashMap<usize, String>>>> = Lazy::new(|| {
    let mut m = HashMap::new();
    m.insert(0, "null".to_string());
    Arc::new(Mutex::new(m))
});

enum SeqType {
    Fix,
    Anchor,
    Index,
}

#[derive(Clone, Debug)]
pub struct BlockInfo {
    pub idx: String,
    pub seq_type: String,
    pub seqs: HashMap<String, Vec<u8>>,
    pub max_mismatch: usize,
    pub query_start: Option<usize>,
    pub query_end: Option<usize>,
    pub seq_len: Option<usize>, // pub aligner: Box
    pub aligner: Option<BAligner>,
    pub flag: usize,
}

impl BlockInfo {
    /// Returns the query slice based on the query_start and query_end
    fn get_query_seq(self, query: &[u8]) -> (&[u8], usize) {
        let query_len = query.len();
        let start = if let Some(start) = self.query_start {
            start
        } else {
            0
        };
        let end = if let Some(end) = self.query_end {
            end
        } else if self.seq_len.is_some() {
            start + self.seq_len.unwrap()
        } else {
            query_len
        };

        (&query[start..end.min(query_len)], start)
    }
}

#[derive(Debug, Clone, Default, Deserialize)]
pub enum AlignMethod {
    #[default]
    SW,
    ANT,
}

#[derive(Debug, Clone, Default, Deserialize)]
struct BlockInfoFile {
    idx: String,
    seq_type: String,
    fasta_file: Option<String>,
    max_mismatch: usize,
    query_start: Option<usize>,
    query_end: Option<usize>,
    seq_len: Option<usize>,
    method: AlignMethod,
}

/// read block.info.tsv 获取 blockinfo 信息方便后面以 blockinfo 为基础, 对 read 进行比对
pub fn get_block_info(file_path: &str) -> Vec<BlockInfo> {
    let file = File::open(file_path).unwrap();
    // let mut rdr = csv::Reader::from_reader(file);
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);
    let mut block_info_vec: Vec<BlockInfo> = vec![];
    let mut flag: usize = 1;
    for result in rdr.deserialize() {
        let record: BlockInfoFile = result.unwrap();
        if &record.seq_type == "Fix" {
            let bi = BlockInfo {
                idx: record.idx.clone(),
                seq_type: record.seq_type,
                seqs: HashMap::new(),
                max_mismatch: record.max_mismatch,
                query_start: record.query_start,
                query_end: record.query_end,
                seq_len: None,
                aligner: None,
                flag: 0,
            };
            block_info_vec.push(bi);
            continue;
        }
        let seqs = if let Some(fasta_file) = &record.fasta_file {
            read_fasta(fasta_file).unwrap()
        } else {
            HashMap::new()
        };

        let aligner = BAligner::new(
            record.method,
            &record.fasta_file.unwrap(),
            record.max_mismatch,
        );
        let bi = BlockInfo {
            idx: record.idx.clone(),
            seq_type: record.seq_type,
            seqs,
            max_mismatch: record.max_mismatch,
            query_start: record.query_start,
            query_end: record.query_end,
            seq_len: record.seq_len,
            aligner: Some(aligner),
            flag,
        };
        BLOCKFLAGS.lock().unwrap().insert(flag, record.idx);
        flag *= 2;
        block_info_vec.push(bi);
    }
    block_info_vec
}

#[test]
fn test_get_block_info() {
    let workdir = env::current_dir();
    dbg!(workdir.unwrap());
    let block_info_vec = get_block_info("test/block.info.tsv");
    dbg!(block_info_vec);
    // File::open("./t");
}
