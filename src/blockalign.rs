#![allow(unused)]
use std::collections::HashMap;
use std::io::Read;

use bio::io::fastq::Record;
use rayon::vec;
use serde::{Deserialize, Serialize};
use tracing::{debug, error, info, warn, Level};

use crate::aligner::Alignment;
use crate::blockinfo::{get_block_info_fasta, BlockInfo};
use crate::utils::dna_to_spans;

#[derive(Debug, Clone)]
pub struct BlockAlign {
    pub info: BlockInfo,
    pub align: Option<Alignment>,
    pub n_match: usize,
    pub best_index: String,
}

impl BlockAlign{

    pub fn get_query_start(&self) -> Option<usize> {
        self.align.as_ref().map(|align| align.query_start)
    }

    pub fn get_query_end(&self) -> Option<usize> {
        self.align.as_ref().map(|align| align.query_end)
    }
}


#[derive(Debug, Clone)]
pub struct BlockAlignAbbr {
    pub best_index: String,
    pub index_start: usize,
    pub index_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub n_match: usize,
    pub strand: char,
    pub flag: usize,
}

impl BlockAlignAbbr {
    pub fn to_str(&self) -> String {
        format!(
            "{}:{}:{}:{}:{}:{}:{}",
            self.best_index,
            self.index_start,
            self.index_end,
            self.query_start,
            self.query_end,
            self.n_match,
            self.strand
        )
    }
}

impl BlockAlign {
    pub fn to_abbr(&self) -> Option<BlockAlignAbbr> {
        // let mut hash = HashMap::new();
        self.align.as_ref().map(|aln| BlockAlignAbbr {
            best_index: self.best_index.to_owned(),
            index_start: aln.index_start,
            index_end: aln.index_end,
            query_start: aln.query_start,
            query_end: aln.query_end,
            n_match: self.n_match,
            strand: aln.strand,
            // flag: aln.flag
            flag: self.info.flag,
        })
    }
}

/// read sequence mapping aganist the block sequence
/// return BlockAlign vector
pub fn block_align_read(read: &[u8], block_info_list: &[BlockInfo]) -> Vec<Option<BlockAlign>> {
    // let aligner2 = &aligner;
    let mut block_align_list: Vec<Option<BlockAlign>> = vec![];
    let mut pre_query_end: usize = 0;

    for block_info in block_info_list.iter().filter(|x| x.seq_type != "Fix") {
        let mut block_info = block_info.to_owned();
        // block_info.query_start = Some(pre_query_end.saturating_sub(OFFSET));
        // let ba = block_info.clone().align(read, aligner);
        let align = block_info.aligner.clone().unwrap().align(read);
        let ba = if let Some(align) = align {
            let ba = BlockAlign {
                info: block_info.clone(),
                align: Some(align.clone()),
                n_match: align.clone().n_match,
                best_index: align.best_index,
            };
            // ba.to_abbr()
            Some(ba)
        } else {
            None
        };
        block_align_list.push(ba);
    }
    block_align_list
}

pub fn block_align_read_with_insert(
    read: &[u8],
    block_info_list: &[BlockInfo],
) -> HashMap<String, Option<BlockAlign>> {
    // let aligner2 = &aligner;
    let mut block_align_hash: HashMap<String, Option<BlockAlign>> = HashMap::new();
    let mut pre_query_end: usize = 0;
    let read_len = read.len();

    for block_info in block_info_list.iter().filter(|x| x.seq_type != "Fix") {
        let mut block_info = block_info.clone();
        let idx = (&block_info.idx).clone();
        // block_info.query_start = Some(pre_query_end.saturating_sub(OFFSET));
        // let ba = block_info.clone().align(read, aligner);
        let align = block_info.aligner.clone().unwrap().align(read);
        let ba = if let Some(align) = align {
            let ba = BlockAlign {
                info: block_info.clone(),
                align: Some(align.clone()),
                n_match: align.clone().n_match,
                best_index: align.best_index,
            };
            // ba.to_abbr()
            Some(ba)
        } else {
            None
        };
        // block_align_list.push(ba);
        block_align_hash.insert(idx, ba);
    }

    // 处理可变序列
    let block_info_len = block_info_list.len();
    for block_ii in 0..block_info_len {
        let mut block_info = block_info_list[block_ii].clone();
        let idx = (&block_info.idx).to_string();
        if block_info.seq_type != "Variable" {
            continue;
        }

        let mut query_start: Option<usize> = None;
        // 开头模块
        if block_ii == 0 {
            query_start = Some(0);
        } else {
            let pre_block_info = &block_info_list[block_ii - 1];
            let ba = block_align_hash.get(&pre_block_info.idx).unwrap().clone().unwrap();
            query_start = ba.get_query_end().map(|pos| pos + 1);

        }

        // 末尾模块
        let mut query_end: Option<usize> = None;
        if block_ii == (block_info_len - 1) {
            query_end = Some(read_len);
        }
        let next_block_info = &block_info_list[block_ii + 1];
        let ba = block_align_hash.get(&next_block_info.idx).unwrap().clone().unwrap();
        query_end = ba.get_query_start().map(|pos| pos - 1);
        
        let align = Alignment{
            best_index: "".to_string(),
            index_start: 0,
            index_end: 0,
            query_start: query_start.unwrap() ,
            query_end: query_end.unwrap(),
            n_match: 0,
            strand: '+',
            operations: ba.clone().align.unwrap().operations,
            // operations: vec![AlignmentOperation::]
        };

        let ba = BlockAlign {
            info: block_info.clone(),
            align: Some(align),
            n_match: 0,
            best_index: "".to_string(),
        };

        block_align_hash.insert(idx, Some(ba));
    }
    block_align_hash
}

static OFFSET: usize = 4;

/// save the fastq record and its block align information
#[derive(Debug, Clone)]
pub struct ReadBlockAlign {
    pub record: Record,
    pub block_align: Vec<Option<BlockAlign>>,
}

impl ReadBlockAlign {
    /// mapping read againt the block sequence
    /// return ReadBlockAlign struct
    pub fn read_block_info(record: &Record, block_info_list: &[BlockInfo]) -> Self {
        let read_name = record.id();
        let read_seq = record.seq();
        let ba_list = block_align_read(read_seq, block_info_list);
        Self {
            record: record.clone(),
            block_align: ba_list,
        }
    }

    /// get the block flag of the read
    pub fn get_block_flag(&self) -> usize {
        let mut flag = 0;
        for ba in self.block_align.clone().iter() {
            if let Some(ba) = ba {
                flag |= ba.info.flag;
            } else {
            }
        }
        flag
    }

    /// tostring
    pub fn get_block_str(&self) -> String {
        let mut block_str_list: Vec<String> = vec![];
        for ba in self.block_align.clone().iter() {
            if let Some(ba) = ba {
                block_str_list.push(ba.to_abbr().unwrap().to_str().to_owned());
            } else {
            }
        }

        block_str_list.join(";")
    }

    /// generate the new fastq record that the sequcne only include the export_block
    /// which is in the blockinfo
    pub fn get_new_record(
        &self,
        block_info_list: &[BlockInfo],
        export_block: &[String],
    ) -> Option<Record> {
        let record = &self.record;
        let seq_len = record.seq().len();
        // get the position of every block of new record
        let mut block_hash = HashMap::new();
        block_info_list
            .iter()
            .filter(|&x| x.seq_type == "Fix")
            .for_each(|x| {
                block_hash.insert(
                    x.idx.to_owned(),
                    (x.query_start.unwrap_or(0), x.query_end.unwrap_or(0)),
                );
            });

        //
        self.block_align.iter().for_each(|x| {
            if let Some(x) = x {
                //   let ba=  x.to_abbr().unwrap();
                if x.align.is_some() {
                    let idx = x.info.idx.to_owned();
                    let query_start = x.align.as_ref().unwrap().clone().query_start;
                    let query_end = x.align.as_ref().unwrap().clone().query_end;
                    block_hash.insert(idx, (query_start, query_end));
                }
            }
        });

        let mut new_seq_vec = Vec::new();
        let mut new_qual_vec = Vec::new();
        for block in export_block.iter() {
            if !block_hash.contains_key(block) {
                return None;
            }
            let (start, end) = block_hash.get(block).unwrap();
            let new_seq = &record.seq().to_vec()[*start.min(&seq_len)..*end.min(&seq_len)];
            new_seq_vec.extend(new_seq);
            let new_qual = &record.qual().to_vec()[*start.min(&seq_len)..*end.min(&seq_len)];
            new_qual_vec.extend(new_qual);
        }

        let new_record =
            Record::with_attrs(record.id(), record.desc(), &new_seq_vec, &new_qual_vec);
        Some(new_record)
    }

    /// blockinfo
    pub fn get_seq_hashmap(
        &self,
        block_info_list: &[BlockInfo],
    ) -> Option<HashMap<String, SeqOut>> {
        let record = &self.record;
        let seq_len = record.seq().len();
        let mut seq_hash: HashMap<String, SeqOut> = HashMap::new();
        // get the position of every block of new record
        let mut block_hash = HashMap::new();

        // init the block_hash
        // block_info_list 是最初输入的
        block_info_list
            .iter()
            // .filter(|&x| x.seq_type == "Fix")
            .for_each(|x| {
                block_hash.insert(
                    x.idx.to_owned(),
                    (x.query_start.unwrap_or(0), x.query_end.unwrap_or(0)),
                );
                // seq_hash.insert(x.idx.to_owned(), );
            });

        // todo:
        let mut block_align_count = 0;
        self.block_align.iter().for_each(|x| {
            if let Some(x) = x {
                //   let ba=  x.to_abbr().unwrap();
                let idx = x.info.idx.to_owned();
                if x.align.is_some() {
                    block_align_count += 1;
                    let query_start = x.align.as_ref().unwrap().clone().query_start;
                    let query_end = x.align.as_ref().unwrap().clone().query_end;
                    block_hash.insert(idx, (query_start, query_end));
                }
            }
        });
        block_hash.iter().for_each(|(name, (start, end))| {
            let seq = &record.seq().to_vec()[*start.min(&seq_len)..*end.min(&seq_len)];
            let qual = &record.qual().to_vec()[*start.min(&seq_len)..*end.min(&seq_len)];
            let seq = Some(String::from_utf8(seq.to_vec()).unwrap());
            let qual = Some(String::from_utf8(qual.to_vec()).unwrap());
            let seqout = SeqOut::new(name.to_string(), None, seq, qual);
            seq_hash.insert(name.to_owned(), seqout);
        });

        // add fastq tag
        seq_hash.insert(
            "read".to_string(),
            SeqOut::new(
                (&record.id()).to_string(),
                (&record.desc().map(|s| s.to_string())).to_owned(),
                None,
                None,
            ),
        );
        if block_align_count > 0 {
            Some(seq_hash)
        } else {
            None
        }
    }

    /// blockinfo
    /// 主要是将插入片段加进去
    pub fn get_seq_hashmap2(
        &self,
        block_info_list: &[BlockInfo],
    ) -> Option<HashMap<String, SeqOut>> {
        let record = &self.record;
        let seq_len = record.seq().len();
        let mut seq_hash: HashMap<String, SeqOut> = HashMap::new();
        // get the position of every block of new record
        let mut block_hash = HashMap::new();

        // init the block_hash
        // block_info_list 是最初输入的
        block_info_list
            .iter()
            // .filter(|&x| x.seq_type == "Fix")
            .for_each(|x| {
                block_hash.insert(
                    x.idx.to_owned(),
                    (x.query_start.unwrap_or(0), x.query_end.unwrap_or(0)),
                );
                // seq_hash.insert(x.idx.to_owned(), );
            });

        let mut block_align_count = 0;
        let pre_variable_seq_bool = false;
        let mut block_list: Vec<(String, usize, usize)> = Vec::new();

        for ii in 0..self.block_align.len() {
            if let Some(x) = &self.block_align[ii] {
                // if let Some()
                match x.info.seq_type.as_str() {
                    "Fix" => {
                        let idx = x.info.idx.to_owned();
                        if x.align.is_some() {
                            block_align_count += 1;
                            let query_start = x.align.as_ref().unwrap().clone().query_start;
                            let query_end = x.align.as_ref().unwrap().clone().query_end;
                            // block_hash.insert(idx, (query_start, query_end));
                            // block_list.push((idx, query_start, query_end));
                            block_list[ii] = (idx, query_start, query_end);
                            if pre_variable_seq_bool {
                                block_list[ii - 1].2 = query_end;
                            }
                        }
                    }
                    "Variabe" => {}
                    _ => {
                        todo!()
                    }
                }
            }
        }
        // todo:
        let mut block_align_count = 0;
        self.block_align.iter().for_each(|x| {
            if let Some(x) = x {
                //   let ba=  x.to_abbr().unwrap();
                let idx = x.info.idx.to_owned();
                if x.align.is_some() {
                    block_align_count += 1;
                    let query_start = x.align.as_ref().unwrap().clone().query_start;
                    let query_end = x.align.as_ref().unwrap().clone().query_end;
                    block_hash.insert(idx, (query_start, query_end));
                }
            }
        });
        block_hash.iter().for_each(|(name, (start, end))| {
            let seq = &record.seq().to_vec()[*start.min(&seq_len)..*end.min(&seq_len)];
            let qual = &record.qual().to_vec()[*start.min(&seq_len)..*end.min(&seq_len)];
            let seq = Some(String::from_utf8(seq.to_vec()).unwrap());
            let qual = Some(String::from_utf8(qual.to_vec()).unwrap());
            let seqout = SeqOut::new(name.to_string(), None, seq, qual);
            seq_hash.insert(name.to_owned(), seqout);
        });

        // add fastq tag
        seq_hash.insert(
            "read".to_string(),
            SeqOut::new(
                (&record.id()).to_string(),
                (&record.desc().map(|s| s.to_string())).to_owned(),
                None,
                None,
            ),
        );
        if block_align_count > 0 {
            Some(seq_hash)
        } else {
            None
        }
    }

    /// tojson for wasm
    ///
    pub fn to_pretty(&self) -> ReadBlockAlignPretty {
        // todo!();
        let seq = std::str::from_utf8(self.record.seq()).unwrap();
        let read_name = self.record.id().to_string();
        let mut range_vec = vec![];
        self.block_align.iter().for_each(|block_align| {
            if let Some(block_align) = block_align {
                let block_name = &block_align.info.idx;
                let block_align_abbr = block_align.to_abbr();
                if let Some(block_align_abbr) = block_align_abbr {
                    let class_name = block_align_abbr.best_index; //.clone();

                    range_vec.push((
                        block_align_abbr.query_start..block_align_abbr.query_end,
                        block_name.clone(),
                    ));
                }
            }
        });
        let html = dna_to_spans(seq, &range_vec, "other");
        ReadBlockAlignPretty { read_name, html }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadBlockAlignPretty {
    read_name: String,
    html: String,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct SeqOut {
    name: String,
    desc: String,
    seq: String,
    qual: String,
}

pub fn to_str(s: Option<String>) -> String {
    if let Some(s) = s {
        s
    } else {
        "".to_string()
    }
}

impl SeqOut {
    pub fn new(
        name: String,
        desc: Option<String>,
        seq: Option<String>,
        qual: Option<String>,
    ) -> Self {
        SeqOut {
            name,
            desc: to_str(desc),
            seq: to_str(seq),
            qual: to_str(qual),
        }
    }
}

#[test]
fn test_block_align_read_with_insert() {
    let blockinfo_str = "idx	seq_type	fasta_seq_id	max_mismatch	query_start	query_end	seq_len	method
aa	Anchor	aa1,aa2	2				ANT
bb	Variable	bb1,bb2	2				ANT
cc	Anchor	cc1,cc2	2				ANT";
    // println!("{}", blockinfo_str);
    let fasta_file = ">aa1
ATCGATCGTAAAAA
>aa2
ATCCTAAATTACCA
>bb1
ATCGATAGTAAAAA
>bb2
ATCGCTCGTAAAAA
>cc1
ATCGATCTTAAAAA
>cc2
ATCGATCGTACAAA";
    let blockinfo_vec = get_block_info_fasta(blockinfo_str, fasta_file).unwrap();
    let read = b"CTCGATCGATCGTAAAAACGCCTCGCTATATCGTATCGATCGTACAAA";
    let block_align = block_align_read_with_insert(read, &blockinfo_vec);
    
    for (k, v) in block_align{
        // println!("{k}");
        println!("{k}: {}", v.unwrap().to_abbr().unwrap().to_str());
    }
    // dbg!(block_align);
}
