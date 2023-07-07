#![allow(unused)]
use std::collections::HashMap;

use bio::io::fastq::Record;
use rayon::vec;
use serde::{Deserialize, Serialize};
use tracing::{debug, error, info, warn, Level};

use crate::aligner::Alignment;
use crate::blockinfo::BlockInfo;
use crate::utils::dna_to_spans;

#[derive(Debug, Clone)]
pub struct BlockAlign {
    pub info: BlockInfo,
    pub align: Option<Alignment>,
    pub n_match: usize,
    pub best_index: String,
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

    /// generate the new fastq record that the sequcne only include the export_block which is in the blockinfo
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
