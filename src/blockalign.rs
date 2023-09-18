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

impl BlockAlign {
    pub fn new(info: &BlockInfo, align: &Alignment) -> BlockAlign {
        BlockAlign {
            info: info.clone(),
            align: Some(align.clone()),
            n_match: align.n_match,
            best_index: align.best_index.clone(),
        }
    }

    pub fn get_query_start(&self) -> Option<usize> {
        self.align.as_ref().map(|align| align.query_start)
    }

    pub fn get_query_end(&self) -> Option<usize> {
        self.align.as_ref().map(|align| align.query_end)
    }

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

static OFFSET: usize = 4;
