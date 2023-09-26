#![allow(unused)]
use crate::{
    blockinfo::{AlignMethod, BlockInfo},
    utils::Strand,
};
use antlib::{
    align::{align_read as ant_align, AlignOpts as AntAlignOpts},
    index::Index as AntIndex,
};
use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::AlignmentOperation::{self, *};
use std::{
    collections::HashMap,
    fmt::Debug,
    path::{Path, PathBuf},
};

use crate::utils::{read_fasta, revcomp};

#[derive(Clone)]
pub struct AntAligner {
    index: AntIndex,
    opts: AntAlignOpts,
    max_mismatch: usize,
}

impl Debug for AntAligner {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct BandedAligner {
    name: String,
    gap_open: i32,
    gap_extend: i32,
    // match_fn: F,
    mmatch: i32,
    mismatch: i32,
    k: usize,
    w: usize,
    seq: Vec<u8>,
    max_mismatch: usize,
}

#[derive(Debug, Clone)]
pub struct Alignment {
    pub best_index: String,
    pub index_start: usize,
    pub index_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub n_match: usize,
    pub strand: Strand,
    pub operations: Option<Vec<AlignmentOperation>>,
}

impl Alignment {
    pub fn to_str(&self) -> String {
        format!(
            "{}:{}:{}:{}:{}:{}:{:#?}",
            self.best_index,
            self.index_start,
            self.index_end,
            self.query_start,
            self.query_end,
            self.n_match,
            self.strand,
        )
    }
}

pub trait Align {
    fn align(&self, seq: &[u8]) -> Option<Alignment>;
}

impl Align for AntAligner {
    fn align(&self, seq: &[u8]) -> Option<Alignment> {
        // self.index.align(x, y)
        // let index = Index::create_from_files("test/index.fa", 1, 1).unwrap();
        let binding = ant_align(&self.index, seq, &self.opts);
        let align = binding.get(0);
        if let Some(align) = align {
            let gx_aln = align.gx_aln.clone();
            let n_match = align
                .gx_aln
                .operations
                .iter()
                .filter(|&&x| x == Match)
                .count();

            // dbg!(self.max_mismatch);
            if (gx_aln.ylen - n_match) > self.max_mismatch {
                return None;
            }
            // if n_match < self.max
            Some(Alignment {
                best_index: align.ref_name.clone(),
                index_start: align.gx_aln.ystart,
                index_end: align.gx_aln.yend,
                query_start: align.gx_aln.xstart,
                query_end: align.gx_aln.xend,
                n_match,
                strand: if align.strand {
                    Strand::Plus
                } else {
                    Strand::Minus
                },
                operations: Some(gx_aln.operations),
            })
        } else {
            None
        }
    }
}

impl Align for BandedAligner {
    fn align(&self, seq: &[u8]) -> Option<Alignment> {
        let score = |a: u8, b: u8| if a == b { self.mmatch } else { self.mismatch };
        let mut aligner = Aligner::new(self.gap_open, self.gap_extend, score, self.k, self.w);
        let forward_aln = aligner.local(&self.seq, seq);
        // let n_match = aln.
        let forward_n_match = forward_aln
            .operations
            .iter()
            .filter(|&&x| x == Match)
            .count();
        let seq_len = forward_aln.xlen;
        if forward_n_match == seq_len {
            let aln = forward_aln;
            return Some(Alignment {
                best_index: (self.name).to_owned(),
                index_start: aln.xstart,
                index_end: aln.xend,
                query_start: aln.ystart,
                query_end: aln.yend,
                n_match: forward_n_match,
                strand: Strand::Plus,
                operations: Some(aln.operations),
            });
        }
        let seq = revcomp(seq);
        let revcom_aln = aligner.local(&self.seq, &seq);
        let revcom_n_match = revcom_aln
            .operations
            .iter()
            .filter(|&&x| x == Match)
            .count();

        let (aln, n_match, strand) = if (revcom_n_match > forward_n_match) {
            (revcom_aln, revcom_n_match, Strand::Minus)
        } else {
            (forward_aln, forward_n_match, Strand::Plus)
        };

        if (seq_len - n_match) > self.max_mismatch {
            None
        } else {
            Some(Alignment {
                best_index: (self.name).to_owned(),
                index_start: aln.xstart,
                index_end: aln.xend,
                query_start: aln.ystart,
                query_end: aln.yend,
                n_match,
                strand,
                operations: Some(aln.operations),
            })
        }
    }
}

#[derive(Clone, Debug)]
pub enum BAligner {
    BandedAligner(BandedAligner),
    AntAligner(AntAligner),
}

impl BAligner {
    // impl BAligner {
    pub fn new(
        method: AlignMethod,
        seq_hash: &HashMap<String, Vec<u8>>,
        max_mismatch: usize,
    ) -> BAligner {
        // let seq_hash = read_fasta(fasta_file).unwrap();

        let (seq_name, seq) = seq_hash.iter().next().unwrap();
        match method {
            AlignMethod::SW => {
                let k = 3; // kmer match length
                let w = 5; // Window size for creating the band
                           // let mut aligner = Aligner::new(-5, -1, score, k, w);
                let bandedaligner = BandedAligner {
                    gap_open: -2,
                    gap_extend: -1,
                    k,
                    w,
                    seq: seq.clone(),
                    name: seq_name.to_owned(),
                    mmatch: 1,
                    mismatch: -1,
                    max_mismatch,
                };
                // Box::new(bandedaligner)
                BAligner::BandedAligner(bandedaligner)
            }
            AlignMethod::ANT => {
                // let seq_hash = read_fasta(fasta_file).unwrap();
                // let index = AntIndex::create_from_files(fasta_file, 1, 1).unwrap();
                let index = AntIndex::create_from_hashmap(seq_hash, 1, 1).unwrap();
                let align_opts = AntAlignOpts {
                    min_seed_len: 11,
                    min_match_counts_percent: 0.0,
                    min_aln_score: 0,
                    multimap_score_range: 20,
                };
                // let align = ant_align(&index, read, &align_opts);
                let antaligner = AntAligner {
                    index: index,
                    opts: align_opts,
                    max_mismatch,
                };
                BAligner::AntAligner(antaligner)
            }
        }
    }

    pub fn align(&self, seq: &[u8]) -> Option<Alignment> {
        // self.aligner.align(seq)
        match self {
            BAligner::BandedAligner(aligner) => aligner.align(seq),
            BAligner::AntAligner(aligner) => aligner.align(seq),
        }
    }
}

#[test]
fn test_baligner() {
    use std::time::{Duration, Instant};
    let now = Instant::now();
    let seq_hash = read_fasta("test/test.fasta").expect("read fasta error");
    let alingner = BAligner::new(AlignMethod::ANT, &seq_hash, 5);
    // aligner
    // let aligner = get_aligner(AlignMethod::ANT, "test/index.fa");
    let read = b"CACAAAGACAAAAAAAAAAAACCAACAACTACTT";
    // let read = b"CACAAAGACACCAACAACTACTT";
    let aln = alingner.align(read);
    // for _ in 0..10000{
    //     let aln = aligner.align(read);
    // }
    println!("ANT: spend {}", now.elapsed().as_secs());

    dbg!(aln);

    // let now = Instant::now();
    // let aligner = BAligner::new(AlignMethod::SW, "test/index.fa", 2);
    // let read =
    //     b"ATCGATCGATCATCATCTCCCTATATATCGATCGATCATCATCGATCGATCATCATCGATCGATCATCATCGATCGATCATC";
    // let aln = aligner.align(read);
    // // for _ in 0..10000{
    // //     let aln = aligner.align(read);
    // // }
    // println!("SW: spend {}", now.elapsed().as_secs());
    // dbg!(aln);
}

#[test]
fn test_bio_banded_align() {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let k = 4; // kmer match length
    let w = 5; // Window size for creating the band
    let mut aligner = Aligner::new(-5, -1, score, k, w);
    let alignment = aligner.local(b"ATCGATCTG", b"ATCGATCTG");
    dbg!(alignment);
}

#[test]
fn test_ant() {
    let read = b"AATGTCTATGTACATACTTGACTGGTTTCATCTGCTAATGATTGCAGCAACCACAAGATCTACACCACAAAGACACCAACAACTACTTCACTCTTTCCCTACACAAGCACTTCTTAAGATGTGTGAGTACAGGTTTCATCAATAATCATTTCTTATATGAGTGCCTCATTACATGCAGTATTTATACTAAGCATTTACCATCTTAGCTTCTATCAAAATTATGGTATATCACTCACACCTCATGTCCTCCCCTTTACTATGCCTGAAGGAATAATACTATCAGTGTTCCCATTATAGCTACTCTCATGACCCTAGACACCCACTTTTCCCTCTTAGCCAATATTGTGCCTATTCCCATTACTAAGTCTTTTGCACACTAGAAAGCAGCAGTTGGCCTACCCTCTAGTCTCAATCTTCCAACACATGGCCTCGACAGATGAGAAAGAGCACACATCTGAACTTCCAGTTCACATATTTTCATCAGAATGAATCCTTGTCTCGTATGCCATTCTTCTTGCAGCCAATCATC";
    let align_opts = AntAlignOpts {
        min_seed_len: 11,
        min_match_counts_percent: 0.0,
        min_aln_score: 0,
        multimap_score_range: 0,
        // intron_mode: false,
    };
    let index = AntIndex::create_from_files("/Users/pidong/tmp/ont_p5.index.fa", 1, 1).unwrap();
    let align = ant_align(&index, read, &align_opts);
    dbg!(align);
}
