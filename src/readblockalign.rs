use std::collections::HashMap;

use bio::io::fastq::Record;
use bio_types::sequence::SequenceRead;
use minijinja::value::Value;
use minijinja::Template;
use serde::{Deserialize, Serialize};

use crate::aligner::Alignment;
use crate::blockalign;
use crate::blockinfo::get_block_info_fasta;
use crate::utils::Strand;
use crate::utils::{check_vec_equal, dna_to_spans};
use crate::{blockalign::BlockAlign, blockinfo::BlockInfo};

/// save the fastq record and its block align information
#[derive(Debug, Clone)]
pub struct ReadBlockAlign {
    pub block_idx_list: Vec<String>,
    pub record: Record,
    pub block_align: HashMap<String, Option<BlockAlign>>, // blockname: block_align
    strand: Strand,
}

impl ReadBlockAlign {
    pub fn new(
        block_idx_list: &[String],
        record: &Record,
        block_align: &HashMap<String, Option<BlockAlign>>,
        strand: Strand,
    ) -> Self {
        ReadBlockAlign {
            block_idx_list: block_idx_list.to_vec(),
            record: record.clone(),
            block_align: block_align.clone(),
            strand: strand.clone(),
        }
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


impl ReadBlockAlign {
    /// mapping read againt the block sequence
    /// return ReadBlockAlign struct
    // pub fn read_block_info(record: &Record, block_info_list: &[BlockInfo]) -> Self {
    //     let read_name = record.id();
    //     let read_seq = record.seq();
    //     let block_align = block_align_read_with_insert(&record, block_info_list);
    //     let block_idx_list: Vec<String> =
    //         block_info_list.iter().map(|x| x.idx.to_string()).collect();

    //     Self {
    //         record: record.clone(),
    //         block_align,
    //         block_idx_list,
    //     }
    // }

    /// 获取 blockinfo 各个 block 的比对情况，包括 fix 和 variable
    pub fn read_block_info(record: &Record, block_info_list: &[BlockInfo]) -> Self {
        let read = record.seq();
        // let aligner2 = &aligner;
        let mut block_align_hash: HashMap<String, Option<BlockAlign>> = HashMap::new();
        let mut pre_query_end: usize = 0;
        let read_len = read.len();

        let mut strand_vec: Vec<Strand> = vec![];
        // 处理 anchor/fix 序列
        for block_info in block_info_list.iter().filter(|x| x.seq_type == "Fix") {
            let block_info = block_info.clone();
            let idx = block_info.idx.clone();
            let align: Option<Alignment> = block_info.aligner.clone().and_then(|x| x.align(read));
            let ba = align.map(|x| BlockAlign::new(&block_info, &x));
            block_align_hash.insert(idx, ba.clone());
            if let Some(ba) = ba {
                strand_vec.push(ba.get_query_strand().unwrap());
            } else {
                strand_vec.push(Strand::Ambiguous);
            }
        }
        // 判断正负链
        let strand_vec_equal = check_vec_equal(&strand_vec);
        let strand = if strand_vec_equal {
            strand_vec.first().unwrap_or(&Strand::Ambiguous)
        } else {
            &Strand::Ambiguous
        };

        // 处理可变序列
        let block_info_len = block_info_list.len();
        for block_ii in 0..block_info_len {
            let mut block_info = block_info_list[block_ii].clone();
            let idx = block_info.idx.to_string();
            if block_info.seq_type != "Variable" {
                continue;
            }

            let mut query_start: Option<usize> = None;
            // 开头模块
            if block_ii == 0 {
                // query_start = Some(0);
                query_start = if strand.is_reverse().unwrap(){ Some(read_len)} else {Some(0)};
            } else {
                let pre_block_info = &block_info_list[block_ii - 1];
                let ba = block_align_hash.get(&pre_block_info.idx).unwrap();
                if let Some(ba) = ba {
                    query_start = ba.get_query_end().map(|pos| pos + 1);
                } else {
                    block_align_hash.insert(idx, None);
                    continue;
                }
            }

            // 末尾模块
            let mut query_end: Option<usize> = None;
            if block_ii == (block_info_len - 1) {
                // query_end = Some(read_len);
                query_end = if strand.is_reverse().unwrap(){ Some(0)} else {Some(read_len)};
            }
            let next_block_info = &block_info_list[block_ii + 1];
            let ba = block_align_hash.get(&next_block_info.idx).unwrap().clone();

            if let Some(ba) = ba {
                // query_start = ba.get_query_end().map(|pos| pos + 1);
                query_end = ba.get_query_start().map(|pos| pos - 1);
            } else {
                block_align_hash.insert(idx, None);
                continue;
            }
            if strand.is_reverse().unwrap(){
                (query_start, query_end) = (query_end, query_start);
            }
            // 存在异常
            // if query_start > query_end {
            //     block_align_hash.insert(idx, None);
            //     continue;
            // }
            // let
            let align = Alignment {
                best_index: "".to_string(),
                index_start: 0,
                index_end: 0,
                query_start: query_start.unwrap(),
                query_end: query_end.unwrap(),
                n_match: 0,
                strand: strand.clone(),
                operations: None,
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

        let block_idx_list: Vec<String> =
            block_info_list.iter().map(|x| x.idx.to_string()).collect();

        Self {
            record: record.clone(),
            block_align: block_align_hash,
            block_idx_list,
            strand: strand.clone(),
        }
    }

    /// get the block flag of the read
    pub fn get_block_flag(&self) -> usize {
        let mut flag = 0;
        for idx in self.block_idx_list.iter() {
            let ba = self
                .block_align
                .get(idx)
                .unwrap_or_else(|| panic!("not found the idx: {}", idx));
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
        // for ba in self.block_align.clone().iter() {
        for idx in self.block_idx_list.iter() {
            let ba = self
                .block_align
                .get(idx)
                .expect(&format!("not found the idx: {}", idx));
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
        self.block_align
            .iter()
            .filter(|(x, y)| export_block.contains(x))
            .for_each(|(_, x)| {
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
    /// 基于 readblockalign 生成 hashmap 用于对应模版
    pub fn get_seq_hashmap(&self) -> Option<HashMap<String, JinjaSeq>> {
        if self.strand == Strand::Ambiguous {
            return None;
        }
        let record = &self.record;
        let mut seq_hash: HashMap<String, JinjaSeq> = HashMap::new();
        let read_name = String::from_utf8(self.record.name().to_vec()).unwrap_or("".to_string());
        // todo:
        let mut block_align_count = 0;
        self.block_align
            .iter()
            .for_each(|(block_name, block_align)| {
                if let Some(block_align) = block_align {
                    //   let ba=  x.to_abbr().unwrap();
                    let idx = block_align.info.idx.to_owned();
                    let best_index = &block_align.best_index;
                    if block_align.align.is_some() {
                        block_align_count += 1;

                        let query_start = block_align.get_query_start().unwrap_or(0);
                        let query_end = block_align.get_query_end().unwrap_or(0);
                        let strand = block_align.get_query_strand().unwrap();
                        // if query_end < query_start {
                        //     // dbg!(&block_align);
                        //     // dbg!(self.record.id());
                        //     return None
                        // }
                        // 如果 variable 两端的 fix 是颠倒的，会有混乱，需要修改
                        let seqout =
                            JinjaSeq::new(best_index, record, query_start, query_end, &strand);
                        seq_hash.insert(idx, seqout);
                    }
                } else {
                    // dbg!(format!("not found {block_name}"));
                    seq_hash.insert(block_name.to_string(), JinjaSeq::default());
                }
            });
        // add fastq tag
        seq_hash.insert(
            "read".to_string(),
            JinjaSeq {
                name: String::from_utf8(self.record.name().to_vec()).unwrap_or("".to_string()),
                ..Default::default()
            },
        );
        if block_align_count > 0 {
            Some(seq_hash)
        } else {
            None
        }
    }

    pub fn template_str(
        &self,
        template: &Template<'_, '_>,
        block_info_list: &[BlockInfo],
    ) -> Option<String> {
        let seq_hash = self.get_seq_hashmap();
        let out = seq_hash.as_ref().map(|x| {
            let ctx = Value::from_serializable(x);
            template.render(ctx).expect("无法渲染模板")
        });
        out
    }

    /// tojson for wasm
    ///
    pub fn to_pretty(&self) -> ReadBlockAlignPretty {
        // todo!();
        let seq = std::str::from_utf8(self.record.seq()).unwrap();
        let read_name = self.record.id().to_string();
        let mut range_vec = vec![];
        self.block_align.iter().for_each(|(_, block_align)| {
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

#[derive(Debug, Deserialize, Serialize, Default)]
pub struct JinjaSeq {
    name: String,
    desc: String,
    seq: String,
    qual: String,
    start: usize,
    end: usize,
    strand: char,
}

pub fn to_str(s: Option<String>) -> String {
    if let Some(s) = s {
        s
    } else {
        "".to_string()
    }
}

impl JinjaSeq {
    pub fn new(name: &str, record: &Record, start: usize, end: usize, strand: &Strand) -> Self {
        let mut new_end = end;
        let mut new_start = start;
        // if strand.is_reverse().is_some(){
        //     new_end = start;
        //     new_start = end;
        // }else{
        //     new_end = end;
        //     new_start = start;
        // }
        
        let seq_len = record.seq().len();
        let seq = &record.seq().to_vec()[new_start.min(seq_len)..new_end.min(seq_len)];
        let qual = &record.qual().to_vec()[new_start.min(seq_len)..new_end.min(seq_len)];
        let desc = &record.desc();
        let seq = Some(String::from_utf8(seq.to_vec()).unwrap());
        let qual = Some(String::from_utf8(qual.to_vec()).unwrap());
        // let seqout = SeqOut::new(name.to_string(), None, seq, qual);
        JinjaSeq {
            name: name.to_string(),
            desc: to_str(desc.map(|x| x.to_string())),
            seq: to_str(seq),
            qual: to_str(qual),
            start,
            end,
            strand: strand.to_char(),
        }
    }
}

#[test]
fn test_block_align_read_with_insert() {
    let blockinfo_str = "idx	seq_type	fasta_seq_id	max_mismatch	query_start	query_end	seq_len	method
aa	Fix	aa1,aa2	2				ANT
bb	Variable	bb1,bb2	2				ANT
cc	Fix	cc1,cc2	2				ANT";
    // println!("{}", blockinfo_str);
    let fasta_file = ">aa1
AAAAAAAAAAAAA
>aa2
ATCCTAAATTACCA
>bb1
TTTTTTTTTTTTTT
>bb2
ATCGCTCGTAAAAA
>cc1
ATCGATCTTAAAAA
>cc2
CCCCCCCCCCCCC";
    let blockinfo_vec = get_block_info_fasta(blockinfo_str, fasta_file).unwrap();
    let read =
        b"CTGGGGGGGGGGGGGGCGATCGATCGTAAACCCCCCCCCCCCCCAACGCTTTTTTTTTTTTTTCTCGCTATATCGTATCGATGTAC";
    let record = Record::with_attrs("read01_rev", None, read, read);
    let read_block_align = ReadBlockAlign::read_block_info(&record, &blockinfo_vec);
    let block_align = &read_block_align.block_align;
    for (k, v) in block_align {
        println!("{k}");
        if let Some(v) = v {
            println!("{k}: {}", v.to_abbr().unwrap().to_str());
        }
    }
    let seq_hashmap = read_block_align.get_seq_hashmap().unwrap();
    dbg!(seq_hashmap);
    // dbg!(block_align);
}
