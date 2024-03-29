use crate::blockinfo::get_block_info_fasta;
use crate::readblockalign::ReadBlockAlign;
use bio::io::fastq::{self, Record};
use wasm_bindgen::prelude::*;
use wasm_bindgen_test::*;

#[wasm_bindgen]
pub fn read_align(seq_info: &str, blockinfo_str: &str, fasta_file: &str) -> JsValue {
    // let seq_info: Vec<&str> = seq_info.split('|').collect();

    let block_info_vec = get_block_info_fasta(blockinfo_str, fasta_file).expect("blockinfo error");
    console_log!("the length of block_info_vec is {}", block_info_vec.len());
    let records = fastq::Reader::new(seq_info.as_bytes()).records();
    let mut read_count = 0;
    let mut read_block_align_vec = vec![];
    records.into_iter().for_each(|record| {
        let record = record.unwrap();
        let read_block_align = ReadBlockAlign::read_block_info(&record, &block_info_vec);
        read_count += 1;
        read_block_align_vec.push(read_block_align.to_pretty());
    });
    console_log!("the length of fastq record {}", read_count);
    serde_wasm_bindgen::to_value(&read_block_align_vec).unwrap()
}

// #[wasm_bindgen_test]
#[test]
fn test_read_align() {
    let seq_info = "@HISEQ_HU01:89:H7YRLADXX:1:1101:1573:2113 1:N:0:ATCACG
ATCGATCGTAAAAA
+
ATCGATCGTAAAAA";
    let blockinfo = "idx	seq_type	fasta_seq_id	max_mismatch	query_start	query_end	seq_len	method
aa	Anchor	aa1,aa2	2				ANT
bb	Anchor	bb1,bb2	2				ANT
cc	Anchor	cc1,cc2	2				ANT";
    let fasta_info = ">aa1
ATCGATCGTAAAAA
>aa2
ATCGATCGTAAAAA
>bb1
ATCGATCGTAAAAA
>bb2
ATCGATCGTAAAAA
>cc1
ATCGATCGTAAAAA
>cc2
ATCGATCGTAAAAA";
    let out = read_align(seq_info, blockinfo, fasta_info);
    // println!("{}", &out);
    dbg!(out);
}

#[wasm_bindgen]
pub struct FqRecord {
    records: Vec<Record>,
}

#[wasm_bindgen]
impl FqRecord {
    #[wasm_bindgen(constructor)]
    pub fn new(seq_info: &str) -> FqRecord {
        let records = fastq::Reader::new(seq_info.as_bytes()).records();
        let record_vec: Vec<Record> = records.into_iter().map(|x| x.unwrap()).collect();
        FqRecord {
            records: record_vec,
        }
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }
}
