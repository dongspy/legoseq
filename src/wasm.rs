use wasm_bindgen::prelude::*;
use crate::blockinfo::get_block_info_fasta;
use crate::blockalign::ReadBlockAlign;
use bio::io::fastq;

#[wasm_bindgen]
pub fn add(a: i32, b: i32) -> i32 {
    a + b
}

#[wasm_bindgen]
pub fn test(seq_info: &str, blockinfo_str: &str, fasta_file: &str) -> String {

    "abcd".to_string()
}

pub fn add_text(a: &str, b: &str) -> String {
    format!("{} and {}", a, b)
}

#[wasm_bindgen]
pub fn read_align(seq_info: &str, blockinfo_str: &str, fasta_file: &str) -> String {
    // let seq_info: Vec<&str> = seq_info.split('|').collect();
    let block_info_vec = get_block_info_fasta(blockinfo_str, fasta_file);    
    let records = fastq::Reader::new(seq_info.as_bytes()).records();
    
    let read_block_align =
                    ReadBlockAlign::read_block_info(&records.into_iter().next().unwrap().unwrap(), &block_info_vec.unwrap());
                // let block_align_list = read_block_align.block_align;
    // "abcd".to_string()
    read_block_align.get_block_str()
}