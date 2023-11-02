use bio::io::{fastq, fasta};
use bio_types::sequence::SequenceRead;
use clap::command;
use clap::Parser;
use dashmap::DashMap;
use minijinja::{Environment, Template};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};
use tracing::info;

use super::blockinfo::{get_block_info_fasta_from_file, BLOCKFLAGS};
use super::readblockalign::ReadBlockAlign;
use super::utils::get_reader;
use super::blockinfo::BlockInfo;

pub trait Record {
    fn id(&self) -> &str;
    fn seq(&self) -> &[u8];
    fn desc(&self) -> Option<&str>;
    fn qual(&self) -> &[u8];
}

impl Record for fasta::Record {
    fn id(&self) -> &str {
        self.id()
        
    }
    fn seq(&self) -> &[u8] {
        self.seq()
    }
    fn qual(&self) -> &[u8] {
        &[]
    }

    fn desc(&self) -> Option<&str> {
        self.desc()
    }
}

impl Record for fastq::Record {
    fn id(&self) -> &str {
        self.id()
    }
    fn seq(&self) -> &[u8] {
        self.seq()
    }
    fn qual(&self) -> &[u8] {
        self.qual()
    }


    fn desc(&self) -> Option<&str> {
        self.desc()
    }
}

#[derive(Clone)]
pub struct FastxRecord<R>{
    pub record1: R,
    pub record2: Option<R>
}



pub fn process_record<R:Record + Clone>(
    record: FastxRecord<R>, 
    block_info_list:&[BlockInfo], 
    prefix:&str,
    outdir: &Path,
    barcode_handle_hash: DashMap<String, Arc<Mutex<Vec<File>>>>,
    template: Template<'_, '_>,
    out_fq_handle_vec: Arc<Mutex<Vec<File>>>,
    ud_fq_handle_vec: Arc<Mutex<Vec<File>>>,
    read_info_handle: Arc<Mutex<File>>,
    flag_stat_hash: Arc<Mutex<HashMap<usize, usize>>>
) {
    // let _ = record;

    let record_r1 = record.record1;
    let record_r2 = record.record2.unwrap();
    let read_name = record_r1.id();
    let read_block_align =
        ReadBlockAlign::read_block_info(&record_r1, &block_info_list);
    let flag = read_block_align.get_block_flag();
    let best_index_vec = read_block_align.get_best_index();
    // let outdir = outdir.to_string();
    // if barcode index existed, demultiplex
    if !best_index_vec.is_empty() {
        let best_index_str = best_index_vec.join("_");
        let barcode_file =
            outdir.join(format!("{}.{}.{}", prefix, best_index_str, "r1.fastq"));
        let barcode_file2 =
            outdir.join(format!("{}.{}.{}", prefix, best_index_str, "r2.fastq"));
        
        let barcode_handle = barcode_handle_hash
            .entry(barcode_file.to_str().unwrap().to_string())
            .or_insert_with(|| {
                Arc::new(Mutex::new(vec![
                    File::create(barcode_file).unwrap(),
                    File::create(barcode_file2).unwrap(),
                ]))
            });
        let template_str = read_block_align.template_str(&template);
        if let Some(template_str) = template_str {
            writeln!(
                barcode_handle.lock().unwrap().get(0).unwrap(),
                "{}",
                template_str
            )
            .unwrap();
            writeln!(
                barcode_handle.lock().unwrap().get(1).unwrap(),
                "@{} {}\n{}\n+\n{}",
                record_r2.id(),
                // record_r2.desc().expect(""),
                "",
                String::from_utf8_lossy(record_r2.seq()),
                String::from_utf8_lossy(record_r2.qual()) //    String::from
            )
            .unwrap();
        } else {
            let out_fq_handle = ud_fq_handle_vec.lock().unwrap();
            writeln!(
                out_fq_handle.get(0).unwrap(),
                "@{} {}\n{}\n+\n{}",
                record_r1.id(),
                record_r1.desc().unwrap_or(""),
                String::from_utf8_lossy(record_r1.seq()),
                String::from_utf8_lossy(record_r1.qual())
            )
            .unwrap();
            writeln!(
                out_fq_handle.get(1).unwrap(),
                "@{} {}\n{}\n+\n{}",
                record_r2.id(),
                record_r2.desc().unwrap_or(""),
                // "",
                String::from_utf8_lossy(record_r2.seq()),
                String::from_utf8_lossy(record_r2.qual())
            )
            .unwrap();
        }
    } else {
        // 没有 barcode 拆分
        // export to file based on the jinja template
        let template_str = read_block_align.template_str(&template);

        if let Some(template_str) = template_str {
            let out_fq_handle = out_fq_handle_vec.lock().unwrap();
            // println!("{}", &template_str);
            writeln!(out_fq_handle.get(0).unwrap(), "{}", template_str).unwrap();
            writeln!(
                out_fq_handle.get(1).unwrap(),
                "@{} {}\n{}\n+\n{}",
                record_r2.id(),
                record_r2.desc().unwrap(),
                String::from_utf8_lossy(record_r2.seq()),
                String::from_utf8_lossy(record_r2.qual())
            )
            .unwrap();
        } else {
            let out_fq_handle = ud_fq_handle_vec.lock().unwrap();
            writeln!(
                out_fq_handle.get(0).unwrap(),
                "@{} {}\n{}\n+\n{}",
                record_r1.id(),
                record_r1.desc().unwrap_or(""),
                String::from_utf8_lossy(record_r1.seq()),
                String::from_utf8_lossy(record_r1.qual())
            )
            .unwrap();
            writeln!(
                out_fq_handle.get(1).unwrap(),
                "@{} {}\n{}\n+\n{}",
                record_r2.id(),
                record_r2.desc().unwrap_or(""),
                // "",
                String::from_utf8_lossy(record_r2.seq()),
                String::from_utf8_lossy(record_r2.qual())
            )
            .unwrap();
        }
    }

    let output_merge_str = read_block_align.get_block_str();

    writeln!(
        read_info_handle.lock().unwrap(),
        "{}\t{}\t{}",
        read_name,
        flag,
        output_merge_str
    )
    .unwrap();
    *flag_stat_hash.lock().unwrap().entry(flag).or_insert(0) += 1;
}
