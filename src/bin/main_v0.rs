#![allow(dead_code)]
#![allow(unused)]
use bio::alignment::pairwise::{banded::Aligner, MatchFunc};
use bio::io::fastq::{self, Record};
// use bio_types::alignment::Alignment;
use clap::command;
use clap::Parser;
use crossbeam_channel::{unbounded, Receiver, Sender};
use once_cell::sync::Lazy;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};
use tracing::{debug, error, info, warn, Level};
use tracing_subscriber::fmt::format;

use legoseq::aligner::Alignment;
use legoseq::blockalign::ReadBlockAlign;
use legoseq::blockalign::{block_align_read, BlockAlign, BlockAlignAbbr};
use legoseq::blockinfo::{get_block_info, BlockInfo, BLOCKFLAGS};
use legoseq::utils::get_reader;

static CLI: Lazy<Cli> = Lazy::new(Cli::parse);

#[derive(Parser)]
#[command(version, author, about, long_about = None)]
struct Cli {
    /// Optional name to operate on
    /// Sets a custom block information file
    #[arg(long, value_name = "BLOCKINFO")]
    block_info: String,
    /// fastq file
    #[arg(long, value_name = "FILE")]
    fq1: String,
    /// fastq file, optional
    #[arg(long, value_name = "FILE")]
    fq2: Option<String>,
    /// threads
    #[arg(long, value_name = "FILE")]
    threads: usize,
    /// output directory
    #[arg(long, value_name = "PATH")]
    outdir: String,
    /// the prefix of output file
    #[arg(long, value_name = "PATH")]
    prefix: String,

    /// export block
    #[arg(long, value_name = "String")]
    export_blocks: Option<String>,
}

fn main() {
    let threads = &CLI.threads.to_owned();
    let outdir = &CLI.outdir;
    let r1_file = &CLI.fq1;
    let r2_file = &CLI.fq2;
    let block_info_file = &CLI.block_info;
    let prefix = &CLI.prefix;
    let export_blocks = &CLI.export_blocks;
    tracing_subscriber::fmt::init();
    info!("Start");

    rayon::ThreadPoolBuilder::new()
        .num_threads(*threads)
        .build_global()
        .unwrap();

    // output read info  statistics
    let outdir = Path::new(outdir);
    let read_info_file = outdir.join(format!("{}.{}", prefix, "read_info.stat.tsv"));
    info!(
        "write read information to file: {}",
        read_info_file.display()
    );
    let mut read_info_handle = Arc::new(Mutex::new(File::create(read_info_file).unwrap()));
    BLOCKFLAGS.lock().unwrap().iter().for_each(|(k, v)| {
        write!(read_info_handle.lock().unwrap(), "#idx:flag={}:{}\n", k, v);
    });

    let block_info_list = get_block_info(block_info_file);

    // export block
    // 命令行中指定了 --export_blocks 才会执行，同时会根据 --fq2 是否指定来导出fastq 文件
    let mut export_block_list: Option<Vec<String>> = None;
    let mut out_fq1 = None;
    let mut out_fq2 = None;

    if let Some(export_blocks) = export_blocks {
        export_block_list = Some(export_blocks.split('-').map(|x| x.to_string()).collect());
        let out_fq1_file = outdir.join(format!("{}.{}", prefix, "blocks.r1.fastq"));
        let out_fq2_file = outdir.join(format!("{}.{}", prefix, "blocks.r2.fastq"));
        info!(
            "write block seq information to fastq file: {}",
            out_fq1_file.display()
        );
        out_fq1 = Some(Arc::new(Mutex::new(
            fastq::Writer::to_file(out_fq1_file).unwrap(),
        )));
        out_fq2 = Some(Arc::new(Mutex::new(
            fastq::Writer::to_file(out_fq2_file).unwrap(),
        )));
    } else {
        // None
    };
    // let export_block_list: Option<Vec<String>> = export_blocks
    //     .as_ref()
    //     .map(|export_blocks| export_blocks.split('-').map(|x| x.to_string()).collect());
    // 统计所有 flag 的数目
    let flag_stat_hash: Arc<Mutex<HashMap<usize, usize>>> = Arc::new(Mutex::new(HashMap::new()));

    let record_r1 = fastq::Reader::new(get_reader(r1_file)).records();
    if let Some(r2_file) = r2_file {
        let record_r2 = fastq::Reader::new(get_reader(r2_file)).records();

        record_r1
            .into_iter()
            .zip(record_r2.into_iter())
            .par_bridge()
            .for_each(|(record_r1, record_r2)| {
                let record_r1 = record_r1.unwrap();
                let record_r2 = record_r2.unwrap();
                let read_name = record_r1.id();
                let read_seq = record_r1.seq();
                let read_block_align =
                    ReadBlockAlign::read_block_info(&record_r1, &block_info_list);
                // let block_align_list = read_block_align.block_align;

                let flag = read_block_align.get_block_flag();
                *flag_stat_hash.lock().unwrap().entry(flag).or_insert(0) += 1;

                let output_merge_str = read_block_align.get_block_str();

                writeln!(
                    read_info_handle.lock().unwrap(),
                    "{}\t{}\t{}",
                    read_name,
                    flag,
                    output_merge_str
                )
                .unwrap();

                // demultiplexed the fastq file
            });
    } else {
        record_r1.into_iter().par_bridge().for_each(|(record_r1)| {
            let record_r1 = record_r1.unwrap();
            let read_name = record_r1.id();
            let read_seq = record_r1.seq();
            let seq_len = read_seq.len();
            let read_block_align = ReadBlockAlign::read_block_info(&record_r1, &block_info_list);

            let flag = read_block_align.get_block_flag();
            *flag_stat_hash.lock().unwrap().entry(flag).or_insert(0) += 1;
            let output_merge_str = read_block_align.get_block_str();
            if let Some(export_block_list) = &export_block_list {
                let new_record =
                    read_block_align.get_new_record(&block_info_list, export_block_list);
                if new_record.is_some() {
                    let out_fq1 = out_fq1.clone();
                    out_fq1
                        .unwrap()
                        .lock()
                        .unwrap()
                        .write_record(&new_record.unwrap())
                        .unwrap();

                    // read2
                    let new_seq = &record_r1.seq().to_vec()[seq_len - 150..];
                    let new_qual = &record_r1.qual().to_vec()[seq_len - 150..];
                    let new_record2 =
                        Record::with_attrs(record_r1.id(), record_r1.desc(), new_seq, new_qual);
                    let out_fq2 = out_fq2.clone();
                    out_fq2
                        .unwrap()
                        .lock()
                        .unwrap()
                        .write_record(&new_record2)
                        .unwrap();
                }
            }

            write!(
                read_info_handle.lock().unwrap(),
                "{}\t{}\t{}\n",
                read_name,
                flag,
                output_merge_str
            )
            .unwrap();
            *flag_stat_hash.lock().unwrap().entry(flag).or_insert(0) += 1;
            // demultiplexed the fastq file
        })
    }

    let flag_stat_file = outdir.join(format!("{}.{}", prefix, "block_flag.stat.tsv"));
    info!(
        "write block flag stat to file: {}",
        flag_stat_file.display()
    );
    let mut flag_stat_handle = File::create(flag_stat_file).unwrap();

    BLOCKFLAGS.lock().unwrap().iter().for_each(|(k, v)| {
        write!(flag_stat_handle, "#idx:flag={}:{}\n", k, v);
    });
    flag_stat_hash.lock().unwrap().iter().for_each(|(k, v)| {
        write!(flag_stat_handle, "{}\t{}\n", k, v);
    });

    info!("End");
}
