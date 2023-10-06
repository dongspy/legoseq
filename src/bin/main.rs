// #![allow(dead_code)]
// #![allow(unused)]
use bio::alignment::pairwise::{banded::Aligner, MatchFunc};
use bio::io::fastq::{self, Record};
// use bio_types::alignment::Alignment;
use clap::command;
use clap::Parser;
use crossbeam_channel::{unbounded, Receiver, Sender};
use minijinja::{context, value::Value, Environment, Template};
use once_cell::sync::Lazy;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};
use tracing::{debug, error, info, warn, Level};
use tracing_subscriber::fmt::format;

use legoseq::aligner::Alignment;
// use legoseq::blockalign::{BlockAlign, BlockAlignAbbr};
use legoseq::blockinfo::{get_block_info_fasta_from_file, BlockInfo, BLOCKFLAGS};
use legoseq::readblockalign::ReadBlockAlign;
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
    /// fasta file
    #[arg(long, value_name = "FILE")]
    fasta: String,
    /// threads
    #[arg(long, value_name = "FILE")]
    threads: usize,
    /// output directory
    #[arg(long, value_name = "PATH")]
    outdir: String,
    /// the prefix of output file
    #[arg(long, value_name = "PATH")]
    prefix: String,

    /// the template file
    #[arg(long, value_name = "PATH")]
    template: Option<String>,
}

fn main() {
    let threads = &CLI.threads.to_owned();
    let outdir = &CLI.outdir;
    let r1_file = &CLI.fq1;
    let r2_file = &CLI.fq2;
    let fasta_file = &CLI.fasta;
    let block_info_file = &CLI.block_info;
    let prefix = &CLI.prefix;
    // let export_blocks = &CLI.export_blocks;
    let template: &Option<String> = &CLI.template;
    tracing_subscriber::fmt::init();
    info!("Start");

    rayon::ThreadPoolBuilder::new()
        .num_threads(*threads)
        .build_global()
        .unwrap();

    // output read info  statistics
    let outdir = Path::new(outdir);
    let read_info_file = outdir.join(format!("{}.{}", prefix, "read_info.stat.tsv"));
    let out_fq_file = outdir.join(format!("{}.{}", prefix, "template.out.fastq"));
    info!(
        "write read information to file: {}",
        read_info_file.display()
    );

    let block_info_list = get_block_info_fasta_from_file(block_info_file, fasta_file).unwrap();

    let mut read_info_handle = Arc::new(Mutex::new(File::create(read_info_file).unwrap()));
    let mut out_fq_handle = Arc::new(Mutex::new(File::create(out_fq_file.clone()).unwrap()));

    BLOCKFLAGS.lock().unwrap().iter().for_each(|(k, v)| {
        writeln!(read_info_handle.lock().unwrap(), "#idx:flag={}:{}", k, v);
    });


    //minijinja
    let template_string = fs::read_to_string(template.clone().unwrap()).expect("无法读取模板文件");
    // 创建一个新的 MiniJinja 环境
    let env = Environment::new();
    // 从字符串创建一个模板
    let template: Template<'_, '_> = env
        .template_from_str(&template_string)
        .expect("无法从字符串创建模板");

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
                    "{}\t{}\t{}1",
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
            let read_block_align = ReadBlockAlign::read_block_info(&record_r1, &block_info_list);
            let flag = read_block_align.get_block_flag();
            let output_merge_str = read_block_align.get_block_str();

            // export to file based on the jinja template
            let template_str = read_block_align.template_str(&template);
            if let Some(template_str) = template_str {
                writeln!(out_fq_handle.lock().unwrap(), "{}", template_str).unwrap();
            }

            writeln!(
                read_info_handle.lock().unwrap(),
                "{}\t{}\t{}",
                read_name,
                flag,
                output_merge_str
            )
            .unwrap();
            *flag_stat_hash.lock().unwrap().entry(flag).or_insert(0) += 1;
        })
    }

    let flag_stat_file = outdir.join(format!("{}.{}", prefix, "block_flag.stat.tsv"));
    info!(
        "write block flag stat to file: {}",
        flag_stat_file.display()
    );
    info!(
        "write template information to file: {}",
        out_fq_file.display()
    );
    let mut flag_stat_handle = File::create(flag_stat_file).unwrap();

    BLOCKFLAGS.lock().unwrap().iter().for_each(|(k, v)| {
        writeln!(flag_stat_handle, "#idx:flag={}:{}", k, v);
    });
    flag_stat_hash.lock().unwrap().iter().for_each(|(k, v)| {
        writeln!(flag_stat_handle, "{}\t{}", k, v);
    });

    info!("End");
}
