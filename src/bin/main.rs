#![allow(dead_code)]
#![allow(unused)]
use bio::io::fasta;
use bio::io::fastq;
use clap::command;
use clap::Parser;
use dashmap::DashMap;
use legoseq::record::process_record_pair;
use legoseq::record::process_record_single;
use minijinja::{Environment, Template};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};
use tracing::info;

use legoseq::blockinfo::{get_block_info_fasta_from_file, BLOCKFLAGS};
use legoseq::utils::get_reader;

#[derive(Parser)]
#[command(version, author, about, long_about = None)]
struct Cli {
    /// fastq file
    #[arg(long, value_name = "FILE")]
    in1: String,
    /// fastq file, optional
    #[arg(long, value_name = "FILE")]
    in2: Option<String>,
    /// input type, fasta or fastq, default is fastq
    #[arg(long, value_name = "INPUT_TYPE", default_value = "fastq")]
    input_type: Option<String>,
    /// fasta file
    #[arg(long, value_name = "FILE")]
    fasta: String,
    /// custom block information file
    #[arg(long, value_name = "BLOCKINFO")]
    block_info: String,
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
    /// the output file extension, the default is input_type
    #[arg(long, value_name = "EXT")]
    ext: Option<String>,
}

fn main() {
    let cli = Cli::parse();
    let threads = &cli.threads.to_owned();
    let outdir = &cli.outdir;
    let r1_file = &cli.in1;
    let r2_file = &cli.in2;
    let input_type = &cli.input_type.unwrap();
    let fasta_file = &cli.fasta;
    let block_info_file = &cli.block_info;
    let prefix = &cli.prefix;
    let template: &Option<String> = &cli.template;
    let ext = &cli.ext.unwrap_or(input_type.to_string());
    tracing_subscriber::fmt::init();
    info!("Start");

    rayon::ThreadPoolBuilder::new()
        .num_threads(*threads)
        .build_global()
        .unwrap();

    //minijinja
    let template_string = fs::read_to_string(template.clone().unwrap()).expect("无法读取模板文件");
    // 创建一个新的 MiniJinja 环境
    let env = Environment::new();
    // 从字符串创建一个模板
    let template: Template<'_, '_> = env
        .template_from_str(&template_string)
        .expect("无法从字符串创建模板");

    let outdir = Path::new(outdir);
    let read_info_file = outdir.join(format!("{}.{}", prefix, "read_info.stat.tsv"));
    let read_info_handle: Arc<Mutex<File>> =
        Arc::new(Mutex::new(File::create(read_info_file.clone()).unwrap()));
    
    // 统计所有 flag 的数目
    let flag_stat_hash: Arc<Mutex<HashMap<usize, usize>>> = Arc::new(Mutex::new(HashMap::new()));
    let block_info_list = get_block_info_fasta_from_file(block_info_file, fasta_file).unwrap();
    let out_fq_handle_vec: Arc<Mutex<Vec<File>>>;
    let ud_fq_handle_vec: Arc<Mutex<Vec<File>>>;
    let barcode_handle_hash: Arc<Mutex<DashMap<String, Vec<File>>>> = Default::default();
    BLOCKFLAGS.lock().unwrap().iter().for_each(|(k, v)| {
        writeln!(read_info_handle.lock().unwrap(), "#idx:flag={}:{}", k, v).unwrap();
    });
    if r2_file.is_some() {
        // output read info  statistics
        let out_fq_file_r1 = outdir.join(format!("{}.{}.{}", prefix, "template.r1", ext));
        let out_fq_file_r2 = outdir.join(format!("{}.{}.{}", prefix, "template.r2", ext));
        let ud_fq_file_r1 = outdir.join(format!("{}.{}.{}", prefix, "undetermined.r1", input_type));
        let ud_fq_file_r2 = outdir.join(format!("{}.{}.{}", prefix, "undetermined.r2", input_type));
        out_fq_handle_vec = Arc::new(Mutex::new(vec![
            File::create(out_fq_file_r1.clone()).unwrap(),
            File::create(out_fq_file_r2.clone()).unwrap(),
        ]));
        ud_fq_handle_vec = Arc::new(Mutex::new(vec![
            File::create(ud_fq_file_r1.clone()).unwrap(),
            File::create(ud_fq_file_r2.clone()).unwrap(),
        ]));
    } else {
        let out_fq_file_r1 = outdir.join(format!("{}.{}.{}", prefix, "template", ext));
        let ud_fq_file_r1 = outdir.join(format!("{}.{}.{}", prefix, "undetermined", input_type));
        out_fq_handle_vec = Arc::new(Mutex::new(vec![
            File::create(out_fq_file_r1.clone()).unwrap()
        ]));
        ud_fq_handle_vec = Arc::new(Mutex::new(vec![
            File::create(ud_fq_file_r1.clone()).unwrap()
        ]));
    }
    if input_type == "fastq" {
        let record_r1 = fastq::Reader::new(get_reader(r1_file)).records();
        if let Some(r2_file) = r2_file {
            let record_r2 = fastq::Reader::new(get_reader(r2_file)).records();
            record_r1
                .into_iter()
                .zip(record_r2.into_iter())
                .par_bridge()
                .for_each(|(record_r1, record_r2)| {
                    let mut barcode_handle_hash = barcode_handle_hash.clone();
                    process_record_pair(
                        record_r1.unwrap(),
                        record_r2.unwrap(),
                        ext,
                        &block_info_list,
                        prefix,
                        outdir,
                        barcode_handle_hash,
                        template.clone(),
                        out_fq_handle_vec.clone(),
                        ud_fq_handle_vec.clone(),
                        read_info_handle.clone(),
                        flag_stat_hash.clone(),
                    )
                })
        } else {
            record_r1.into_iter().par_bridge().for_each(|record_r1| {
                let barcode_handle_hash = barcode_handle_hash.clone();
                process_record_single(
                    record_r1.unwrap(),
                    ext,
                    &block_info_list,
                    prefix,
                    outdir,
                    barcode_handle_hash,
                    template.clone(),
                    out_fq_handle_vec.clone(),
                    ud_fq_handle_vec.clone(),
                    read_info_handle.clone(),
                    flag_stat_hash.clone(),
                )
            })
        }
    }else if input_type == "fasta" {
        let record_r1 = fasta::Reader::new(get_reader(r1_file)).records();
        if let Some(r2_file) = r2_file {
            let record_r2 = fasta::Reader::new(get_reader(r2_file)).records();
            record_r1
                .into_iter()
                .zip(record_r2.into_iter())
                .par_bridge()
                .for_each(|(record_r1, record_r2)| {
                    let barcode_handle_hash = barcode_handle_hash.clone();
                    process_record_pair(
                        record_r1.unwrap(),
                        record_r2.unwrap(),
                        ext,
                        &block_info_list,
                        prefix,
                        outdir,
                        barcode_handle_hash,
                        template.clone(),
                        out_fq_handle_vec.clone(),
                        ud_fq_handle_vec.clone(),
                        read_info_handle.clone(),
                        flag_stat_hash.clone(),
                    )
                })
        } else {
            record_r1.into_iter().par_bridge().for_each(|record_r1| {
                let barcode_handle_hash = barcode_handle_hash.clone();
                process_record_single(
                    record_r1.unwrap(),
                    ext,
                    &block_info_list,
                    prefix,
                    outdir,
                    barcode_handle_hash,
                    template.clone(),
                    out_fq_handle_vec.clone(),
                    ud_fq_handle_vec.clone(),
                    read_info_handle.clone(),
                    flag_stat_hash.clone(),
                )
            })
        }
    }

    // write read flag stat file
    let flag_stat_file = outdir.join(format!("{}.{}", prefix, "block_flag.stat.tsv"));
    let mut flag_stat_handle = File::create(flag_stat_file).unwrap();
    BLOCKFLAGS.lock().unwrap().iter().for_each(|(k, v)| {
        writeln!(flag_stat_handle, "#idx:flag={}:{}", k, v).unwrap();
    });
    flag_stat_hash.lock().unwrap().iter().for_each(|(k, v)| {
        write!(flag_stat_handle, "{}\t{}\n", k, v);
    });
    info!("End");

}