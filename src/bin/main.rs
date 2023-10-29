// #![allow(dead_code)]
// #![allow(unused)]
use bio::io::fastq;
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

use legoseq::blockinfo::{get_block_info_fasta_from_file, BLOCKFLAGS};
use legoseq::readblockalign::ReadBlockAlign;
use legoseq::utils::get_reader;


#[derive(Parser)]
#[command(version, author, about, long_about = None)]
struct Cli {
    /// fastq file
    #[arg(long, value_name = "FILE")]
    fq1: String,
    /// fastq file, optional
    #[arg(long, value_name = "FILE")]
    fq2: Option<String>,
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
}

fn main() {
    let cli = Cli::parse();
    let threads = &cli.threads.to_owned();
    let outdir = &cli.outdir;
    let r1_file = &cli.fq1;
    let r2_file = &cli.fq2;
    let fasta_file = &cli.fasta;
    let block_info_file = &cli.block_info;
    // let block_info_file2 = &cli.block_info2;
    let prefix = &cli.prefix;
    // let export_blocks = &CLI.export_blocks;
    let template: &Option<String> = &cli.template;
    // let template2: &Option<String> = &cli.template2;
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
    let read_info_handle = Arc::new(Mutex::new(File::create(read_info_file.clone()).unwrap()));
    BLOCKFLAGS.lock().unwrap().iter().for_each(|(k, v)| {
        writeln!(read_info_handle.lock().unwrap(), "#idx:flag={}:{}", k, v).unwrap();
    });
    // 统计所有 flag 的数目
    let flag_stat_hash: Arc<Mutex<HashMap<usize, usize>>> = Arc::new(Mutex::new(HashMap::new()));

    let record_r1 = fastq::Reader::new(get_reader(r1_file)).records();
    if let Some(r2_file) = r2_file {
        // output read info  statistics

        let out_fq_file_r1 = outdir.join(format!("{}.{}", prefix, "template.r1.fastq"));
        let out_fq_file_r2 = outdir.join(format!("{}.{}", prefix, "template.r2.fastq"));
        let ud_fq_file_r1 = outdir.join(format!("{}.{}", prefix, "undetermined.r1.fastq"));
        let ud_fq_file_r2 = outdir.join(format!("{}.{}", prefix, "undetermined.r2.fastq"));

        let block_info_list = get_block_info_fasta_from_file(block_info_file, fasta_file).unwrap();
        let out_fq_handle_vec = Arc::new(Mutex::new(vec![
            File::create(out_fq_file_r1.clone()).unwrap(),
            File::create(out_fq_file_r2.clone()).unwrap(),
        ]));
        let ud_fq_handle_vec = Arc::new(Mutex::new(vec![
            File::create(ud_fq_file_r1.clone()).unwrap(),
            File::create(ud_fq_file_r2.clone()).unwrap(),
        ]));
        

        let barcode_handle_hash: DashMap<String, Arc<Mutex<Vec<File>>>> = DashMap::new();

        let record_r2 = fastq::Reader::new(get_reader(r2_file)).records();
        record_r1
            .into_iter()
            .zip(record_r2.into_iter())
            .par_bridge()
            .for_each(|(record_r1, record_r2)| {
                let record_r1 = record_r1.unwrap();
                let record_r2 = record_r2.unwrap();
                let read_name = record_r1.id();
                let read_block_align =
                    ReadBlockAlign::read_block_info(&record_r1, &block_info_list);
                let flag = read_block_align.get_block_flag();
                let best_index_vec = read_block_align.get_best_index();

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
            })
    } else {
        let out_fq_file = outdir.join(format!("{}.{}", prefix, "template.fastq"));
        let ud_fq_file = outdir.join(format!("{}.{}", prefix, "undetermined.fastq"));

        let block_info_list = get_block_info_fasta_from_file(block_info_file, fasta_file).unwrap();
        let read_info_handle = Arc::new(Mutex::new(File::create(read_info_file.clone()).unwrap()));
        let out_fq_handle: Arc<Mutex<File>> =
            Arc::new(Mutex::new(File::create(out_fq_file.clone()).unwrap()));
        let ud_fq_handle: Arc<Mutex<File>> =
            Arc::new(Mutex::new(File::create(ud_fq_file.clone()).unwrap()));
        let barcode_handle_hash: DashMap<String, Arc<Mutex<Vec<File>>>> = DashMap::new();

        BLOCKFLAGS.lock().unwrap().iter().for_each(|(k, v)| {
            writeln!(read_info_handle.lock().unwrap(), "#idx:flag={}:{}", k, v).unwrap();
        });

        record_r1.into_iter().par_bridge().for_each(|record_r1| {
            let record_r1 = record_r1.unwrap();
            let read_name = record_r1.id();
            let read_block_align = ReadBlockAlign::read_block_info(&record_r1, &block_info_list);
            let flag = read_block_align.get_block_flag();
            let best_index_vec = read_block_align.get_best_index();

            // if barcode index existed, demultiplex
            if !best_index_vec.is_empty() {
                let best_index_str = best_index_vec.join("_");
                let barcode_file =
                    outdir.join(format!("{}.{}.{}", prefix, best_index_str, "fastq"));
                let barcode_handle = barcode_handle_hash
                    .entry(barcode_file.to_str().unwrap().to_string())
                    .or_insert_with(|| {
                        Arc::new(Mutex::new(vec![File::create(barcode_file).unwrap()]))
                    });

                let template_str = read_block_align.template_str(&template);
                if let Some(template_str) = template_str {
                    writeln!(
                        barcode_handle.lock().unwrap().first().unwrap(),
                        "{}",
                        template_str
                    )
                    .unwrap();
                }else{
                    // let out_fq_handle = out_fq_handle.lock().unwrap();
                        writeln!(
                            ud_fq_handle.lock().unwrap(),
                            "@{} {}\n+\n{}\n{}",
                            record_r1.id(),
                            record_r1.desc().unwrap_or(""),
                            String::from_utf8_lossy(record_r1.seq()),
                            String::from_utf8_lossy(record_r1.qual())
                        )
                        .unwrap();
                }
            } else {
                // export to file based on the jinja template
                let template_str = read_block_align.template_str(&template);
                if let Some(template_str) = template_str {
                    writeln!(out_fq_handle.lock().unwrap(), "{}", template_str).unwrap();
                }else{
                    // let out_fq_handle = out_fq_handle.lock().unwrap();
                        writeln!(
                            ud_fq_handle.lock().unwrap(),
                            "@{} {}\n+\n{}\n{}",
                            record_r1.id(),
                            record_r1.desc().unwrap_or(""),
                            String::from_utf8_lossy(record_r1.seq()),
                            String::from_utf8_lossy(record_r1.qual())
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
        })
    }

    let flag_stat_file = outdir.join(format!("{}.{}", prefix, "block_flag.stat.tsv"));
    info!(
        "write read information to file: {}",
        read_info_file.display()
    );
    info!(
        "write block flag stat to file: {}",
        flag_stat_file.display()
    );
    // info!(
    //     "write template information to file: {}",
    //     out_fq_file.display()
    // );
    let mut flag_stat_handle = File::create(flag_stat_file).unwrap();

    BLOCKFLAGS.lock().unwrap().iter().for_each(|(k, v)| {
        writeln!(flag_stat_handle, "#idx:flag={}:{}", k, v).unwrap();
    });
    flag_stat_hash.lock().unwrap().iter().for_each(|(k, v)| {
        writeln!(flag_stat_handle, "{}\t{}", k, v).unwrap();
    });

    info!("End");
}
