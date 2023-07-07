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

use crate::aligner::Alignment;
use crate::blockalign::ReadBlockAlign;
use crate::blockalign::{block_align_read, BlockAlign, BlockAlignAbbr};
use crate::blockinfo::{get_block_info, BlockInfo, BLOCKFLAGS};
use crate::utils::get_reader;
