# legoseq


Manipulating sequencing data is like playing with LEGO bricks.


## Features

* The DNA sequence is segmented into blocks based on the Fix (Anchor) sequence, with the flexibility of freely combining these blocks using jinja2 templates.

* Features support for the division of sample indices.

* Compatible with a variety of sequencing platforms, including Illumina, Pacbio, and Nanopore.

* Accepts inputs in fasta/fastq formats.

* Supports single end and pair end inputs; for paired-end sequencing, currently, only read1 is available for free combination.

* Facilitates parallel computing, ensuring rapid processing speeds.



## Install

### Source

If the [Rust](https://www.rust-lang.org/tools/install) compiler and associated [Cargo](https://github.com/rust-lang/cargo/) are installed, legoseq may be installed via

```
git clone git@github.com:dongspy/legoseq.git
cd legoseq
cargo build --release
target/release/legoseq --help

```

## Quick Start

### test

```
cd test

sh  run_pair_end.sh #  pair-end fastq as input

sh run_pair_end_fa.sh # pair-end fasta as input

sh run_single_end.sh # single-end fastq as input
```