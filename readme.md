# legoseq


Manipulating sequencing data is like playing with LEGO bricks.


## Features

* The DNA sequence is segmented into blocks based on the Fix (Anchor) sequence, with the flexibility of freely combining these blocks using [jinja2](https://github.com/mitsuhiko/minijinja/blob/main/COMPATIBILITY.md#expressions) templates.

* Supports for the demultiplex of sequence based on the index(barcode).

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

### blockinfo

The software supports block types of Fix, Variable, Index.

* Fix refers to a given sequence of bases, which needs to be provided in advance, and legoseq uses an alignment algorithm to determine the position of the Fix sequence in the read;

* Variable refers to the variable sequence, which needs to be determined according to the upstream and downstream Fix sequences;

* Index refers to the barcode sequence of a sample, which can be used for demultiplexing;

For details, see test/data/blockinfo.tsv
 

### template

Generate the corresponding sequence file according to the jinja2 template.

```
>{{read.name}} {{Fix_2.seq[:6]}} len={{read.seq|length}}
ATCG{{read.seq}}AAA{{Variable_3.seq}}
```

### test

```
cd test

sh  run_pair_end.sh #  pair-end fastq as input

sh run_pair_end_fa.sh # pair-end fasta as input

sh run_single_end.sh # single-end fastq as input
```