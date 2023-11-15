../target/release/legoseq \
    --block-info data/blockinfo.tsv \
    --in1 data/test_input.fastq --in2 data/test_input.fastq \
    --fasta data/blockinfo.fasta --outdir output --prefix pair_end \
    --template data/template.txt --threads 4
