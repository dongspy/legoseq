../target/release/legoseq \
    --block-info data/blockinfo.tsv \
    --in1 data/test_input.fastq --fasta data/blockinfo.fasta \
    --threads 4 --outdir output --prefix single_end \
    --template data/template.txt
