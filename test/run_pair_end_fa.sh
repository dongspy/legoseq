../target/debug/legoseq \
    --input-type fasta \
    --block-info data/blockinfo.tsv \
    --in1 data/test_input.fasta --in2 data/test_input.fasta \
    --fasta data/test.fasta --outdir output --prefix pair_end_fa \
    --template data/template.txt --threads 4 
