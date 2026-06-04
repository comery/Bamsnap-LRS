#!/bin/bash
# ====================== auto check FASTA file ======================
FASTA_FILE="hg38.fasta"

if [ -f "$FASTA_FILE" ]; then
    echo "reference exist $FASTA_FILE, run command with --fa"
    FASTA_ARG="--fa $FASTA_FILE"
else
    echo "no reference $FASTA_FILE, run command no --fa"
    echo "reference can download from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/"
    FASTA_ARG=""
fi

python  ../../bin/bamsnap-lrs highlight  --highlight-vcf region_chr21.vcf \
  --bam sample1.hignlignt.bam --bam sample2.hignlignt.bam \
  --pos chr21:20462000-20467000 --out highlight.svg --max-reads 20 $FASTA_ARG
