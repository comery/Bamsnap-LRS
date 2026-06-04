#!/bin/bash

# ====================== auto check FASTA file ======================
FASTA_FILE="chm13v2.0.fasta"

if [ -f "$FASTA_FILE" ]; then
    echo "reference exist $FASTA_FILE, run command with --fa"
    FASTA_ARG="--fa $FASTA_FILE"
else
    echo "no reference $FASTA_FILE, run command no --fa"
    echo "reference can download from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
    FASTA_ARG=""
fi

## insertion
python  ../../bin/bamsnap-lrs dna \
  --bam insertion.Cyclone.bam --bam insertion.HIFI.bam --bam insertion.ONT.bam \
  --regions insertion.vcf  --out-prefix insertion.svg  --show-axis --show-coverage --padding 50  \
  $FASTA_ARG 
python  ../../bin/bamsnap-lrs dna \
  --bam insertion.HIFI.bam --regions insertion.vcf  --out-prefix insertion.onebam.svg  --show-axis --show-coverage --padding 50  \
  $FASTA_ARG

## deletion
python  ../../bin/bamsnap-lrs dna \
  --bam deletion.Cyclone.bam --bam deletion.HIFI.bam --bam deletion.ONT.bam \
  --regions deletion.vcf  --out-prefix deletion.svg  --show-axis --show-coverage \
  $FASTA_ARG
  
python  ../../bin/bamsnap-lrs dna \
  --bam deletion.HIFI.bam --regions deletion.vcf  --out-prefix deletion.onebam.svg  --show-axis --show-coverage \
  $FASTA_ARG

## duplication
python  ../../bin/bamsnap-lrs dna \
  --bam duplication.HIFI.bam --regions duplication.vcf --out-prefix duplication.svg --show-axis --show-coverage --show-supp \
  $FASTA_ARG

## inversion
python  ../../bin/bamsnap-lrs dna \
  --bam inversion.HIFI.bam --regions inversion.vcf --out-prefix inversion.svg --show-axis --show-coverage --show-supp \
  $FASTA_ARG

