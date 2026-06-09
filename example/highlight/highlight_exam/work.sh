samtools view -bh example_reads.sam > example_reads.bam
samtools index example_reads.bam

python  ../../../bin/bamsnap-lrs highlight  --highlight-vcf example.vcf --bam example_reads.bam --pos scaffold1:0-500 --out highlight.one.bam.svg --fa example_ref.fa --width 800 --max-reads 10
