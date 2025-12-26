~/Desktop/Work/数字化地球所/BioEarthDigital/CycSim/target/release/cycsim sim -d 30 -M 10000 -m 3000 ~/Desktop/Work/数字化地球所/BioEarthDigital/CycSim/model/hifi_model.v1.1.cy chm13v2.chrM.modi.fasta > chm13v2.chrM.modi.sim.bam
samtools sort -o chm13v2.chrM.modi.sim.sort.bam chm13v2.chrM.modi.sim.bam
samtools fastq -F 0x900 chm13v2.chrM.modi.sim.sort.bam  > all_reads.fq
minimap2 -ax map-hifi  -t 4  -c chm13v2.chrM.fasta all_reads.fq|samtools view -bS - > sim.mapping.bam
samtools sort -o sim.mapping.sort.bam sim.mapping.bam
samtools index sim.mapping.sort.bam

# sim NGS and mapping
~/Desktop/Work/Software/biosoft/wgsim/wgsim -N 100000 -1 100 -2 100 chm13v2.chrM.modi.fasta chm13v2.chrM.ngs.read1.fq chm13v2.chrM.ngs.read2.fq
bwa index chm13v2.chrM.fasta
bwa mem chm13v2.chrM.fasta chm13v2.chrM.ngs.read1.fq chm13v2.chrM.ngs.read2.fq |samtools view -bS - > chm13v2.chrM.ngs.mapping.bam
samtools sort -o chm13v2.chrM.ngs.mapping.sort.bam chm13v2.chrM.ngs.mapping.bam
samtools index chm13v2.chrM.ngs.mapping.sort.bam
