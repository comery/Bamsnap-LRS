################################################################# hawkeye
#################################################################
## insertion
(time hawkeye.py sv_browse -i insertion.Cyclone.bam,insertion.HIFI.bam,insertion.ONT.bam  -b insertion.vcf  -f vcf -t 1 -r chm13v2.0.fasta  -g chm13 -o ins.example ) 2> SVhewkeye.ins.runtime.log
## deletion
(time hawkeye.py sv_browse -i deletion.Cyclone.bam,deletion.HIFI.bam,deletion.ONT.bam -b deletion.vcf  -f vcf -t 1 -r chm13v2.0.fasta  -g chm13 -o del.example ) 2> SVhewkeye.del.runtime.log
## duplication
(time hawkeye.py sv_browse -i duplication.HIFI.bam -b duplication.vcf -f vcf -t 1 -r chm13v2.0.fasta  -g chm13 -o dup.example ) 2> SVhewkeye.dup.runtime.log
## inversion
(time hawkeye.py sv_browse -i inversion.HIFI.bam -b inversion.vcf -f vcf -t 1 -r chm13v2.0.fasta  -g chm13 -o inv.example ) 2> SVhewkeye.inv.runtime.log


################################################################# bamsnap_lrs
#################################################################
time ( python ../../bin/bamsnap-lrs dna --bam insertion.Cyclone.bam --bam insertion.HIFI.bam --bam insertion.ONT.bam  --regions insertion.vcf  --out-prefix insertion.png --fa chm13v2.0.fasta  --show-axis --show-coverage --width 1200  ) 2>bamsnap_ins.runtime.log

time ( python ../../bin/bamsnap-lrs dna --bam deletion.Cyclone.bam --bam deletion.HIFI.bam --bam deletion.ONT.bam  --regions deletion.vcf  --out-prefix deletion.png  --show-axis --show-coverage --width 1200 ) 2> bamsnap_del.runtime.log

time ( python ../../bin/bamsnap-lrs dna --bam duplication.HIFI.bam --regions duplication.vcf  --out-prefix duplication.png --fa chm13v2.0.fasta  --show-axis --show-coverage --width 1200 --show-supp  ) 2>bamsnap_dup.runtime.log

time ( python ../../bin/bamsnap-lrs dna --bam inversion.HIFI.bam --regions inversion.vcf  --out-prefix inversion.png --fa chm13v2.0.fasta  --show-axis --show-coverage --width 1200 --show-supp  ) 2>bamsnap_inv.runtime.log

