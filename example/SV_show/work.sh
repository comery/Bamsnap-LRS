## insertion
python  /data/work/01.bamsnap_lrs/07.modify.v7/Bamsnap-LRS-main/bin/bamsnap-lrs dna --bam insertion.Cyclone.bam --bam insertion.HIFI.bam --bam insertion.ONT.bam --regions insertion.vcf  --out-prefix insertion.svg  --show-axis --show-coverage 
python  /data/work/01.bamsnap_lrs/07.modify.v7/Bamsnap-LRS-main/bin/bamsnap-lrs dna --bam insertion.HIFI.bam --regions insertion.vcf  --out-prefix insertion.onebam.svg  --show-axis --show-coverage 


## deletion
python  /data/work/01.bamsnap_lrs/07.modify.v7/Bamsnap-LRS-main/bin/bamsnap-lrs dna --bam deletion.Cyclone.bam --bam deletion.HIFI.bam --bam deletion.ONT.bam --regions deletion.vcf  --out-prefix deletion.svg  --show-axis --show-coverage 
python  /data/work/01.bamsnap_lrs/07.modify.v7/Bamsnap-LRS-main/bin/bamsnap-lrs dna --bam deletion.HIFI.bam --regions deletion.vcf  --out-prefix deletion.onebam.svg  --show-axis --show-coverage 


## duplication
python  /data/work/01.bamsnap_lrs/07.modify.v7/Bamsnap-LRS-main/bin/bamsnap-lrs dna --bam duplication.HIFI.bam --regions duplication.vcf --out-prefix duplication.svg --show-axis --show-coverage --show-supp


## inversion
python  /data/work/01.bamsnap_lrs/07.modify.v7/Bamsnap-LRS-main/bin/bamsnap-lrs dna --bam inversion.HIFI.bam --regions inversion.vcf --out-prefix inversion.svg --show-axis --show-coverage --show-supp

