# random_add_ps.py: add PS information to phased variant file

# HG002_CHM13v2.0_v5.0q_smvar.het.randomPS.vcf.gz.gz
The raw small variant set downloaded from https://42basepairs.com/download/web/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_CHM13v2.0_v5.0q_smvar.vcf.gz
bcftools view -m2 -M2 -v snps -g het -Oz -o HG002_CHM13v2.0_v5.0q_smvar.het.vcf.gz HG002_CHM13v2.0_v5.0q_smvar.vcf.gz
python random_add_ps.py -i HG002_CHM13v2.0_v5.0q_smvar.het.vcf.gz -o HG002_CHM13v2.0_v5.0q_smvar.het.randomPS.vcf \
    --min-snps 3 --max-snps 8 --seed 208
bgzip HG002_CHM13v2.0_v5.0q_smvar.het.randomPS.vcf
