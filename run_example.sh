#!/bin/bash
# 测试脚本 - 使用示例文件生成不同格式的输出

# 设置 PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:$(pwd)/src"

# 测试区域
CHROM="chrM"
START=1000
END=3000
BAM="example/sim.mapping.sort.bam"
FASTA="example/chm13v2.chrM.fasta"

echo "测试 Bamsnap-LRS 输出功能..."
echo "使用文件: $BAM"
echo "参考序列: $FASTA"
echo "区域: ${CHROM}:${START}-${END}"
echo ""


# SVG
echo "1. single pos, output SVG"
bin/bamsnap-lrs dna \
    --bam "$BAM" \
    --pos "${CHROM}:1-10000" \
    --out example/chrM_1_10000.bamsnapLRS.svg \
    --fa "$FASTA" \
	-g example/chrM.mitos2.gff \
    --show-axis \
    --show-coverage \
    --track-title "Test Reads - SVG" \
    --width 1500

if [ $? -eq 0 ]; then
    echo "✓ SVG output success: example/chrM_1000_3000.bamsnapLRS.svg"
else
    echo "✗ SVG output failed"
fi

# PDF
echo "2. single pos, output PDF"
bin/bamsnap-lrs dna \
    --bam "$BAM" \
    --pos "${CHROM}:${START}-${END}" \
    --out example/chrM_1000_3000.bamsnapLRS.pdf \
    --fa "$FASTA" \
	--bed example/chrM.annot.bed \
    --show-axis \
    --show-coverage \
    --track-title "Test Reads - PDF" \
    --width 1500 \
	--detail high

if [ $? -eq 0 ]; then
    echo "✓ PDF output success: example/chrM_1000_3000.bamsnapLRS.pdf"
else
    echo "✗ PDF output failed"
fi

region.bed

# PDF
echo "2.1, multiple pos, output PDF"
bin/bamsnap-lrs dna \
    --bam "$BAM" \
    --regions example/region.bed \
	--out-prefix example/batch \
    --fa "$FASTA" \
	--bed example/chrM.annot.bed \
    --show-axis \
    --show-coverage \
    --track-title "Test Reads - PDF" \
    --width 1500 \
	--detail high

if [ $? -eq 0 ]; then
    echo "✓ PDF output success: example/batch_*"
else
    echo "✗ PDF output failed"
fi


<<EOF
# PDF
echo "3. Multiple Bam..."
bin/bamsnap-lrs dna \
	--bam "$BAM" \
	--bam "$BAM" \
    --pos "${CHROM}:${START}-${END}" \
    --out example/chrM_1000_3000.double.bamsnapLRS.pdf \
    --fa "$FASTA" \
    --show-axis \
    --show-coverage \
    --track-title "Test Reads - PDF" \
    --width 1500

if [ $? -eq 0 ]; then
    echo "✓ PDF output success: example/chrM_1000_3000.double.bamsnapLRS.pdf"
else
    echo "✗ PDF output failed"
fi


# RNA
echo "4. RNA 比对可视化..."
bin/bamsnap-lrs rna \
	--bam example/rna/test.rna.bam \
	--pos Cdec_SDR_X:200000-300000 \
	--out example/rna/test.rna.svg \
	--fa example/rna/ref.fasta \
	--show-axis \
	--width 4000

if [ $? -eq 0 ]; then
    echo "✓ SVG output success: example/rna/testbam.look/output.svg"
else
    echo "✗ SVG output failed"
fi
EOF
