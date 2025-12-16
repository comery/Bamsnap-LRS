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
echo "1. 生成 SVG 格式..."
python3 -m bamsnap_lrs snap \
    --bam "$BAM" \
    --pos "${CHROM}:${START}-${END}" \
    --out test_output.svg \
    --use-fa \
    --fa "$FASTA" \
    --show-axis \
    --show-coverage \
    --track-title "Test Reads - SVG" \
    --width 1200

if [ $? -eq 0 ]; then
    echo "✓ SVG output success: test_output.svg"
else
    echo "✗ SVG output failed"
fi

# PDF
echo "2. 生成 PDF 格式..."
python3 -m bamsnap_lrs snap \
    --bam "$BAM" \
    --pos "${CHROM}:${START}-${END}" \
    --out test_output.pdf \
    --use-fa \
    --fa "$FASTA" \
    --show-axis \
    --show-coverage \
    --track-title "Test Reads - PDF" \
    --width 1200

if [ $? -eq 0 ]; then
    echo "✓ PDF output success: test_output.pdf"
else
    echo "✗ PDF output failed"
fi

