# 测试说明

## 安装依赖

在运行测试之前，需要安装必要的 Python 依赖：

```bash
# 使用虚拟环境（推荐）
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# 或者使用 --user 标志（如果系统允许）
pip3 install --user -r requirements.txt
```

## 运行测试

### 方法 1: 使用测试脚本

```bash
./test_example.sh
```

### 方法 2: 手动运行

```bash
# 设置 PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:$(pwd)/src"

# 生成 PNG 格式
python3 -m bamsnap_lrs dna \
    --bam example/sim.mapping.sort.bam \
    --pos chrM:1000-2000 \
    --out test_output.png \
    --fa example/chm13v2.chrM.fasta \
    --show-axis \
    --show-coverage \
    --track-title "Test Reads" \
    --width 1200

# 生成 SVG 格式
python3 -m bamsnap_lrs dna \
    --bam example/sim.mapping.sort.bam \
    --pos chrM:1000-2000 \
    --out test_output.svg \
    --fa example/chm13v2.chrM.fasta \
    --show-axis \
    --show-coverage \
    --track-title "Test Reads" \
    --width 1200

# 生成 PDF 格式
python3 -m bamsnap_lrs dna \
    --bam example/sim.mapping.sort.bam \
    --pos chrM:1000-2000 \
    --out test_output.pdf \
    --fa example/chm13v2.chrM.fasta \
    --show-axis \
    --show-coverage \
    --track-title "Test Reads" \
    --width 1200
```

## 测试文件

- BAM 文件: `example/sim.mapping.sort.bam`
- 参考序列: `example/chm13v2.chrM.fasta`
- 测试区域: `chrM:1000-2000`

## 输出格式

程序会根据输出文件扩展名自动选择格式：
- `.png` → PNG 格式（位图）
- `.svg` → SVG 格式（矢量图）
- `.pdf` → PDF 格式（文档）

## 功能特性

测试将验证以下功能：
1. ✅ Coverage 堆积图显示
2. ✅ 根据参考序列显示变异（灰色=无变异，彩色=有变异）
3. ✅ 插入/缺失标注
4. ✅ 多格式输出（PNG/SVG/PDF）
5. ✅ Track 系统（标题栏、分隔线）

