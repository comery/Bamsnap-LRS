#!/usr/bin/env python3
"""
测试脚本：解析 CIGAR 字符串

用法：
    python tests/test_parse_cigar.py --cigar "100M1D50M" [--md MD_TAG] [--cs CS_TAG]
"""

import argparse
import json
import sys
import os

# 添加 src 目录到 Python 路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from bamsnap_lrs.cigar import from_cigar_md_cs


def main():
    p = argparse.ArgumentParser(description="解析 CIGAR 字符串并输出 JSON 格式的 segments")
    p.add_argument("--cigar", required=True, help="CIGAR 字符串，如 '100M1D50M'")
    p.add_argument("--md", help="MD 标签（可选）")
    p.add_argument("--cs", help="cs 标签（可选，minimap2 差异字符串）")
    args = p.parse_args()
    
    segs = from_cigar_md_cs(args.cigar, args.md, args.cs)
    print(json.dumps([s.__dict__ for s in segs], indent=2))


if __name__ == "__main__":
    main()

