#这个脚本（ref.py）的作用是：从参考基因组FASTA文件中提取指定染色体和区间的序列片段。
#主要功能包括：
#get_ref_subseq：优先用pyfaidx库高效提取指定区间的序列，如果pyfaidx不可用则自动切换到_naive_get_subseq的纯Python实现。
#_naive_get_subseq：手动遍历FASTA文件，找到目标染色体后拼接序列并返回指定区间。
#简言之：这是Bamsnap-LRS的“参考序列提取”模块，支持高效或兼容地获取基因组任意区间的碱基序列。

def get_ref_subseq(fa_path: str, chrom: str, start: int, end: int):
    try:
        from pyfaidx import Fasta
        fa = Fasta(fa_path, as_raw=True, sequence_always_upper=True)
        return str(fa[chrom][start:end])
    except Exception:
        return _naive_get_subseq(fa_path, chrom, start, end)


def _naive_get_subseq(fa_path: str, chrom: str, start: int, end: int):
    seq = []
    cur = None
    s = start
    e = end
    with open(fa_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                cur = name
                continue
            if cur != chrom:
                continue
            seq.append(line.strip())
            if len("".join(seq)) >= e:
                break
    if not seq:
        return None
    joined = "".join(seq).upper()
    return joined[s:e]