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