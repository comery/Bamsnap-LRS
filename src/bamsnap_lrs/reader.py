from dataclasses import dataclass
from typing import List, Optional

from .cigar import from_cigar_md_cs, from_cigar_with_ref, Segment
from .ref import get_ref_subseq


@dataclass
class Read:
    qname: str
    chrom: str
    start: int
    end: int
    reverse: bool
    mapq: int
    primary: bool
    supplementary: bool
    secondary: bool
    segments: List[Segment]
    seq: Optional[str]


def fetch_reads(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    max_reads: int = 500,
    mapq_min: int = 0,
    show_supp: bool = True,
    show_secondary: bool = False,
    downsample_strategy: str = "mapq",
    use_md: bool = False,
    use_cs: bool = True,
    use_ref: bool = False,
    fa_path: Optional[str] = None,
) -> List[Read]:
    import pysam

    reads: List[Read] = []
    with pysam.AlignmentFile(bam_path, "rb") as af:
        for r in af.fetch(chrom, start, end):
            if r.mapping_quality < mapq_min:
                continue
            supp = bool(r.flag & 0x800)
            sec = bool(r.flag & 0x100)
            prim = not supp and not sec
            if not show_supp and supp:
                continue
            if not show_secondary and sec:
                continue
            md = None
            cs = None
            try:
                md = r.get_tag("MD") if use_md else None
            except Exception:
                md = None
            try:
                cs = r.get_tag("cs") if use_cs else None
            except Exception:
                cs = None
            ref_seq = None
            if use_ref and fa_path:
                ref_seq = get_ref_subseq(fa_path, chrom, r.reference_start, r.reference_end)
            if use_ref and ref_seq:
                segs = from_cigar_with_ref(r.cigarstring, r.query_sequence, ref_seq)
            else:
                segs = from_cigar_md_cs(r.cigarstring, md, cs)
            reads.append(
                Read(
                    qname=r.query_name,
                    chrom=chrom,
                    start=r.reference_start,
                    end=r.reference_end,
                    reverse=r.is_reverse,
                    mapq=r.mapping_quality,
                    primary=prim,
                    supplementary=supp,
                    secondary=sec,
                    segments=segs,
                    seq=r.query_sequence,
                )
            )
            if len(reads) >= max_reads * 3:
                break
    if len(reads) > max_reads:
        if downsample_strategy == "mapq":
            reads.sort(key=lambda x: (not x.primary, -x.mapq, -(x.end - x.start)))
            reads = reads[:max_reads]
        else:
            reads = reads[:max_reads]
    return reads