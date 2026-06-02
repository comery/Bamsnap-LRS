#从BAM/CRAM文件中读取指定基因组区间的测序比对数据，并将其解析为标准化的Read对象列表，供后续可视化和分析使用。
#主要功能包括：
#Read数据结构：封装了每条read的基本信息（名称、区间、链向、MAPQ、主/辅/次比对、片段列表、序列等）。
#fetch_reads：读取DNA测序数据，支持多种过滤（MAPQ、主/辅/次比对）、下采样、MD/CS/参考序列解析等，输出标准化Read对象。
#fetch_rna_reads：专为RNA测序设计，能自动收集同一转录本的所有拼接片段（包括supplementary），便于展示剪接结构。
#简言之：比对数据读取与标准化模块，把原始BAM/CRAM比对结果转为统一的结构化对象，方便后续可视化和统计。


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
    # Alignment span on the original read, 0-based half-open.
    # For reverse-strand alignments, pysam query_alignment_start/end are
    # converted back to the original read orientation.
    query_start: int = 0
    query_end: int = 0
    query_length: int = 0

def _query_span_on_original_read(r):
    """
    Return alignment span on the original read as 0-based half-open coordinates.

    pysam query_alignment_start/query_alignment_end follow the stored SEQ/CIGAR
    orientation. For reverse-strand alignments, convert them back to the
    original read direction.

    Example:
        original read length = 11545
        reverse alignment CIGAR: 4578S6967M
        stored query span = 4578-11545
        original read span = 0-6967
    """
    q_len = r.infer_query_length(always=True) or r.query_length or len(r.query_sequence or "")

    qa_start = r.query_alignment_start
    qa_end = r.query_alignment_end

    if qa_start is None:
        qa_start = 0
    if qa_end is None:
        qa_end = q_len

    if r.is_reverse:
        query_start = q_len - qa_end
        query_end = q_len - qa_start
    else:
        query_start = qa_start
        query_end = qa_end

    return query_start, query_end, q_len
    
    
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
    import os

    reads: List[Read] = []
    
    # Detect file format and prepare open parameters
    is_cram = bam_path.lower().endswith('.cram')
    open_mode = "rc" if is_cram else "rb"
    open_kwargs = {}
    
    # For CRAM files, provide reference genome if available
    if is_cram and fa_path and os.path.exists(fa_path):
        open_kwargs['reference_filename'] = fa_path
    
    with pysam.AlignmentFile(bam_path, open_mode, **open_kwargs) as af:
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
                
            q_start, q_end, q_len = _query_span_on_original_read(r)
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
                    query_start=q_start,
                    query_end=q_end,
                    query_length=q_len,
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


def fetch_rna_reads(
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
    """Fetch RNA reads, handling spliced alignments.
    
    For RNA-seq data, a single transcript may be split into multiple alignment segments
    (exons) due to alternative splicing. These segments appear as separate records in BAM/CRAM,
    typically marked as supplementary alignments. This function collects all segments
    belonging to the same transcript (same qname) and returns them as separate Read objects.
    The renderer will connect these segments with lines to show the splicing structure.
    """
    import pysam
    import os
    
    reads: List[Read] = []
    # Track transcripts by qname to ensure we get all segments
    transcript_segments: dict[str, List[Read]] = {}
    
    # Detect file format and prepare open parameters
    is_cram = bam_path.lower().endswith('.cram')
    open_mode = "rc" if is_cram else "rb"
    open_kwargs = {}
    
    # For CRAM files, provide reference genome if available
    if is_cram and fa_path and os.path.exists(fa_path):
        open_kwargs['reference_filename'] = fa_path
    
    with pysam.AlignmentFile(bam_path, open_mode, **open_kwargs) as af:
        for r in af.fetch(chrom, start, end):
            if r.mapping_quality < mapq_min:
                continue
            supp = bool(r.flag & 0x800)
            sec = bool(r.flag & 0x100)
            prim = not supp and not sec
            
            # For RNA, we want to include supplementary alignments (spliced segments)
            # but filter secondary if requested
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
            
            read = Read(
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
            
            # Group by transcript name (qname)
            if read.qname not in transcript_segments:
                transcript_segments[read.qname] = []
            transcript_segments[read.qname].append(read)
    
    # Flatten and sort segments by transcript, then by position
    all_reads: List[Read] = []
    for qname, segments in transcript_segments.items():
        # Sort segments by start position (for forward) or end position (for reverse)
        if segments[0].reverse:
            segments.sort(key=lambda x: -x.end)  # Reverse strand: right to left
        else:
            segments.sort(key=lambda x: x.start)  # Forward strand: left to right
        all_reads.extend(segments)
    
    # Downsample if needed
    if len(all_reads) > max_reads:
        if downsample_strategy == "mapq":
            # Sort by primary status, then mapq, then length
            all_reads.sort(key=lambda x: (not x.primary, -x.mapq, -(x.end - x.start)))
            all_reads = all_reads[:max_reads]
        else:
            all_reads = all_reads[:max_reads]
    
    return all_reads