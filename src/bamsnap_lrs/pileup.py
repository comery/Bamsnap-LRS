from typing import Dict, List

from .reader import Read


def base_pileup(reads: List[Read], start: int, end: int) -> List[Dict[str, int]]:
    n = end - start
    out: List[Dict[str, int]] = [{"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "depth": 0} for _ in range(n)]
    for r in reads:
        if not r.seq:
            continue
        rc = 0
        qc = 0
        for s in r.segments:
            if s.ref_consumed == 0 and s.read_consumed == 0:
                continue
            if s.ref_consumed > 0 and s.read_consumed > 0:
                for i in range(s.ref_consumed):
                    pos = r.start + rc + i
                    if pos < start or pos >= end:
                        continue
                    bi = qc + i
                    if bi < 0 or bi >= len(r.seq):
                        continue
                    b = r.seq[bi].upper()
                    if b not in ("A", "C", "G", "T"):
                        b = "N"
                    idx = pos - start
                    out[idx][b] += 1
                    out[idx]["depth"] += 1
                rc += s.ref_consumed
                qc += s.read_consumed
            elif s.ref_consumed > 0 and s.read_consumed == 0:
                rc += s.ref_consumed
            elif s.read_consumed > 0 and s.ref_consumed == 0:
                qc += s.read_consumed
    return out


def pileup_to_pixels(pile: List[Dict[str, int]], width: int) -> List[Dict[str, int]]:
    """Aggregate base positions to pixel level, taking average values"""
    n = len(pile)
    out: List[Dict[str, int]] = []
    for x in range(width):
        s = int(x * n / width)
        e = int((x + 1) * n / width)
        if e <= s:
            e = s + 1
        # Calculate average base distribution within this pixel range
        agg = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "depth": 0}
        count = 0
        for i in range(s, min(e, n)):
            for k in agg.keys():
                agg[k] += pile[i][k]
            count += 1
        # Take average (rounded)
        if count > 1:
            for k in agg.keys():
                agg[k] = round(agg[k] / count)
        out.append(agg)
    return out


def pileup_to_pixels_with_ref(pile: List[Dict[str, int]], width: int, ref_seq: str = None) -> List[Dict[str, int]]:
    """Aggregate base positions to pixel level, calculating reference match and variant portions
    
    Note: To preserve variant information, variant bases use ceil (round up) to ensure
    even small variants are displayed
    """
    import math
    n = len(pile)
    out: List[Dict[str, int]] = []
    for x in range(width):
        s = int(x * n / width)
        e = int((x + 1) * n / width)
        if e <= s:
            e = s + 1
        # Calculate base distribution within this pixel range
        # Separately count: bases matching reference (ref_match), bases not matching reference
        ref_match = 0
        variant_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
        total_depth = 0
        
        for i in range(s, min(e, n)):
            # Get reference base at this position
            ref_base = None
            if ref_seq and i < len(ref_seq):
                ref_base = ref_seq[i].upper()
            
            for base in ["A", "C", "G", "T", "N"]:
                count = pile[i].get(base, 0)
                if count > 0:
                    if ref_base and base == ref_base:
                        ref_match += count
                    else:
                        variant_counts[base] += count
            total_depth += pile[i].get("depth", 0)
        
        # Calculate average
        num_positions = e - s
        if num_positions > 1:
            # Reference bases use rounding
            ref_match = round(ref_match / num_positions)
            # Variant bases use ceil (round up) to ensure variants are not lost due to averaging
            for base in variant_counts:
                if variant_counts[base] > 0:
                    variant_counts[base] = max(1, math.ceil(variant_counts[base] / num_positions))
            total_depth = round(total_depth / num_positions)
        
        agg = {
            "ref_match": ref_match,
            "A": variant_counts["A"],
            "C": variant_counts["C"],
            "G": variant_counts["G"],
            "T": variant_counts["T"],
            "N": variant_counts["N"],
            "depth": total_depth
        }
        out.append(agg)
    return out
