from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass
class Segment:
    type: str
    op: str
    length: int
    ref_consumed: int
    read_consumed: int
    read_seq: Optional[str] = None

# Parse CIGAR string. Input: 100M5I10D, Output: [(M, 100), (I, 5), (D, 10)]
def parse_cigar_string(cigar: str) -> List[Tuple[str, int]]:
    ops = []
    num = []
    for ch in cigar:
        if ch.isdigit():
            num.append(ch)
        else:
            if not num:
                raise ValueError("invalid cigar")
            ops.append((ch, int("".join(num))))
            num = []
    if num:
        raise ValueError("invalid cigar")
    return ops

# Merge adjacent segments with the same type and op to reduce redundancy. e.g. [(match, =, 10), (match, =, 5)] -> [(match, =, 15)]
def merge_segments(segments: List[Segment]) -> List[Segment]:
    if not segments:
        return segments
    out = [segments[0]]
    for s in segments[1:]:
        last = out[-1]
        if s.type == last.type and s.op == last.op:
            out[-1] = Segment(
                type=last.type,
                op=last.op,
                length=last.length + s.length,
                ref_consumed=last.ref_consumed + s.ref_consumed,
                read_consumed=last.read_consumed + s.read_consumed,
            )
        else:
            out.append(s)
    return out

# Parse MD tag (alignment info from SAM/BAM describing mismatches and deletions).
# Returns a list of (op, length) tuples.
# Input: 10A5^AC10, Output: [("=", 10), ("X", 1), ("=", 5), ("D", 2), ("=", 10)]
def parse_md(md: str) -> List[Tuple[str, int]]:
    res = []
    i = 0
    n = len(md)
    while i < n:
        if md[i].isdigit():
            j = i
            while j < n and md[j].isdigit():
                j += 1
            res.append(("=", int(md[i:j])))
            i = j
        elif md[i] == "^":
            j = i + 1
            while j < n and md[j].isalpha():
                j += 1
            res.append(("D", j - (i + 1)))
            i = j
        else:
            res.append(("X", 1))
            i += 1
    return res

# Parse CS tag (detailed alignment from tools like minimap2), supporting match,
# mismatch, insertion, deletion, splicing, etc.
# Input: 100*ag+ct-gg, Output: [("=", 100, None), ("X", 1, "g"), ("I", 2, "ct"), ("D", 2, None)]
def parse_cs(cs: str) -> List[Tuple[str, int, Optional[str]]]:
    res = []
    i = 0
    n = len(cs)
    while i < n:
        c = cs[i]
        if c == ":":
            j = i + 1
            while j < n and cs[j].isdigit():
                j += 1
            res.append(("=", int(cs[i + 1:j]), None))
            i = j
        elif c == "*":
            res.append(("X", 1, cs[i+2])) # *ag (ref a, read g)
            i += 3
        elif c == "+":
            j = i + 1
            while j < n and cs[j].isalpha():
                j += 1
            ins_seq = cs[i+1:j]
            res.append(("I", j - (i + 1), ins_seq))
            i = j
        elif c == "-":
            j = i + 1
            while j < n and cs[j].isalpha():
                j += 1
            res.append(("D", j - (i + 1), None))
            i = j
        elif c == "~":
            # Intron: ~[a-z][a-z][0-9]+[a-z][a-z]
            # Skip 2 chars (donor)
            j = i + 3
            # Read length
            k = j
            while k < n and cs[k].isdigit():
                k += 1
            length = int(cs[j:k])
            # Skip 2 chars (acceptor)
            i = k + 2
            res.append(("N", length, None))
        else:
            i += 1
    return res

# Combine CIGAR, MD, and CS info to produce a Segment list. Prefer CS tag,
# then CIGAR+MD, otherwise fall back to CIGAR alone.
def from_cigar_md_cs(cigar: str, md: Optional[str] = None, cs: Optional[str] = None) -> List[Segment]:
    ops = parse_cigar_string(cigar)
    out: List[Segment] = []

    # CS path: use cs for detailed alignment states, but preserve CIGAR-only
    # operations such as soft/hard clipping and padding.
    if cs:
        detail = parse_cs(cs)
        d_idx = 0
        d_off = 0

        def take_cs(allowed_ops, need):
            nonlocal d_idx, d_off
            while need > 0:
                if d_idx >= len(detail):
                    raise ValueError("cs tag is shorter than CIGAR")

                d_op, d_len, d_seq = detail[d_idx]
                remain = d_len - d_off
                if remain <= 0:
                    d_idx += 1
                    d_off = 0
                    continue

                if d_op not in allowed_ops:
                    raise ValueError(f"cs op {d_op} is inconsistent with CIGAR")

                take = min(need, remain)

                if d_op == "=":
                    out.append(Segment("match", "=", take, take, take))
                elif d_op == "X":
                    out.append(Segment("mismatch", "X", take, take, take))
                elif d_op == "I":
                    seq = d_seq[d_off:d_off + take] if d_seq is not None else None
                    out.append(Segment("ins", "I", take, 0, take, read_seq=seq))
                elif d_op == "D":
                    out.append(Segment("del", "D", take, take, 0))
                elif d_op == "N":
                    out.append(Segment("ref_skip", "N", take, take, 0))

                need -= take
                d_off += take
                if d_off == d_len:
                    d_idx += 1
                    d_off = 0

        for op, l in ops:
            if op == "M":
                take_cs({"=", "X"}, l)
            elif op == "=":
                take_cs({"="}, l)
            elif op == "X":
                take_cs({"X"}, l)
            elif op == "I":
                take_cs({"I"}, l)
            elif op == "D":
                take_cs({"D"}, l)
            elif op == "N":
                take_cs({"N"}, l)
            elif op == "S":
                out.append(Segment("soft", "S", l, 0, l))
            elif op == "H":
                out.append(Segment("hard", "H", l, 0, 0))
            elif op == "P":
                out.append(Segment("pad", "P", l, 0, 0))
            else:
                raise ValueError("unsupported op")

        return merge_segments(out)

    # MD path: keep one MD cursor across the entire CIGAR so that MD state is
    # not restarted for each M block.
    detail = parse_md(md) if md else None
    md_idx = 0
    md_off = 0

    def skip_empty_md():
        nonlocal md_idx, md_off
        while detail and md_idx < len(detail) and detail[md_idx][1] - md_off == 0:
            md_idx += 1
            md_off = 0

    def take_md_aligned(need, emit=True):
        nonlocal md_idx, md_off
        while need > 0:
            skip_empty_md()
            if not detail or md_idx >= len(detail):
                raise ValueError("MD tag is shorter than CIGAR")

            d_op, d_len = detail[md_idx]
            if d_op == "D":
                raise ValueError("unexpected MD deletion inside aligned CIGAR block")

            remain = d_len - md_off
            take = min(need, remain)

            if emit:
                if d_op == "=":
                    out.append(Segment("match", "=", take, take, take))
                elif d_op == "X":
                    out.append(Segment("mismatch", "X", take, take, take))

            need -= take
            md_off += take
            if md_off == d_len:
                md_idx += 1
                md_off = 0

    def take_md_deletion(need):
        nonlocal md_idx, md_off
        while need > 0:
            skip_empty_md()
            if not detail or md_idx >= len(detail):
                raise ValueError("MD tag is shorter than CIGAR deletion")

            d_op, d_len = detail[md_idx]
            if d_op != "D":
                raise ValueError("CIGAR deletion is inconsistent with MD tag")

            remain = d_len - md_off
            take = min(need, remain)
            need -= take
            md_off += take
            if md_off == d_len:
                md_idx += 1
                md_off = 0

    for op, l in ops:
        if op == "M":
            if detail:
                take_md_aligned(l, emit=True)
            else:
                out.append(Segment("match", "M", l, l, l))

        elif op in ("=", "X"):
            if detail:
                # Keep the MD cursor synchronized even though =/X are already
                # explicit in the CIGAR.
                take_md_aligned(l, emit=False)
            t = "match" if op == "=" else "mismatch"
            out.append(Segment(t, op, l, l, l))

        elif op == "I":
            out.append(Segment("ins", "I", l, 0, l))

        elif op == "D":
            if detail:
                take_md_deletion(l)
            out.append(Segment("del", "D", l, l, 0))

        elif op == "N":
            # N is represented by CIGAR and is not consumed from MD.
            out.append(Segment("ref_skip", "N", l, l, 0))

        elif op == "S":
            out.append(Segment("soft", "S", l, 0, l))

        elif op == "H":
            out.append(Segment("hard", "H", l, 0, 0))

        elif op == "P":
            out.append(Segment("pad", "P", l, 0, 0))

        else:
            raise ValueError("unsupported op")

    return merge_segments(out)

# Base-by-base alignment using CIGAR string, read sequence, and reference
# sequence to produce a detailed Segment list (distinguishing match and mismatch).
def from_cigar_with_ref(cigar: str, read_seq: str, ref_seq: str) -> List[Segment]:
    ops = parse_cigar_string(cigar)
    out: List[Segment] = []
    ri = 0
    fi = 0
    for op, l in ops:
        if op == "M":
            k = 0
            while k < l:
                m = 0
                while m < l - k and ri + m < len(read_seq) and fi + m < len(ref_seq) and read_seq[ri + m].upper() == ref_seq[fi + m].upper():
                    m += 1
                if m > 0:
                    out.append(Segment("match", "=", m, m, m))
                    k += m
                    ri += m
                    fi += m
                    if k >= l:
                        break
                x = 0
                while x < l - k and ri + x < len(read_seq) and fi + x < len(ref_seq) and read_seq[ri + x].upper() != ref_seq[fi + x].upper():
                    x += 1
                if x > 0:
                    out.append(Segment("mismatch", "X", x, x, x))
                    k += x
                    ri += x
                    fi += x
        elif op == "=":
            out.append(Segment("match", "=", l, l, l))
            ri += l
            fi += l
        elif op == "X":
            out.append(Segment("mismatch", "X", l, l, l))
            ri += l
            fi += l
        elif op == "I":
            ins_seq = read_seq[ri:ri+l] if read_seq and ri+l <= len(read_seq) else None
            out.append(Segment("ins", "I", l, 0, l, read_seq=ins_seq))
            ri += l
        elif op == "D":
            out.append(Segment("del", "D", l, l, 0))
            fi += l
        elif op == "N":
            out.append(Segment("ref_skip", "N", l, l, 0))
            fi += l
        elif op == "S":
            out.append(Segment("soft", "S", l, 0, l))
            ri += l
        elif op == "H":
            out.append(Segment("hard", "H", l, 0, 0))
        elif op == "P":
            out.append(Segment("pad", "P", l, 0, 0))
        else:
            raise ValueError("unsupported op")
    return merge_segments(out)
