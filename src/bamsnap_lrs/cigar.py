from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass
class Segment:
    type: str
    op: str
    length: int
    ref_consumed: int
    read_consumed: int


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


def parse_cs(cs: str) -> List[Tuple[str, int]]:
    res = []
    i = 0
    n = len(cs)
    while i < n:
        c = cs[i]
        if c == ":":
            j = i + 1
            while j < n and cs[j].isdigit():
                j += 1
            res.append(("=", int(cs[i + 1:j])))
            i = j
        elif c == "*":
            res.append(("X", 1))
            i += 3
        elif c == "+":
            j = i + 1
            while j < n and cs[j].isalpha():
                j += 1
            res.append(("I", j - (i + 1)))
            i = j
        elif c == "-":
            j = i + 1
            while j < n and cs[j].isalpha():
                j += 1
            res.append(("D", j - (i + 1)))
            i = j
        else:
            i += 1
    return res


def from_cigar_md_cs(cigar: str, md: Optional[str] = None, cs: Optional[str] = None) -> List[Segment]:
    ops = parse_cigar_string(cigar)
    out: List[Segment] = []
    if cs:
        detail = parse_cs(cs)
        for d_op, d_len in detail:
            if d_op == "=":
                out.append(Segment("match", "=", d_len, d_len, d_len))
            elif d_op == "X":
                out.append(Segment("mismatch", "X", 1, 1, 1))
            elif d_op == "I":
                out.append(Segment("ins", "I", d_len, 0, d_len))
            elif d_op == "D":
                out.append(Segment("del", "D", d_len, d_len, 0))
        return merge_segments(out)
    detail = parse_md(md) if md else None
    for op, l in ops:
        if op in ("=", "X"):
            t = "match" if op == "=" else "mismatch"
            out.append(Segment(t, op, l, l, l))
        elif op == "M":
            if detail:
                k = 0
                rem = l
                while rem > 0 and k < len(detail):
                    d_op, d_len = detail[k]
                    take = min(rem, d_len)
                    if d_op == "=":
                        out.append(Segment("match", "=", take, take, take))
                    elif d_op == "X":
                        out.append(Segment("mismatch", "X", take, take, take))
                    elif d_op == "D":
                        out.append(Segment("del", "D", take, take, 0))
                    rem -= take
                    if take == d_len:
                        k += 1
                    else:
                        detail[k] = (d_op, d_len - take)
            else:
                out.append(Segment("match", "M", l, l, l))
        elif op == "I":
            out.append(Segment("ins", "I", l, 0, l))
        elif op == "D":
            out.append(Segment("del", "D", l, l, 0))
        elif op == "N":
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
            out.append(Segment("ins", "I", l, 0, l))
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