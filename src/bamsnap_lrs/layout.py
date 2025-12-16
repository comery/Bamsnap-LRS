from dataclasses import dataclass
from typing import List, Tuple

from .cigar import Segment


@dataclass
class DrawRect:
    x0: int
    x1: int
    y: int
    h: int
    type: str


def segments_to_pixels(segments: List[Segment], read_start: int, region_start: int, bp_per_px: float, detail: str = "mid") -> List[Tuple[str, int, int]]:
    out: List[Tuple[str, int, int]] = []
    ref_cursor = read_start
    last_x1 = None  # Track the x1 of the previous segment to avoid overlap
    for s in segments:
        if s.ref_consumed == 0:
            if s.type == "ins":
                x = int((ref_cursor - region_start) / bp_per_px)
                out.append(("ins", x, x))
            continue
        x0 = int((ref_cursor - region_start) / bp_per_px)
        # If the previous segment's x1 is greater than current x0, use x1 as starting point to avoid overlap
        if last_x1 is not None and last_x1 > x0:
            x0 = last_x1
        x1 = int((ref_cursor + s.ref_consumed - region_start) / bp_per_px)
        x1 = max(x1, x0 + 1)
        t = s.type
        if detail == "low" and t == "mismatch":
            t = "match"
        out.append((t, x0, x1))
        last_x1 = x1
        ref_cursor += s.ref_consumed
    return out


def assign_stacks(read_spans: List[Tuple[int, int]], max_stack: int) -> List[int]:
    stacks: List[List[Tuple[int, int]]] = [[] for _ in range(max_stack)]
    res: List[int] = []
    for s, e in read_spans:
        placed = -1
        for i in range(max_stack):
            if not stacks[i] or stacks[i][-1][1] <= s:
                stacks[i].append((s, e))
                placed = i
                break
        if placed == -1:
            res.append(max_stack - 1)
        else:
            res.append(placed)
    return res
