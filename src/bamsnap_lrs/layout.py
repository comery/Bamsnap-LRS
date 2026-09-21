from dataclasses import dataclass
from typing import List, Tuple, Dict, Any, Optional

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
    for s in segments:
        if s.ref_consumed == 0:
            if s.type == "ins":
                x = int((ref_cursor - region_start) / bp_per_px)
                out.append(("ins", x, x))
            continue
        x0 = int((ref_cursor - region_start) / bp_per_px)
        x1 = int((ref_cursor + s.ref_consumed - region_start) / bp_per_px)
        x1 = max(x1, x0 + 1)
        t = s.type
        if detail == "low" and t == "mismatch":
            t = "match"
        out.append((t, x0, x1))
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


def assign_packed_stacks(read_spans: List[Tuple[int, int]]) -> List[int]:
    """Greedily pack non-overlapping read spans while preserving input priority.

    Unlike ``assign_stacks()``, this helper does not assume that spans arrive in
    genomic order.  That matters in highlight mode because reads are first
    sorted by haplotype signature rather than by start coordinate.  The input
    order therefore determines which reads get first choice of the upper rows,
    while any later read may reuse an existing row when its genomic interval
    does not overlap intervals already placed on that row.

    Intervals are treated as half-open [start, end), so touching boundaries do
    not count as overlap.
    """
    rows: List[List[Tuple[int, int]]] = []
    result: List[int] = []

    for start, end in read_spans:
        if end < start:
            start, end = end, start

        placed_row = None
        for row_idx, occupied in enumerate(rows):
            overlaps = any(start < occ_end and occ_start < end
                           for occ_start, occ_end in occupied)
            if not overlaps:
                occupied.append((start, end))
                occupied.sort(key=lambda x: (x[0], x[1]))
                placed_row = row_idx
                break

        if placed_row is None:
            rows.append([(start, end)])
            placed_row = len(rows) - 1

        result.append(placed_row)

    return result


def assign_stacks_grouped(
    read_spans: List[Tuple[int, int]],
    group_keys: List,
) -> List[int]:
    """Assign stack rows while keeping alignments from the same qname compact.

    Rules:
      1. Alignments with the same qname are treated as one visual group.
      2. Non-overlapping pieces from the same qname share one row.
      3. Overlapping pieces from the same qname are placed on adjacent rows.
      4. Different qname groups are packed globally to save vertical space.
    """
    n = len(read_spans)
    if n == 0:
        return []

    # Build qname groups, preserving first appearance.
    group_to_members: Dict = {}
    for i, key in enumerate(group_keys):
        group_to_members.setdefault(key, []).append(i)

    # Sort groups by genomic left boundary.
    ordered_groups = sorted(
        group_to_members.values(),
        key=lambda members: min(read_spans[m][0] for m in members)
    )

    res = [0] * n

    # row_free_at[row] = rightmost genomic coordinate already occupied on this row.
    row_free_at: List[int] = []

    def ensure_rows(k: int):
        while len(row_free_at) < k:
            row_free_at.append(0)

    for members in ordered_groups:
        members_sorted = sorted(
            members,
            key=lambda m: (read_spans[m][0], read_spans[m][1])
        )

        # Local layout within one qname group.
        # Non-overlapping pieces share the same local row.
        local_rows: List[List[int]] = []
        local_row_end: List[int] = []

        for m in members_sorted:
            s, e = read_spans[m]
            placed = False

            for j in range(len(local_rows)):
                if local_row_end[j] <= s:
                    local_rows[j].append(m)
                    local_row_end[j] = e
                    placed = True
                    break

            if not placed:
                local_rows.append([m])
                local_row_end.append(e)

        k = len(local_rows)

        local_min_start = [
            min(read_spans[m][0] for m in row)
            for row in local_rows
        ]
        local_max_end = [
            max(read_spans[m][1] for m in row)
            for row in local_rows
        ]

        # Place this qname group as a compact block into global rows.
        base_row = 0
        while True:
            ensure_rows(base_row + k)
            ok = True

            for j in range(k):
                if row_free_at[base_row + j] > local_min_start[j]:
                    ok = False
                    break

            if ok:
                break

            base_row += 1

        for j, row in enumerate(local_rows):
            global_row = base_row + j
            for m in row:
                res[m] = global_row
            row_free_at[global_row] = max(row_free_at[global_row], local_max_end[j])

    return res

def _query_span_from_segments(read: Any) -> Tuple[int, int]:
    """Return the alignment span on the original read using terminal CIGAR clips.

    This is used as a robust fallback for split/supplementary alignments. In
    particular, reverse-strand supplementary records with hard clipping can
    have misleading query_start/query_end values if those values were derived
    only from query_alignment_start/query_alignment_end.
    """
    segments = getattr(read, "segments", None) or []
    q_len = int(getattr(read, "query_length", 0) or 0)
    if not segments or q_len <= 0:
        return int(getattr(read, "query_start", 0) or 0), int(getattr(read, "query_end", 0) or 0)

    left_clip = 0
    for s in segments:
        if s.type not in ("soft", "hard"):
            break
        left_clip += int(getattr(s, "length", 0) or 0)

    right_clip = 0
    for s in reversed(segments):
        if s.type not in ("soft", "hard"):
            break
        right_clip += int(getattr(s, "length", 0) or 0)

    if getattr(read, "reverse", False):
        q_start, q_end = right_clip, q_len - left_clip
    else:
        q_start, q_end = left_clip, q_len - right_clip

    if q_end < q_start:
        return int(getattr(read, "query_start", 0) or 0), int(getattr(read, "query_end", 0) or 0)
    return q_start, q_end


def query_order_key(read: Any):
    """Sort split alignments by their physical order on the original read."""
    q_start, q_end = _query_span_from_segments(read)
    return (q_start, q_end, read.start, read.end)


def ref_pos_for_query_side(read: Any, side: str) -> int:
    """
    Convert one side of the original-read alignment span to a reference position.

    Forward alignment:
        query_start -> ref_start
        query_end   -> ref_end

    Reverse alignment:
        query_start -> ref_end
        query_end   -> ref_start
    """
    if side == "start":
        return read.end if read.reverse else read.start
    if side == "end":
        return read.start if read.reverse else read.end
    raise ValueError("side must be 'start' or 'end'")


def _norm_interval(a: int, b: int) -> Tuple[int, int]:
    """Return a half-open interval with increasing coordinates."""
    return (a, b) if a <= b else (b, a)


def _clip_interval(interval: Tuple[int, int], region_start: int, region_end: int) -> Tuple[int, int]:
    s, e = interval
    return max(region_start, s), min(region_end, e)


def _interval_overlaps(a: Tuple[int, int], b: Tuple[int, int]) -> bool:
    """Half-open interval overlap. Touching boundaries are allowed."""
    return a[0] < b[1] and b[0] < a[1]


def _interval_is_clear(interval: Tuple[int, int], occupied: List[Tuple[int, int]]) -> bool:
    """Return True if interval does not overlap any occupied interval."""
    s, e = interval
    if e <= s:
        return True
    for occ in occupied:
        if _interval_overlaps(interval, occ):
            return False
    return True


def _all_intervals_are_clear(
    intervals: List[Tuple[int, int]],
    occupied: List[Tuple[int, int]]
) -> bool:
    for interval in intervals:
        if not _interval_is_clear(interval, occupied):
            return False
    return True


def assign_split_read_stacks(
    reads: List[Any],
    region_start: int,
    region_end: int,
) -> List[int]:
    """
    Assign stack rows for DNA supplementary/split alignments.

    v7 policy: supplementary/split reads are treated as read-groups, and rows
    used by those groups are kept separate from ordinary single-alignment reads.

    Why this layout is used:
      1. Alignments from the same qname are ordered by their original-read
         query coordinates, so the top-to-bottom order still carries the order
         of segments on the original read.
      2. Each primary/supplementary piece of a multi-alignment read gets its
         own local row.  Pieces are not collapsed into one row.
      3. The whole qname group is reserved as an adjacent global block.  This
         prevents unrelated reads from being inserted into empty genomic gaps
         inside a split-read group.
      4. Rows are typed as either "supp" or "normal".  Ordinary reads are
         placed only on normal rows, while supplementary/split groups are
         placed only on supp rows.  This avoids the previous issue where
         compact spacing for supplementary rows also affected ordinary gray
         reads sharing the same row.
    """
    n = len(reads)
    if n == 0:
        return []

    spans: List[Tuple[int, int]] = []
    for r in reads:
        s = max(r.start, region_start)
        e = min(r.end, region_end)
        if e <= s:
            s, e = r.start, r.end
        spans.append((s, e))

    group_to_members: Dict[Any, List[int]] = {}
    for i, r in enumerate(reads):
        group_to_members.setdefault(r.qname, []).append(i)

    ordered_groups = sorted(
        group_to_members.values(),
        key=lambda members: min(spans[m][0] for m in members)
    )

    res = [0] * n
    global_rows: List[List[Tuple[int, int]]] = []
    # None = unused row, "normal" = ordinary reads only, "supp" = supplementary/split groups only.
    row_kinds: List[Optional[str]] = []

    def ensure_global_rows(k: int):
        while len(global_rows) < k:
            global_rows.append([])
            row_kinds.append(None)

    def row_is_clear(intervals: List[Tuple[int, int]], global_row: int) -> bool:
        for interval in intervals:
            if not _interval_is_clear(interval, global_rows[global_row]):
                return False
        return True

    def is_supp_group(members: List[int]) -> bool:
        return len(members) > 1 or any(getattr(reads[m], "supplementary", False) for m in members)

    def row_can_hold_kind(row: int, kind: str) -> bool:
        return row_kinds[row] is None or row_kinds[row] == kind

    for members in ordered_groups:
        if is_supp_group(members):
            # Supplementary/split group: preserve query order as adjacent rows.
            members_sorted = sorted(members, key=lambda m: query_order_key(reads[m]))
            local_rows: List[List[int]] = [[m] for m in members_sorted]
            k = len(local_rows)
            if k == 0:
                continue

            # Keep one qname as one visual group by reserving the full visible
            # span on every local row.  Ordinary reads are not allowed on these
            # rows at all; other supplementary groups may reuse them only when
            # their reserved spans do not overlap.
            group_left = min(spans[m][0] for m in members_sorted)
            group_right = max(spans[m][1] for m in members_sorted)
            group_span = _clip_interval((group_left, group_right), region_start, region_end)
            if group_span[1] <= group_span[0]:
                group_span = (group_left, group_right)

            global_reserve: List[List[Tuple[int, int]]] = [[group_span] for _ in range(k)]

            base_row = 0
            while True:
                ensure_global_rows(base_row + k)
                ok = True
                for j in range(k):
                    row = base_row + j
                    if not row_can_hold_kind(row, "supp") or not row_is_clear(global_reserve[j], row):
                        ok = False
                        break
                if ok:
                    break
                base_row += 1

            for j, row_members in enumerate(local_rows):
                global_row = base_row + j
                row_kinds[global_row] = "supp"
                for m in row_members:
                    res[m] = global_row
                global_rows[global_row].extend(global_reserve[j])
            continue

        # Ordinary single-alignment reads: place only on normal rows.  They are
        # no longer allowed to share rows with supplementary/split read groups.
        m = members[0]
        block_interval = spans[m]
        row = 0
        while True:
            ensure_global_rows(row + 1)
            if row_can_hold_kind(row, "normal") and _interval_is_clear(block_interval, global_rows[row]):
                break
            row += 1
        row_kinds[row] = "normal"
        res[m] = row
        global_rows[row].append(block_interval)

    return res

def assign_bed_stacks(feature_spans: List[Tuple[int, int]], min_distance: int = 10, max_stack: int = None) -> List[int]:
    """
    Assign stack levels for BED features with distance threshold.
    
    If two consecutive features are closer than min_distance, they will be placed
    on different stack levels to avoid visual overlap.
    
    Args:
        feature_spans: List of (start, end) tuples for features
        min_distance: Minimum distance (in bp) between features to allow same stack level
        max_stack: Maximum number of stack levels (default: len(feature_spans))
    
    Returns:
        List of stack indices for each feature
    """
    if not feature_spans:
        return []
    
    if max_stack is None:
        max_stack = len(feature_spans)
    
    stacks: List[List[Tuple[int, int]]] = [[] for _ in range(max_stack)]
    res: List[int] = []
    
    for idx, (s, e) in enumerate(feature_spans):
        placed = -1
        
        # Try to find a stack level where this feature can be placed
        for i in range(max_stack):
            if not stacks[i]:
                # Empty stack, can place here
                stacks[i].append((s, e))
                placed = i
                break
            else:
                # Check if this feature can be placed after the last feature in this stack
                last_end = stacks[i][-1][1]
                distance = s - last_end
                
                # If features don't overlap and distance is >= min_distance, can place here
                if last_end <= s and distance >= min_distance:
                    stacks[i].append((s, e))
                    placed = i
                    break
                # If features overlap or too close, try next stack
        
        # If couldn't place, use the last available stack
        if placed == -1:
            placed = min(max_stack - 1, idx)
            stacks[placed].append((s, e))
        
        res.append(placed)
    
    return res
