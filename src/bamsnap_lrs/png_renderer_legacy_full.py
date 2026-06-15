from typing import List, Dict, Optional, Tuple, Any

# from PIL import Image, ImageDraw, ImageFont

from .layout import (
    segments_to_pixels,
    assign_stacks,
    assign_bed_stacks,
    assign_stacks_grouped,
    assign_split_read_stacks,
    query_order_key,
    ref_pos_for_query_side,
)
from .reader import Read
from .pileup import base_pileup, pileup_to_pixels
from .styles import color_for_type, color_for_base, shade_by_mapq, STRAND_COLORS, MISMATCH_COLORS
from .highlight import HighlightSite, SampleCall


def _read_base_at_pos(read: Read, pos: int) -> Optional[str]:
    """Walk read.segments to find which read base is aligned to reference
    position `pos` (0-based). Returns the uppercase base, or None if the
    read does not cover that position with a match/mismatch segment (e.g.
    deletion or outside read range).
    """
    if not read.seq:
        return None
    if pos < read.start or pos >= read.end:
        return None
    ref_cur = read.start
    read_cur = 0
    for s in read.segments:
        if s.ref_consumed > 0 and s.read_consumed > 0:
            # match or mismatch — read covers this stretch
            if ref_cur <= pos < ref_cur + s.ref_consumed:
                k = pos - ref_cur
                idx = read_cur + k
                if 0 <= idx < len(read.seq):
                    return read.seq[idx].upper()
                return None
        elif s.ref_consumed > 0 and s.read_consumed == 0:
            # deletion in read — position is "absent" on the read
            if ref_cur <= pos < ref_cur + s.ref_consumed:
                return None
        ref_cur += s.ref_consumed
        read_cur += s.read_consumed
    return None


def _hap_signature(read: Read, site_positions: List[int]) -> str:
    """Build a per-read signature string from the read's base at each VCF
    site. '?' means the read does not cover that site. Used to cluster reads
    by haplotype: reads with the same signature share the same observed
    haplotype across the displayed SNPs.
    """
    out = []
    for p in site_positions:
        b = _read_base_at_pos(read, p)
        out.append(b if b in ("A", "C", "G", "T") else "?")
    return "".join(out)

def _focus_region_to_pixels(
    focus_region: Optional[Tuple[int, int]],
    start: int,
    end: int,
    bp_per_px: float,
    margin: int,
    width: int,
) -> Optional[Tuple[int, int]]:
    """
    Convert a genomic focus region to drawable x-pixel coordinates.
    """
    if not focus_region:
        return None

    focus_start, focus_end = focus_region

    x0 = margin + int((focus_start - start) / bp_per_px)
    x1 = margin + int((focus_end - start) / bp_per_px)

    x0 = max(margin, min(width - margin, x0))
    x1 = max(margin, min(width - margin, x1))

    # Make sure tiny regions are still visible
    if x1 <= x0:
        x1 = min(width - margin, x0 + 2)

    return x0, x1


# Pale gray used for mismatches that fall OUTSIDE the VCF highlight sites.
# Slightly lighter than the match color so the mismatch shape is still
# perceptible, but visually almost merges with the rest of the read.
MUTED_MISMATCH_COLOR = (200, 200, 200)


# Supplementary/split-read fill colors.  Reads without supplementary/split
# alignments keep the original light-gray body.  For reads with visible
# supplementary pieces, each qname gets one color and pieces along the original
# read get gradually lower opacity according to query_start/query_end order.
SUPPLEMENTARY_READ_COLORS = [
    (80, 170, 120),   # green
    (70, 150, 210),   # blue
    (240, 170, 70),   # orange
    (165, 120, 205),  # purple
    (220, 110, 110),  # red/salmon
    (60, 190, 180),   # cyan
    (190, 160, 80),   # olive/gold
    (120, 170, 210),  # light blue
    (150, 190, 95),   # yellow-green
    (205, 130, 170),  # pink
]
SUPPLEMENTARY_OPACITY_STEP = 0.20       # legacy compatibility; opacity now uses piece-count-based spacing
SUPPLEMENTARY_MAX_TRANSPARENCY = 0.80   # user-adjustable cap, now applied directly


def _is_supplementary_read_group(reads: List[Read], idxs: List[int]) -> bool:
    """Only treat a qname as a supplementary/split group when multiple
    alignment pieces of that qname are visible in the current window.

    A single visible piece, even if that piece itself has supplementary=True,
    is rendered like an ordinary read because the rest of the group is not
    present in the current view.
    """
    return len(idxs) > 1


def _build_supplementary_read_styles(
    reads: List[Read],
    groups: Dict[str, List[int]],
    stacks: Optional[List[int]] = None,
) -> Dict[int, Tuple[Tuple[int, int, int], int]]:
    """Return {read_index: (rgb, alpha)} for supplementary/split reads.

    The color is constant within one qname.  Alpha follows the alignment order
    on the original read and mirrors svg_renderer.py: two-piece groups use
    0.00/0.40 transparency, while groups with three or more pieces are evenly
    distributed from 0.00 to SUPPLEMENTARY_MAX_TRANSPARENCY.  Colors are
    assigned greedily so nearby read groups avoid reusing the same color when
    the palette allows it.
    """
    candidates = []
    for qname, idxs in groups.items():
        if not _is_supplementary_read_group(reads, idxs):
            continue
        if stacks is not None and idxs:
            group_rows = [stacks[i] for i in idxs]
            row_min = min(group_rows)
            row_max = max(group_rows)
            center = (row_min + row_max) / 2.0
        else:
            row_min = row_max = min(idxs) if idxs else 0
            center = float(row_min)
        left = min(getattr(reads[i], "start", 0) for i in idxs)
        candidates.append((center, row_min, row_max, left, qname, idxs))

    candidates.sort(key=lambda x: (x[0], x[3], x[4]))
    styles: Dict[int, Tuple[Tuple[int, int, int], int]] = {}
    placed: List[Tuple[float, int, int]] = []  # center, row range, color index

    for group_order, (center, row_min, row_max, _left, _qname, idxs) in enumerate(candidates):
        nearby_colors = {
            color_idx
            for prev_center, _prev_row_min, _prev_row_max, color_idx in placed
            if abs(prev_center - center) <= 4
        }
        color_idx = group_order % len(SUPPLEMENTARY_READ_COLORS)
        for offset in range(len(SUPPLEMENTARY_READ_COLORS)):
            trial = (group_order + offset) % len(SUPPLEMENTARY_READ_COLORS)
            if trial not in nearby_colors:
                color_idx = trial
                break

        color = SUPPLEMENTARY_READ_COLORS[color_idx]
        idxs_sorted = sorted(idxs, key=lambda i: query_order_key(reads[i]))
        n_pieces = len(idxs_sorted)

        # Keep PNG opacity exactly consistent with svg_renderer.py:
        #   1 piece  -> 0.00 transparency, though single-piece groups are not styled
        #   2 pieces -> 0.00, 0.40
        #   >=3      -> evenly distributed from 0.00 to SUPPLEMENTARY_MAX_TRANSPARENCY
        if n_pieces <= 1:
            transparency_values = [0.0]
        elif n_pieces == 2:
            transparency_values = [0.0, 0.4]
        else:
            transparency_values = [
                SUPPLEMENTARY_MAX_TRANSPARENCY * piece_order / (n_pieces - 1)
                for piece_order in range(n_pieces)
            ]

        min_alpha = max(0, min(255, int(round(255 * (1.0 - SUPPLEMENTARY_MAX_TRANSPARENCY)))))
        for idx, transparency in zip(idxs_sorted, transparency_values):
            transparency = max(0.0, min(SUPPLEMENTARY_MAX_TRANSPARENCY, transparency))
            alpha = max(min_alpha, min(255, int(round(255 * (1.0 - transparency)))))
            styles[idx] = (color, alpha)
        placed.append((center, row_min, row_max, color_idx))

    return styles


def _draw_rect_with_alpha(
    dr: ImageDraw.ImageDraw,
    dr_alpha: ImageDraw.ImageDraw,
    box,
    color: Tuple[int, int, int],
    alpha: int,
):
    if alpha >= 255:
        dr.rectangle(box, fill=color)
    else:
        dr_alpha.rectangle(box, fill=(color[0], color[1], color[2], alpha))


def _draw_polygon_with_alpha(
    dr: ImageDraw.ImageDraw,
    dr_alpha: ImageDraw.ImageDraw,
    points,
    color: Tuple[int, int, int],
    alpha: int,
):
    if alpha >= 255:
        dr.polygon(points, fill=color)
    else:
        dr_alpha.polygon(points, fill=(color[0], color[1], color[2], alpha))


# Vertical layout tuning for read rows.
# - Alignments with the same qname (primary + supplementary pieces of one read)
#   stay compact when they occupy adjacent stack rows.
# - Adjacent rows belonging to different reads get a larger visual gap.
def _stack_row_offsets(
    stacks: List[int],
    reads: List[Read],
    read_height: int,
    same_read_gap: int = 1,
    different_read_gap: int = 4,
    bottom_padding: int = 5,
) -> Tuple[List[int], int]:
    """Return per-stack y offsets for read rows.

    v6 policy: supplementary/split alignments are allowed to be compact only
    when the two neighbouring stack rows are both *pure rows* of the same
    supplementary qname.  If either row also contains ordinary reads, the row
    spacing stays normal.  This keeps the compact split-read appearance without
    pulling unrelated gray reads on the same stack rows closer together.
    """
    num_rows = (max(stacks) + 1) if stacks else 1
    row_qnames = [set() for _ in range(num_rows)]
    row_intervals = [[] for _ in range(num_rows)]

    groups: Dict[str, List[int]] = {}
    for i, r in enumerate(reads):
        groups.setdefault(r.qname, []).append(i)
    supplementary_qnames = {
        qname for qname, idxs in groups.items()
        if _is_supplementary_read_group(reads, idxs)
    }

    for r, stack in zip(reads, stacks):
        if 0 <= stack < num_rows:
            row_qnames[stack].add(r.qname)
            s = getattr(r, "start", 0)
            e = getattr(r, "end", s)
            if e > s:
                row_intervals[stack].append((s, e, r.qname))

    def intervals_overlap(a, b) -> bool:
        return a[0] < b[1] and b[0] < a[1]

    def shared_qname_blocks_overlap(row_a: int, row_b: int, shared_name: str) -> bool:
        a_blocks = [x for x in row_intervals[row_a] if x[2] == shared_name]
        b_blocks = [x for x in row_intervals[row_b] if x[2] == shared_name]
        for a in a_blocks:
            for b in b_blocks:
                if intervals_overlap(a, b):
                    return True
        return False

    def is_pure_same_supp_pair(row_a: int, row_b: int):
        prev_names = row_qnames[row_a]
        cur_names = row_qnames[row_b]
        shared = (prev_names & cur_names) & supplementary_qnames
        if len(shared) != 1:
            return False, None
        qname = next(iter(shared))
        # Important: if there are ordinary reads, or another qname, on either
        # row, do not compress the row itself.  Otherwise unrelated reads on
        # the same stack row would move together with the split-read group.
        if prev_names != {qname} or cur_names != {qname}:
            return False, None
        return True, qname

    nonoverlap_same_read_gap = 0

    offsets = [0] * num_rows
    for row in range(1, num_rows):
        pure_pair, qname = is_pure_same_supp_pair(row - 1, row)
        if pure_pair and qname is not None:
            if shared_qname_blocks_overlap(row - 1, row, qname):
                gap = same_read_gap
            else:
                gap = nonoverlap_same_read_gap
        else:
            gap = different_read_gap

        offsets[row] = offsets[row - 1] + read_height + gap

    reads_area_h = offsets[-1] + read_height + bottom_padding if offsets else read_height + bottom_padding
    return offsets, reads_area_h

def _draw_per_track_coverage_png(
    img: Image.Image,
    dr: ImageDraw.ImageDraw,
    reads: List[Read],
    title: str,
    current_y: int,
    width: int,
    content_width: int,
    coverage_height: int,
    start: int,
    end: int,
    ref_seq: Optional[str],
    coverage_max_depth: Optional[int],
    margin: int,
    is_rna: bool,
    detail: str,
    bp_per_px: float,
) -> Tuple[int, Image.Image, ImageDraw.ImageDraw]:
    """Draw one coverage track for one BAM/track and return updated y/img/draw."""
    coverage_title = f"Coverage: {title}" if title else "Coverage"
    current_y += draw_track_header(dr, coverage_title, current_y, width, coverage_height)

    # Leave vertical room for RNA splice/split-read arcs above the depth bars.
    coverage_track_top = current_y + coverage_height
    current_y = coverage_track_top

    base_distribution = calculate_base_distribution(reads, start, end, ref_seq)
    draw_coverage_track(
        dr,
        base_distribution,
        current_y,
        content_width,
        coverage_height,
        start,
        end,
        ref_seq,
        coverage_max_depth,
        margin,
        is_rna,
        detail,
    )

    if is_rna:
        arc_layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
        dr_arc = ImageDraw.Draw(arc_layer)
        coverage_track_bottom = current_y + coverage_height
        arc_color = (253, 209, 211, 150)
        arc_anchor_y = coverage_track_bottom - coverage_height // 2

        # Arcs for ref_skip operations within reads.
        for r in reads:
            ref_cursor = r.start
            for seg in r.segments:
                if seg.type == "ref_skip":
                    seg_start = ref_cursor
                    seg_end = ref_cursor + seg.ref_consumed
                    xa = margin + int((seg_start - start) / bp_per_px)
                    xb = margin + int((seg_end - start) / bp_per_px)
                    xa = max(margin, min(width - margin - 1, xa))
                    xb = max(margin, min(width - margin - 1, xb))
                    if xb - xa >= 2:
                        w = xb - xa
                        h = min(w // 2, coverage_height)
                        bbox = [xa, arc_anchor_y - h, xb, arc_anchor_y + h]
                        dr_arc.arc(bbox, start=180, end=0, fill=arc_color, width=1)
                ref_cursor += seg.ref_consumed

        # Arcs for split alignments belonging to the same qname.
        groups_temp: Dict[str, List[int]] = {}
        for idx_r, r in enumerate(reads):
            groups_temp.setdefault(r.qname, []).append(idx_r)
        for _qname, idxs in groups_temp.items():
            if len(idxs) > 1:
                idxs_sorted = sorted(idxs, key=lambda i: reads[i].start)
                for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                    if reads[a].end < reads[b].start:
                        xa = margin + int((reads[a].end - start) / bp_per_px)
                        xb = margin + int((reads[b].start - start) / bp_per_px)
                        xa = max(margin, min(width - margin - 1, xa))
                        xb = max(margin, min(width - margin - 1, xb))
                        if xb - xa >= 2:
                            w = xb - xa
                            h = min(w // 2, coverage_height)
                            bbox = [xa, arc_anchor_y - h, xb, arc_anchor_y + h]
                            dr_arc.arc(bbox, start=180, end=0, fill=arc_color, width=1)

        img = img.convert("RGBA")
        img = Image.alpha_composite(img, arc_layer)
        img = img.convert("RGB")
        dr = ImageDraw.Draw(img)

    current_y += coverage_height + 15
    return current_y, img, dr


def render_snapshot(
    tracks: List[Any],  # Can be List[Read] or List[Dict]
    chrom: str,
    start: int,
    end: int,
    width: int = 1200,
    read_height: int = 6,
    detail: str = "mid",
    show_axis: bool = True,
    show_composition: bool = False,
    composition_height: Optional[int] = None,
    comp_max_depth: Optional[int] = None,
    style: str = "jbrowse",
    color_by: str = "type",
    ref_seq: Optional[str] = None,
    show_coverage: bool = True,
    coverage_height: int = 15,  # Default reduced to 15
    track_title: str = "Reads",
    show_insertion_labels: bool = True,
    coverage_max_depth: Optional[int] = None,
    is_rna: bool = False,
    gff_genes: Optional[List[Any]] = None,
    bed_features: Optional[List[Any]] = None,
    highlight_sites: Optional[List[HighlightSite]] = None,
    highlight_samples: Optional[List[str]] = None,
    no_hap_sort: bool = False,
    no_hap_filter: bool = False,
    focus_region: Optional[Tuple[int, int]] = None,
):
    # Handle backward compatibility: if tracks is list of Read, wrap it
    if tracks and not isinstance(tracks[0], dict):
        tracks = [{'reads': tracks, 'title': track_title}]

    # If using jbrowse style and coverage enabled, use new rendering function
    if style == "jbrowse" and show_coverage:
        return render_jbrowse_style(
            tracks=tracks,
            chrom=chrom,
            start=start,
            end=end,
            width=width,
            read_height=read_height,
            is_rna=is_rna,
            detail=detail,
            show_axis=show_axis,
            show_coverage=show_coverage,
            coverage_height=coverage_height,
            style=style,
            color_by=color_by,
            ref_seq=ref_seq,
            show_insertion_labels=show_insertion_labels,
            coverage_max_depth=coverage_max_depth,
            gff_genes=gff_genes,
            bed_features=bed_features,
            highlight_sites=highlight_sites,
            highlight_samples=highlight_samples,
            no_hap_sort=no_hap_sort,
            no_hap_filter=no_hap_filter,
            focus_region=focus_region,
        )

    # For default style, flatten reads from all tracks
    reads = []
    for t in tracks:
        reads.extend(t['reads'])

    bp_per_px = (end - start) / float(width)
    if is_rna:
        spans = [(max(r.start, start), min(r.end, end)) for r in reads]
        keys = [r.qname for r in reads]
        stacks = assign_stacks_grouped(spans, keys)
    else:
        stacks = assign_split_read_stacks(reads, start, end)
    row_offsets, reads_area_h = _stack_row_offsets(stacks, reads, read_height)
    top = 0
    if show_axis:
        top += 20
    comp_h = composition_height if composition_height is not None else (read_height * 2)
    if show_composition:
        top += comp_h + 10
    height = top + reads_area_h
    img = Image.new("RGB", (width, height), (255, 255, 255))
    dr = ImageDraw.Draw(img)
    groups: Dict[str, List[int]] = {}
    for i, r in enumerate(reads):
        groups.setdefault(r.qname, []).append(i)
    supplementary_read_styles = _build_supplementary_read_styles(reads, groups, stacks)
    dr_alpha = ImageDraw.Draw(img, "RGBA")
    if show_axis:
        dr.line([(0, 10), (width - 1, 10)], fill=(0, 0, 0))
        span = end - start
        step = max(1, span // 10)
        for pos in range(start, end + 1, step):
            x = int((pos - start) / bp_per_px)
            dr.line([(x, 7), (x, 13)], fill=(0, 0, 0))
            dr.text((x + 2, 1), str(pos), fill=(0, 0, 0))
    
    focus_px = _focus_region_to_pixels(
        focus_region,
        start,
        end,
        bp_per_px,
        0,      # default style 这里没有左右 margin
        width,
    )

    if focus_px is not None:
        fx0, fx1 = focus_px
        shade_top = 20 if show_axis else 0
        dr.rectangle(
            [(fx0, shade_top), (fx1, height)],
            fill=(235, 242, 255)
        )
    
    if show_composition:
        pile = base_pileup(reads, start, end)
        bins = pileup_to_pixels(pile, width)
        off = 10 + (20 if show_axis else 0)
        maxd = comp_max_depth if comp_max_depth is not None else max((b["depth"] for b in bins), default=1)
        if maxd <= 0:
            maxd = 1
        for x, agg in enumerate(bins):
            d = max(1, agg["depth"])
            h = comp_h
            hx = int(h * min(d / maxd, 1.0))
            if hx <= 0:
                continue
            a = int(hx * agg["A"] / d)
            c = int(hx * agg["C"] / d)
            g = int(hx * agg["G"] / d)
            t = int(hx * agg["T"] / d)
            used = a + c + g + t
            n = max(0, hx - used)
            y = off + h
            if a > 0:
                dr.rectangle([(x, y - a), (x + 1, y)], fill=(80, 160, 80))
                y -= a
            if c > 0:
                dr.rectangle([(x, y - c), (x + 1, y)], fill=(80, 120, 200))
                y -= c
            if g > 0:
                dr.rectangle([(x, y - g), (x + 1, y)], fill=(240, 200, 70))
                y -= g
            if t > 0:
                dr.rectangle([(x, y - t), (x + 1, y)], fill=(220, 90, 90))
                y -= t
            if n > 0:
                dr.rectangle([(x, y - n), (x + 1, y)], fill=(160, 160, 160))
    for idx, r in enumerate(reads):
        y = top + row_offsets[stacks[idx]]
        rects = segments_to_pixels(r.segments, r.start, start, bp_per_px, detail=detail)
        supp_style = supplementary_read_styles.get(idx)
        for t, x0, x1 in rects:
            if color_by == "type":
                color = color_for_type(t)
            elif color_by == "strand":
                color = STRAND_COLORS["-" if r.reverse else "+"]
            elif color_by == "mapq":
                base_color = color_for_type(t)
                color = shade_by_mapq(base_color, r.mapq)
            else:
                color = color_for_type(t)
            if t == "ins":
                ins_line_width = max(2, read_height // 3)
                dr.line([(x0, y), (x0, y + read_height)], fill=color, width=ins_line_width)
            elif t == "ref_skip":
                # Intron: draw as line
                y_center = y + read_height // 2
                dr.line([(x0, y_center), (x1, y_center)], fill=(176, 196, 222), width=1)
            elif t == "del":
                dr.rectangle([(x0, y), (x1, y + read_height)], outline=(120, 120, 120), fill=color)
            else:
                if supp_style and t in ("match", "soft", "hard"):
                    fill_color, fill_alpha = supp_style
                    _draw_rect_with_alpha(dr, dr_alpha, [(x0, y), (x1, y + read_height)], fill_color, fill_alpha)
                else:
                    dr.rectangle([(x0, y), (x1, y + read_height)], fill=color)
        if color_by == "base" and detail == "high" and bp_per_px <= 1 and r.seq:
            ref_cursor = r.start
            read_cursor = 0
            for s in r.segments:
                if s.ref_consumed > 0 and s.read_consumed > 0:
                    for i in range(s.ref_consumed):
                        pos = ref_cursor + i
                        if pos < start or pos >= end:
                            continue
                        bi = read_cursor + i
                        if bi < 0 or bi >= len(r.seq):
                            continue
                        x = int((pos - start) / bp_per_px)
                        b = r.seq[bi].upper()
                        if ref_seq is not None:
                            rb = ref_seq[pos - start].upper()
                            bc = color_for_type("match") if b == rb else color_for_base(b)
                        else:
                            bc = color_for_base(b)
                        dr.rectangle([(x, y), (x + 1, y + read_height)], fill=bc)
                    ref_cursor += s.ref_consumed
                    read_cursor += s.read_consumed
                elif s.ref_consumed > 0:
                    ref_cursor += s.ref_consumed
                elif s.read_consumed > 0:
                    read_cursor += s.read_consumed
        if style == "jbrowse":
            w = x1 - x0
            if w >= 6:
                head = max(4, min(max(read_height + 1, 5), 8, max(1, w - 1)))
                arrow_color, arrow_alpha = ((255, 255, 255), 255) if supp_style else ((100, 100, 100), 255)
                if r.reverse:
                    _draw_polygon_with_alpha(dr, dr_alpha, [(x0, y), (x0 + head, y + read_height // 2), (x0, y + read_height)], arrow_color, arrow_alpha)
                else:
                    _draw_polygon_with_alpha(dr, dr_alpha, [(x1, y), (x1 - head, y + read_height // 2), (x1, y + read_height)], arrow_color, arrow_alpha)
    for qname, idxs in groups.items():
        if len(idxs) > 1:
            if is_rna:
                # Keep RNA connection semantics unchanged: reference-based splicing-style links.
                idxs_sorted = sorted(idxs, key=lambda i: (reads[i].start, reads[i].end))
                for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                    ya = top + row_offsets[stacks[a]]
                    yb = top + row_offsets[stacks[b]]
                    xa = int((reads[a].end - start) / bp_per_px)
                    xb = int((reads[b].start - start) / bp_per_px)
                    xa = max(0, min(width - 1, xa))
                    xb = max(0, min(width - 1, xb))
                    dr.line(
                        [(xa, ya + read_height // 2), (xb, yb + read_height // 2)],
                        fill=(80, 160, 80)
                    )
            else:
                # DNA supplementary/split-read connector lines are intentionally
                # disabled.  The relationship is now encoded by qname color plus
                # query-order opacity, so no extra linking line is drawn.
                pass
    return img


def calculate_coverage(reads: List[Read], start: int, end: int, width: int) -> List[int]:
    """Calculate coverage at each pixel position"""
    coverage = [0] * width
    bp_per_px = (end - start) / float(width)
    for r in reads:
        r_start_px = int((max(r.start, start) - start) / bp_per_px)
        r_end_px = int((min(r.end, end) - start) / bp_per_px)
        for x in range(max(0, r_start_px), min(width, r_end_px + 1)):
            coverage[x] += 1
    return coverage


def calculate_base_distribution(
    reads: List[Read], start: int, end: int, ref_seq: Optional[str] = None
) -> List[Dict[str, int]]:
    """Calculate base distribution at each position (no aggregation, original resolution)

    Returns base distribution at each position, including ref_match and variant base counts
    """
    # Calculate base distribution at each position
    pile = base_pileup(reads, start, end)

    # Calculate ref_match and variant bases for each position
    result = []
    for i, p in enumerate(pile):
        ref_base = None
        if ref_seq and i < len(ref_seq):
            ref_base = ref_seq[i].upper()

        ref_match = 0
        variant_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}

        for base in ["A", "C", "G", "T", "N"]:
            count = p.get(base, 0)
            if count > 0:
                if ref_base and base == ref_base:
                    ref_match += count
                else:
                    variant_counts[base] += count

        result.append({
            "ref_match": ref_match,
            "A": variant_counts["A"],
            "C": variant_counts["C"],
            "G": variant_counts["G"],
            "T": variant_counts["T"],
            "N": variant_counts["N"],
            "depth": p.get("depth", 0)
        })

    return result


def draw_track_header(
    dr: ImageDraw.ImageDraw,
    title: str,
    y: int,
    width: int,
    track_height: int,
    font_size: int = 11,
) -> int:
    """Draw track header bar, returns header height"""
    header_height = 15  # Reduced height
    # Background color (light gray)
    dr.rectangle([(0, y), (width, y + header_height)], fill=(245, 245, 245))
    # Bottom border
    dr.line([(0, y + header_height - 1), (width, y + header_height - 1)], fill=(200, 200, 200), width=1)
    # Title text
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", font_size)
    except:
        try:
            font = ImageFont.truetype("arial.ttf", font_size)
        except:
            font = ImageFont.load_default()
    dr.text((5, y + 1), title, fill=(0, 0, 0), font=font)
    return header_height


def draw_coverage_track(
    dr: ImageDraw.ImageDraw,
    base_distribution: List[Dict[str, int]],
    y: int,
    width: int,
    height: int,
    start: int,
    end: int,
    ref_seq: Optional[str] = None,
    max_depth: Optional[int] = None,
    margin: int = 0,  # Left margin
    is_rna: bool = False,
    detail: str = "mid",
) -> int:
    """Draw coverage stacked bar chart, showing variants based on reference (JBrowse style)

    base_distribution is distribution at each base position (original resolution)
    Drawing maps pixel positions to base positions
    - First draw variant bases (A→C→G→T→N)
    - Finally draw reference base (gray)
    """
    from .styles import MISMATCH_COLORS  # Use same colors as read mismatch

    if not base_distribution:
        return height

    num_bases = len(base_distribution)  # Number of base positions

    # Calculate maximum depth
    max_cov = max_depth
    if max_cov is None:
        max_cov = max((b["depth"] for b in base_distribution), default=1)
    if max_cov <= 0:
        max_cov = 1

    # Draw y-axis scale
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 8)
    except:
        try:
            font = ImageFont.truetype("arial.ttf", 8)
        except:
            font = ImageFont.load_default()

    # Draw y-axis line (within left margin)
    axis_x = margin - 2
    bar_bottom = y + height
    bar_top = y
    dr.line([(axis_x, bar_top), (axis_x, bar_bottom)], fill=(100, 100, 100), width=1)

    # Draw ticks (0, middle, maximum)
    # Maximum tick (top)
    dr.line([(axis_x - 3, bar_top), (axis_x, bar_top)], fill=(100, 100, 100), width=1)
    dr.text((2, bar_top - 2), str(max_cov), fill=(80, 80, 80), font=font)

    # Middle tick
    mid_y = (bar_top + bar_bottom) // 2
    mid_val = max_cov // 2
    dr.line([(axis_x - 3, mid_y), (axis_x, mid_y)], fill=(100, 100, 100), width=1)
    dr.text((2, mid_y - 4), str(mid_val), fill=(80, 80, 80), font=font)

    # Bottom tick (0)
    dr.line([(axis_x - 3, bar_bottom), (axis_x, bar_bottom)], fill=(100, 100, 100), width=1)
    dr.text((2, bar_bottom - 8), "0", fill=(80, 80, 80), font=font)

    # Draw by pixel position, each pixel aggregates all base positions it covers
    import math
    for px in range(width):
        # Calculate the range of base positions this pixel covers
        base_start_idx = int(px * num_bases / width)
        base_end_idx = int((px + 1) * num_bases / width)
        if base_end_idx <= base_start_idx:
            base_end_idx = base_start_idx + 1
        if base_end_idx > num_bases:
            base_end_idx = num_bases

        # Aggregate all base positions in this range
        agg_ref_match = 0
        agg_variants = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
        agg_depth = 0

        for base_idx in range(base_start_idx, base_end_idx):
            dist = base_distribution[base_idx]
            agg_ref_match += dist.get("ref_match", 0)
            for base in ["A", "C", "G", "T", "N"]:
                agg_variants[base] += dist.get(base, 0)
            agg_depth += dist.get("depth", 0)

        # Average over the number of positions (but use ceil for variants to preserve rare variants)
        num_positions = base_end_idx - base_start_idx
        if num_positions > 1:
            agg_ref_match = round(agg_ref_match / num_positions)
            for base in agg_variants:
                if agg_variants[base] > 0:
                    # Use ceil to ensure variants are not lost
                    agg_variants[base] = max(1, math.ceil(agg_variants[base] / num_positions))
            agg_depth = round(agg_depth / num_positions)

        if agg_depth == 0:
            continue

        # Actual drawing position (add margin)
        draw_x = margin + px

        # Calculate total bar height for this position (proportional to depth)
        bar_height = int(height * min(agg_depth / max_cov, 1.0))
        if bar_height <= 0:
            continue

        # Calculate height for each part
        base_heights = {}
        total_count = agg_ref_match
        for base in ["A", "C", "G", "T", "N"]:
            total_count += agg_variants[base]

        if total_count == 0:
            continue

        # Reference match portion
        if agg_ref_match > 0:
            base_heights["ref"] = int(bar_height * agg_ref_match / total_count)

        # Variant portion - ensure at least 1 pixel height if variant exists
        for base in ["A", "C", "G", "T", "N"]:
            count = agg_variants[base]
            if count > 0:
                h = int(bar_height * count / total_count)
                # Ensure at least 1 pixel for visible variants
                base_heights[base] = max(1, h) if h == 0 else h

        # Adjust heights to ensure full fill
        total_height_used = sum(base_heights.values())
        if total_height_used < bar_height and base_heights:
            max_base = max(base_heights.items(), key=lambda x: x[1])[0]
            base_heights[max_base] += (bar_height - total_height_used)
        elif total_height_used > bar_height:
            # If total exceeds bar_height (due to ensuring min 1px), reduce ref portion
            excess = total_height_used - bar_height
            if "ref" in base_heights and base_heights["ref"] > excess:
                base_heights["ref"] -= excess

        # Draw coverage bar (same for RNA and DNA mode)
        bar_bottom = y + height
        current_y = bar_bottom

        # If detail is low, just draw a simple gray depth bar
        if detail == "low":
            dr.rectangle([(draw_x, bar_bottom - bar_height), (draw_x + 1, bar_bottom)], fill=(180, 180, 180))
            continue

        # First draw variant portion (in fixed order A→C→G→T→N, with corresponding colors)
        for base in ["A", "C", "G", "T", "N"]:
            if base not in base_heights or base_heights[base] <= 0:
                continue
            h = base_heights[base]
            color = MISMATCH_COLORS.get(base, (160, 160, 160))
            dr.rectangle([(draw_x, current_y - h), (draw_x + 1, current_y)], fill=color)
            current_y -= h

        # Finally draw reference match portion (gray, on top)
        if "ref" in base_heights and base_heights["ref"] > 0:
            h = base_heights["ref"]
            dr.rectangle([(draw_x, current_y - h), (draw_x + 1, current_y)], fill=(180, 180, 180))

    return height


def draw_gene_track(
    dr: ImageDraw.ImageDraw,
    genes: List[Any],
    y: int,
    width: int,
    start: int,
    end: int,
    bp_per_px: float,
    margin: int,
    stacks: List[int],
) -> int:
    """Draw gene annotation track"""
    header_h = draw_track_header(dr, "Gene Annotation", y, width + 2 * margin, 15)
    current_y = y + header_h + 5

    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 10)
    except:
        try:
            font = ImageFont.truetype("arial.ttf", 10)
        except:
            font = ImageFont.load_default()

    for i, gene in enumerate(genes):
        stack = stacks[i]
        gene_y = current_y + stack * 20

        gx0 = margin + int((max(gene.start, start) - start) / bp_per_px)
        gx1 = margin + int((min(gene.end, end) - start) / bp_per_px)

        if gx1 <= gx0:
            continue

        # Colors for gene features
        color_utr = (100, 149, 237)  # Cornflower Blue
        color_cds = (200, 160, 50)   # Brownish Yellow
        feature_height = 10
        mid_y = gene_y + 5

        # Draw gene line (intron)
        dr.line([(gx0, mid_y), (gx1, mid_y)], fill=(0, 0, 0), width=1)

        # Strand arrow - solid black arrow outside the feature with connecting line
        head_size = 6
        arrow_width = 4
        connector_length = 3  # Short line connecting feature to arrow
        if gene.strand == '+':
            # Forward strand: arrow pointing right, outside the feature
            arrow_x = gx1 + connector_length
            # Draw connecting line
            dr.line([(gx1, mid_y), (arrow_x, mid_y)], fill=(0, 0, 0), width=1)
            # Draw arrow triangle (pointing right)
            dr.polygon([
                (arrow_x, mid_y - arrow_width),
                (arrow_x + head_size, mid_y),
                (arrow_x, mid_y + arrow_width)
            ], fill=(0, 0, 0), outline=(0, 0, 0))
        elif gene.strand == '-':
            # Reverse strand: arrow pointing left, outside the feature
            arrow_x = gx0 - connector_length
            # Draw connecting line
            dr.line([(gx0, mid_y), (arrow_x, mid_y)], fill=(0, 0, 0), width=1)
            # Draw arrow triangle (pointing left)
            dr.polygon([
                (arrow_x, mid_y - arrow_width),
                (arrow_x - head_size, mid_y),
                (arrow_x, mid_y + arrow_width)
            ], fill=(0, 0, 0), outline=(0, 0, 0))

        # Draw exons (UTR color)
        for exon in gene.exons:
            ex0 = margin + int((max(exon.start, start) - start) / bp_per_px)
            ex1 = margin + int((min(exon.end, end) - start) / bp_per_px)
            if ex1 > ex0:
                dr.rectangle([(ex0, mid_y - feature_height // 2), (ex1, mid_y + feature_height // 2)], fill=color_utr)

        # Draw CDS (CDS color, same height)
        for cds in gene.cds:
            cx0 = margin + int((max(cds.start, start) - start) / bp_per_px)
            cx1 = margin + int((min(cds.end, end) - start) / bp_per_px)
            if cx1 > cx0:
                dr.rectangle([(cx0, mid_y - feature_height // 2), (cx1, mid_y + feature_height // 2)], fill=color_cds)

        # Draw gene name
        dr.text((gx0, mid_y + 6), gene.name, fill=(0, 0, 0), font=font)

    num_stacks = max(stacks) + 1 if stacks else 0
    return header_h + num_stacks * 20 + 10


def draw_bed_track(
    dr: ImageDraw.ImageDraw,
    features: List[Any],
    y: int,
    width: int,
    start: int,
    end: int,
    bp_per_px: float,
    margin: int,
    stacks: List[int],
) -> int:
    """Draw BED annotation track"""
    header_h = draw_track_header(dr, "BED Annotation", y, width + 2 * margin, 15)
    current_y = y + header_h + 5

    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 10)
    except:
        try:
            font = ImageFont.truetype("arial.ttf", 10)
        except:
            font = ImageFont.load_default()

    for i, feat in enumerate(features):
        stack = stacks[i]
        feat_y = current_y + stack * 20
        mid_y = feat_y + 5

        fx0 = margin + int((max(feat.start, start) - start) / bp_per_px)
        fx1 = margin + int((min(feat.end, end) - start) / bp_per_px)

        if fx1 <= fx0:
            continue

        # Default colors
        default_color = (100, 149, 237)  # Cornflower Blue
        thick_color = (200, 160, 50)      # Brownish Yellow

        # Use itemRgb if available
        if feat.item_rgb:
            color = feat.item_rgb
        else:
            color = default_color

        feature_height = 10

        # Draw main feature rectangle
        dr.rectangle([(fx0, mid_y - feature_height // 2), (fx1, mid_y + feature_height // 2)], fill=color, outline=(0, 0, 0))

        # Draw thickStart/thickEnd if specified (CDS-like region)
        if feat.thick_start is not None and feat.thick_end is not None:
            thick_x0 = margin + int((max(feat.thick_start, start) - start) / bp_per_px)
            thick_x1 = margin + int((min(feat.thick_end, end) - start) / bp_per_px)
            if thick_x1 > thick_x0:
                thick_color_rgb = feat.item_rgb if feat.item_rgb else thick_color
                dr.rectangle([(thick_x0, mid_y - feature_height // 2), (thick_x1, mid_y + feature_height // 2)], fill=thick_color_rgb, outline=(0, 0, 0))

        # Draw blocks (exons) if available
        if feat.blocks:
            for block in feat.blocks:
                bx0 = margin + int((max(block.start, start) - start) / bp_per_px)
                bx1 = margin + int((min(block.end, end) - start) / bp_per_px)
                if bx1 > bx0:
                    dr.rectangle([(bx0, mid_y - feature_height // 2), (bx1, mid_y + feature_height // 2)], fill=color, outline=(0, 0, 0))

            # Draw intron lines between blocks
            if len(feat.blocks) > 1:
                sorted_blocks = sorted(feat.blocks, key=lambda b: b.start)
                for j in range(len(sorted_blocks) - 1):
                    block_end = sorted_blocks[j].end
                    next_block_start = sorted_blocks[j + 1].start
                    if block_end < next_block_start:
                        line_x0 = margin + int((max(block_end, start) - start) / bp_per_px)
                        line_x1 = margin + int((min(next_block_start, end) - start) / bp_per_px)
                        if line_x1 > line_x0:
                            dr.line([(line_x0, mid_y), (line_x1, mid_y)], fill=(0, 0, 0), width=1)
        else:
            # If no blocks, draw intron line for the whole feature
            dr.line([(fx0, mid_y), (fx1, mid_y)], fill=(0, 0, 0), width=1)

        # Strand arrow - solid black arrow outside the feature with connecting line
        head_size = 6
        arrow_width = 4
        connector_length = 3  # Short line connecting feature to arrow
        if feat.strand == '+':
            # Forward strand: arrow pointing right, outside the feature
            arrow_x = fx1 + connector_length
            # Draw connecting line
            dr.line([(fx1, mid_y), (arrow_x, mid_y)], fill=(0, 0, 0), width=1)
            # Draw arrow triangle (pointing right)
            dr.polygon([
                (arrow_x, mid_y - arrow_width),
                (arrow_x + head_size, mid_y),
                (arrow_x, mid_y + arrow_width)
            ], fill=(0, 0, 0), outline=(0, 0, 0))
        elif feat.strand == '-':
            # Reverse strand: arrow pointing left, outside the feature
            arrow_x = fx0 - connector_length
            # Draw connecting line
            dr.line([(fx0, mid_y), (arrow_x, mid_y)], fill=(0, 0, 0), width=1)
            # Draw arrow triangle (pointing left)
            dr.polygon([
                (arrow_x, mid_y - arrow_width),
                (arrow_x - head_size, mid_y),
                (arrow_x, mid_y + arrow_width)
            ], fill=(0, 0, 0), outline=(0, 0, 0))

        # Feature name
        if feat.name:
            dr.text((fx0, mid_y + 6), feat.name, fill=(0, 0, 0), font=font)

    num_stacks = max(stacks) + 1 if stacks else 0
    return header_h + num_stacks * 20 + 10


def _safe_font(size: int = 10):
    """Try Helvetica -> Arial -> default; used across several helpers here."""
    try:
        return ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", size)
    except Exception:
        try:
            return ImageFont.truetype("arial.ttf", size)
        except Exception:
            return ImageFont.load_default()


# ----------------------------------------------------------------------------
# VCF highlight track
# ----------------------------------------------------------------------------
# Layout (top to bottom):
#   [ header bar "VCF Highlights" ]
#   [ SNP marker row ] -- thin row showing where each variant site sits
#   For each sample s:
#       [ "<sname> REF" row ]  -- colored cell per site, colored by REF base
#       [ "<sname> ALT" row ]  -- colored cell per site, colored by ALT base
#                                (empty/gray when the sample is hom-REF)
#
# This matches the IGV "variant track" look: colored tick marks at each SNP,
# followed by per-sample REF/ALT color strips.
#
# Layout notes:
#   - We reserve a left-side gutter (VCF_GUTTER) inside the track for per-row
#     labels, so the colored cells never overlap the text. The gutter sits to
#     the right of the image margin and to the left of the plotting area.
#   - The coverage/reads tracks are unaffected; they still use the full width.

VCF_MARKER_ROW_H = 10     # height of the top "where are the SNPs" row
VCF_BASE_ROW_H = 11       # height of each REF or ALT row (tall enough that
                          # the 9-pt label is clearly legible)
VCF_ROW_GAP = 2           # vertical gap between rows
VCF_GUTTER = 90           # left-side gutter inside the track for row labels


def _sample_ploidy(sites, sname: str) -> int:
    """Return the maximum ploidy seen for this sample across all sites
    in the region. Falls back to 2 if no phased call is found (the default
    diploid REF/ALT layout)."""
    p = 0
    for s in sites:
        call = s.sample_calls.get(sname)
        if call is None:
            continue
        if call.is_phased and call.hap_bases:
            p = max(p, len(call.hap_bases))
    return p if p > 0 else 2


def _vcf_highlight_track_height_rows(rows_per_sample: List[int]) -> int:
    """Total pixel height when each sample contributes a variable number of
    rows (depending on its observed ploidy in the region)."""
    header_h = 15
    total_rows = 1 + sum(rows_per_sample)
    if total_rows <= 0:
        return header_h + 30
    rows_px = (
        VCF_MARKER_ROW_H
        + sum(rows_per_sample) * VCF_BASE_ROW_H
        + (total_rows - 1) * VCF_ROW_GAP
    )
    return header_h + 4 + rows_px + 6


def _vcf_highlight_track_height(num_samples: int) -> int:
    """Legacy fixed-2-rows-per-sample height (kept for the height pre-pass
    if we don't have the actual sites in hand). Defers to the per-sample
    variant when possible."""
    return _vcf_highlight_track_height_rows([2] * num_samples)


def _phase_block_spans(sites: "List[HighlightSite]", sname: str):
    """Return a list of (ps_id, [site_indices...]) for one sample.

    Block-coverage rule (per VCF semantics):
      For each PS value PS_k that appears in this sample's calls, compute
      the FIRST and LAST site index whose call carries PS=PS_k. EVERY site
      in the range [first, last] is then considered part of PS_k -- even if
      that intermediate site is homozygous (and thus has no PS tag of its
      own), or even if it's missing/het-without-PS. The reasoning: PS is
      the phasing tool's declaration "these variants belong to one
      block"; biologically any site lying inside the spatial extent of
      that block is covered by it.

    If two PS spans overlap (which shouldn't happen for a well-formed VCF
    but might if two PS values are interleaved), the LATER PS in the input
    "wins" the overlap; we just pick whichever PS is most recently seen
    when assigning indices. Standalone runs (sites with no PS at all) keep
    going into a None bucket.
    """
    n = len(sites)

    # Step 1: compute first/last index for each PS in this sample
    ps_first: Dict[str, int] = {}
    ps_last: Dict[str, int] = {}
    for i, s in enumerate(sites):
        call = s.sample_calls.get(sname)
        if call is None or not call.is_phased or call.phase_set is None:
            continue
        ps = call.phase_set
        if ps not in ps_first:
            ps_first[ps] = i
        ps_last[ps] = i

    # Step 2: build eff_ps[i] -- the PS that covers site i (or None)
    eff_ps: List[Optional[str]] = [None] * n
    # For determinism + overlapping-block safety, process PSs in the order
    # they first appear so later-starting PSs overwrite earlier ones in
    # any overlap region.
    ordered_ps = sorted(ps_first.keys(), key=lambda k: ps_first[k])
    for ps in ordered_ps:
        for i in range(ps_first[ps], ps_last[ps] + 1):
            eff_ps[i] = ps

    # Step 3: collapse consecutive eff_ps runs into spans
    out = []
    cur_ps = "__SENTINEL__"
    cur_indices = []
    for i in range(n):
        ps = eff_ps[i]
        if ps != cur_ps:
            if cur_indices:
                out.append((cur_ps if cur_ps != "__SENTINEL__" else None, cur_indices))
            cur_ps = ps
            cur_indices = [i]
        else:
            cur_indices.append(i)
    if cur_indices:
        out.append((cur_ps if cur_ps != "__SENTINEL__" else None, cur_indices))
    return out


def draw_vcf_highlight_track(
    dr: ImageDraw.ImageDraw,
    sites: List[HighlightSite],
    samples: List[str],
    y: int,
    width: int,
    start: int,
    end: int,
    bp_per_px: float,
    margin: int,
) -> int:
    """Draw the VCF highlight track and return its total height.

    Per sample we draw 2 rows. For each SNP position we choose what to put
    in those rows based on whether the GT is phased:
      * Phased (GT uses '|' AND PS tag present):
          row 1 (top)    -> hap1 allele color
          row 2 (bottom) -> hap2 allele color
          A semi-transparent dark-gray rectangle is drawn behind the cells
          spanning all SNPs that share the same PS, so you can SEE which
          variants are co-phased into the same block.
      * Unphased (0/0, 0/1 with '/', 1/1, ./. ...):
          row 1 -> REF base color iff sample has_ref
          row 2 -> ALT base color iff sample has alt_base
        i.e. the previous REF/ALT semantics still apply when phasing is
        unavailable. Falls back gracefully on legacy VCFs.
    """
    header_h = draw_track_header(dr, "VCF Highlights", y, width, 15)
    cursor_y = y + header_h + 4

    label_font = _safe_font(9)

    plot_x0 = margin + VCF_GUTTER
    plot_x1 = margin + width

    def site_to_x(pos: int) -> int:
        return margin + int((pos - start) / bp_per_px)

    cell_w = max(2, int(round(1.0 / bp_per_px)) if bp_per_px < 1 else 2)

    # --- SNP marker row (top) ------------------------------------------------
    marker_y0 = cursor_y
    marker_y1 = cursor_y + VCF_MARKER_ROW_H
    dr.rectangle([(margin, marker_y0), (margin + width, marker_y1)],
                 fill=(250, 250, 250))
    for site in sites:
        x = site_to_x(site.pos)
        if x < margin or x >= margin + width:
            continue
        any_ref = any(site.sample_calls.get(s, SampleCall()).has_ref for s in samples)
        alt_base_shown = None
        for sname in samples:
            call = site.sample_calls.get(sname, SampleCall())
            if call.alt_base is not None:
                alt_base_shown = call.alt_base
                break
        mid = (marker_y0 + marker_y1) // 2
        if any_ref:
            ref_color = MISMATCH_COLORS.get(site.ref, (120, 120, 120)) if site.is_snv else (120, 120, 120)
            dr.rectangle([(x, marker_y0), (x + cell_w, mid)], fill=ref_color)
        if alt_base_shown is not None:
            alt_color = MISMATCH_COLORS.get(alt_base_shown, (120, 120, 120)) if site.is_snv else (120, 120, 120)
            dr.rectangle([(x, mid), (x + cell_w, marker_y1)], fill=alt_color)
    dr.rectangle([(margin, marker_y0), (plot_x0, marker_y1)], fill=(240, 240, 240))
    dr.text((margin + 2, marker_y0 - 1), "SNP sites", fill=(70, 70, 70), font=label_font)
    cursor_y = marker_y1 + VCF_ROW_GAP

    # --- Per-sample rows (N rows per sample, N = observed ploidy) -----------
    # Color used for phase-block backgrounds
    PHASE_BG = (216, 216, 216)  # warm gray; clearly visible but not overwhelming
    UNPHASED_BG = (248, 248, 248)

    for sname in samples:
        is_single_synthetic = len(samples) == 1 and sname == "VCF"
        short_name = sname if len(sname) <= 12 else sname[:11] + "…"

        # How many rows does this sample need?
        # When phased calls exist, ploidy = max number of '|'-separated
        # alleles seen. When NOT phased, we keep 2 rows (REF on top, ALT on
        # bottom) — matches v4 unphased semantics.
        ploidy = _sample_ploidy(sites, sname)
        any_phased = any(
            sites[i].sample_calls.get(sname, SampleCall()).is_phased
            for i in range(len(sites))
        )
        n_rows = ploidy if any_phased else 2  # REF / ALT default when unphased

        # Compute row y-coordinates
        row_y0 = [cursor_y + r * (VCF_BASE_ROW_H + VCF_ROW_GAP) for r in range(n_rows)]
        row_y1 = [y0 + VCF_BASE_ROW_H for y0 in row_y0]
        sample_top = row_y0[0]
        sample_bottom = row_y1[-1]

        # Compute phase-block spans for this sample
        blocks = _phase_block_spans(sites, sname)

        # Plain row backgrounds
        for y0, y1 in zip(row_y0, row_y1):
            dr.rectangle([(margin, y0), (margin + width, y1)], fill=UNPHASED_BG)

        # Phase-block backdrops: span all N rows of this sample.
        for ps_id, idx_list in blocks:
            if ps_id is None:
                continue
            first_x = site_to_x(sites[idx_list[0]].pos)
            last_x = site_to_x(sites[idx_list[-1]].pos) + cell_w
            first_x = max(plot_x0, first_x)
            last_x = min(plot_x1, last_x)
            if last_x <= first_x:
                continue
            dr.rectangle([(first_x, sample_top), (last_x, sample_bottom)], fill=PHASE_BG)

        # Set of site indices folded into some named PS block
        sites_in_ps = set()
        for ps_id, idx_list in blocks:
            if ps_id is not None:
                sites_in_ps.update(idx_list)

        # Draw cells
        for site_idx, site in enumerate(sites):
            x = site_to_x(site.pos)
            if x < margin or x >= margin + width:
                continue
            call = site.sample_calls.get(sname, SampleCall())

            if call.is_phased and call.phase_set is not None and call.hap_bases:
                # Directly phased site: paint one cell per hap row.
                # If this site has lower ploidy than the sample-row count,
                # the extra (top) rows for this site are left blank.
                for r_idx, hap_b in enumerate(call.hap_bases):
                    if r_idx >= n_rows or hap_b is None:
                        continue
                    if site.is_snv:
                        c = MISMATCH_COLORS.get(hap_b, (150, 150, 150))
                    else:
                        c = (150, 150, 150)
                    dr.rectangle([(x, row_y0[r_idx]), (x + cell_w, row_y1[r_idx])], fill=c)
            elif site_idx in sites_in_ps and any_phased:
                # Hom site folded into a PS block. Paint the same base on
                # every row (all haps carry it).
                hom_base = None
                if call.has_ref and call.alt_base is None:
                    hom_base = site.ref       # 0/0
                elif (not call.has_ref) and call.alt_base is not None:
                    hom_base = call.alt_base  # 1/1 (or 1/2)
                if hom_base is not None:
                    c = MISMATCH_COLORS.get(hom_base, (150, 150, 150)) if site.is_snv else (150, 150, 150)
                    for r_idx in range(n_rows):
                        dr.rectangle([(x, row_y0[r_idx]), (x + cell_w, row_y1[r_idx])], fill=c)
            else:
                # Unphased site outside any PS block: top row REF, bottom row ALT
                if call.has_ref:
                    c = MISMATCH_COLORS.get(site.ref, (150, 150, 150)) if site.is_snv else (150, 150, 150)
                    dr.rectangle([(x, row_y0[0]), (x + cell_w, row_y1[0])], fill=c)
                if call.alt_base is not None and n_rows >= 2:
                    c = MISMATCH_COLORS.get(call.alt_base, (150, 150, 150)) if site.is_snv else (150, 150, 150)
                    dr.rectangle([(x, row_y0[-1]), (x + cell_w, row_y1[-1])], fill=c)

        # Labels in left gutter (one per row)
        if is_single_synthetic:
            labels = ["REF", "ALT"]
        elif any_phased:
            labels = [f"{short_name} h{r + 1}" for r in range(n_rows)]
        else:
            labels = [f"{short_name} REF", f"{short_name} ALT"]
        for r_idx in range(n_rows):
            dr.rectangle([(margin, row_y0[r_idx]), (plot_x0, row_y1[r_idx])], fill=(240, 240, 240))
            if r_idx < len(labels):
                dr.text((margin + 2, row_y0[r_idx]), labels[r_idx],
                        fill=(70, 70, 70), font=label_font)
        cursor_y = sample_bottom + VCF_ROW_GAP

    # Thin separator line at the bottom of the track for visual anchoring
    dr.line([(margin, cursor_y), (margin + width, cursor_y)], fill=(220, 220, 220), width=1)
    return (cursor_y - y) + 6


def render_jbrowse_style(
    tracks: List[Any],
    chrom: str,
    start: int,
    end: int,
    width: int = 1200,
    read_height: int = 6,
    detail: str = "mid",
    show_axis: bool = True,
    show_coverage: bool = True,
    coverage_height: int = 15,  # Default reduced to 15
    track_title: str = "Reads",
    style: str = "jbrowse",
    color_by: str = "type",
    ref_seq: Optional[str] = None,
    show_insertion_labels: bool = True,
    coverage_max_depth: Optional[int] = None,
    margin: int = 20,  # Left and right margins
    is_rna: bool = False,
    gff_genes: Optional[List[Any]] = None,
    bed_features: Optional[List[Any]] = None,
    highlight_sites: Optional[List[HighlightSite]] = None,
    highlight_samples: Optional[List[str]] = None,
    no_hap_sort: bool = False,
    no_hap_filter: bool = False,
    focus_region: Optional[Tuple[int, int]] = None,
) -> Image.Image:
    """Render JBrowse-style snapshot with track system and coverage chart"""
    # Handle backward compatibility: if tracks is list of Read, wrap it
    if tracks and not isinstance(tracks[0], dict):
        tracks = [{'reads': tracks, 'title': track_title}]

    # Actual drawing area width (excluding left and right margins)
    content_width = width - 2 * margin
    bp_per_px = (end - start) / float(content_width)

    # When a highlight VCF is active, sort reads in each track by their
    # observed VCF-site haplotype signature so reads carrying the same hap
    # cluster together visually. Sort is stable so reads with identical
    # signatures keep their relative order. We do this BEFORE assign_stacks
    # When a highlight VCF is active, two things happen to reads per-track:
    #   (1) Filter: drop reads that don't overlap ANY VCF site in this region
    #       (unless --no-hap-filter). These reads carry no SNP-linkage info
    #       and just consume vertical space, so the default is to hide them.
    #   (2) Sort: reorder by observed hap signature so same-hap reads cluster
    #       (unless --no-hap-sort).
    do_hap_sort = bool(highlight_sites) and not no_hap_sort
    do_hap_filter = bool(highlight_sites) and not no_hap_filter
    if highlight_sites:
        site_positions = sorted(s.pos for s in highlight_sites)
        if site_positions:
            first_p, last_p = site_positions[0], site_positions[-1]
            site_set = set(site_positions)

            def _read_covers_any_site(r):
                # Quick window prune first
                if r.end <= first_p or r.start > last_p:
                    return False
                # Then exact membership check
                for p in site_positions:
                    if r.start <= p < r.end:
                        return True
                return False

            for track in tracks:
                if do_hap_filter:
                    track['reads'] = [r for r in track['reads'] if _read_covers_any_site(r)]
                if do_hap_sort:
                    track['reads'] = sorted(
                        track['reads'],
                        key=lambda r: (_hap_signature(r, site_positions), r.start),
                    )

    # Pre-calculate stacks and row offsets for each track. Row offsets use a
    # larger gap between different reads but preserve a compact gap between
    # supplementary pieces of the same read/qname.
    track_stacks = []
    track_row_offsets = []
    track_reads_area_heights = []
    for track in tracks:
        reads = track['reads']
        if do_hap_sort:
            # 1-read-per-row so the hap-sorted order is preserved on screen
            stacks = list(range(len(reads)))
        else:
            if is_rna:
                spans = [(max(r.start, start), min(r.end, end)) for r in reads]
                keys = [r.qname for r in reads]
                stacks = assign_stacks_grouped(spans, keys)
            else:
                stacks = assign_split_read_stacks(reads, start, end)
        row_offsets, reads_area_h = _stack_row_offsets(stacks, reads, read_height)
        track_stacks.append(stacks)
        track_row_offsets.append(row_offsets)
        track_reads_area_heights.append(reads_area_h)

    # Calculate total height
    top = 0
    if show_axis:
        top += 20  # Axis area (reduced)

    total_height = top

    # Each BAM/track gets its own coverage track. This avoids mixing depth
    # from multiple BAM files into one aggregate profile.
    per_track_coverage_h = 15 + coverage_height + coverage_height + 15 if show_coverage else 0

    # Calculate annotation track height if enabled
    annotation_track_h = 0
    gene_stacks = []
    bed_stacks = []
    if gff_genes:
        gene_spans = [(g.start, g.end) for g in gff_genes]
        gene_stacks = assign_stacks(gene_spans, max_stack=len(gff_genes))
        # Header (15) + Genes area
        num_stacks = max(gene_stacks) + 1 if gene_stacks else 0
        annotation_track_h = 15 + num_stacks * 20 + 10
        total_height += annotation_track_h
    elif bed_features:
        bed_spans = [(f.start, f.end) for f in bed_features]
        bed_stacks = assign_bed_stacks(bed_spans, min_distance=10, max_stack=len(bed_features))
        # Header (15) + Features area
        num_stacks = max(bed_stacks) + 1 if bed_stacks else 0
        annotation_track_h = 15 + num_stacks * 20 + 10
        total_height += annotation_track_h

    # Calculate VCF highlight track height if enabled.
    # Each sample takes max(ploidy_seen, 2) rows; samples without any phased
    # calls fall back to 2 (REF/ALT semantics).
    vcf_track_h = 0
    if highlight_sites and highlight_samples:
        rows_per_sample = [_sample_ploidy(highlight_sites, s) for s in highlight_samples]
        vcf_track_h = _vcf_highlight_track_height_rows(rows_per_sample)
        total_height += vcf_track_h

    # Calculate height for each track
    for i, track in enumerate(tracks):
        track_header_height = 15
        reads_area_height = track_reads_area_heights[i]
        total_height += per_track_coverage_h + track_header_height + reads_area_height + 5

    img = Image.new("RGB", (width, total_height), (255, 255, 255))
    dr = ImageDraw.Draw(img)

    current_y = 0

    # Collect tick positions for drawing guide lines
    tick_positions = []

    # Draw coordinate axis
    if show_axis:
        # Axis line position adjusted to bottom
        axis_y = current_y + 18
        dr.line([(margin, axis_y), (width - margin - 1, axis_y)], fill=(0, 0, 0), width=1)
        span = end - start
        step = max(1, span // 10)
        try:
            font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 9)
        except:
            try:
                font = ImageFont.truetype("arial.ttf", 9)
            except:
                font = ImageFont.load_default()
        for pos in range(start, end + 1, step):
            x = margin + int((pos - start) / bp_per_px)
            tick_positions.append(x)
            # Tick mark
            dr.line([(x, axis_y - 3), (x, axis_y + 3)], fill=(0, 0, 0), width=1)
            # Text moved up, separated from axis
            pos_str = f"{pos:,}"
            dr.text((x + 2, current_y + 2), pos_str, fill=(0, 0, 0), font=font)
        current_y += 22
    
    focus_px = _focus_region_to_pixels(
        focus_region,
        start,
        end,
        bp_per_px,
        margin,
        width,
    )

    if focus_px is not None:
        fx0, fx1 = focus_px

        # Very light background for the original target region
        shade_top = current_y
        shade_bottom = total_height

        # very subtle blue-gray
        focus_bg = (235, 242, 255)

        dr.rectangle(
            [(fx0, shade_top), (fx1, shade_bottom)],
            fill=focus_bg
        )
    
    
    # Draw vertical guide lines (dashed) - from axis to bottom of image
    for tick_x in tick_positions:
        # Draw dashed line: draw 2px line every 4 pixels
        for dash_y in range(current_y, total_height, 6):
            dash_end = min(dash_y + 3, total_height)
            dr.line([(tick_x, dash_y), (tick_x, dash_end)], fill=(220, 220, 220), width=1)

    # Per-BAM coverage tracks are drawn immediately above each read track.

    # Draw annotation track if enabled
    if gff_genes:
        current_y += draw_gene_track(
            dr,
            gff_genes,
            current_y,
            content_width,
            start,
            end,
            bp_per_px,
            margin,
            gene_stacks
        )
    elif bed_features:
        current_y += draw_bed_track(
            dr,
            bed_features,
            current_y,
            content_width,
            start,
            end,
            bp_per_px,
            margin,
            bed_stacks
        )

    # Draw VCF highlight track if enabled. Sits right above the reads so the
    # user can visually line up the colored REF/ALT cells with each read.
    if highlight_sites and highlight_samples:
        current_y += draw_vcf_highlight_track(
            dr,
            highlight_sites,
            highlight_samples,
            current_y,
            content_width,
            start,
            end,
            bp_per_px,
            margin,
        )

    # Build a pos -> HighlightSite index once; used inside the reads loop to
    # decide whether a mismatch base should be colored or rendered in pale
    # gray. Empty dict means "no highlight active" which falls back to the
    # original full-color behavior.
    highlight_index = {}
    if highlight_sites:
        highlight_index = {s.pos: s for s in highlight_sites}

    # Iterate over tracks
    for i, track in enumerate(tracks):
        reads = track['reads']
        current_track_title = track.get('title', track_title)
        stacks = track_stacks[i]
        row_offsets = track_row_offsets[i]
        reads_area_height = track_reads_area_heights[i]

        if show_coverage:
            current_y, img, dr = _draw_per_track_coverage_png(
                img, dr, reads, current_track_title, current_y, width,
                content_width, coverage_height, start, end, ref_seq,
                coverage_max_depth, margin, is_rna, detail, bp_per_px
            )

        # Read track header
        current_y += draw_track_header(dr, current_track_title, current_y, width, reads_area_height)

        # Draw reads
        # (Rest of the read drawing logic remains same)
        reads_start_y = current_y
        groups: Dict[str, List[int]] = {}
        for idx_r, r in enumerate(reads):
            groups.setdefault(r.qname, []).append(idx_r)
        supplementary_read_styles = _build_supplementary_read_styles(reads, groups, stacks)
        dr_alpha = ImageDraw.Draw(img, "RGBA")

        for idx, r in enumerate(reads):
            y = reads_start_y + row_offsets[stacks[idx]]
            rects = segments_to_pixels(r.segments, r.start, start, bp_per_px, detail=detail)
            supp_style = supplementary_read_styles.get(idx)

            # Create rect index to segment mapping, for getting insertion length
            rect_to_seg = {}
            ref_cursor = r.start
            rect_idx = 0
            for seg_idx, seg in enumerate(r.segments):
                if seg.ref_consumed == 0:
                    if seg.type == "ins":
                        x = int((ref_cursor - start) / bp_per_px)
                        if rect_idx < len(rects) and rects[rect_idx][0] == "ins":
                            rect_to_seg[rect_idx] = seg_idx
                            rect_idx += 1
                    continue
                if rect_idx < len(rects):
                    rect_to_seg[rect_idx] = seg_idx
                    rect_idx += 1
                ref_cursor += seg.ref_consumed

            for rect_idx, (t, x0, x1) in enumerate(rects):
                # Add margin offset
                x0_draw = margin + x0
                x1_draw = margin + x1

                # For insertion, check if original position is within visible area
                if t == "ins":
                    if x0 < 0 or x0 > (width - 2 * margin):
                        continue

                # Limit to visible area
                x0_draw = max(margin, min(width - margin, x0_draw))
                x1_draw = max(margin, min(width - margin, x1_draw))

                if x1_draw <= x0_draw and t != "ins":
                    continue

                if color_by == "type":
                    color = color_for_type(t)
                elif color_by == "mapq":
                    base_color = color_for_type(t)
                    color = shade_by_mapq(base_color, r.mapq)
                else:
                    color = color_for_type(t)

                if t == "ins":
                    # In highlight mode the purple ins tick and the 'I(..)' label
                    # are both distracting — the user is hunting SNP linkage on
                    # the reads, and long-read data has ins noise everywhere.
                    # Skip drawing ins entirely when highlight is active so
                    # only VCF sites show up colored.
                    if highlight_index:
                        pass  # suppressed
                    else:
                        ins_color = (128, 0, 128)  # Purple
                        ins_line_width = max(2, read_height // 3)
                        dr.line([(x0_draw, y), (x0_draw, y + read_height)], fill=ins_color, width=ins_line_width)
                        if detail == "high" and not is_rna and show_insertion_labels and rect_idx in rect_to_seg:
                            seg = r.segments[rect_to_seg[rect_idx]]
                            if seg and seg.length > 0:
                                try:
                                    font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 8)
                                except:
                                    try:
                                        font = ImageFont.truetype("arial.ttf", 8)
                                    except:
                                        font = ImageFont.load_default()
                                label = f"I({seg.length})"
                                dr.text((x0_draw + 2, y - 12), label, fill=ins_color, font=font)
                elif t == "ref_skip":
                    intron_color = (176, 196, 222)
                    y_center = y + read_height // 2
                    dr.line([(x0_draw, y_center), (x1_draw, y_center)], fill=intron_color, width=1)
                elif t == "del":
                    # In highlight mode: draw del in the SAME pale muted gray
                    # used for non-highlighted mismatches — indistinguishable
                    # from a plain match, so it doesn't draw the eye.
                    if highlight_index:
                        del_color = MUTED_MISMATCH_COLOR
                    else:
                        del_color = (115, 115, 115)
                    dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=del_color)
                    # Skip the del-length label too when highlight is active
                    if not highlight_index and detail == "high" and not is_rna and show_insertion_labels and rect_idx in rect_to_seg:
                        seg = r.segments[rect_to_seg[rect_idx]]
                        if seg and seg.length > 0:
                            try:
                                font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 8)
                            except:
                                try:
                                    font = ImageFont.truetype("arial.ttf", 8)
                                except:
                                    font = ImageFont.load_default()
                            label = str(seg.length)
                            del_width = x1_draw - x0_draw
                            if del_width > 4:
                                dr.text(((x0_draw + x1_draw) // 2 - 5, y + 1), label, fill=(255, 255, 255), font=font)
                            else:
                                dr.text((x0_draw, y - 10), label, fill=(0, 0, 0), font=font)
                elif t == "mismatch":
                    # Default path: no highlight VCF active -> original behavior
                    # (color mismatch by the read base).
                    #
                    # Highlight path: if a highlight VCF is active, we instead
                    # walk the mismatch segment base-by-base. At each ref
                    # position that corresponds to a VCF site, we color with
                    # the read base (this is the site the user wants to see).
                    # Positions NOT in the VCF get painted in MUTED gray so
                    # they recede visually.
                    #
                    # Note: a single "mismatch" rect may represent several
                    # consecutive X ops that were merged; we therefore cannot
                    # just use one color for the whole rect when highlight is
                    # on.
                    if detail == "low":
                        fill_color, fill_alpha = supp_style if supp_style else ((195, 195, 195), 255)
                        _draw_rect_with_alpha(dr, dr_alpha, [(x0_draw, y), (x1_draw, y + read_height)], fill_color, fill_alpha)
                    elif not highlight_index:
                        # Original (un-highlighted) path.
                        if r.seq and rect_idx in rect_to_seg:
                            seg_idx = rect_to_seg[rect_idx]
                            read_cursor = sum(s.read_consumed for s in r.segments[:seg_idx])
                            if read_cursor < len(r.seq):
                                base = r.seq[read_cursor].upper()
                                mismatch_color = MISMATCH_COLORS.get(base, (200, 60, 60))
                                dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=mismatch_color)
                            else:
                                fill_color, fill_alpha = supp_style if supp_style else ((195, 195, 195), 255)
                                _draw_rect_with_alpha(dr, dr_alpha, [(x0_draw, y), (x1_draw, y + read_height)], fill_color, fill_alpha)
                        else:
                            fill_color, fill_alpha = supp_style if supp_style else ((195, 195, 195), 255)
                            _draw_rect_with_alpha(dr, dr_alpha, [(x0_draw, y), (x1_draw, y + read_height)], fill_color, fill_alpha)
                    else:
                        # Highlight-aware mismatch rendering.
                        # First paint the whole rect muted gray, then overdraw
                        # per-base cells for the ref positions that are VCF
                        # sites.
                        fill_color, fill_alpha = supp_style if supp_style else (MUTED_MISMATCH_COLOR, 255)
                        _draw_rect_with_alpha(dr, dr_alpha, [(x0_draw, y), (x1_draw, y + read_height)], fill_color, fill_alpha)
                        if r.seq and rect_idx in rect_to_seg:
                            seg_idx = rect_to_seg[rect_idx]
                            seg = r.segments[seg_idx]
                            # ref position of the 1st base covered by this seg
                            ref_pos0 = r.start + sum(
                                s.ref_consumed for s in r.segments[:seg_idx]
                            )
                            # read index of the 1st base covered by this seg
                            read_pos0 = sum(s.read_consumed for s in r.segments[:seg_idx])
                            # For mismatch segments, ref_consumed == read_consumed == length.
                            for k in range(seg.length):
                                ref_p = ref_pos0 + k
                                if ref_p not in highlight_index:
                                    continue
                                # This base is at a VCF-highlighted site.
                                rk = read_pos0 + k
                                if rk >= len(r.seq):
                                    continue
                                base = r.seq[rk].upper()
                                color = MISMATCH_COLORS.get(base, (200, 60, 60))
                                bx0 = margin + int((ref_p - start) / bp_per_px)
                                bx1 = margin + int((ref_p + 1 - start) / bp_per_px)
                                bx1 = max(bx1, bx0 + 1)
                                bx0 = max(margin, min(width - margin, bx0))
                                bx1 = max(margin, min(width - margin, bx1))
                                dr.rectangle([(bx0, y), (bx1, y + read_height)], fill=color)
                else:
                    if t == "match":
                        # Highlight mode: at VCF positions inside this match
                        # segment, the read agrees with REF — so paint that
                        # one base with the REF color (consistent with the
                        # VCF Highlights track above). Everything else stays
                        # the regular match gray. This lets every read be
                        # read as a "colored barcode" across all VCF sites,
                        # regardless of match/mismatch.
                        fill_color, fill_alpha = supp_style if supp_style else ((195, 195, 195), 255)
                        _draw_rect_with_alpha(dr, dr_alpha, [(x0_draw, y), (x1_draw, y + read_height)], fill_color, fill_alpha)
                        if highlight_index and r.seq and rect_idx in rect_to_seg:
                            seg_idx = rect_to_seg[rect_idx]
                            seg = r.segments[seg_idx]
                            ref_pos0 = r.start + sum(
                                s.ref_consumed for s in r.segments[:seg_idx]
                            )
                            read_pos0 = sum(s.read_consumed for s in r.segments[:seg_idx])
                            for k in range(seg.length):
                                ref_p = ref_pos0 + k
                                if ref_p not in highlight_index:
                                    continue
                                rk = read_pos0 + k
                                if rk >= len(r.seq):
                                    continue
                                base = r.seq[rk].upper()
                                color = MISMATCH_COLORS.get(base, (150, 150, 150))
                                bx0 = margin + int((ref_p - start) / bp_per_px)
                                bx1 = margin + int((ref_p + 1 - start) / bp_per_px)
                                bx1 = max(bx1, bx0 + 1)
                                bx0 = max(margin, min(width - margin, bx0))
                                bx1 = max(margin, min(width - margin, bx1))
                                dr.rectangle([(bx0, y), (bx1, y + read_height)], fill=color)
                    else:
                        if supp_style and t in ("match", "soft", "hard"):
                            fill_color, fill_alpha = supp_style
                            _draw_rect_with_alpha(dr, dr_alpha, [(x0_draw, y), (x1_draw, y + read_height)], fill_color, fill_alpha)
                        else:
                            dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=color)

            # Read direction arrow
            if style == "jbrowse" and rects:
                read_end_px = margin + int((min(r.end, end) - start) / bp_per_px)
                read_start_px = margin + int((max(r.start, start) - start) / bp_per_px)
                read_end_px = max(margin, min(width - margin - 1, read_end_px))
                read_start_px = max(margin, min(width - margin - 1, read_start_px))
                read_w = max(1, read_end_px - read_start_px)
                head = max(4, min(max(read_height + 1, 5), 8, max(1, read_w - 1)))
                arrow_fill, arrow_alpha = ((255, 255, 255), 255) if supp_style else ((211, 211, 211), 255)
                arrow_truncated = (200, 200, 200)

                if r.reverse:
                    arrow_x = read_start_px
                    is_truncated = r.start < start
                    points = [(arrow_x, y), (arrow_x - head, y + read_height // 2), (arrow_x, y + read_height)]
                    if is_truncated:
                        dr.polygon(points, outline=arrow_truncated, fill=(255, 255, 255))
                    else:
                        _draw_polygon_with_alpha(dr, dr_alpha, points, arrow_fill, arrow_alpha)
                else:
                    arrow_x = read_end_px
                    is_truncated = r.end > end
                    points = [(arrow_x, y), (arrow_x + head, y + read_height // 2), (arrow_x, y + read_height)]
                    if is_truncated:
                        dr.polygon(points, outline=arrow_truncated, fill=(255, 255, 255))
                    else:
                        _draw_polygon_with_alpha(dr, dr_alpha, points, arrow_fill, arrow_alpha)

        # Connection lines
        if is_rna:
            for qname, idxs in groups.items():
                if len(idxs) > 1:
                    idxs_sorted = sorted(idxs, key=lambda i: reads[i].start)
                    for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                        if reads[a].end < reads[b].start:
                            ya_center = reads_start_y + row_offsets[stacks[a]] + read_height // 2
                            yb_center = reads_start_y + row_offsets[stacks[b]] + read_height // 2
                            xa = margin + int((reads[a].end - start) / bp_per_px)
                            xb = margin + int((reads[b].start - start) / bp_per_px)
                            xa = max(margin, min(width - margin - 1, xa))
                            xb = max(margin, min(width - margin - 1, xb))
                            dr.line([(xa, ya_center), (xb, yb_center)], fill=(176, 196, 222), width=1)
        else:
            # DNA supplementary/split-read connector lines are intentionally
            # disabled.  The relationship is now encoded by qname color plus
            # query-order opacity, so no extra linking line is drawn.
            pass

        current_y = reads_start_y + reads_area_height + 5

    return img
