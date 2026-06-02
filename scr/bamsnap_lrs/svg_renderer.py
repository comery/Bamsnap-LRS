'''
作用：负责将同样的信息渲染为SVG格式的矢量图。
SVG支持缩放不失真，适合网页交互、精细展示和后期编辑。
支持丰富的可视化细节（如基因注释、BED注释、reads分层、变异高亮、坐标轴、箭头等），并能灵活调整样式和内容。
'''
"""SVG renderer for bamsnap visualization"""
from typing import List, Dict, Optional, Any, Tuple
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom

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

# Matches PNG renderer: mid-gray for mismatches outside the VCF highlight,
# darker than the match color (#c3c3c3) so the mismatch shape stays visible.
# Pale gray for non-VCF mismatches; lighter than match color (#c3c3c3) so
# the mismatch shape stays visible but doesn't visually compete.
MUTED_MISMATCH_HEX = "#c8c8c8"


# Supplementary/split-read fill colors.  SVG uses fill-opacity for the same
# query-order gradient as png_renderer.py.
SUPPLEMENTARY_READ_COLORS = [
    (80, 170, 120),
    (70, 150, 210),
    (240, 170, 70),
    (165, 120, 205),
    (220, 110, 110),
    (60, 190, 180),
    (190, 160, 80),
    (120, 170, 210),
    (150, 190, 95),
    (205, 130, 170),
]
# Kept as a legacy constant, but the current opacity gradient is calculated
# from the number of visible pieces in each qname group.
SUPPLEMENTARY_OPACITY_STEP = 0.20
SUPPLEMENTARY_MAX_TRANSPARENCY = 0.80


def _is_supplementary_read_group(reads: List[Read], idxs: List[int]) -> bool:
    """Return True only when a qname has multiple visible alignments.

    A single visible alignment marked supplementary is not styled as a
    supplementary/split-read group, because the other pieces of that read are
    outside the current plotting window and there is no within-window group to
    display.
    """
    return len(idxs) > 1


def _build_supplementary_read_styles(
    reads: List[Read],
    groups: Dict[str, List[int]],
    stacks: Optional[List[int]] = None,
) -> Dict[int, Tuple[Tuple[int, int, int], float]]:
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
    styles: Dict[int, Tuple[Tuple[int, int, int], float]] = {}
    placed: List[Tuple[float, int, int, int]] = []

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

        # Distribute transparency evenly across the visible pieces in this
        # qname group: first piece is fully opaque, last piece reaches
        # SUPPLEMENTARY_MAX_TRANSPARENCY, and middle pieces are interpolated.
        # Example with max=0.80:
        #   2 pieces -> 0.00, 0.40
        #   3 pieces -> 0.00, 0.40, 0.80
        #   5 pieces -> 0.00, 0.20, 0.40, 0.60, 0.80
        if n_pieces <= 1:
            transparency_values = [0.0]
        elif n_pieces==2:
            transparency_values = [0.0,0.4]
        else:
            transparency_values = [
                SUPPLEMENTARY_MAX_TRANSPARENCY * piece_order / (n_pieces - 1)
                for piece_order in range(n_pieces)
            ]

        min_opacity = max(0.0, 1.0 - SUPPLEMENTARY_MAX_TRANSPARENCY)
        for idx, transparency in zip(idxs_sorted, transparency_values):
            transparency = min(SUPPLEMENTARY_MAX_TRANSPARENCY, max(0.0, transparency))
            opacity = max(min_opacity, 1.0 - transparency)
            styles[idx] = (color, opacity)
        placed.append((center, row_min, row_max, color_idx))

    return styles


def _svg_rect_attrs_with_optional_opacity(
    x: int,
    y: int,
    width: int,
    height: int,
    fill: str,
    opacity: Optional[float] = None,
) -> Dict[str, str]:
    attrs = {
        "x": str(x), "y": str(y),
        "width": str(width), "height": str(height),
        "fill": fill, "stroke": "none",
    }
    if opacity is not None and opacity < 0.999:
        attrs["fill-opacity"] = f"{opacity:.3f}"
    return attrs


def rgb_to_hex(rgb: tuple) -> str:
    """Convert RGB tuple to hex color"""
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
    
def _read_base_at_pos(read: Read, pos: int) -> Optional[str]:
    """Return the read base aligned to reference position `pos` (0-based).

    This local copy removes svg_renderer.py's previous dependency on
    png_renderer.py for haplotype sorting.
    """
    if not read.seq:
        return None
    if pos < read.start or pos >= read.end:
        return None
    ref_cur = read.start
    read_cur = 0
    for s in read.segments:
        if s.ref_consumed > 0 and s.read_consumed > 0:
            if ref_cur <= pos < ref_cur + s.ref_consumed:
                idx = read_cur + (pos - ref_cur)
                if 0 <= idx < len(read.seq):
                    return read.seq[idx].upper()
                return None
        elif s.ref_consumed > 0 and s.read_consumed == 0:
            if ref_cur <= pos < ref_cur + s.ref_consumed:
                return None
        ref_cur += s.ref_consumed
        read_cur += s.read_consumed
    return None


def _hap_signature(read: Read, site_positions: List[int]) -> str:
    """Build a per-read haplotype signature over highlighted VCF sites."""
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

    if x1 <= x0:
        x1 = min(width - margin, x0 + 2)

    return x0, x1

def draw_svg_track_header(svg, title, y, width):
    """Draw track header in SVG"""
    header_h = 15
    SubElement(svg, "rect", {
        "x": "0", "y": str(y),
        "width": str(width), "height": str(header_h),
        "fill": "#f5f5f5"
    })
    SubElement(svg, "line", {
        "x1": "0", "y1": str(y + header_h - 1),
        "x2": str(width), "y2": str(y + header_h - 1),
        "stroke": "#c8c8c8", "stroke-width": "1"
    })
    SubElement(svg, "text", {
        "x": "5", "y": str(y + 11),
        "font-size": "11", "fill": "black"
    }).text = title
    return header_h


# Same vertical spacing policy as png_renderer.py.
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

def _svg_base_distribution(reads: List[Read], start: int, end: int, ref_seq: Optional[str]):
    pile = base_pileup(reads, start, end)
    base_distribution = []
    for i_base, p in enumerate(pile):
        ref_base = None
        if ref_seq and i_base < len(ref_seq):
            ref_base = ref_seq[i_base].upper()

        ref_match = 0
        variant_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
        for base in ["A", "C", "G", "T", "N"]:
            count = p.get(base, 0)
            if count > 0:
                if ref_base and base == ref_base:
                    ref_match += count
                else:
                    variant_counts[base] += count

        base_distribution.append({
            "ref_match": ref_match,
            "A": variant_counts["A"],
            "C": variant_counts["C"],
            "G": variant_counts["G"],
            "T": variant_counts["T"],
            "N": variant_counts["N"],
            "depth": p.get("depth", 0),
        })
    return base_distribution


def draw_svg_base_composition_track(
    svg,
    reads: List[Read],
    y: int,
    width: int,
    content_width: int,
    start: int,
    end: int,
    composition_height: Optional[int],
    comp_max_depth: Optional[int],
    margin: int,
    read_height: int,
) -> int:
    """Draw the old PNG-style base composition chart in SVG.

    This is intentionally separate from the coverage track: coverage shows
    REF-match versus variant bases, while composition simply stacks A/C/G/T/N
    counts for all reads across the displayed window.
    """
    if not reads:
        return 0

    comp_h = composition_height if composition_height is not None else (read_height * 2)
    comp_h = max(1, int(comp_h))
    top_padding = 5
    bottom_padding = 5
    chart_top = y + top_padding
    chart_bottom = chart_top + comp_h

    pile = base_pileup(reads, start, end)
    bins = pileup_to_pixels(pile, content_width)

    max_depth = comp_max_depth if comp_max_depth is not None else max((b.get("depth", 0) for b in bins), default=1)
    if max_depth <= 0:
        max_depth = 1

    base_colors = {
        "A": "#50a050",
        "C": "#5078c8",
        "G": "#f0c846",
        "T": "#dc5a5a",
        "N": "#a0a0a0",
    }

    # A very light baseline/frame keeps the composition track readable without
    # making it look like a separate heavy coverage track.
    SubElement(svg, "line", {
        "x1": str(margin), "y1": str(chart_bottom),
        "x2": str(margin + content_width), "y2": str(chart_bottom),
        "stroke": "#d0d0d0", "stroke-width": "1",
    })

    for x, agg in enumerate(bins):
        depth = agg.get("depth", 0)
        if depth <= 0:
            continue

        bar_height = int(comp_h * min(depth / max_depth, 1.0))
        if bar_height <= 0:
            continue

        # Scale each base to the visible bar height.  The remainder is assigned
        # to N to preserve the same stacked-height behavior as the old PNG path.
        a_h = int(bar_height * agg.get("A", 0) / depth)
        c_h = int(bar_height * agg.get("C", 0) / depth)
        g_h = int(bar_height * agg.get("G", 0) / depth)
        t_h = int(bar_height * agg.get("T", 0) / depth)
        used = a_h + c_h + g_h + t_h
        n_h = max(0, bar_height - used)

        draw_x = margin + x
        cursor_y = chart_bottom
        for base, h in (("A", a_h), ("C", c_h), ("G", g_h), ("T", t_h), ("N", n_h)):
            if h <= 0:
                continue
            SubElement(svg, "rect", {
                "x": str(draw_x),
                "y": str(cursor_y - h),
                "width": "1",
                "height": str(h),
                "fill": base_colors[base],
                "stroke": "none",
            })
            cursor_y -= h

    return top_padding + comp_h + bottom_padding


def draw_svg_per_track_coverage(
    svg,
    reads: List[Read],
    title: str,
    y: int,
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
) -> int:
    """Draw one SVG coverage track for one BAM/track."""
    track_y_start = y
    current_y = y + draw_svg_track_header(svg, f"Coverage: {title}" if title else "Coverage", y, width)

    # Reserve space for RNA arcs above the coverage bars.
    current_y += coverage_height

    base_distribution = _svg_base_distribution(reads, start, end, ref_seq)
    num_bases = len(base_distribution)
    max_cov = coverage_max_depth
    if max_cov is None:
        max_cov = max((b["depth"] for b in base_distribution), default=1)
    if max_cov <= 0:
        max_cov = 1

    bar_bottom = current_y + coverage_height
    bar_top = current_y
    axis_x = margin - 2

    SubElement(svg, "line", {
        "x1": str(axis_x), "y1": str(bar_top),
        "x2": str(axis_x), "y2": str(bar_bottom),
        "stroke": "#646464", "stroke-width": "1"
    })
    for tick_y, tick_label in [(bar_top, str(max_cov)), ((bar_top + bar_bottom) // 2, str(max_cov // 2)), (bar_bottom, "0")]:
        SubElement(svg, "line", {
            "x1": str(axis_x - 3), "y1": str(tick_y),
            "x2": str(axis_x), "y2": str(tick_y),
            "stroke": "#646464", "stroke-width": "1"
        })
        SubElement(svg, "text", {
            "x": "2", "y": str(tick_y + (6 if tick_y == bar_top else 3)),
            "font-size": "8", "fill": "#505050"
        }).text = tick_label

    import math
    if num_bases > 0:
        for px in range(content_width):
            base_start_idx = int(px * num_bases / content_width)
            base_end_idx = int((px + 1) * num_bases / content_width)
            if base_end_idx <= base_start_idx:
                base_end_idx = base_start_idx + 1
            if base_end_idx > num_bases:
                base_end_idx = num_bases

            agg_ref_match = 0
            agg_variants = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
            agg_depth = 0
            for base_idx in range(base_start_idx, base_end_idx):
                dist = base_distribution[base_idx]
                agg_ref_match += dist.get("ref_match", 0)
                for base in ["A", "C", "G", "T", "N"]:
                    agg_variants[base] += dist.get(base, 0)
                agg_depth += dist.get("depth", 0)

            num_positions = base_end_idx - base_start_idx
            if num_positions > 1:
                agg_ref_match = round(agg_ref_match / num_positions)
                for base in agg_variants:
                    if agg_variants[base] > 0:
                        agg_variants[base] = max(1, math.ceil(agg_variants[base] / num_positions))
                agg_depth = round(agg_depth / num_positions)
            if agg_depth == 0:
                continue

            draw_x = margin + px
            bar_height = int(coverage_height * min(agg_depth / max_cov, 1.0))
            if bar_height <= 0:
                continue

            total_count = agg_ref_match + sum(agg_variants.values())
            if total_count == 0:
                continue

            if detail == "low":
                SubElement(svg, "rect", {
                    "x": str(draw_x), "y": str(bar_bottom - bar_height),
                    "width": "1", "height": str(bar_height), "fill": "#b4b4b4"
                })
                continue

            base_heights = {}
            if agg_ref_match > 0:
                base_heights["ref"] = int(bar_height * agg_ref_match / total_count)
            for base in ["A", "C", "G", "T", "N"]:
                count = agg_variants[base]
                if count > 0:
                    h = int(bar_height * count / total_count)
                    base_heights[base] = max(1, h) if h == 0 else h

            total_height_used = sum(base_heights.values())
            if total_height_used < bar_height and base_heights:
                max_base = max(base_heights.items(), key=lambda x: x[1])[0]
                base_heights[max_base] += (bar_height - total_height_used)
            elif total_height_used > bar_height:
                excess = total_height_used - bar_height
                if "ref" in base_heights and base_heights["ref"] > excess:
                    base_heights["ref"] -= excess

            current_stack_y = bar_bottom
            for base in ["A", "C", "G", "T", "N"]:
                h = base_heights.get(base, 0)
                if h <= 0:
                    continue
                SubElement(svg, "rect", {
                    "x": str(draw_x), "y": str(current_stack_y - h),
                    "width": "1", "height": str(h),
                    "fill": rgb_to_hex(MISMATCH_COLORS.get(base, (160, 160, 160)))
                })
                current_stack_y -= h
            h = base_heights.get("ref", 0)
            if h > 0:
                SubElement(svg, "rect", {
                    "x": str(draw_x), "y": str(current_stack_y - h),
                    "width": "1", "height": str(h), "fill": "#b4b4b4"
                })

    if is_rna:
        arc_color_hex = "#fdd1d3"
        arc_anchor_y = track_y_start + 15 + coverage_height + coverage_height // 2
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
                        mid_x = (xa + xb) / 2
                        ctrl_y = arc_anchor_y - 2 * h
                        SubElement(svg, "path", {
                            "d": f"M {xa} {arc_anchor_y} Q {mid_x} {ctrl_y} {xb} {arc_anchor_y}",
                            "fill": "none", "stroke": arc_color_hex,
                            "stroke-opacity": "0.6", "stroke-width": "0.5"
                        })
                ref_cursor += seg.ref_consumed

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
                            mid_x = (xa + xb) / 2
                            ctrl_y = arc_anchor_y - 2 * h
                            SubElement(svg, "path", {
                                "d": f"M {xa} {arc_anchor_y} Q {mid_x} {ctrl_y} {xb} {arc_anchor_y}",
                                "fill": "none", "stroke": arc_color_hex,
                                "stroke-opacity": "0.6", "stroke-width": "1"
                            })

    return track_y_start + 15 + coverage_height + coverage_height + 15


def draw_svg_gene_track(svg, genes, y, width, start, end, bp_per_px, margin, stacks):
    """Draw gene track in SVG"""
    header_h = draw_svg_track_header(svg, "Gene Annotation", y, width + 2 * margin)
    current_y = y + header_h + 5

    for i, gene in enumerate(genes):
        stack = stacks[i]
        gene_y = current_y + stack * 20
        mid_y = gene_y + 5

        gx0 = margin + int((max(gene.start, start) - start) / bp_per_px)
        gx1 = margin + int((min(gene.end, end) - start) / bp_per_px)

        if gx1 <= gx0:
            continue

        # Colors for gene features
        color_utr = "#6495ed"  # Cornflower Blue
        color_cds = "#c8a032"  # Brownish Yellow
        feature_height = 10

        # Intron line
        SubElement(svg, "line", {
            "x1": str(gx0), "y1": str(mid_y),
            "x2": str(gx1), "y2": str(mid_y),
            "stroke": "black", "stroke-width": "1"
        })

        # Strand arrow - solid black arrow outside the feature with connecting line
        head_size = 6
        arrow_width = 4
        connector_length = 3  # Short line connecting feature to arrow
        if gene.strand == '+':
            # Forward strand: arrow pointing right, outside the feature
            arrow_x = gx1 + connector_length
            # Draw connecting line
            SubElement(svg, "line", {
                "x1": str(gx1), "y1": str(mid_y),
                "x2": str(arrow_x), "y2": str(mid_y),
                "stroke": "black", "stroke-width": "1"
            })
            # Draw arrow triangle (pointing right)
            SubElement(svg, "polygon", {
                "points": f"{arrow_x},{mid_y-arrow_width} {arrow_x+head_size},{mid_y} {arrow_x},{mid_y+arrow_width}",
                "fill": "black", "stroke": "black", "stroke-width": "0.5"
            })
        elif gene.strand == '-':
            # Reverse strand: arrow pointing left, outside the feature
            arrow_x = gx0 - connector_length
            # Draw connecting line
            SubElement(svg, "line", {
                "x1": str(gx0), "y1": str(mid_y),
                "x2": str(arrow_x), "y2": str(mid_y),
                "stroke": "black", "stroke-width": "1"
            })
            # Draw arrow triangle (pointing left)
            SubElement(svg, "polygon", {
                "points": f"{arrow_x},{mid_y-arrow_width} {arrow_x-head_size},{mid_y} {arrow_x},{mid_y+arrow_width}",
                "fill": "black", "stroke": "black", "stroke-width": "0.5"
            })

        # Exons (UTR color)
        for exon in gene.exons:
            ex0 = margin + int((max(exon.start, start) - start) / bp_per_px)
            ex1 = margin + int((min(exon.end, end) - start) / bp_per_px)
            if ex1 > ex0:
                SubElement(svg, "rect", {
                    "x": str(ex0), "y": str(mid_y - feature_height / 2),
                    "width": str(ex1 - ex0), "height": str(feature_height),
                    "fill": color_utr
                })

        # CDS (CDS color, same height)
        for cds in gene.cds:
            cx0 = margin + int((max(cds.start, start) - start) / bp_per_px)
            cx1 = margin + int((min(cds.end, end) - start) / bp_per_px)
            if cx1 > cx0:
                SubElement(svg, "rect", {
                    "x": str(cx0), "y": str(mid_y - feature_height / 2),
                    "width": str(cx1 - cx0), "height": str(feature_height),
                    "fill": color_cds
                })

        # Gene name
        SubElement(svg, "text", {
            "x": str(gx0), "y": str(mid_y + 16),
            "font-size": "10", "fill": "black"
        }).text = gene.name

    num_stacks = max(stacks) + 1 if stacks else 0
    return header_h + num_stacks * 20 + 10


def draw_svg_bed_track(svg, features, y, width, start, end, bp_per_px, margin, stacks):
    """Draw BED track in SVG"""
    header_h = draw_svg_track_header(svg, "BED Annotation", y, width + 2 * margin)
    current_y = y + header_h + 5

    for i, feat in enumerate(features):
        stack = stacks[i]
        feat_y = current_y + stack * 20
        mid_y = feat_y + 5

        fx0 = margin + int((max(feat.start, start) - start) / bp_per_px)
        fx1 = margin + int((min(feat.end, end) - start) / bp_per_px)

        if fx1 <= fx0:
            continue

        # Default colors
        default_color = "#6495ed"  # Cornflower Blue
        thick_color = "#c8a032"    # Brownish Yellow

        # Use itemRgb if available
        if feat.item_rgb:
            color_hex = rgb_to_hex(feat.item_rgb)
        else:
            color_hex = default_color

        feature_height = 10

        # Draw main feature rectangle
        SubElement(svg, "rect", {
            "x": str(fx0), "y": str(mid_y - feature_height / 2),
            "width": str(fx1 - fx0), "height": str(feature_height),
            "fill": color_hex, "stroke": "black", "stroke-width": "0.5"
        })

        # Draw thickStart/thickEnd if specified (CDS-like region)
        if feat.thick_start is not None and feat.thick_end is not None:
            thick_x0 = margin + int((max(feat.thick_start, start) - start) / bp_per_px)
            thick_x1 = margin + int((min(feat.thick_end, end) - start) / bp_per_px)
            if thick_x1 > thick_x0:
                thick_color_hex = rgb_to_hex(feat.item_rgb) if feat.item_rgb else thick_color
                SubElement(svg, "rect", {
                    "x": str(thick_x0), "y": str(mid_y - feature_height / 2),
                    "width": str(thick_x1 - thick_x0), "height": str(feature_height),
                    "fill": thick_color_hex, "stroke": "black", "stroke-width": "0.5"
                })

        # Draw blocks (exons) if available
        if feat.blocks:
            for block in feat.blocks:
                bx0 = margin + int((max(block.start, start) - start) / bp_per_px)
                bx1 = margin + int((min(block.end, end) - start) / bp_per_px)
                if bx1 > bx0:
                    SubElement(svg, "rect", {
                        "x": str(bx0), "y": str(mid_y - feature_height / 2),
                        "width": str(bx1 - bx0), "height": str(feature_height),
                        "fill": color_hex, "stroke": "black", "stroke-width": "0.5"
                    })

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
                            SubElement(svg, "line", {
                                "x1": str(line_x0), "y1": str(mid_y),
                                "x2": str(line_x1), "y2": str(mid_y),
                                "stroke": "black", "stroke-width": "1"
                            })
        else:
            # If no blocks, draw intron line for the whole feature
            SubElement(svg, "line", {
                "x1": str(fx0), "y1": str(mid_y),
                "x2": str(fx1), "y2": str(mid_y),
                "stroke": "black", "stroke-width": "1"
            })

        # Strand arrow - solid black arrow outside the feature with connecting line
        head_size = 6
        arrow_width = 4
        connector_length = 3  # Short line connecting feature to arrow
        if feat.strand == '+':
            # Forward strand: arrow pointing right, outside the feature
            arrow_x = fx1 + connector_length
            # Draw connecting line
            SubElement(svg, "line", {
                "x1": str(fx1), "y1": str(mid_y),
                "x2": str(arrow_x), "y2": str(mid_y),
                "stroke": "black", "stroke-width": "1"
            })
            # Draw arrow triangle (pointing right)
            SubElement(svg, "polygon", {
                "points": f"{arrow_x},{mid_y-arrow_width} {arrow_x+head_size},{mid_y} {arrow_x},{mid_y+arrow_width}",
                "fill": "black", "stroke": "black", "stroke-width": "0.5"
            })
        elif feat.strand == '-':
            # Reverse strand: arrow pointing left, outside the feature
            arrow_x = fx0 - connector_length
            # Draw connecting line
            SubElement(svg, "line", {
                "x1": str(fx0), "y1": str(mid_y),
                "x2": str(arrow_x), "y2": str(mid_y),
                "stroke": "black", "stroke-width": "1"
            })
            # Draw arrow triangle (pointing left)
            SubElement(svg, "polygon", {
                "points": f"{arrow_x},{mid_y-arrow_width} {arrow_x-head_size},{mid_y} {arrow_x},{mid_y+arrow_width}",
                "fill": "black", "stroke": "black", "stroke-width": "0.5"
            })

        # Feature name
        if feat.name:
            SubElement(svg, "text", {
                "x": str(fx0), "y": str(mid_y + 16),
                "font-size": "10", "fill": "black"
            }).text = feat.name

    num_stacks = max(stacks) + 1 if stacks else 0
    return header_h + num_stacks * 20 + 10


# ----------------------------------------------------------------------------
# SVG VCF highlight track
# ----------------------------------------------------------------------------
# Mirrors draw_vcf_highlight_track in png_renderer.py — same layout, same
# phasing-block semantics.

SVG_VCF_MARKER_ROW_H = 10
SVG_VCF_BASE_ROW_H = 11
SVG_VCF_ROW_GAP = 2
SVG_VCF_GUTTER = 90


def _svg_sample_ploidy(sites, sname: str) -> int:
    p = 0
    for s in sites:
        call = s.sample_calls.get(sname)
        if call is None:
            continue
        if call.is_phased and call.hap_bases:
            p = max(p, len(call.hap_bases))
    return p if p > 0 else 2


def _svg_vcf_highlight_track_height_rows(rows_per_sample):
    header_h = 15
    total_rows = 1 + sum(rows_per_sample)
    rows_px = (
        SVG_VCF_MARKER_ROW_H
        + sum(rows_per_sample) * SVG_VCF_BASE_ROW_H
        + (total_rows - 1) * SVG_VCF_ROW_GAP
    )
    return header_h + 4 + rows_px + 6


def svg_vcf_highlight_track_height(num_samples: int) -> int:
    """Legacy 2-rows-per-sample height (kept for ABI back-compat)."""
    return _svg_vcf_highlight_track_height_rows([2] * num_samples)


def _svg_phase_block_spans(sites, sname):
    """Compute (ps_id, [site_indices]) runs for one sample.

    Mirrors png_renderer._phase_block_spans: a PS block covers ALL sites
    between the first and last occurrence of that PS in this sample's
    calls, regardless of whether the intermediate sites are hom / het / no
    PS / missing. See PNG version for full rationale.
    """
    n = len(sites)

    ps_first = {}
    ps_last = {}
    for i, s in enumerate(sites):
        call = s.sample_calls.get(sname)
        if call is None or not call.is_phased or call.phase_set is None:
            continue
        ps = call.phase_set
        if ps not in ps_first:
            ps_first[ps] = i
        ps_last[ps] = i

    eff_ps = [None] * n
    ordered_ps = sorted(ps_first.keys(), key=lambda k: ps_first[k])
    for ps in ordered_ps:
        for i in range(ps_first[ps], ps_last[ps] + 1):
            eff_ps[i] = ps

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


def draw_svg_vcf_highlight_track(
    svg,
    sites: List[HighlightSite],
    samples: List[str],
    y: int,
    width: int,
    start: int,
    end: int,
    bp_per_px: float,
    margin: int,
) -> int:
    """SVG counterpart of draw_vcf_highlight_track (PNG).

    Each sample gets 2 rows. Phased sites (GT uses '|' and PS present) are
    rendered as hap1/hap2 colors and a dark-gray backdrop spans all SNPs
    sharing a PS. Unphased sites fall back to REF/ALT semantics.
    """
    header_h = draw_svg_track_header(svg, "VCF Highlights", y, width + 2 * margin)
    cursor_y = y + header_h + 4

    plot_x0 = margin + SVG_VCF_GUTTER
    plot_x1 = margin + width

    def site_to_x(pos: int) -> int:
        return margin + int((pos - start) / bp_per_px)

    cell_w = max(2, int(round(1.0 / bp_per_px)) if bp_per_px < 1 else 2)

    # --- SNP marker row ------------------------------------------------------
    marker_y0 = cursor_y
    marker_y1 = cursor_y + SVG_VCF_MARKER_ROW_H
    SubElement(svg, "rect", {
        "x": str(margin), "y": str(marker_y0),
        "width": str(width), "height": str(SVG_VCF_MARKER_ROW_H),
        "fill": "#fafafa"
    })
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
            ref_rgb = MISMATCH_COLORS.get(site.ref, (120, 120, 120)) if site.is_snv else (120, 120, 120)
            SubElement(svg, "rect", {
                "x": str(x), "y": str(marker_y0),
                "width": str(cell_w), "height": str(mid - marker_y0),
                "fill": rgb_to_hex(ref_rgb)
            })
        if alt_base_shown is not None:
            alt_rgb = MISMATCH_COLORS.get(alt_base_shown, (120, 120, 120)) if site.is_snv else (120, 120, 120)
            SubElement(svg, "rect", {
                "x": str(x), "y": str(mid),
                "width": str(cell_w), "height": str(marker_y1 - mid),
                "fill": rgb_to_hex(alt_rgb)
            })
    SubElement(svg, "rect", {
        "x": str(margin), "y": str(marker_y0),
        "width": str(SVG_VCF_GUTTER), "height": str(SVG_VCF_MARKER_ROW_H),
        "fill": "#f0f0f0"
    })
    SubElement(svg, "text", {
        "x": str(margin + 2), "y": str(marker_y0 + SVG_VCF_MARKER_ROW_H - 2),
        "font-size": "9", "fill": "#464646"
    }).text = "SNP sites"
    cursor_y = marker_y1 + SVG_VCF_ROW_GAP

    # --- Per-sample rows (N rows per sample, N = observed ploidy) -----------
    PHASE_BG = "#d8d8d8"
    UNPHASED_BG = "#f8f8f8"
    for sname in samples:
        is_single_synthetic = len(samples) == 1 and sname == "VCF"
        short_name = sname if len(sname) <= 12 else sname[:11] + "…"

        blocks = _svg_phase_block_spans(sites, sname)

        ploidy = _svg_sample_ploidy(sites, sname)
        any_phased = any(
            sites[i].sample_calls.get(sname, SampleCall()).is_phased
            for i in range(len(sites))
        )
        n_rows = ploidy if any_phased else 2

        row_y0 = [cursor_y + r * (SVG_VCF_BASE_ROW_H + SVG_VCF_ROW_GAP) for r in range(n_rows)]
        sample_top = row_y0[0]
        sample_bottom = row_y0[-1] + SVG_VCF_BASE_ROW_H

        # Plain row backgrounds
        for ry in row_y0:
            SubElement(svg, "rect", {
                "x": str(margin), "y": str(ry),
                "width": str(width), "height": str(SVG_VCF_BASE_ROW_H),
                "fill": UNPHASED_BG
            })

        # Phase-block backdrops (span all N rows)
        for ps_id, idx_list in blocks:
            if ps_id is None:
                continue
            first_x = site_to_x(sites[idx_list[0]].pos)
            last_x = site_to_x(sites[idx_list[-1]].pos) + cell_w
            first_x = max(plot_x0, first_x)
            last_x = min(plot_x1, last_x)
            if last_x <= first_x:
                continue
            SubElement(svg, "rect", {
                "x": str(first_x), "y": str(sample_top),
                "width": str(last_x - first_x),
                "height": str(sample_bottom - sample_top),
                "fill": PHASE_BG
            })

        sites_in_ps = set()
        for ps_id, idx_list in blocks:
            if ps_id is not None:
                sites_in_ps.update(idx_list)

        # Cells
        for site_idx, site in enumerate(sites):
            x = site_to_x(site.pos)
            if x < margin or x >= margin + width:
                continue
            call = site.sample_calls.get(sname, SampleCall())
            if call.is_phased and call.phase_set is not None and call.hap_bases:
                for r_idx, hap_b in enumerate(call.hap_bases):
                    if r_idx >= n_rows or hap_b is None:
                        continue
                    if site.is_snv:
                        rgb = MISMATCH_COLORS.get(hap_b, (150, 150, 150))
                    else:
                        rgb = (150, 150, 150)
                    SubElement(svg, "rect", {
                        "x": str(x), "y": str(row_y0[r_idx]),
                        "width": str(cell_w), "height": str(SVG_VCF_BASE_ROW_H),
                        "fill": rgb_to_hex(rgb)
                    })
            elif site_idx in sites_in_ps and any_phased:
                hom_base = None
                if call.has_ref and call.alt_base is None:
                    hom_base = site.ref
                elif (not call.has_ref) and call.alt_base is not None:
                    hom_base = call.alt_base
                if hom_base is not None:
                    rgb = MISMATCH_COLORS.get(hom_base, (150, 150, 150)) if site.is_snv else (150, 150, 150)
                    for r_idx in range(n_rows):
                        SubElement(svg, "rect", {
                            "x": str(x), "y": str(row_y0[r_idx]),
                            "width": str(cell_w), "height": str(SVG_VCF_BASE_ROW_H),
                            "fill": rgb_to_hex(rgb)
                        })
            else:
                if call.has_ref:
                    rgb = MISMATCH_COLORS.get(site.ref, (150, 150, 150)) if site.is_snv else (150, 150, 150)
                    SubElement(svg, "rect", {
                        "x": str(x), "y": str(row_y0[0]),
                        "width": str(cell_w), "height": str(SVG_VCF_BASE_ROW_H),
                        "fill": rgb_to_hex(rgb)
                    })
                if call.alt_base is not None and n_rows >= 2:
                    rgb = MISMATCH_COLORS.get(call.alt_base, (150, 150, 150)) if site.is_snv else (150, 150, 150)
                    SubElement(svg, "rect", {
                        "x": str(x), "y": str(row_y0[-1]),
                        "width": str(cell_w), "height": str(SVG_VCF_BASE_ROW_H),
                        "fill": rgb_to_hex(rgb)
                    })

        # Gutter overlays + labels
        if is_single_synthetic:
            labels = ["REF", "ALT"]
        elif any_phased:
            labels = [f"{short_name} h{r + 1}" for r in range(n_rows)]
        else:
            labels = [f"{short_name} REF", f"{short_name} ALT"]
        for r_idx in range(n_rows):
            ry = row_y0[r_idx]
            SubElement(svg, "rect", {
                "x": str(margin), "y": str(ry),
                "width": str(SVG_VCF_GUTTER), "height": str(SVG_VCF_BASE_ROW_H),
                "fill": "#f0f0f0"
            })
            if r_idx < len(labels):
                SubElement(svg, "text", {
                    "x": str(margin + 2), "y": str(ry + SVG_VCF_BASE_ROW_H - 2),
                    "font-size": "9", "fill": "#464646"
                }).text = labels[r_idx]
        cursor_y = sample_bottom + SVG_VCF_ROW_GAP

    SubElement(svg, "line", {
        "x1": str(margin), "y1": str(cursor_y),
        "x2": str(margin + width), "y2": str(cursor_y),
        "stroke": "#dcdcdc", "stroke-width": "1"
    })
    return (cursor_y - y) + 6


def render_svg_snapshot(
    tracks: List[Dict[str, Any]],
    chrom: str,
    start: int,
    end: int,
    width: int = 1200,
    read_height: int = 6,
    detail: str = "mid",
    show_axis: bool = True,
    show_coverage: bool = True,
    coverage_height: int = 15,
    show_composition: bool = False,
    composition_height: Optional[int] = None,
    comp_max_depth: Optional[int] = None,
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
    **kwargs
) -> str:
    """Render snapshot as SVG string"""
    # Actual drawing area width (excluding left and right margins)
    content_width = width - 2 * margin
    bp_per_px = (end - start) / float(content_width)

    # Calculate total height
    top = 0
    axis_height = 34  # chromosome label + coordinate labels + axis line
    if show_axis:
        top += axis_height

    total_height = top

    composition_track_h = 0
    if show_composition:
        comp_h = composition_height if composition_height is not None else (read_height * 2)
        composition_track_h = max(1, int(comp_h)) + 10
        total_height += composition_track_h

    # Pre-calculate stacks and height for each track
    track_meta = []

    # Each BAM/track gets its own coverage track. Do not merge multiple BAMs
    # into one aggregate depth profile.
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

    # VCF highlight track height
    vcf_track_h = 0
    if highlight_sites and highlight_samples:
        rows_per_sample = [_svg_sample_ploidy(highlight_sites, s) for s in highlight_samples]
        vcf_track_h = _svg_vcf_highlight_track_height_rows(rows_per_sample)
        total_height += vcf_track_h

    # When highlight is active:
    #   - filter reads that don't cover any VCF site (unless --no-hap-filter)
    #   - sort by observed hap signature (unless --no-hap-sort)
    do_hap_sort = bool(highlight_sites) and not no_hap_sort
    do_hap_filter = bool(highlight_sites) and not no_hap_filter
    if highlight_sites:
        site_positions = sorted(s.pos for s in highlight_sites)
        if site_positions:
            first_p, last_p = site_positions[0], site_positions[-1]

            def _covers_any(r):
                if r.end <= first_p or r.start > last_p:
                    return False
                for p in site_positions:
                    if r.start <= p < r.end:
                        return True
                return False

            for track in tracks:
                if do_hap_filter:
                    track['reads'] = [r for r in track['reads'] if _covers_any(r)]
                if do_hap_sort:
                    track['reads'] = sorted(
                        track['reads'],
                        key=lambda r: (_hap_signature(r, site_positions), r.start),
                    )

    for track in tracks:
        reads = track['reads']
        if do_hap_sort:
            # VCF highlight 模式下，为保持 haplotype 排序，每条 read 一行
            stacks = list(range(len(reads)))
        else:
            if is_rna:
                spans = [(max(r.start, start), min(r.end, end)) for r in reads]
                keys = [r.qname for r in reads]
                stacks = assign_stacks_grouped(spans, keys)
            else:
                stacks = assign_split_read_stacks(reads, start, end)

        row_offsets, reads_area_h = _stack_row_offsets(stacks, reads, read_height)

        track_h = per_track_coverage_h + 15 + reads_area_h + 5
        total_height += track_h

        track_meta.append({
            'stacks': stacks,
            'row_offsets': row_offsets,
            'reads_area_h': reads_area_h
        })

    # Create SVG root element
    svg = Element("svg", {
        "width": str(width),
        "height": str(total_height),
        "xmlns": "http://www.w3.org/2000/svg"
    })

    # Background
    SubElement(svg, "rect", {
        "x": "0", "y": "0",
        "width": str(width),
        "height": str(total_height),
        "fill": "white"
    })

    current_y = 0

    # Collect tick positions for drawing guide lines
    tick_positions = []

    # Draw coordinate axis
    if show_axis:
        # Axis layout: chromosome name on top, coordinate labels below it,
        # and the axis line at the bottom of the reserved axis area.
        chrom_label_y = current_y + 10
        tick_label_y = current_y + 23
        axis_y = current_y + 30

        SubElement(svg, "text", {
            "x": str(margin),
            "y": str(chrom_label_y),
            "font-size": "10",
            "font-weight": "bold",
            "fill": "black",
            "text-anchor": "start"
        }).text = chrom

        SubElement(svg, "line", {
            "x1": str(margin), "y1": str(axis_y),
            "x2": str(width - margin), "y2": str(axis_y),
            "stroke": "black",
            "stroke-width": "1"
        })

        span = max(0, end - start)
        step = max(1, span // 10) if span > 0 else 1
        tick_bp_positions = list(range(start, end + 1, step))
        if not tick_bp_positions or tick_bp_positions[-1] != end:
            tick_bp_positions.append(end)

        for pos in tick_bp_positions:
            x = margin + int((pos - start) / bp_per_px) if bp_per_px else margin
            x = max(margin, min(width - margin, x))
            tick_positions.append(x)

            # Tick mark
            SubElement(svg, "line", {
                "x1": str(x), "y1": str(axis_y - 3),
                "x2": str(x), "y2": str(axis_y + 3),
                "stroke": "black",
                "stroke-width": "1"
            })

            # Keep edge labels inside the SVG canvas.  The rightmost label was
            # previously clipped because it was drawn at x+2 with the default
            # left anchor.
            label_x = x
            text_anchor = "middle"
            if x <= margin + 1:
                label_x = margin
                text_anchor = "start"
            elif x >= width - margin - 1:
                label_x = width - margin
                text_anchor = "end"

            text = SubElement(svg, "text", {
                "x": str(label_x),
                "y": str(tick_label_y),
                "font-size": "9",
                "fill": "black",
                "text-anchor": text_anchor
            })
            text.text = f"{pos:,}"
        current_y += axis_height
    
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
    
        SubElement(svg, "rect", {
            "x": str(fx0),
            "y": str(current_y),
            "width": str(max(1, fx1 - fx0)),
            "height": str(total_height - current_y),
            "fill": "#ebf2ff"
        })
    
    # Draw vertical guide lines (dashed) - from axis to bottom of image
    for tick_x in tick_positions:
        SubElement(svg, "line", {
            "x1": str(tick_x), "y1": str(current_y),
            "x2": str(tick_x), "y2": str(total_height),
            "stroke": "#dcdcdc",
            "stroke-width": "1",
            "stroke-dasharray": "3,3"  # Dashed: 3px line, 3px gap
        })

    # Per-BAM coverage tracks are drawn immediately above each read track.

    if show_composition:
        composition_reads: List[Read] = []
        for track in tracks:
            composition_reads.extend(track.get('reads', []))
        current_y += draw_svg_base_composition_track(
            svg,
            composition_reads,
            current_y,
            width,
            content_width,
            start,
            end,
            composition_height,
            comp_max_depth,
            margin,
            read_height,
        )

    # Draw annotation track if enabled
    if gff_genes:
        current_y += draw_svg_gene_track(
            svg,
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
        current_y += draw_svg_bed_track(
            svg,
            bed_features,
            current_y,
            content_width,
            start,
            end,
            bp_per_px,
            margin,
            bed_stacks
        )

    # Draw VCF highlight track right above the reads
    if highlight_sites and highlight_samples:
        current_y += draw_svg_vcf_highlight_track(
            svg,
            highlight_sites,
            highlight_samples,
            current_y,
            content_width,
            start,
            end,
            bp_per_px,
            margin,
        )

    # Pos -> HighlightSite index; empty when no highlight VCF is active
    # (preserves original full-color mismatch behavior).
    svg_highlight_index = {}
    if highlight_sites:
        svg_highlight_index = {s.pos: s for s in highlight_sites}

    # Iterate over tracks
    for i, track in enumerate(tracks):
        reads = track['reads']
        track_title = track.get('title', 'Reads')
        meta = track_meta[i]
        stacks = meta['stacks']
        row_offsets = meta['row_offsets']
        reads_area_h = meta['reads_area_h']

        if show_coverage:
            current_y = draw_svg_per_track_coverage(
                svg, reads, track_title, current_y, width, content_width,
                coverage_height, start, end, ref_seq, coverage_max_depth,
                margin, is_rna, detail, bp_per_px
            )

        # Draw read track header
        header_y = current_y
        SubElement(svg, "rect", {
            "x": "0", "y": str(header_y),
            "width": str(width),
            "height": "15",
            "fill": "#f5f5f5"
        })
        SubElement(svg, "line", {
            "x1": "0", "y1": str(header_y + 14),
            "x2": str(width), "y2": str(header_y + 14),
            "stroke": "#c8c8c8",
            "stroke-width": "1"
        })
        header_text = SubElement(svg, "text", {
            "x": "5",
            "y": str(header_y + 11),
            "font-size": "11",
            "fill": "black"
        })
        header_text.text = track_title
        current_y += 15

        # Draw reads
        reads_start_y = current_y
        groups: Dict[str, List[int]] = {}
        for idx_r, r in enumerate(reads):
            groups.setdefault(r.qname, []).append(idx_r)
        supplementary_read_styles = _build_supplementary_read_styles(reads, groups, stacks)

        for idx, r in enumerate(reads):
            y = reads_start_y + row_offsets[stacks[idx]]
            rects = segments_to_pixels(r.segments, r.start, start, bp_per_px, detail=detail)
            supp_style = supplementary_read_styles.get(idx)
            # Draw insertion blocks after the read body. If they are drawn in
            # segment order, the following match block can cover part of the
            # purple insertion because insertion consumes no reference pixels.
            # The visual block width is set to the same one-reference-base
            # pixel width used by a single-base mismatch rectangle.
            insertion_overlays = []

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
                x0_draw = margin + x0
                x1_draw = margin + x1

                if t == "ins":
                    if x0 < 0 or x0 > (width - 2 * margin):
                        continue

                x0_draw = max(margin, min(width - margin, x0_draw))
                x1_draw = max(margin, min(width - margin, x1_draw))

                if x1_draw <= x0_draw and t != "ins":
                    continue

                if color_by == "type":
                    color = color_for_type(t)
                elif color_by == "strand":
                    color = STRAND_COLORS["-" if r.reverse else "+"]
                elif color_by == "mapq":
                    base_color = color_for_type(t)
                    color = shade_by_mapq(base_color, r.mapq)
                else:
                    color = color_for_type(t)

                color_hex = rgb_to_hex(color)

                if t == "ins":
                    # In highlight mode, suppress ins ticks/labels entirely
                    # (long-read ins noise distracts from VCF-site analysis).
                    if svg_highlight_index:
                        pass
                    else:
                        label = None
                        if detail == "high" and show_insertion_labels:
                            seg_idx = rect_to_seg.get(rect_idx)
                            if seg_idx is not None:
                                seg = r.segments[seg_idx]
                                if seg and seg.length > 0:
                                    label = f"I({seg.length})"
                        # Render insertion as a virtual 1-bp block so its
                        # visual size follows a single-base mismatch rectangle.
                        # segments_to_pixels stores insertion as x0 == x1
                        # because it consumes no reference bases, so we add the
                        # one-base display width here and draw it later on top.
                        one_base_px = max(1, int(round(1.0 / bp_per_px)))
                        ins_x0 = x0_draw
                        ins_x1 = min(width - margin, ins_x0 + one_base_px)
                        if ins_x1 <= ins_x0:
                            ins_x0 = max(margin, ins_x0 - one_base_px)
                            ins_x1 = min(width - margin, ins_x0 + one_base_px)
                        insertion_overlays.append((ins_x0, ins_x1, y, read_height, color_hex, label))
                elif t == "ref_skip":
                    y_center = y + read_height // 2
                    SubElement(svg, "line", {
                        "x1": str(x0_draw), "y1": str(y_center),
                        "x2": str(x1_draw), "y2": str(y_center),
                        "stroke": "#b0c4de", "stroke-width": "1"
                    })
                elif t == "del":
                    # Highlight mode: draw del in muted-gray so it blends in
                    # with matches and doesn't steal attention from VCF sites.
                    del_fill = MUTED_MISMATCH_HEX if svg_highlight_index else "#808080"
                    SubElement(svg, "rect", {
                        "x": str(x0_draw), "y": str(y),
                        "width": str(x1_draw - x0_draw), "height": str(read_height),
                        "fill": del_fill, "stroke": "none", "stroke-width": "1"
                    })
                    if not svg_highlight_index and detail == "high" and show_insertion_labels:
                        seg_idx = rect_to_seg.get(rect_idx)
                        if seg_idx is not None:
                            seg = r.segments[seg_idx]
                            if seg and seg.length > 0:
                                del_width = x1_draw - x0_draw
                                label = str(seg.length)
                                if del_width > 12:  # Enough space inside
                                    SubElement(svg, "text", {
                                        "x": str((x0_draw + x1_draw) / 2), "y": str(y + read_height - 1),
                                        "font-size": "7", "fill": "white",
                                        "text-anchor": "middle"
                                    }).text = label
                                else:  # Too small, put above
                                    SubElement(svg, "text", {
                                        "x": str((x0_draw + x1_draw) / 2), "y": str(y - 2),
                                        "font-size": "8", "fill": "black",
                                        "text-anchor": "middle"
                                    }).text = label
                elif t == "mismatch":
                    # See png_renderer.py for the rationale. When no highlight
                    # VCF is active we keep the original behavior; otherwise
                    # we paint the rect muted gray then overdraw per-base
                    # colored cells only where the ref position is a VCF site.
                    if detail == "low":
                        if supp_style:
                            fill_rgb, fill_opacity = supp_style
                            SubElement(svg, "rect", _svg_rect_attrs_with_optional_opacity(
                                x0_draw, y, x1_draw - x0_draw, read_height,
                                rgb_to_hex(fill_rgb), fill_opacity
                            ))
                        else:
                            SubElement(svg, "rect", {
                                "x": str(x0_draw), "y": str(y),
                                "width": str(x1_draw - x0_draw), "height": str(read_height),
                                "fill": "#c3c3c3", "stroke": "none"
                            })
                    elif not svg_highlight_index:
                        current_color_hex = color_hex
                        if r.seq and rect_idx in rect_to_seg:
                            seg_idx = rect_to_seg[rect_idx]
                            read_cursor = sum(s.read_consumed for s in r.segments[:seg_idx])
                            if read_cursor < len(r.seq):
                                base = r.seq[read_cursor].upper()
                                mismatch_color = MISMATCH_COLORS.get(base, (200, 60, 60))
                                current_color_hex = rgb_to_hex(mismatch_color)
                            else:
                                current_color_hex = "#c3c3c3"
                        else:
                            current_color_hex = "#c3c3c3"
                        SubElement(svg, "rect", {
                            "x": str(x0_draw), "y": str(y),
                            "width": str(x1_draw - x0_draw), "height": str(read_height),
                            "fill": current_color_hex, "stroke": "none"
                        })
                    else:
                        # Highlight-aware SVG mismatch rendering.
                        if supp_style:
                            fill_rgb, fill_opacity = supp_style
                            SubElement(svg, "rect", _svg_rect_attrs_with_optional_opacity(
                                x0_draw, y, x1_draw - x0_draw, read_height,
                                rgb_to_hex(fill_rgb), fill_opacity
                            ))
                        else:
                            SubElement(svg, "rect", {
                                "x": str(x0_draw), "y": str(y),
                                "width": str(x1_draw - x0_draw), "height": str(read_height),
                                "fill": MUTED_MISMATCH_HEX, "stroke": "none"
                            })
                        if r.seq and rect_idx in rect_to_seg:
                            seg_idx = rect_to_seg[rect_idx]
                            seg = r.segments[seg_idx]
                            ref_pos0 = r.start + sum(
                                s.ref_consumed for s in r.segments[:seg_idx]
                            )
                            read_pos0 = sum(s.read_consumed for s in r.segments[:seg_idx])
                            for k in range(seg.length):
                                ref_p = ref_pos0 + k
                                if ref_p not in svg_highlight_index:
                                    continue
                                rk = read_pos0 + k
                                if rk >= len(r.seq):
                                    continue
                                base = r.seq[rk].upper()
                                col = MISMATCH_COLORS.get(base, (200, 60, 60))
                                bx0 = margin + int((ref_p - start) / bp_per_px)
                                bx1 = margin + int((ref_p + 1 - start) / bp_per_px)
                                if bx1 <= bx0:
                                    bx1 = bx0 + 1
                                bx0 = max(margin, min(width - margin, bx0))
                                bx1 = max(margin, min(width - margin, bx1))
                                SubElement(svg, "rect", {
                                    "x": str(bx0), "y": str(y),
                                    "width": str(bx1 - bx0), "height": str(read_height),
                                    "fill": rgb_to_hex(col), "stroke": "none"
                                })
                else:
                    if t == "match":
                        # First paint the whole match block in muted gray,
                        # then in highlight mode overdraw per-base color
                        # cells at each VCF site using the read's actual base
                        # (matches png_renderer behavior).
                        if supp_style:
                            fill_rgb, fill_opacity = supp_style
                            SubElement(svg, "rect", _svg_rect_attrs_with_optional_opacity(
                                x0_draw, y, x1_draw - x0_draw, read_height,
                                rgb_to_hex(fill_rgb), fill_opacity
                            ))
                        else:
                            SubElement(svg, "rect", {
                                "x": str(x0_draw), "y": str(y),
                                "width": str(x1_draw - x0_draw), "height": str(read_height),
                                "fill": "#c3c3c3", "stroke": "none"
                            })
                        if svg_highlight_index and r.seq and rect_idx in rect_to_seg:
                            seg_idx = rect_to_seg[rect_idx]
                            seg = r.segments[seg_idx]
                            ref_pos0 = r.start + sum(
                                s.ref_consumed for s in r.segments[:seg_idx]
                            )
                            read_pos0 = sum(s.read_consumed for s in r.segments[:seg_idx])
                            for k in range(seg.length):
                                ref_p = ref_pos0 + k
                                if ref_p not in svg_highlight_index:
                                    continue
                                rk = read_pos0 + k
                                if rk >= len(r.seq):
                                    continue
                                base = r.seq[rk].upper()
                                col = MISMATCH_COLORS.get(base, (150, 150, 150))
                                bx0 = margin + int((ref_p - start) / bp_per_px)
                                bx1 = margin + int((ref_p + 1 - start) / bp_per_px)
                                if bx1 <= bx0:
                                    bx1 = bx0 + 1
                                bx0 = max(margin, min(width - margin, bx0))
                                bx1 = max(margin, min(width - margin, bx1))
                                SubElement(svg, "rect", {
                                    "x": str(bx0), "y": str(y),
                                    "width": str(bx1 - bx0), "height": str(read_height),
                                    "fill": rgb_to_hex(col), "stroke": "none"
                                })
                    else:
                        if supp_style and t in ("match", "soft", "hard"):
                            fill_rgb, fill_opacity = supp_style
                            SubElement(svg, "rect", _svg_rect_attrs_with_optional_opacity(
                                x0_draw, y, x1_draw - x0_draw, read_height,
                                rgb_to_hex(fill_rgb), fill_opacity
                            ))
                        else:
                            SubElement(svg, "rect", {
                                "x": str(x0_draw), "y": str(y),
                                "width": str(x1_draw - x0_draw), "height": str(read_height),
                                "fill": color_hex, "stroke": "none"
                            })

            # Draw insertion blocks last so adjacent match blocks cannot cover them.
            for ins_x0, ins_x1, ins_y, ins_h, ins_color, ins_label in insertion_overlays:
                ins_w = max(1, ins_x1 - ins_x0)
                SubElement(svg, "rect", {
                    "x": str(ins_x0), "y": str(ins_y),
                    "width": str(ins_w), "height": str(ins_h),
                    "fill": ins_color, "stroke": "none",
                    "pointer-events": "none",
                })
                if ins_label:
                    SubElement(svg, "text", {
                        "x": str(ins_x0 + ins_w / 2.0), "y": str(ins_y - 2),
                        "font-size": "8", "fill": "#800080",
                        "text-anchor": "middle",
                        "pointer-events": "none",
                    }).text = ins_label

            # # Arrow head
            # if style == "jbrowse":
            #     w = x1_draw - x0_draw
            #     if w >= 6:
            #         head = max(3, min(6, w // 6))
            #         arrow_pts = ""
            #         if r.reverse:
            #             p1 = f"{x0_draw},{y}"
            #             p2 = f"{x0_draw+head},{y+read_height/2}"
            #             p3 = f"{x0_draw},{y+read_height}"
            #             arrow_pts = f"{p1} {p2} {p3}"
            #         else:
            #             p1 = f"{x1_draw},{y}"
            #             p2 = f"{x1_draw-head},{y+read_height/2}"
            #             p3 = f"{x1_draw},{y+read_height}"
            #             arrow_pts = f"{p1} {p2} {p3}"
            #         SubElement(svg, "polygon", {
            #             "points": arrow_pts, "fill": "#646464"
            #         })
            
                    # Draw read direction arrow once per read, not once per segment
            if style == "jbrowse" and not is_rna:
                read_x0 = margin + int((max(r.start, start) - start) / bp_per_px)
                read_x1 = margin + int((min(r.end, end) - start) / bp_per_px)

                # Clamp to plotting area.  Keep the arrow inside the content
                # region, not inside the full SVG width.
                plot_x0 = margin
                plot_x1 = width - margin - 1
                read_x0 = max(plot_x0, min(plot_x1, read_x0))
                read_x1 = max(plot_x0, min(plot_x1, read_x1))

                read_w = read_x1 - read_x0
                if read_w >= 6:
                    head = max(4, min(max(read_height + 1, 5), 8, max(1, read_w - 1)))
                    mid_y = y + read_height / 2

                    if r.reverse:
                        # reverse strand: arrow points left
                        arrow_pts = (
                            f"{read_x0},{mid_y} "
                            f"{read_x0 + head},{y-1.5} "
                            f"{read_x0 + head},{y + read_height+1.5}"
                        )
                    else:
                        # forward strand: arrow points right
                        arrow_pts = (
                            f"{read_x1},{mid_y} "
                            f"{read_x1 - head},{y-1.5} "
                            f"{read_x1 - head},{y + read_height+1.5}"
                        )

                    arrow_attrs = {
                        "points": arrow_pts,
                        "fill": "#ff0000"
                    }
                    if supp_style:
                        arrow_attrs["fill"] = "#ff0000"
                    SubElement(svg, "polygon", arrow_attrs)

        # Connect segments
        if is_rna:
            for qname, idxs in groups.items():
                if len(idxs) > 1:
                    idxs_sorted = sorted(idxs, key=lambda i: reads[i].start)
                    for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                        if reads[a].end < reads[b].start:
                            ya = reads_start_y + row_offsets[stacks[a]]
                            yb = reads_start_y + row_offsets[stacks[b]]
                            xa = margin + int((reads[a].end - start) / bp_per_px)
                            xb = margin + int((reads[b].start - start) / bp_per_px)
                            xa = max(margin, min(width - margin - 1, xa))
                            xb = max(margin, min(width - margin - 1, xb))
                            ya_center = ya + read_height // 2
                            yb_center = yb + read_height // 2
                            SubElement(svg, "line", {
                                "x1": str(xa), "y1": str(ya_center),
                                "x2": str(xb), "y2": str(yb_center),
                                "stroke": "#B0C4DE", "stroke-width": "0.5"
                            })
        else:
            # DNA supplementary/split-read connector lines are intentionally
            # disabled.  The relationship is now encoded by qname color plus
            # query-order opacity, so no extra linking line is drawn.
            pass

        current_y = reads_start_y + reads_area_h + 5

    # Convert to string
    rough_string = tostring(svg, encoding='unicode')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
