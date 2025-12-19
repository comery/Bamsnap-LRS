from typing import List, Dict, Optional, Tuple, Any

from PIL import Image, ImageDraw, ImageFont

from .layout import segments_to_pixels, assign_stacks
from .reader import Read
from .pileup import base_pileup, pileup_to_pixels
from .styles import color_for_type, color_for_base, shade_by_mapq, STRAND_COLORS, MISMATCH_COLORS


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
        )
    
    # For default style, flatten reads from all tracks
    reads = []
    for t in tracks:
        reads.extend(t['reads'])

    bp_per_px = (end - start) / float(width)
    spans = [(max(r.start, start), min(r.end, end)) for r in reads]
    stacks = assign_stacks(spans, max_stack=max(1, len(reads)))
    top = 0
    if show_axis:
        top += 20
    comp_h = composition_height if composition_height is not None else (read_height * 2)
    if show_composition:
        top += comp_h + 10
    height = top + (max(stacks) + 1) * (read_height + 2) + 10
    img = Image.new("RGB", (width, height), (255, 255, 255))
    dr = ImageDraw.Draw(img)
    groups: Dict[str, List[int]] = {}
    for i, r in enumerate(reads):
        groups.setdefault(r.qname, []).append(i)
    if show_axis:
        dr.line([(0, 10), (width - 1, 10)], fill=(0, 0, 0))
        span = end - start
        step = max(1, span // 10)
        for pos in range(start, end + 1, step):
            x = int((pos - start) / bp_per_px)
            dr.line([(x, 7), (x, 13)], fill=(0, 0, 0))
            dr.text((x + 2, 1), str(pos), fill=(0, 0, 0))
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
        y = top + stacks[idx] * (read_height + 2)
        rects = segments_to_pixels(r.segments, r.start, start, bp_per_px, detail=detail)
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
                dr.line([(x0, y), (x0, y + read_height)], fill=color)
            elif t == "ref_skip":
                # Intron: draw as line
                y_center = y + read_height // 2
                dr.line([(x0, y_center), (x1, y_center)], fill=(176, 196, 222), width=1)
            elif t == "del":
                dr.rectangle([(x0, y), (x1, y + read_height)], outline=(120, 120, 120), fill=color)
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
                head = max(3, min(6, w // 6))
                if r.reverse:
                    dr.polygon([(x0, y), (x0 + head, y + read_height // 2), (x0, y + read_height)], fill=(100, 100, 100))
                else:
                    dr.polygon([(x1, y), (x1 - head, y + read_height // 2), (x1, y + read_height)], fill=(100, 100, 100))
    for qname, idxs in groups.items():
        if len(idxs) > 1:
            idxs_sorted = sorted(idxs, key=lambda i: (reads[i].start, reads[i].end))
            for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                ya = top + stacks[a] * (read_height + 2)
                yb = top + stacks[b] * (read_height + 2)
                xa = int((reads[a].end - start) / bp_per_px)
                xb = int((reads[b].start - start) / bp_per_px)
                xa = max(0, min(width - 1, xa))
                xb = max(0, min(width - 1, xb))
                dr.line([(xa, ya + read_height // 2), (xb, yb + read_height // 2)], fill=(80, 160, 80))
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
        
        # Draw strand arrow only at the end of the gene
        head_size = 5
        if gene.strand == '+':
            dr.line([(gx1, mid_y), (gx1 - head_size, mid_y - head_size // 2)], fill=(0, 0, 0), width=1)
            dr.line([(gx1, mid_y), (gx1 - head_size, mid_y + head_size // 2)], fill=(0, 0, 0), width=1)
        elif gene.strand == '-':
            dr.line([(gx0, mid_y), (gx0 + head_size, mid_y - head_size // 2)], fill=(0, 0, 0), width=1)
            dr.line([(gx0, mid_y), (gx0 + head_size, mid_y + head_size // 2)], fill=(0, 0, 0), width=1)

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
) -> Image.Image:
    """Render JBrowse-style snapshot with track system and coverage chart"""
    # Handle backward compatibility: if tracks is list of Read, wrap it
    if tracks and not isinstance(tracks[0], dict):
        tracks = [{'reads': tracks, 'title': track_title}]

    # Actual drawing area width (excluding left and right margins)
    content_width = width - 2 * margin
    bp_per_px = (end - start) / float(content_width)

    # Pre-calculate stacks for each track
    track_stacks = []
    for track in tracks:
        reads = track['reads']
        spans = [(max(r.start, start), min(r.end, end)) for r in reads]
        stacks = assign_stacks(spans, max_stack=max(1, len(reads)))
        track_stacks.append(stacks)

    # Calculate total height
    top = 0
    if show_axis:
        top += 20  # Axis area (reduced)

    total_height = top
    
    # Calculate aggregate coverage height if enabled
    aggregate_coverage_h = 0
    if show_coverage:
        # Header (15) + Padding for arcs (coverage_height) + Coverage track (coverage_height) + Bottom margin (15)
        aggregate_coverage_h = 15 + coverage_height + coverage_height + 15
        total_height += aggregate_coverage_h
    
    # Calculate gene track height if enabled
    gene_track_h = 0
    gene_stacks = []
    if gff_genes:
        gene_spans = [(g.start, g.end) for g in gff_genes]
        gene_stacks = assign_stacks(gene_spans, max_stack=len(gff_genes))
        # Header (15) + Genes area
        num_stacks = max(gene_stacks) + 1 if gene_stacks else 0
        gene_track_h = 15 + num_stacks * 20 + 10
        total_height += gene_track_h

    # Calculate height for each track
    for i, track in enumerate(tracks):
        stacks = track_stacks[i]
        
        # Read track (reduced spacing, total height halved)
        track_header_height = 15  # Reduced header height
        reads_area_height = (max(stacks) + 1) * (read_height + 1) + 5  # Reduced read spacing
        
        total_height += track_header_height + reads_area_height + 5
    
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
    
    # Draw vertical guide lines (dashed) - from axis to bottom of image
    for tick_x in tick_positions:
        # Draw dashed line: draw 2px line every 4 pixels
        for dash_y in range(current_y, total_height, 6):
            dash_end = min(dash_y + 3, total_height)
            dr.line([(tick_x, dash_y), (tick_x, dash_end)], fill=(220, 220, 220), width=1)
    
    # Draw aggregate coverage track if enabled
    if show_coverage:
        all_reads = []
        for track in tracks:
            all_reads.extend(track['reads'])
            
        current_y += draw_track_header(dr, "Aggregate Coverage", current_y, width, coverage_height)
        
        # Add padding for arcs
        padding_for_arcs = coverage_height
        current_y += padding_for_arcs
        
        # Calculate base distribution for all reads
        base_distribution = calculate_base_distribution(all_reads, start, end, ref_seq)
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
            detail
        )
        
        # Draw pink connection lines in coverage track for RNA mode
        if is_rna:
            # Create a transparent overlay for arcs
            arc_layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
            dr_arc = ImageDraw.Draw(arc_layer)
            
            # Calculate coverage track Y position
            coverage_track_bottom = current_y + coverage_height
            
            # Pink color for arcs: (253, 209, 211)
            arc_color = (253, 209, 211, 150)
            
            # Anchor y position: 1/2 of coverage height
            arc_anchor_y = coverage_track_bottom - coverage_height // 2

            # 1. Draw arcs for ref_skip (introns) within each read
            for r in all_reads:
                ref_cursor = r.start
                for seg in r.segments:
                    if seg.type == "ref_skip":
                        seg_start = ref_cursor
                        seg_end = ref_cursor + seg.ref_consumed
                        
                        xa = margin + int((seg_start - start) / bp_per_px)
                        xb = margin + int((seg_end - start) / bp_per_px)
                        xa = max(margin, min(width - margin - 1, xa))
                        xb = max(margin, min(width - margin - 1, xb))
                        
                        if xb - xa < 2:
                            ref_cursor += seg.ref_consumed
                            continue
                            
                        # Draw arc
                        w = xb - xa
                        h = min(w // 2, coverage_height)
                        bbox = [xa, arc_anchor_y - h, xb, arc_anchor_y + h]
                        dr_arc.arc(bbox, start=180, end=0, fill=arc_color, width=1)
                        
                    ref_cursor += seg.ref_consumed

            # 2. Draw arcs for split reads
            groups_temp: Dict[str, List[int]] = {}
            for idx_r, r in enumerate(all_reads):
                groups_temp.setdefault(r.qname, []).append(idx_r)
            
            for qname, idxs in groups_temp.items():
                if len(idxs) > 1:
                    # Sort segments by genomic position
                    idxs_sorted = sorted(idxs, key=lambda i: all_reads[i].start)
                    for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                        if all_reads[a].end < all_reads[b].start:
                            xa = margin + int((all_reads[a].end - start) / bp_per_px)
                            xb = margin + int((all_reads[b].start - start) / bp_per_px)
                            xa = max(margin, min(width - margin - 1, xa))
                            xb = max(margin, min(width - margin - 1, xb))
                            
                            if xb - xa < 2:
                                continue
                                
                            w = xb - xa
                            h = min(w // 2, coverage_height)
                            bbox = [xa, arc_anchor_y - h, xb, arc_anchor_y + h]
                            dr_arc.arc(bbox, start=180, end=0, fill=arc_color, width=1)
            
            # Composite arc layer
            img = img.convert("RGBA")
            img = Image.alpha_composite(img, arc_layer)
            img = img.convert("RGB")
            dr = ImageDraw.Draw(img)
        
        current_y += coverage_height + 15  # Padding after coverage

    # Draw gene track if enabled
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

    # Iterate over tracks
    for i, track in enumerate(tracks):
        reads = track['reads']
        current_track_title = track.get('title', track_title)
        stacks = track_stacks[i]
        
        # Read track header
        reads_area_height = (max(stacks) + 1) * (read_height + 1) + 5
        current_y += draw_track_header(dr, current_track_title, current_y, width, reads_area_height)
        
        # Draw reads
        # (Rest of the read drawing logic remains same)
        reads_start_y = current_y
        groups: Dict[str, List[int]] = {}
        for idx_r, r in enumerate(reads):
            groups.setdefault(r.qname, []).append(idx_r)
        
        for idx, r in enumerate(reads):
            y = reads_start_y + stacks[idx] * (read_height + 1)
            rects = segments_to_pixels(r.segments, r.start, start, bp_per_px, detail=detail)
            
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
                    ins_color = (128, 0, 128)  # Purple
                    dr.line([(x0_draw, y), (x0_draw, y + read_height)], fill=ins_color, width=1)
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
                    del_color = (115, 115, 115)
                    dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=del_color)
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
                            label = str(seg.length)
                            del_width = x1_draw - x0_draw
                            if del_width > 4:
                                dr.text(((x0_draw + x1_draw) // 2 - 5, y + 1), label, fill=(255, 255, 255), font=font)
                            else:
                                dr.text((x0_draw, y - 10), label, fill=(0, 0, 0), font=font)
                elif t == "mismatch":
                    if detail != "low" and r.seq and rect_idx in rect_to_seg:
                        seg_idx = rect_to_seg[rect_idx]
                        read_cursor = sum(s.read_consumed for s in r.segments[:seg_idx])
                        if read_cursor < len(r.seq):
                            base = r.seq[read_cursor].upper()
                            mismatch_color = MISMATCH_COLORS.get(base, (200, 60, 60))
                            dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=mismatch_color)
                        else:
                            dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=(195, 195, 195))
                    else:
                        dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=(195, 195, 195))
                else:
                    if t == "match":
                        color = (195, 195, 195)
                    dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=color)
            
            # Read direction arrow
            if style == "jbrowse" and rects:
                read_end_px = margin + int((min(r.end, end) - start) / bp_per_px)
                read_start_px = margin + int((max(r.start, start) - start) / bp_per_px)
                read_end_px = max(margin, min(width - margin - 1, read_end_px))
                read_start_px = max(margin, min(width - margin - 1, read_start_px))
                head = max(3, min(read_height // 2, 5))
                arrow_fill = (211, 211, 211)
                arrow_truncated = (200, 200, 200)
                
                if r.reverse:
                    arrow_x = read_start_px
                    is_truncated = r.start < start
                    points = [(arrow_x, y), (arrow_x - head, y + read_height // 2), (arrow_x, y + read_height)]
                    if is_truncated:
                        dr.polygon(points, outline=arrow_truncated, fill=(255, 255, 255))
                    else:
                        dr.polygon(points, fill=arrow_fill)
                else:
                    arrow_x = read_end_px
                    is_truncated = r.end > end
                    points = [(arrow_x, y), (arrow_x + head, y + read_height // 2), (arrow_x, y + read_height)]
                    if is_truncated:
                        dr.polygon(points, outline=arrow_truncated, fill=(255, 255, 255))
                    else:
                        dr.polygon(points, fill=arrow_fill)
        
        # Connection lines
        if is_rna:
            for qname, idxs in groups.items():
                if len(idxs) > 1:
                    idxs_sorted = sorted(idxs, key=lambda i: reads[i].start)
                    for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                        if reads[a].end < reads[b].start:
                            ya_center = reads_start_y + stacks[a] * (read_height + 1) + read_height // 2
                            yb_center = reads_start_y + stacks[b] * (read_height + 1) + read_height // 2
                            xa = margin + int((reads[a].end - start) / bp_per_px)
                            xb = margin + int((reads[b].start - start) / bp_per_px)
                            xa = max(margin, min(width - margin - 1, xa))
                            xb = max(margin, min(width - margin - 1, xb))
                            dr.line([(xa, ya_center), (xb, yb_center)], fill=(176, 196, 222), width=1)
        else:
            for qname, idxs in groups.items():
                if len(idxs) > 1:
                    idxs_sorted = sorted(idxs, key=lambda i: (reads[i].start, reads[i].end))
                    for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                        ya_center = reads_start_y + stacks[a] * (read_height + 1) + read_height // 2
                        yb_center = reads_start_y + stacks[b] * (read_height + 1) + read_height // 2
                        xa = margin + int((reads[a].end - start) / bp_per_px)
                        xb = margin + int((reads[b].start - start) / bp_per_px)
                        xa = max(margin, min(width - margin - 1, xa))
                        xb = max(margin, min(width - margin - 1, xb))
                        dr.line([(xa, ya_center), (xb, yb_center)], fill=(80, 160, 80), width=1)

        current_y = reads_start_y + reads_area_height + 5
    
    return img
