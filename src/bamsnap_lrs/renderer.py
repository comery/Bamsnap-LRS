from typing import List, Dict, Optional, Tuple

from PIL import Image, ImageDraw, ImageFont

from .layout import segments_to_pixels, assign_stacks
from .reader import Read
from .pileup import base_pileup, pileup_to_pixels
from .styles import color_for_type, color_for_base, shade_by_mapq, STRAND_COLORS, MISMATCH_COLORS


def render_snapshot(
    reads: List[Read],
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
):
    # If using jbrowse style and coverage enabled, use new rendering function
    if style == "jbrowse" and show_coverage:
        return render_jbrowse_style(
            reads=reads,
            chrom=chrom,
            start=start,
            end=end,
            width=width,
            read_height=read_height,
            detail=detail,
            show_axis=show_axis,
            show_coverage=show_coverage,
            coverage_height=coverage_height,
            track_title=track_title,
            style=style,
            color_by=color_by,
            ref_seq=ref_seq,
            show_insertion_labels=show_insertion_labels,
            coverage_max_depth=coverage_max_depth,
        )
    
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
                dr.line([(x0, y - 1), (x0, y + read_height + 1)], fill=color)
            elif t == "del" or t == "ref_skip":
                dr.rectangle([(x0, y), (x1, y + read_height)], outline=(120, 120, 120), fill=None)
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
    
    # Draw by pixel position, each pixel maps to corresponding base position
    for px in range(width):
        # Calculate base position index for this pixel
        base_idx = int(px * num_bases / width)
        if base_idx >= num_bases:
            base_idx = num_bases - 1
        
        dist = base_distribution[base_idx]
        depth = dist["depth"]
        if depth == 0:
            continue
        
        # Actual drawing position (add margin)
        draw_x = margin + px
        
        # Calculate total bar height for this position (proportional to depth)
        bar_height = int(height * min(depth / max_cov, 1.0))
        if bar_height <= 0:
            continue
        
        # Check if using data with reference information
        has_ref_match = "ref_match" in dist
        
        if has_ref_match:
            # Use pre-calculated reference match and variant data
            ref_match = dist.get("ref_match", 0)
            
            # Calculate height for each part
            base_heights = {}
            total_count = ref_match
            for base in ["A", "C", "G", "T", "N"]:
                total_count += dist.get(base, 0)
            
            if total_count == 0:
                continue
            
            # Reference match portion
            if ref_match > 0:
                base_heights["ref"] = int(bar_height * ref_match / total_count)
            
            # Variant portion
            for base in ["A", "C", "G", "T", "N"]:
                count = dist.get(base, 0)
                if count > 0:
                    base_heights[base] = int(bar_height * count / total_count)
            
            # Adjust heights to ensure full fill
            total_height_used = sum(base_heights.values())
            if total_height_used < bar_height and base_heights:
                max_base = max(base_heights.items(), key=lambda x: x[1])[0]
                base_heights[max_base] += (bar_height - total_height_used)
            
            # Draw from bottom to top
            bar_bottom = y + height
            current_y = bar_bottom
            
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
                current_y -= h
        else:
            # Old logic (used when no reference sequence)
            bp_per_px = (end - start) / float(width)
            ref_pos = start + int((x + 0.5) * bp_per_px)
            ref_base = None
            if ref_seq and ref_pos >= start and ref_pos < end:
                ref_idx = ref_pos - start
                if ref_idx < len(ref_seq):
                    ref_base = ref_seq[ref_idx].upper()
            
            # Calculate height for each base (proportional)
            base_heights = {}
            total_height_used = 0
            for base in ["A", "C", "G", "T", "N"]:
                count = dist.get(base, 0)
                if count > 0:
                    h = int(bar_height * count / depth)
                    base_heights[base] = h
                    total_height_used += h
            
            # Adjust heights to ensure full fill
            if total_height_used < bar_height and base_heights:
                max_base = max(base_heights.items(), key=lambda x: x[1])[0]
                base_heights[max_base] += (bar_height - total_height_used)
            
            # Bottom position
            bar_bottom = y + height
            current_y = bar_bottom
            
            # First draw variant bases (non-reference bases)
            for base in ["A", "C", "G", "T", "N"]:
                if base not in base_heights:
                    continue
                h = base_heights[base]
                if h <= 0:
                    continue
                
                # Skip reference base, draw later
                if ref_base and base == ref_base:
                    continue
                
                color = MISMATCH_COLORS.get(base, (160, 160, 160))
                dr.rectangle([(draw_x, current_y - h), (draw_x + 1, current_y)], fill=color)
                current_y -= h
            
            # Finally draw reference base (gray, on top)
            if ref_base and ref_base in base_heights and base_heights[ref_base] > 0:
                h = base_heights[ref_base]
                dr.rectangle([(draw_x, current_y - h), (draw_x + 1, current_y)], fill=(180, 180, 180))
    
    return height


def render_jbrowse_style(
    reads: List[Read],
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
) -> Image.Image:
    """Render JBrowse-style snapshot with track system and coverage chart"""
    # Actual drawing area width (excluding left and right margins)
    content_width = width - 2 * margin
    bp_per_px = (end - start) / float(content_width)
    spans = [(max(r.start, start), min(r.end, end)) for r in reads]
    stacks = assign_stacks(spans, max_stack=max(1, len(reads)))
    
    # Calculate total height (halved)
    top = 0
    if show_axis:
        top += 20  # Axis area (reduced)
    
    # Coverage track (reduced height)
    coverage_track_height = 0
    if show_coverage:
        coverage_track_height = coverage_height + 15  # 15 is header height (reduced)
    
    # Read track (reduced spacing, total height halved)
    track_header_height = 15  # Reduced header height
    reads_area_height = (max(stacks) + 1) * (read_height + 1) + 5  # Reduced read spacing
    
    total_height = top + coverage_track_height + track_header_height + reads_area_height + 5
    
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
    
    # Draw coverage track
    if show_coverage:
        current_y += draw_track_header(dr, f"{track_title} - Coverage", current_y, width, coverage_height)
        # Calculate base distribution (original resolution, independent for each base position)
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
            margin  # Pass margin
        )
        current_y += coverage_height + 3  # Reduced spacing
    
    # Draw read track header
    reads_y_start = current_y
    current_y += draw_track_header(dr, track_title, current_y, width, reads_area_height)
    
    # Draw reads
    groups: Dict[str, List[int]] = {}
    for i, r in enumerate(reads):
        groups.setdefault(r.qname, []).append(i)
    
    # Get read main color (for arrows)
    read_main_color = (180, 180, 180)  # Default gray
    
    for idx, r in enumerate(reads):
        y = current_y + stacks[idx] * (read_height + 1)  # Reduced spacing
        rects = segments_to_pixels(r.segments, r.start, start, bp_per_px, detail=detail)
        
        # Create rect index to segment mapping, for getting insertion length
        rect_to_seg = {}
        ref_cursor = r.start
        rect_idx = 0
        for seg_idx, seg in enumerate(r.segments):
            if seg.ref_consumed == 0:
                if seg.type == "ins":
                    x = int((ref_cursor - start) / bp_per_px)
                    if rect_idx < len(rects) and rects[rect_idx][0] == "ins" and rects[rect_idx][1] == x:
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
            # If insertion is outside the area, skip (don't clip to boundary)
            if t == "ins":
                if x0 < 0 or x0 > (width - 2 * margin):
                    continue
            
            # Limit to visible area (don't exceed x-axis)
            x0_draw = max(margin, min(width - margin, x0_draw))
            x1_draw = max(margin, min(width - margin, x1_draw))
            
            # If clipped width is 0, skip (except for insertion, which is a vertical line)
            if x1_draw <= x0_draw and t != "ins":
                continue
            
            # Unified type colors regardless of strand
            if color_by == "type":
                color = color_for_type(t)
            elif color_by == "mapq":
                base_color = color_for_type(t)
                color = shade_by_mapq(base_color, r.mapq)
            else:
                # Default use type colors (regardless of strand)
                color = color_for_type(t)
            
            if t == "ins":
                # Draw insertion vertical line
                ins_color = (128, 0, 128)  # Purple
                dr.line([(x0_draw, y - 1), (x0_draw, y + read_height + 1)], fill=ins_color, width=2)
                # Label insertion length
                if show_insertion_labels and rect_idx in rect_to_seg:
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
                        # Label length above insertion position (purple)
                        dr.text((x0_draw + 2, y - 12), label, fill=ins_color, font=font)
            elif t == "del" or t == "ref_skip":
                # Deletion shown as dark gray rectangle (same width as normal read)
                del_color = (80, 80, 80)  # Dark gray
                dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=del_color)
                # Label deletion length
                if show_insertion_labels and rect_idx in rect_to_seg:
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
                        # Determine text color based on deletion width
                        del_width = x1_draw - x0_draw
                        if del_width > 4:
                            text_color = (255, 255, 255)  # White
                            text_x = (x0_draw + x1_draw) // 2
                            dr.text((text_x - 5, y + 1), label, fill=text_color, font=font)
                        else:
                            text_color = (0, 0, 0)  # Black
                            # When width is small, text goes above deletion
                            dr.text((x0_draw, y - 10), label, fill=text_color, font=font)
            elif t == "mismatch":
                # For mismatches, use more vibrant colors
                if r.seq and detail == "high" and rect_idx in rect_to_seg:
                    seg_idx = rect_to_seg[rect_idx]
                    seg = r.segments[seg_idx]
                    ref_cursor = r.start
                    read_cursor = 0
                    for s_idx, s in enumerate(r.segments):
                        if s_idx < seg_idx:
                            ref_cursor += s.ref_consumed
                            read_cursor += s.read_consumed
                        elif s_idx == seg_idx:
                            if s.read_consumed > 0 and read_cursor < len(r.seq):
                                base = r.seq[read_cursor].upper()
                                mismatch_color = MISMATCH_COLORS.get(base, (200, 60, 60))
                                dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=mismatch_color)
                            break
                else:
                    dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=color)
            else:
                # Positions without variants (match) shown in gray
                if t == "match":
                    color = (180, 180, 180)  # Gray
                dr.rectangle([(x0_draw, y), (x1_draw, y + read_height)], fill=color)
        
        # Draw read direction arrow (only at read end)
        if style == "jbrowse" and rects:
            # Find read's actual end within visible area
            read_end_px = margin + int((min(r.end, end) - start) / bp_per_px)
            read_start_px = margin + int((max(r.start, start) - start) / bp_per_px)
            
            # Ensure within visible area
            read_end_px = max(margin, min(width - margin - 1, read_end_px))
            read_start_px = max(margin, min(width - margin - 1, read_start_px))
            
            # Arrow size (not exceeding half of read height)
            head = max(3, min(read_height // 2, 5))
            
            # Arrow colors
            arrow_fill = (211, 211, 211)      # True end: solid fill
            arrow_truncated = (200, 200, 200)  # Truncated: border color
            
            if r.reverse:
                # Reverse read, arrow at start position (left), pointing left
                arrow_x = read_start_px
                # Check if truncated (read actual start position outside visible area)
                is_truncated = r.start < start
                
                points = [
                    (arrow_x, y), 
                    (arrow_x - head, y + read_height // 2), 
                    (arrow_x, y + read_height)
                ]
                if is_truncated:
                    # Truncated: hollow triangle (border only)
                    dr.polygon(points, outline=arrow_truncated, fill=(255, 255, 255))
                else:
                    # True end: solid fill
                    dr.polygon(points, fill=arrow_fill)
            else:
                # Forward read, arrow at end position (right), pointing right
                arrow_x = read_end_px
                # Check if truncated (read actual end position outside visible area)
                is_truncated = r.end > end
                
                points = [
                    (arrow_x, y), 
                    (arrow_x + head, y + read_height // 2), 
                    (arrow_x, y + read_height)
                ]
                if is_truncated:
                    # Truncated: hollow triangle (border only)
                    dr.polygon(points, outline=arrow_truncated, fill=(255, 255, 255))
                else:
                    # True end: solid fill
                    dr.polygon(points, fill=arrow_fill)
    
    # Draw paired read connection lines
    for qname, idxs in groups.items():
        if len(idxs) > 1:
            idxs_sorted = sorted(idxs, key=lambda i: (reads[i].start, reads[i].end))
            for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                ya = current_y + stacks[a] * (read_height + 1)
                yb = current_y + stacks[b] * (read_height + 1)
                xa = margin + int((reads[a].end - start) / bp_per_px)
                xb = margin + int((reads[b].start - start) / bp_per_px)
                xa = max(margin, min(width - margin - 1, xa))
                xb = max(margin, min(width - margin - 1, xb))
                dr.line([(xa, ya + read_height // 2), (xb, yb + read_height // 2)], fill=(80, 160, 80), width=1)
    
    return img
