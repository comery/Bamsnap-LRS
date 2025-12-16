"""SVG renderer for bamsnap visualization"""
from typing import List, Dict, Optional
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom

from .layout import segments_to_pixels, assign_stacks
from .reader import Read
from .pileup import base_pileup, pileup_to_pixels
from .styles import color_for_type, color_for_base, shade_by_mapq, STRAND_COLORS, MISMATCH_COLORS


def rgb_to_hex(rgb: tuple) -> str:
    """Convert RGB tuple to hex color"""
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def render_svg_snapshot(
    reads: List[Read],
    chrom: str,
    start: int,
    end: int,
    width: int = 1200,
    read_height: int = 6,
    detail: str = "mid",
    show_axis: bool = True,
    show_coverage: bool = True,
    coverage_height: int = 15,
    track_title: str = "Reads",
    style: str = "jbrowse",
    color_by: str = "type",
    ref_seq: Optional[str] = None,
    show_insertion_labels: bool = True,
    coverage_max_depth: Optional[int] = None,
    margin: int = 20,  # Left and right margins
) -> str:
    """Render snapshot as SVG string"""
    # Actual drawing area width (excluding left and right margins)
    content_width = width - 2 * margin
    bp_per_px = (end - start) / float(content_width)
    spans = [(max(r.start, start), min(r.end, end)) for r in reads]
    stacks = assign_stacks(spans, max_stack=max(1, len(reads)))
    
    # Calculate total height
    top = 0
    if show_axis:
        top += 20
    
    coverage_track_height = 0
    if show_coverage:
        coverage_track_height = coverage_height + 15
    
    track_header_height = 15
    reads_area_height = (max(stacks) + 1) * (read_height + 1) + 5
    
    total_height = top + coverage_track_height + track_header_height + reads_area_height + 5
    
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
        # Axis line position adjusted to bottom
        axis_y = current_y + 18
        axis_line = SubElement(svg, "line", {
            "x1": str(margin), "y1": str(axis_y),
            "x2": str(width - margin), "y2": str(axis_y),
            "stroke": "black",
            "stroke-width": "1"
        })
        span = end - start
        step = max(1, span // 10)
        for pos in range(start, end + 1, step):
            x = margin + int((pos - start) / bp_per_px)
            tick_positions.append(x)
            # Tick mark
            tick = SubElement(svg, "line", {
                "x1": str(x), "y1": str(axis_y - 3),
                "x2": str(x), "y2": str(axis_y + 3),
                "stroke": "black",
                "stroke-width": "1"
            })
            # Text moved up, separated from axis
            text = SubElement(svg, "text", {
                "x": str(x + 2),
                "y": str(current_y + 10),
                "font-size": "9",
                "fill": "black"
            })
            text.text = f"{pos:,}"
        current_y += 22
    
    # Draw vertical guide lines (dashed) - from axis to bottom of image
    for tick_x in tick_positions:
        guide_line = SubElement(svg, "line", {
            "x1": str(tick_x), "y1": str(current_y),
            "x2": str(tick_x), "y2": str(total_height),
            "stroke": "#dcdcdc",
            "stroke-width": "1",
            "stroke-dasharray": "3,3"  # Dashed: 3px line, 3px gap
        })
    
    # Draw coverage track
    if show_coverage:
        # Header bar
        header_rect = SubElement(svg, "rect", {
            "x": "0", "y": str(current_y),
            "width": str(width),
            "height": "15",
            "fill": "#f5f5f5"
        })
        header_line = SubElement(svg, "line", {
            "x1": "0", "y1": str(current_y + 14),
            "x2": str(width), "y2": str(current_y + 14),
            "stroke": "#c8c8c8",
            "stroke-width": "1"
        })
        header_text = SubElement(svg, "text", {
            "x": "5",
            "y": str(current_y + 11),
            "font-size": "11",
            "fill": "black"
        })
        header_text.text = f"{track_title} - Coverage"
        current_y += 15
        
        # Calculate base distribution (original resolution, independent for each base position)
        pile = base_pileup(reads, start, end)
        num_bases = len(pile)
        
        # Calculate ref_match and variant bases for each position
        base_distribution = []
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
            
            base_distribution.append({
                "ref_match": ref_match,
                "A": variant_counts["A"],
                "C": variant_counts["C"],
                "G": variant_counts["G"],
                "T": variant_counts["T"],
                "N": variant_counts["N"],
                "depth": p.get("depth", 0)
            })
        
        # Draw coverage stacked bar chart
        max_cov = coverage_max_depth
        if max_cov is None:
            max_cov = max((b["depth"] for b in base_distribution), default=1)
        if max_cov <= 0:
            max_cov = 1
        
        # Draw y-axis scale
        bar_bottom = current_y + coverage_height
        bar_top = current_y
        axis_x = margin - 2
        
        # Y-axis line
        SubElement(svg, "line", {
            "x1": str(axis_x), "y1": str(bar_top),
            "x2": str(axis_x), "y2": str(bar_bottom),
            "stroke": "#646464",
            "stroke-width": "1"
        })
        
        # Maximum value tick (top)
        SubElement(svg, "line", {
            "x1": str(axis_x - 3), "y1": str(bar_top),
            "x2": str(axis_x), "y2": str(bar_top),
            "stroke": "#646464",
            "stroke-width": "1"
        })
        max_text = SubElement(svg, "text", {
            "x": "2",
            "y": str(bar_top + 6),
            "font-size": "8",
            "fill": "#505050"
        })
        max_text.text = str(max_cov)
        
        # Middle tick
        mid_y = (bar_top + bar_bottom) // 2
        mid_val = max_cov // 2
        SubElement(svg, "line", {
            "x1": str(axis_x - 3), "y1": str(mid_y),
            "x2": str(axis_x), "y2": str(mid_y),
            "stroke": "#646464",
            "stroke-width": "1"
        })
        mid_text = SubElement(svg, "text", {
            "x": "2",
            "y": str(mid_y + 3),
            "font-size": "8",
            "fill": "#505050"
        })
        mid_text.text = str(mid_val)
        
        # Bottom tick (0)
        SubElement(svg, "line", {
            "x1": str(axis_x - 3), "y1": str(bar_bottom),
            "x2": str(axis_x), "y2": str(bar_bottom),
            "stroke": "#646464",
            "stroke-width": "1"
        })
        zero_text = SubElement(svg, "text", {
            "x": "2",
            "y": str(bar_bottom),
            "font-size": "8",
            "fill": "#505050"
        })
        zero_text.text = "0"
        
        # Draw by pixel position, each pixel maps to corresponding base position
        for px in range(content_width):
            # Calculate base position index for this pixel
            base_idx = int(px * num_bases / content_width)
            if base_idx >= num_bases:
                base_idx = num_bases - 1
            
            dist = base_distribution[base_idx]
            depth = dist["depth"]
            if depth == 0:
                continue
            
            # Actual drawing position (add margin)
            draw_x = margin + px
            
            # Calculate total bar height for this position (proportional to depth)
            bar_height = int(coverage_height * min(depth / max_cov, 1.0))
            if bar_height <= 0:
                continue
            
            # Check if using aggregated data with reference information
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
                bar_bottom = current_y + coverage_height
                current_stack_y = bar_bottom
                
                # First draw variant portion (in fixed order A→C→G→T→N, with corresponding colors)
                for base in ["A", "C", "G", "T", "N"]:
                    if base not in base_heights or base_heights[base] <= 0:
                        continue
                    h = base_heights[base]
                    color = rgb_to_hex(MISMATCH_COLORS.get(base, (160, 160, 160)))
                    rect = SubElement(svg, "rect", {
                        "x": str(draw_x),
                        "y": str(current_stack_y - h),
                        "width": "1",
                        "height": str(h),
                        "fill": color
                    })
                    current_stack_y -= h
                
                # Finally draw reference match portion (gray, on top)
                if "ref" in base_heights and base_heights["ref"] > 0:
                    h = base_heights["ref"]
                    rect = SubElement(svg, "rect", {
                        "x": str(draw_x),
                        "y": str(current_stack_y - h),
                        "width": "1",
                        "height": str(h),
                        "fill": "#b4b4b4"
                    })
                    current_stack_y -= h
            else:
                # Old logic (used when no reference sequence)
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
                bar_bottom = current_y + coverage_height
                current_stack_y = bar_bottom
                
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
                    
                    color = rgb_to_hex(MISMATCH_COLORS.get(base, (160, 160, 160)))
                    rect = SubElement(svg, "rect", {
                        "x": str(draw_x),
                        "y": str(current_stack_y - h),
                        "width": "1",
                        "height": str(h),
                        "fill": color
                    })
                    current_stack_y -= h
                
                # Finally draw reference base (gray, on top)
                if ref_base and ref_base in base_heights and base_heights[ref_base] > 0:
                    h = base_heights[ref_base]
                    rect = SubElement(svg, "rect", {
                        "x": str(draw_x),
                        "y": str(current_stack_y - h),
                        "width": "1",
                        "height": str(h),
                        "fill": "#b4b4b4"
                    })
        
        current_y += coverage_height + 3
    
    # Draw read track header
    header_rect = SubElement(svg, "rect", {
        "x": "0", "y": str(current_y),
        "width": str(width),
        "height": "15",
        "fill": "#f5f5f5"
    })
    header_line = SubElement(svg, "line", {
        "x1": "0", "y1": str(current_y + 14),
        "x2": str(width), "y2": str(current_y + 14),
        "stroke": "#c8c8c8",
        "stroke-width": "1"
    })
    header_text = SubElement(svg, "text", {
        "x": "5",
        "y": str(current_y + 11),
        "font-size": "11",
        "fill": "black"
    })
    header_text.text = track_title
    current_y += 15
    
    # Draw reads
    groups: Dict[str, List[int]] = {}
    for i, r in enumerate(reads):
        groups.setdefault(r.qname, []).append(i)
    
    # Read main color (for arrows)
    read_main_color = (180, 180, 180)
    
    for idx, r in enumerate(reads):
        y = current_y + stacks[idx] * (read_height + 1)
        rects = segments_to_pixels(r.segments, r.start, start, bp_per_px, detail=detail)
        
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
            
            if color_by == "type":
                color = color_for_type(t)
            elif color_by == "mapq":
                base_color = color_for_type(t)
                color = shade_by_mapq(base_color, r.mapq)
            else:
                color = color_for_type(t)
            
            if t == "ins":
                ins_color = (128, 0, 128)
                line = SubElement(svg, "line", {
                    "x1": str(x0_draw), "y1": str(y - 1),
                    "x2": str(x0_draw), "y2": str(y + read_height + 1),
                    "stroke": rgb_to_hex(ins_color),
                    "stroke-width": "2"
                })
                if show_insertion_labels and rect_idx in rect_to_seg:
                    seg = r.segments[rect_to_seg[rect_idx]]
                    if seg and seg.length > 0:
                        text = SubElement(svg, "text", {
                            "x": str(x0_draw + 2),
                            "y": str(y - 2),
                            "font-size": "8",
                            "fill": rgb_to_hex(ins_color)
                        })
                        text.text = f"I({seg.length})"
            elif t == "del" or t == "ref_skip":
                # Deletion shown as dark gray rectangle (same width as normal read)
                del_color = (80, 80, 80)  # Dark gray
                rect = SubElement(svg, "rect", {
                    "x": str(x0_draw), "y": str(y),
                    "width": str(x1_draw - x0_draw),
                    "height": str(read_height),
                    "fill": rgb_to_hex(del_color)
                })
                # Label deletion length, color based on width
                if show_insertion_labels and rect_idx in rect_to_seg:
                    seg = r.segments[rect_to_seg[rect_idx]]
                    if seg and seg.length > 0:
                        label = str(seg.length)
                        del_width = x1_draw - x0_draw
                        if del_width > 4:
                            # Width > 4, white text inside rectangle
                            text_x = (x0_draw + x1_draw) // 2
                            text = SubElement(svg, "text", {
                                "x": str(text_x),
                                "y": str(y + read_height - 1),
                                "font-size": "8",
                                "fill": "white",
                                "text-anchor": "middle"
                            })
                            text.text = label
                        else:
                            # Width <= 4, black text above rectangle
                            text = SubElement(svg, "text", {
                                "x": str(x0_draw),
                                "y": str(y - 2),
                                "font-size": "8",
                                "fill": "black"
                            })
                            text.text = label
            elif t == "mismatch":
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
                                rect = SubElement(svg, "rect", {
                                    "x": str(x0_draw), "y": str(y),
                                    "width": str(x1_draw - x0_draw),
                                    "height": str(read_height),
                                    "fill": rgb_to_hex(mismatch_color)
                                })
                            break
                else:
                    rect = SubElement(svg, "rect", {
                        "x": str(x0_draw), "y": str(y),
                        "width": str(x1_draw - x0_draw),
                        "height": str(read_height),
                        "fill": rgb_to_hex(color)
                    })
            else:
                if t == "match":
                    color = (180, 180, 180)
                rect = SubElement(svg, "rect", {
                    "x": str(x0_draw), "y": str(y),
                    "width": str(x1_draw - x0_draw),
                    "height": str(read_height),
                    "fill": rgb_to_hex(color)
                })
        
        # Draw direction arrow (only at read end)
        if style == "jbrowse":
            # Find read's actual end within visible area
            read_end_px = margin + int((min(r.end, end) - start) / bp_per_px)
            read_start_px = margin + int((max(r.start, start) - start) / bp_per_px)
            
            # Ensure within visible area
            read_end_px = max(margin, min(width - margin - 1, read_end_px))
            read_start_px = max(margin, min(width - margin - 1, read_start_px))
            
            # Arrow size (not exceeding half of read height)
            head = max(3, min(read_height // 2, 5))
            
            # Arrow colors
            arrow_fill = "#d3d3d3"       # True end: solid fill
            arrow_stroke = "#c8c8c8"     # Truncated: border color (lighter)
            
            if r.reverse:
                # Reverse read, arrow at start position (left), pointing left
                arrow_x = read_start_px
                # Check if truncated (read actual start position outside visible area)
                is_truncated = r.start < start
                
                points = f"{arrow_x},{y} {arrow_x - head},{y + read_height // 2} {arrow_x},{y + read_height}"
                if is_truncated:
                    # Truncated: hollow triangle (white fill + thin gray border)
                    polygon = SubElement(svg, "polygon", {
                        "points": points,
                        "fill": "white",
                        "stroke": arrow_stroke,
                        "stroke-width": "0.5"
                    })
                else:
                    # True end: solid fill
                    polygon = SubElement(svg, "polygon", {
                        "points": points,
                        "fill": arrow_fill
                    })
            else:
                # Forward read, arrow at end position (right), pointing right
                arrow_x = read_end_px
                # Check if truncated (read actual end position outside visible area)
                is_truncated = r.end > end
                
                points = f"{arrow_x},{y} {arrow_x + head},{y + read_height // 2} {arrow_x},{y + read_height}"
                if is_truncated:
                    # Truncated: hollow triangle (white fill + thin gray border)
                    polygon = SubElement(svg, "polygon", {
                        "points": points,
                        "fill": "white",
                        "stroke": arrow_stroke,
                        "stroke-width": "0.5"
                    })
                else:
                    # True end: solid fill
                    polygon = SubElement(svg, "polygon", {
                        "points": points,
                        "fill": arrow_fill
                    })
    
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
                line = SubElement(svg, "line", {
                    "x1": str(xa), "y1": str(ya + read_height // 2),
                    "x2": str(xb), "y2": str(yb + read_height // 2),
                    "stroke": "#50a050",
                    "stroke-width": "1"
                })
    
    # Convert to string
    rough_string = tostring(svg, encoding='unicode')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
