"""SVG renderer for bamsnap visualization"""
from typing import List, Dict, Optional, Any
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom

from .layout import segments_to_pixels, assign_stacks
from .reader import Read
from .pileup import base_pileup
from .styles import color_for_type, color_for_base, shade_by_mapq, STRAND_COLORS, MISMATCH_COLORS


def rgb_to_hex(rgb: tuple) -> str:
    """Convert RGB tuple to hex color"""
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


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
    style: str = "jbrowse",
    color_by: str = "type",
    ref_seq: Optional[str] = None,
    show_insertion_labels: bool = True,
    coverage_max_depth: Optional[int] = None,
    margin: int = 20,  # Left and right margins
    is_rna: bool = False,
    **kwargs
) -> str:
    """Render snapshot as SVG string"""
    # Actual drawing area width (excluding left and right margins)
    content_width = width - 2 * margin
    bp_per_px = (end - start) / float(content_width)
    
    # Calculate total height
    top = 0
    if show_axis:
        top += 22  # Axis area height
    
    total_height = top
    
    # Pre-calculate stacks and height for each track
    track_meta = []
    
    # Calculate aggregate coverage height if needed
    aggregate_coverage_h = 0
    if show_coverage:
        # Header (15) + Padding for arcs (coverage_height) + Coverage track (coverage_height) + Bottom margin (15)
        aggregate_coverage_h = 15 + coverage_height + coverage_height + 15
        total_height += aggregate_coverage_h
    
    for track in tracks:
        reads = track['reads']
        spans = [(max(r.start, start), min(r.end, end)) for r in reads]
        stacks = assign_stacks(spans, max_stack=max(1, len(reads)))
        
        track_header_h = 15
        reads_area_h = (max(stacks) + 1) * (read_height + 1) + 5
        
        track_h = track_header_h + reads_area_h + 5
        total_height += track_h
        
        track_meta.append({
            'stacks': stacks,
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
        # Axis line position adjusted to bottom
        axis_y = current_y + 18
        SubElement(svg, "line", {
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
            SubElement(svg, "line", {
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
        SubElement(svg, "line", {
            "x1": str(tick_x), "y1": str(current_y),
            "x2": str(tick_x), "y2": str(total_height),
            "stroke": "#dcdcdc",
            "stroke-width": "1",
            "stroke-dasharray": "3,3"  # Dashed: 3px line, 3px gap
        })
    
    # Draw aggregate coverage track if enabled
    if show_coverage:
        track_y_start = current_y
        all_reads = []
        for track in tracks:
            all_reads.extend(track['reads'])
            
        # Header bar
        SubElement(svg, "rect", {
            "x": "0", "y": str(current_y),
            "width": str(width),
            "height": "15",
            "fill": "#f5f5f5"
        })
        SubElement(svg, "line", {
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
        header_text.text = "Aggregate Coverage"
        current_y += 15
        
        # Padding for arcs
        padding_for_arcs = coverage_height
        current_y += padding_for_arcs
        
        # Calculate base distribution for all reads
        pile = base_pileup(all_reads, start, end)
        num_bases = len(pile)
        
        # Calculate ref_match and variant bases for each position
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
        
        # Draw by pixel position
        import math
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
            
            base_heights = {}
            total_count = agg_ref_match
            for base in ["A", "C", "G", "T", "N"]:
                total_count += agg_variants[base]
            
            if total_count == 0:
                continue
            
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
            
            # If detail is low, just draw a simple gray depth bar
            if detail == "low":
                SubElement(svg, "rect", {
                    "x": str(draw_x), "y": str(current_stack_y - bar_height),
                    "width": "1", "height": str(bar_height), "fill": "#b4b4b4"
                })
                continue

            for base in ["A", "C", "G", "T", "N"]:
                if base not in base_heights or base_heights[base] <= 0:
                    continue
                h = base_heights[base]
                color = rgb_to_hex(MISMATCH_COLORS.get(base, (160, 160, 160)))
                SubElement(svg, "rect", {
                    "x": str(draw_x), "y": str(current_stack_y - h),
                    "width": "1", "height": str(h), "fill": color
                })
                current_stack_y -= h
            
            if "ref" in base_heights and base_heights["ref"] > 0:
                h = base_heights["ref"]
                SubElement(svg, "rect", {
                    "x": str(draw_x), "y": str(current_stack_y - h),
                    "width": "1", "height": str(h), "fill": "#b4b4b4"
                })
        
        # Arcs for RNA mode in aggregate coverage
        if is_rna:
            arc_color_hex = "#fdd1d3"
            arc_anchor_y = track_y_start + 15 + coverage_height + coverage_height // 2
            
            # 1. Internal arcs (ref_skip) from all reads
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
                        
                        w = xb - xa
                        h = min(w // 2, coverage_height)
                        mid_x = (xa + xb) / 2
                        ctrl_y = arc_anchor_y - 2 * h
                        path_d = f"M {xa} {arc_anchor_y} Q {mid_x} {ctrl_y} {xb} {arc_anchor_y}"
                        SubElement(svg, "path", {
                            "d": path_d, "fill": "none", "stroke": arc_color_hex,
                            "stroke-opacity": "0.6", "stroke-width": "0.5"
                        })
                    ref_cursor += seg.ref_consumed
            
            # 2. Split read arcs from all reads
            groups_temp: Dict[str, List[int]] = {}
            for idx_r, r in enumerate(all_reads):
                groups_temp.setdefault(r.qname, []).append(idx_r)
            
            for qname, idxs in groups_temp.items():
                if len(idxs) > 1:
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
                            mid_x = (xa + xb) / 2
                            ctrl_y = arc_anchor_y - 2 * h
                            path_d = f"M {xa} {arc_anchor_y} Q {mid_x} {ctrl_y} {xb} {arc_anchor_y}"
                            SubElement(svg, "path", {
                                "d": path_d, "fill": "none", "stroke": arc_color_hex,
                                "stroke-opacity": "0.6", "stroke-width": "1"
                            })
        
        current_y = track_y_start + aggregate_coverage_h
    
    # Iterate over tracks
    for i, track in enumerate(tracks):
        reads = track['reads']
        track_title = track.get('title', 'Reads')
        meta = track_meta[i]
        stacks = meta['stacks']
        reads_area_h = meta['reads_area_h']
        
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
        
        for idx, r in enumerate(reads):
            y = reads_start_y + stacks[idx] * (read_height + 1)
            rects = segments_to_pixels(r.segments, r.start, start, bp_per_px, detail=detail)
            
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
                    SubElement(svg, "line", {
                        "x1": str(x0_draw), "y1": str(y),
                        "x2": str(x0_draw), "y2": str(y + read_height),
                        "stroke": color_hex, "stroke-width": "1"
                    })
                    if detail == "high" and show_insertion_labels:
                        seg_idx = rect_to_seg.get(rect_idx)
                        if seg_idx is not None:
                            seg = r.segments[seg_idx]
                            if seg and seg.length > 0:
                                SubElement(svg, "text", {
                                    "x": str(x0_draw), "y": str(y - 2),
                                    "font-size": "8", "fill": "#800080",
                                    "text-anchor": "middle"
                                }).text = f"I({seg.length})"
                elif t == "ref_skip":
                    y_center = y + read_height // 2
                    SubElement(svg, "line", {
                        "x1": str(x0_draw), "y1": str(y_center),
                        "x2": str(x1_draw), "y2": str(y_center),
                        "stroke": "#b0c4de", "stroke-width": "1"
                    })
                elif t == "del":
                    SubElement(svg, "rect", {
                        "x": str(x0_draw), "y": str(y),
                        "width": str(x1_draw - x0_draw), "height": str(read_height),
                        "fill": "#808080", "stroke": "none", "stroke-width": "1"
                    })
                    if detail == "high" and show_insertion_labels:
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
                    current_color_hex = color_hex
                    if detail == "low":
                        current_color_hex = "#c3c3c3"
                    elif r.seq and rect_idx in rect_to_seg:
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
                    if t == "match":
                        color_hex = "#c3c3c3"
                    SubElement(svg, "rect", {
                        "x": str(x0_draw), "y": str(y),
                        "width": str(x1_draw - x0_draw), "height": str(read_height),
                        "fill": color_hex, "stroke": "none"
                    })

            # Arrow head
            if style == "jbrowse":
                w = x1_draw - x0_draw
                if w >= 6:
                    head = max(3, min(6, w // 6))
                    arrow_pts = ""
                    if r.reverse:
                        p1 = f"{x0_draw},{y}"
                        p2 = f"{x0_draw+head},{y+read_height/2}"
                        p3 = f"{x0_draw},{y+read_height}"
                        arrow_pts = f"{p1} {p2} {p3}"
                    else:
                        p1 = f"{x1_draw},{y}"
                        p2 = f"{x1_draw-head},{y+read_height/2}"
                        p3 = f"{x1_draw},{y+read_height}"
                        arrow_pts = f"{p1} {p2} {p3}"
                    SubElement(svg, "polygon", {
                        "points": arrow_pts, "fill": "#646464"
                    })

        # Connect segments
        if is_rna:
            for qname, idxs in groups.items():
                if len(idxs) > 1:
                    idxs_sorted = sorted(idxs, key=lambda i: reads[i].start)
                    for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                        if reads[a].end < reads[b].start:
                            ya = reads_start_y + stacks[a] * (read_height + 1)
                            yb = reads_start_y + stacks[b] * (read_height + 1)
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
            for qname, idxs in groups.items():
                if len(idxs) > 1:
                    idxs_sorted = sorted(idxs, key=lambda i: (reads[i].start, reads[i].end))
                    for a, b in zip(idxs_sorted, idxs_sorted[1:]):
                        ya = reads_start_y + stacks[a] * (read_height + 1)
                        yb = reads_start_y + stacks[b] * (read_height + 1)
                        xa = margin + int((reads[a].end - start) / bp_per_px)
                        xb = margin + int((reads[b].start - start) / bp_per_px)
                        xa = max(margin, min(width - margin - 1, xa))
                        xb = max(margin, min(width - margin - 1, xb))
                        SubElement(svg, "line", {
                            "x1": str(xa), "y1": str(ya + read_height // 2),
                            "x2": str(xb), "y2": str(yb + read_height // 2),
                            "stroke": "#50a050", "stroke-width": "1"
                        })

        current_y = reads_start_y + reads_area_h + 5

    # Convert to string
    rough_string = tostring(svg, encoding='unicode')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
