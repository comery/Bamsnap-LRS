"""PNG renderer compatibility wrapper.
The CLI now renders PNG by converting the SVG output with CairoSVG.  This file
is kept only for backward compatibility with old imports/tests that call
`render_snapshot(...)` directly.
"""
from io import BytesIO
from typing import Any, Dict, List, Optional, Tuple

from PIL import Image

from .highlight import HighlightSite
from .reader import Read


def _read_base_at_pos(read: Read, pos: int) -> Optional[str]:
    """Return the read base aligned to reference position `pos` (0-based)."""
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
    """Build a per-read signature over highlighted VCF sites."""
    out = []
    for p in site_positions:
        b = _read_base_at_pos(read, p)
        out.append(b if b in ("A", "C", "G", "T") else "?")
    return "".join(out)


def render_snapshot(
    tracks: List[Any],
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
    coverage_height: int = 15,
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
) -> Image.Image:
    """Return a PIL Image generated from the SVG renderer.

    The PNG image is produced from the SVG renderer, so SVG and PNG share the
    same drawing logic, including the base-composition track.
    """
    if tracks and not isinstance(tracks[0], dict):
        tracks = [{"reads": tracks, "title": track_title}]

    from .svg_renderer import render_svg_snapshot

    svg_content = render_svg_snapshot(
        tracks=tracks,
        chrom=chrom,
        start=start,
        end=end,
        width=width,
        read_height=read_height,
        detail=detail,
        show_axis=show_axis,
        show_coverage=show_coverage,
        coverage_height=coverage_height,
        show_composition=show_composition,
        composition_height=composition_height,
        comp_max_depth=comp_max_depth,
        style=style,
        color_by=color_by,
        ref_seq=ref_seq,
        show_insertion_labels=show_insertion_labels,
        coverage_max_depth=coverage_max_depth,
        is_rna=is_rna,
        gff_genes=gff_genes,
        bed_features=bed_features,
        highlight_sites=highlight_sites,
        highlight_samples=highlight_samples,
        no_hap_sort=no_hap_sort,
        no_hap_filter=no_hap_filter,
        focus_region=focus_region,
    )

    try:
        import cairosvg
    except ImportError as exc:
        raise ImportError(
            "PNG output through the SVG renderer requires cairosvg. "
            "Install it with `pip install cairosvg` or `conda install -c conda-forge cairosvg`."
        ) from exc

    png_bytes = cairosvg.svg2png(bytestring=svg_content.encode("utf-8"))
    return Image.open(BytesIO(png_bytes)).convert("RGB")
