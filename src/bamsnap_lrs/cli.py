import argparse
import os


def add_common_args(parser):
    """Add common arguments to both dna and rna parsers"""
    parser.add_argument("--bam", required=True, nargs='+', action='extend', help="BAM file path(s)")
    parser.add_argument("--pos", required=True, help="Genomic position, format: chr:start-end or chr:pos")
    parser.add_argument("--out", required=True, help="Output file path (supports .png, .svg, .pdf)")
    parser.add_argument("--max-reads", type=int, default=300, help="Maximum number of reads to display, [300]")
    parser.add_argument("--mapq", type=int, default=0, help="Minimum MAPQ value, [0]")
    parser.add_argument("--show-supp", action="store_true", help="Show supplementary alignments")
    parser.add_argument("--show-secondary", action="store_true", help="Show secondary alignments")
    parser.add_argument("--width", type=int, default=1200, help="Image width (pixels), [1200]")
    parser.add_argument("--read-height", type=int, default=6, help="Height of each read (pixels)")
    parser.add_argument("--detail", choices=["low", "mid", "high"], default="mid", help="Detail level, [mid]")
    parser.add_argument("--downsample-strategy", choices=["mapq", "first"], default="mapq", help="Downsampling strategy, [mapq]")
    parser.add_argument("--use-md", action="store_true", help="Use MD tag to detect mismatches")
    parser.add_argument("--use-cs", action="store_true", help="Use cs tag to detect mismatches")
    parser.add_argument("--fa", help="Reference genome FASTA file path")
    parser.add_argument("--show-axis", action="store_true", help="Show coordinate axis")
    parser.add_argument("--show-composition", action="store_true", help="Show base composition chart")
    parser.add_argument("--comp-height", type=int, default=None, help="Base composition chart height, [None]")
    parser.add_argument("--style", choices=["default", "jbrowse"], default="jbrowse", help="Rendering style, [jbrowse]")
    parser.add_argument("--color-by", choices=["type", "base", "strand", "mapq"], default="type", help="Coloring method, [type]")
    parser.add_argument("--comp-max-depth", type=int, help="Base composition chart maximum depth")
    parser.add_argument("--show-coverage", action="store_true", default=True, help="Show coverage track, [True]")
    parser.add_argument("--no-coverage", dest="show_coverage", action="store_false", help="Hide coverage track")
    parser.add_argument("--coverage-height", type=int, default=15, help="Coverage track height, [15]")
    parser.add_argument("--track-title", type=str, nargs='+', action='extend', help="Track title(s)")
    parser.add_argument("--show-insertion-labels", action="store_true", default=True, help="Show insertion labels, [True]")
    parser.add_argument("--no-insertion-labels", dest="show_insertion_labels", action="store_false", help="Hide insertion labels")
    parser.add_argument("--coverage-max-depth", type=int, help="Coverage track maximum depth")


def render_output(tracks, args, chrom, start, end, ref_seq, is_rna=False):
    """Common rendering logic for both DNA and RNA"""
    # Determine format based on output file extension
    output_format = None
    if args.out.lower().endswith('.svg'):
        output_format = 'svg'
    elif args.out.lower().endswith('.pdf'):
        output_format = 'pdf'
    
    if output_format == 'svg':
        from .svg_renderer import render_svg_snapshot
        svg_content = render_svg_snapshot(
            tracks,
            chrom,
            start,
            end,
            width=args.width,
            read_height=args.read_height,
            detail=args.detail,
            show_axis=args.show_axis,
            show_coverage=args.show_coverage,
            coverage_height=args.coverage_height,
            style=args.style,
            color_by=args.color_by,
            ref_seq=ref_seq,
            show_insertion_labels=args.show_insertion_labels,
            coverage_max_depth=args.coverage_max_depth,
            is_rna=is_rna,
        )
        with open(args.out, 'w') as f:
            f.write(svg_content)
    elif output_format == 'pdf':
        # Use cairosvg to convert SVG to PDF
        try:
            import cairosvg
        except ImportError:
            print("Error: PDF output requires cairosvg package. Please install it first:")
            print("  pip install cairosvg")
            print("  or")
            print("  conda install -c conda-forge cairosvg")
            return
        
        # First generate SVG content
        from .svg_renderer import render_svg_snapshot
        svg_content = render_svg_snapshot(
            tracks,
            chrom,
            start,
            end,
            width=args.width,
            read_height=args.read_height,
            detail=args.detail,
            show_axis=args.show_axis,
            show_coverage=args.show_coverage,
            coverage_height=args.coverage_height,
            style=args.style,
            color_by=args.color_by,
            ref_seq=ref_seq,
            show_insertion_labels=args.show_insertion_labels,
            coverage_max_depth=args.coverage_max_depth,
            is_rna=is_rna,
        )
        # Convert to PDF
        cairosvg.svg2pdf(bytestring=svg_content.encode('utf-8'), write_to=args.out)
    else:
        # Default PNG or other formats
        from .png_renderer import render_snapshot
        # For PNG, we currently assume single track or need to update renderer.py
        # Assuming renderer.py will be updated to take tracks
        img = render_snapshot(
            tracks,
            chrom,
            start,
            end,
            width=args.width,
            read_height=args.read_height,
            detail=args.detail,
            show_axis=args.show_axis,
            show_composition=args.show_composition,
            composition_height=args.comp_height,
            style=args.style,
            color_by=args.color_by,
            ref_seq=ref_seq,
            comp_max_depth=args.comp_max_depth,
            show_coverage=args.show_coverage,
            coverage_height=args.coverage_height,
            show_insertion_labels=args.show_insertion_labels,
            coverage_max_depth=args.coverage_max_depth,
            is_rna=is_rna,
        )
        img.save(args.out)


def main():
    p = argparse.ArgumentParser(description="Bamsnap-LRS: Long-read sequencing data visualization tool")
    sub = p.add_subparsers(dest="cmd")
    
    # dna command - generate genomic region snapshot for DNA data
    dna_parser = sub.add_parser("dna", help="Generate read pileup snapshot for DNA sequencing data")
    add_common_args(dna_parser)
    
    # rna command - generate genomic region snapshot for RNA data
    rna_parser = sub.add_parser("rna", help="Generate read pileup snapshot for RNA sequencing data (handles spliced alignments)")
    add_common_args(rna_parser)
    
    args = p.parse_args()
    
    if args.cmd is None:
        p.print_help()
        return
    
    # Parse position
    chrom, coords = args.pos.split(":")
    if "-" in coords:
        s, e = coords.split("-")
        start = int(s)
        end = int(e)
    else:
        pos = int(coords)
        start = max(0, pos - 250)
        end = pos + 250
    
    # Fetch reads based on command type
    tracks = []
    titles = args.track_title if args.track_title else []
    if len(titles) < len(args.bam):
        # Extend titles if not enough provided
        for i in range(len(titles), len(args.bam)):
            bam_path = args.bam[i]
            # Use filename without extension as default title
            default_title = os.path.splitext(os.path.basename(bam_path))[0]
            titles.append(default_title)

    if args.cmd == "dna":
        from .reader import fetch_reads
        from .ref import get_ref_subseq
        
        for i, bam_path in enumerate(args.bam):
            reads = fetch_reads(
                bam_path,
                chrom,
                start,
                end,
                max_reads=args.max_reads,
                mapq_min=args.mapq,
                show_supp=args.show_supp,
                show_secondary=args.show_secondary,
                downsample_strategy=args.downsample_strategy,
                use_md=args.use_md,
                use_cs=args.use_cs,
                use_ref=bool(args.fa),
                fa_path=args.fa,
            )
            tracks.append({"reads": reads, "title": titles[i]})
            
        ref_seq = None
        if args.fa:
            ref_seq = get_ref_subseq(args.fa, chrom, start, end)
        render_output(tracks, args, chrom, start, end, ref_seq, is_rna=False)
    
    elif args.cmd == "rna":
        from .reader import fetch_rna_reads
        from .ref import get_ref_subseq
        
        for i, bam_path in enumerate(args.bam):
            reads = fetch_rna_reads(
                bam_path,
                chrom,
                start,
                end,
                max_reads=args.max_reads,
                mapq_min=args.mapq,
                show_supp=args.show_supp,
                show_secondary=args.show_secondary,
                downsample_strategy=args.downsample_strategy,
                use_md=args.use_md,
                use_cs=args.use_cs,
                use_ref=bool(args.fa),
                fa_path=args.fa,
            )
            tracks.append({"reads": reads, "title": titles[i]})
            
        ref_seq = None
        if args.fa:
            ref_seq = get_ref_subseq(args.fa, chrom, start, end)
        render_output(tracks, args, chrom, start, end, ref_seq, is_rna=True)


if __name__ == "__main__":
    main()
