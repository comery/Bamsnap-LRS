import argparse
import os
import sys


def add_common_args(parser):
    """Add common arguments to both dna and rna parsers"""
    parser.add_argument("--bam", required=True, nargs='+', action='extend', help="BAM/CRAM file path(s)")
    
    # Position input: either --pos (single region) or --regions (batch mode)
    position_group = parser.add_mutually_exclusive_group(required=True)
    position_group.add_argument("--pos", help="Genomic position, format: chr:start-end or chr:pos")
    position_group.add_argument("--regions", help="BED or VCF file containing regions to process (batch mode). VCF files support structural variants with END tag.")
    
    parser.add_argument("--out", help="Output file path (supports .png, .svg, .pdf). Required for single region mode.")
    parser.add_argument("--out-prefix", help="Output prefix for batch mode: outdir/outprefix_ (e.g., results/sample_.svg). Output files will be named: outprefix_chr_start_end.svg")
    parser.add_argument("--padding", type=int, default=None, help="Manually set padding added to both sides of each region from --regions. If omitted, padding is inferred from each BED/VCF record.")
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
    parser.add_argument("--fa", help="Reference genome FASTA file path (required for CRAM files)")
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
    annotation_group = parser.add_mutually_exclusive_group()
    annotation_group.add_argument("-g", "--gff", help="GFF/GTF file path for gene annotation")
    annotation_group.add_argument("--bed", help="BED file path for feature annotation")
    parser.add_argument("--hide-insertions",action="store_true",help="De-emphasize insertions in the read pileup. Insertions are rendered like normal matches and insertion labels are hidden")
    parser.add_argument("--show-insertion-labels", action="store_true", default=True, help="Show insertion labels, [True]")
    parser.add_argument("--no-insertion-labels", dest="show_insertion_labels", action="store_false", help="Hide insertion labels")
    parser.add_argument("--coverage-max-depth", type=int, help="Coverage track maximum depth")


def _render_svg_content(tracks, args, chrom, start, end, ref_seq, is_rna=False, gff_genes=None, bed_features=None,
                        highlight_sites=None, highlight_samples=None):
    """Render once with svg_renderer.py.

    This is the single source of truth for SVG / PDF / PNG output, so PNG will
    not drift away from SVG when drawing details are adjusted.
    """
    from .svg_renderer import render_svg_snapshot


    return render_svg_snapshot(
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
        show_composition=getattr(args, "show_composition", False),
        composition_height=getattr(args, "comp_height", None),
        comp_max_depth=getattr(args, "comp_max_depth", None),
        style=args.style,
        color_by=args.color_by,
        ref_seq=ref_seq,
        show_insertion_labels=args.show_insertion_labels,
        hide_indels=getattr(args, "hide_insertions", False),
        coverage_max_depth=args.coverage_max_depth,
        is_rna=is_rna,
        gff_genes=gff_genes,
        bed_features=bed_features,
        highlight_sites=highlight_sites,
        highlight_samples=highlight_samples,
        no_hap_sort=getattr(args, "no_hap_sort", False),
        no_hap_filter=getattr(args, "no_hap_filter", False),
        focus_region=getattr(args, "focus_region", None),
    )


def _write_svg_based_output(svg_content: str, output_path: str):
    """Write SVG directly, or convert the same SVG to PNG/PDF with CairoSVG."""
    ext = os.path.splitext(output_path)[1].lower()

    if ext == ".svg":
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(svg_content)
        return

    if ext in {".png", ".pdf"}:
        try:
            import cairosvg
        except ImportError:
            print(
                f"Error: {ext.upper()[1:]} output requires cairosvg package. Please install it first:\n"
                "  pip install cairosvg\n"
                "  or\n"
                "  conda install -c conda-forge cairosvg",
                file=sys.stderr,
            )
            return

        svg_bytes = svg_content.encode("utf-8")
        if ext == ".png":
            cairosvg.svg2png(bytestring=svg_bytes, write_to=output_path)
        else:
            cairosvg.svg2pdf(bytestring=svg_bytes, write_to=output_path)
        return

    raise ValueError(f"Unsupported output extension: {ext or '<none>'}. Use .svg, .png or .pdf")


def render_output(tracks, args, chrom, start, end, ref_seq, is_rna=False, gff_genes=None, bed_features=None,
                  highlight_sites=None, highlight_samples=None):
    """Common rendering logic for DNA/RNA/highlight modes.

    PNG and PDF are now exported from the SVG renderer rather than using a
    separate hand-written PNG renderer. This keeps rendering details consistent.
    """
    svg_content = _render_svg_content(
        tracks,
        args,
        chrom,
        start,
        end,
        ref_seq,
        is_rna=is_rna,
        gff_genes=gff_genes,
        bed_features=bed_features,
        highlight_sites=highlight_sites,
        highlight_samples=highlight_samples,
    )
    _write_svg_based_output(svg_content, args.out)

def add_highlight_args(parser):
    """Add arguments specific to the highlight subcommand.

    The highlight command is a focused variant of `dna`/`rna` that makes
    --highlight-vcf mandatory and sets defaults optimised for SNP-linkage
    analysis (detail=high, show-axis on, hap-sort on, hap-filter on).
    All other common arguments are still available.
    """
    # --- Required: VCF with target sites -------------------------------------
    parser.add_argument(
        "--highlight-vcf", required=True,
        help="VCF file (.vcf or .vcf.gz) with the target SNP / small-variant "
             "sites to highlight.  Reads are coloured at these positions by "
             "the base actually observed, and clustered by their inferred "
             "haplotype signature across the displayed sites.",
    )
    parser.add_argument(
        "--highlight-samples", nargs="+",
        help="Subset of sample names from the VCF header to show in the "
             "highlight track (default: all samples).",
    )

    # --- Sequencing mode (DNA vs RNA) ----------------------------------------
    parser.add_argument(
        "--mode", choices=["dna", "rna"], default="dna",
        help="Sequencing data type: 'dna' (default) or 'rna' (handles spliced "
             "alignments and shows arc connections for introns/split reads).",
    )

    # --- Hap-sort / hap-filter overrides (defined here, not in add_common_args) -
    parser.add_argument(
        "--no-hap-sort", action="store_true",
        help="Do NOT sort/cluster reads by their observed haplotype signature. "
             "Keeps the default coordinate-based stacking.",
    )
    parser.add_argument(
        "--no-hap-filter", action="store_true",
        help="Do NOT hide reads that don't overlap any VCF site. By default "
             "such reads are hidden to reduce visual clutter.",
    )

    # --- Common arguments (shared with dna/rna, skips the highlight-specific ones) -
    # Mark this parser so add_common_args skips --highlight-vcf etc.
    parser._is_highlight_parser = True
    add_common_args(parser)

    # --- Highlight-friendly defaults -----------------------------------------
    parser.set_defaults(
        detail="high",
        show_axis=True,
        no_hap_sort=False,
        no_hap_filter=False,
    )


def main():
    p = argparse.ArgumentParser(description="Bamsnap-LRS: Long-read sequencing data visualization tool")
    sub = p.add_subparsers(dest="cmd")

    # dna command - generate genomic region snapshot for DNA data
    dna_parser = sub.add_parser("dna", help="Generate read pileup snapshot for DNA sequencing data")
    add_common_args(dna_parser)

    # rna command - generate genomic region snapshot for RNA data
    rna_parser = sub.add_parser("rna", help="Generate read pileup snapshot for RNA sequencing data (handles spliced alignments)")
    add_common_args(rna_parser)

    # highlight command - VCF-site-focused SNP linkage visualisation
    highlight_parser = sub.add_parser(
        "highlight",
        help="Highlight VCF target sites on reads for SNP-linkage analysis. "
             "Reads are coloured only at VCF positions and clustered by "
             "observed haplotype; non-VCF mismatches are dimmed to gray.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "highlight mode\n"
            "==============\n"
            "Visualise how individual reads span a set of target SNP / small-\n"
            "variant sites from a VCF file.  Compared to the plain dna/rna\n"
            "commands, highlight mode:\n"
            "  • Requires --highlight-vcf (the VCF is the primary input).\n"
            "  • Adds a 'VCF Highlights' track showing per-sample REF/ALT\n"
            "    alleles as colour strips (IGV-style).\n"
            "  • On each read, only mismatches at VCF positions are coloured\n"
            "    by base; all other mismatches are dimmed to pale gray so\n"
            "    SNP-linkage patterns stand out.\n"
            "  • Reads that do not overlap any VCF site are hidden by default\n"
            "    (use --no-hap-filter to keep them).\n"
            "  • Reads are sorted by their observed haplotype signature\n"
            "    (use --no-hap-sort to disable).\n"
            "  • Defaults to detail=high and --show-axis.\n"
            "  • Supports both DNA (default) and RNA (--mode rna) data.\n"
        ),
    )
    add_highlight_args(highlight_parser)

    args = p.parse_args()

    if args.cmd is None:
        p.print_help()
        return

    # 'highlight' is a specialised wrapper around the dna/rna pipeline.
    # Translate it into the appropriate cmd so all downstream logic is reused.
    if args.cmd == "highlight":
        args.cmd = getattr(args, "mode", "dna")
        # --highlight-vcf is already set as args.highlight_vcf by argparse.
        # --no-hap-sort / --no-hap-filter are already present from add_highlight_args.

    # Determine if batch mode or single region mode
    is_batch_mode = args.regions is not None

    if is_batch_mode:
        # Batch mode: process multiple regions from file
        if not args.out_prefix:
            print("Error: --out-prefix is required when using --regions (batch mode)", file=sys.stderr)
            sys.exit(1)

        # Parse output prefix to extract directory and prefix
        if '/' in args.out_prefix or '\\' in args.out_prefix:
            out_dir = os.path.dirname(args.out_prefix)
            out_prefix = os.path.basename(args.out_prefix)
        else:
            out_dir = '.'
            out_prefix = args.out_prefix

        # Ensure output directory exists
        os.makedirs(out_dir, exist_ok=True)

        # Parse regions file
        from .regions import parse_regions_file
        try:
            # If --padding is omitted, pass None so regions.py can infer per-record padding.
            regions = parse_regions_file(args.regions, padding=getattr(args, 'padding', None))
        except Exception as e:
            print(f"Error: Failed to parse regions file '{args.regions}': {e}", file=sys.stderr)
            sys.exit(1)

        if not regions:
            print("Warning: No regions found in file", file=sys.stderr)
            return

        print(f"Processing {len(regions)} regions...", file=sys.stderr)

        # Determine output format from prefix (or default to svg)
        output_ext = 'svg'  # Default
        if args.out_prefix.lower().endswith('.png'):
            output_ext = 'png'
        elif args.out_prefix.lower().endswith('.pdf'):
            output_ext = 'pdf'
        elif args.out_prefix.lower().endswith('.svg'):
            output_ext = 'svg'

        # Remove extension from prefix if present
        prefix_base = out_prefix
        for ext in ['.png', '.pdf', '.svg']:
            if prefix_base.lower().endswith(ext):
                prefix_base = prefix_base[:-len(ext)]
                break

        # Process each region
        success_count = 0
        failed_regions = []

        for idx, (chrom, start, end, target_start, target_end) in enumerate(regions, 1):
            print(f"[{idx}/{len(regions)}] Processing {chrom}:{start}-{end}...", file=sys.stderr)

            # Generate output filename: outprefix_chr_start_end.ext
            output_filename = f"{prefix_base}_{chrom}_{start}_{end}.{output_ext}"
            output_path = os.path.join(out_dir, output_filename)

            # Create a modified args object for this region
            region_args = argparse.Namespace(**vars(args))
            region_args.out = output_path
            region_args.pos = f"{chrom}:{start}-{end}"
            region_args.focus_region = (target_start, target_end)

            # Process this region (catch errors to continue processing)
            try:
                process_single_region(region_args, chrom, start, end)
                success_count += 1
                print(f"  ✓ Saved to {output_path}", file=sys.stderr)
            except (SystemExit, ValueError) as e:
                # If process_single_region exits due to no reads, continue
                error_msg = str(e) if e else "No reads found"
                failed_regions.append((chrom, start, end, error_msg))
                print(f"  ✗ Skipped: {error_msg}", file=sys.stderr)
            except Exception as e:
                failed_regions.append((chrom, start, end, str(e)))
                print(f"  ✗ Error: {e}", file=sys.stderr)

        print(f"\nCompleted processing {len(regions)} regions.", file=sys.stderr)
        print(f"  Successfully processed: {success_count}", file=sys.stderr)
        if failed_regions:
            print(f"  Failed/Skipped: {len(failed_regions)}", file=sys.stderr)
            for chrom, start, end, reason in failed_regions:
                print(f"    - {chrom}:{start}-{end}: {reason}", file=sys.stderr)
        return

    # Single region mode
    if not args.out:
        print("Error: --out is required when using --pos (single region mode)", file=sys.stderr)
        sys.exit(1)

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

    process_single_region(args, chrom, start, end)


def process_single_region(args, chrom, start, end):
    """Process a single genomic region"""

    # Fetch annotation data if provided
    gff_genes = None
    bed_features = None
    if args.gff:
        from .gff import parse_gff
        gff_genes = parse_gff(args.gff, chrom, start, end)
    elif args.bed:
        from .bed import parse_bed
        bed_features = parse_bed(args.bed, chrom, start, end)

    # Parse highlight VCF if provided. We do this before fetching reads so we
    # can fail fast on malformed VCFs without wasting BAM I/O.
    highlight_sites = None
    highlight_samples = None
    if getattr(args, "highlight_vcf", None):
        from .highlight import parse_highlight_vcf
        try:
            highlight_sites, highlight_samples = parse_highlight_vcf(
                args.highlight_vcf, chrom, start, end
            )
        except Exception as e:
            print(f"Error: Failed to parse highlight VCF '{args.highlight_vcf}': {e}", file=sys.stderr)
            sys.exit(1)
        # Optional sample subset: keep order as user specified
        if getattr(args, "highlight_samples", None):
            requested = args.highlight_samples
            missing = [s for s in requested if s not in highlight_samples]
            if missing:
                print(f"Warning: highlight samples not found in VCF: {missing}", file=sys.stderr)
            highlight_samples = [s for s in requested if s in highlight_samples]
            # Filter the per-site sample_calls dicts down to the requested set
            for site in highlight_sites:
                site.sample_calls = {k: v for k, v in site.sample_calls.items() if k in highlight_samples}
        if not highlight_sites:
            print(f"Note: no highlight VCF sites fall inside {chrom}:{start}-{end}", file=sys.stderr)

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

        # Check if BAM files exist
        missing_bams = []
        for bam_path in args.bam:
            if not os.path.exists(bam_path):
                missing_bams.append(bam_path)

        if missing_bams:
            print("Error: BAM file(s) not found:", file=sys.stderr)
            for bam_path in missing_bams:
                print(f"  - {bam_path}", file=sys.stderr)
            sys.exit(1)

        empty_bams = []
        for i, bam_path in enumerate(args.bam):
            try:
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
            except Exception as e:
                print(f"Error: Failed to read BAM file '{bam_path}': {e}", file=sys.stderr)
                sys.exit(1)

            if not reads:
                empty_bams.append((bam_path, titles[i]))
            else:
                tracks.append({"reads": reads, "title": titles[i]})

        # Check if all BAM/CRAM files have no reads
        if not tracks:
            # In batch mode, raise exception to be caught by caller; in single mode, exit with error
            if hasattr(args, 'regions') and args.regions:
                # Batch mode: raise exception to be caught by caller
                raise ValueError(f"No reads found in region {chrom}:{start}-{end}")
            else:
                # Single mode: print error and exit
                print("Error: No reads found in the specified region.", file=sys.stderr)
                print(f"\nRegion: {chrom}:{start}-{end}", file=sys.stderr)
                print(f"\nBAM/CRAM file(s) checked:", file=sys.stderr)
                for bam_path, title in empty_bams:
                    print(f"  - {bam_path} ({title})", file=sys.stderr)
                print(f"\nPossible reasons:", file=sys.stderr)
                print(f"  1. The region has no aligned reads", file=sys.stderr)
                print(f"  2. MAPQ filter too strict (current: --mapq {args.mapq})", file=sys.stderr)
                if not args.show_supp:
                    print(f"  3. Supplementary alignments are hidden (use --show-supp to include)", file=sys.stderr)
                if not args.show_secondary:
                    print(f"  4. Secondary alignments are hidden (use --show-secondary to include)", file=sys.stderr)
                print(f"  5. Chromosome name mismatch (check if '{chrom}' exists in BAM/CRAM file)", file=sys.stderr)
                sys.exit(1)

        # Warn if some BAM/CRAM files have no reads
        if empty_bams:
            print("Warning: Some BAM/CRAM files have no reads in the specified region:", file=sys.stderr)
            for bam_path, title in empty_bams:
                print(f"  - {bam_path} ({title})", file=sys.stderr)
            print("", file=sys.stderr)

        ref_seq = None
        if args.fa:
            ref_seq = get_ref_subseq(args.fa, chrom, start, end)
        render_output(tracks, args, chrom, start, end, ref_seq, is_rna=False, gff_genes=gff_genes, bed_features=bed_features,
                      highlight_sites=highlight_sites, highlight_samples=highlight_samples)

    elif args.cmd == "rna":
        from .reader import fetch_rna_reads
        from .ref import get_ref_subseq

        # Check if BAM/CRAM files exist
        missing_bams = []
        for bam_path in args.bam:
            if not os.path.exists(bam_path):
                missing_bams.append(bam_path)

        if missing_bams:
            print("Error: BAM/CRAM file(s) not found:", file=sys.stderr)
            for bam_path in missing_bams:
                print(f"  - {bam_path}", file=sys.stderr)
            sys.exit(1)

        # Check if CRAM files have reference genome
        cram_files = [bam_path for bam_path in args.bam if bam_path.lower().endswith('.cram')]
        if cram_files and not args.fa:
            print("Warning: CRAM file(s) detected but no reference genome provided (--fa).", file=sys.stderr)
            print("  pysam will attempt to find the reference from:", file=sys.stderr)
            print("  1. CRAM file header (UR tag)", file=sys.stderr)
            print("  2. REF_PATH or REF_CACHE environment variables", file=sys.stderr)
            print("  If this fails, please provide --fa with the reference FASTA file.", file=sys.stderr)
            print("", file=sys.stderr)

        empty_bams = []
        for i, bam_path in enumerate(args.bam):
            try:
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
            except Exception as e:
                file_type = "CRAM" if bam_path.lower().endswith('.cram') else "BAM"
                print(f"Error: Failed to read {file_type} file '{bam_path}': {e}", file=sys.stderr)
                if bam_path.lower().endswith('.cram') and not args.fa:
                    print("  Hint: CRAM files require a reference genome. Try providing --fa <reference.fasta>", file=sys.stderr)
                sys.exit(1)

            if not reads:
                empty_bams.append((bam_path, titles[i]))
            else:
                tracks.append({"reads": reads, "title": titles[i]})

        # Check if all BAM/CRAM files have no reads
        if not tracks:
            print("Error: No reads found in the specified region.", file=sys.stderr)
            print(f"\nRegion: {chrom}:{start}-{end}", file=sys.stderr)
            print(f"\nBAM/CRAM file(s) checked:", file=sys.stderr)
            for bam_path, title in empty_bams:
                print(f"  - {bam_path} ({title})", file=sys.stderr)
            print(f"\nPossible reasons:", file=sys.stderr)
            print(f"  1. The region has no aligned reads", file=sys.stderr)
            print(f"  2. MAPQ filter too strict (current: --mapq {args.mapq})", file=sys.stderr)
            if not args.show_supp:
                print(f"  3. Supplementary alignments are hidden (use --show-supp to include)", file=sys.stderr)
            if not args.show_secondary:
                print(f"  4. Secondary alignments are hidden (use --show-secondary to include)", file=sys.stderr)
            print(f"  5. Chromosome name mismatch (check if '{chrom}' exists in BAM/CRAM file)", file=sys.stderr)
            print(f"  6. For RNA mode, ensure spliced alignments are properly marked", file=sys.stderr)
            sys.exit(1)

        # Warn if some BAM/CRAM files have no reads
        if empty_bams:
            print("Warning: Some BAM/CRAM files have no reads in the specified region:", file=sys.stderr)
            for bam_path, title in empty_bams:
                print(f"  - {bam_path} ({title})", file=sys.stderr)
            print("", file=sys.stderr)

        ref_seq = None
        if args.fa:
            ref_seq = get_ref_subseq(args.fa, chrom, start, end)
        render_output(tracks, args, chrom, start, end, ref_seq, is_rna=True, gff_genes=gff_genes, bed_features=bed_features,
                      highlight_sites=highlight_sites, highlight_samples=highlight_samples)


if __name__ == "__main__":
    main()
