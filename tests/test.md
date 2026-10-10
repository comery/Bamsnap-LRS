# BamSnap-LRS

**Static visualization of long-read sequencing alignments**

BamSnap-LRS is a command-line tool for generating genome-browser-style snapshots from long-read BAM/CRAM alignments. It supports visualization of DNA alignment views, structural variants, spliced transcript alignments, and SNV linkage, annotation tracks, and batch export to SVG, PDF, or PNG. Figures are generated without a browser or web server.

It is designed for static visualization and manual inspection, **not** for variant calling or computational haplotype phasing.

## Contents

- [Features](#features)
- [Installation](#installation)
- [Input files and indexes](#input-files-and-indexes)
- [Quick start](#quick-start)
- [How the visualization modes work](#how-the-visualization-modes-work)
- [Command-line options](#command-line-options)
- [Examples and gallery](#examples-and-gallery)
- [Limitations and troubleshooting](#limitations-and-troubleshooting)
- [Acknowledgments](#acknowledgments)

## Features

- **DNA and RNA read views:** aligned blocks, mismatches, insertions, deletions, clipping, strand indicators, coverage, and RNA splice junctions.
- **Long-read split alignments:** optional display of supplementary and secondary alignments; visible segments belonging to the same read can be distinguished in the standard read view.
- **SV snapshots:** one or many loci from BED/VCF, automatic or user-defined flanking windows, and shading of the target interval. Intended primarily for INS, DEL, DUP, and INV inspection.
- **VCF Highlight mode:** displays specified VCF sites and genotype/phase-set information, emphasizes alleles on the aligned reads, and groups reads by their observed phase-block signatures.
- **Annotations:** BED3–BED12 features or GFF/GTF annotations.
- **Multiple alignment tracks:** a separate track and coverage profile for each BAM/CRAM file.
- **Consistent exports:** SVG is the rendering source for SVG, PDF, and PNG output.

## Installation

### Requirements

- Python 3.8+
- Packages in [`requirements.txt`](requirements.txt): `pysam`, and `cairosvg`
- Optional: `pyfaidx` for efficient FASTA random access; without it, a plain-text FASTA fallback is used

```bash
git clone https://github.com/comery/Bamsnap-LRS.git
cd Bamsnap-LRS
python -m pip install -r requirements.txt
chmod +x bin/bamsnap-lrs

# Check commands and all available flags
python bin/bamsnap-lrs --help
python bin/bamsnap-lrs dna --help
python bin/bamsnap-lrs rna --help
python bin/bamsnap-lrs highlight --help
```
Run from the cloned repository directory using `python bin/bamsnap-lrs ...` or `./bin/bamsnap-lrs ...`. 

## Input files and indexes

| Input | Supported format | Index requirement | Notes |
|---|---|---|---|
| Alignments (`--bam`) | Coordinate-sorted BAM or CRAM | **Required:** BAM `.bai`/`.csi` or CRAM `.crai`| A reference FASTA is recommended for CRAM decoding. |
| Reference (`--fa`) | Uncompressed FASTA (`.fa`/`.fasta`) | **No pre-existing `.fai` required** | A reference FASTA is recommended for CRAM decoding. |
| Annotation (`--bed`) | Plain or gzip/BGZF-compressed BED | **No `.tbi`/`.csi` required** | BED3–BED12-style features. |
| Annotation (`--gff`) | Plain or gzip/BGZF-compressed GFF/GTF | **No `.tbi`/`.csi` required** | GFF/GTF is read sequentially; use `--gff-display gene|transcript`. |
| Batch loci (`--regions`) | BED, BED.gz, VCF, VCF.gz | **No `.tbi`/`.csi` required** | Compression is recognized by the `.gz` filename suffix; BGZF can be read as gzip. |
| Highlight sites (`--highlight-vcf`) | VCF, VCF.gz | **No `.tbi`/`.csi` required** | Supports multi-sample VCF and `GT`/`PS` for display. |

**Coordinate conventions:** BED uses 0-based half-open intervals; VCF `POS` is 1-based. The batch parser converts VCF records into internal plotting intervals. Chromosome/contig names must match across the alignment, reference, and annotation files (for example, `chr1` and `1` are not automatically interchangeable).

**Multiple BAM/CRAM files:** all tracks use a shared genomic coordinate system and should be aligned to the same reference assembly. Each input file forms its own read track with its own coverage. The default track label is the input filename without its extension; override with `--track-title`. Track titles are **not** automatically inferred from BAM read-group (`RG`) sample tags.

## Quick start

### 1. DNA snapshot

```bash
python bin/bamsnap-lrs dna \
    --bam example/sim.mapping.sort.bam \
    --pos chrM:1000-3000 \
    --fa example/chm13v2.chrM.fasta \
    --out dna.svg \
    --show-axis
```

### 2. RNA splicing and transcript annotations

```bash
python bin/bamsnap-lrs rna \
    --bam aligned_rna.bam \
    --pos chr1:1000000-1005000 \
    --gff annotation.gtf.gz \
    --gff-display transcript \
    --out rna.svg \
    --show-axis
```

Replace the example input paths with your own data. `--gff-display gene` (the default) displays collapsed gene models; `transcript` draws individual transcript models.

### 3. Batch SV snapshots from a VCF

```bash
python bin/bamsnap-lrs dna \
    --bam hifi.bam ont.bam \
    --regions variants.vcf.gz \
    --out-prefix results/sv.svg \
    --fa reference.fa \
    --show-supp \
    --show-axis
```

Batch mode requires `--regions` and `--out-prefix`. For each input record, the program writes a file named `<prefix>_<chrom>_<plot_start>_<plot_end>.svg` (or `.pdf`/`.png` if that extension was supplied). The output directory is created if needed. Failed or empty regions are reported rather than ending the whole DNA batch.

### 4. Highlight phased SNVs in a VCF

```bash
python bin/bamsnap-lrs highlight \
    --bam phased_reads.bam \
    --pos chr21:20000000-20002000 \
    --highlight-vcf phased_sites.vcf.gz \
    --fa reference.fa \
    --hap-layout packed \
    --out highlight.svg
```

Highlight mode defaults to a high-detail view, shows the coordinate axis, hides reads not intersecting highlighted sites, and clusters reads according to observed variant-site patterns. Use `--highlight-samples SAMPLE_A SAMPLE_B` to select samples from a multi-sample VCF.

### 5. Retain base-level differences in a wide window

```bash
python bin/bamsnap-lrs dna \
    --bam aligned_reads.bam \
    --pos chr5:705000-710000 \
    --fa reference.fa \
    --overview-detail show \
    --out overview.svg
```

By default, the renderer simplifies base-level differences when the scale is broader than one base per pixel. `--overview-detail show` requests that mismatches and small indel details remain visible at these scales. Dense details can overlap because multiple bases project onto a single pixel.

## How the visualization modes work

### DNA and RNA

The DNA view renders reads along the reference with optional mismatch, insertion, deletion, clipping, orientation, and supplementary-alignment information. The RNA view additionally shows splice-aware alignment features and junction arcs. Read inclusion is controlled by `--mapq`, `--show-supp`, and `--show-secondary`.

To identify mismatches within generic CIGAR `M` blocks, the program can compare aligned read bases against `--fa`, or use `--use-cs` and/or `--use-md`. If a reference sequence is available, reference-based decoding takes precedence; otherwise `cs` takes precedence over `MD` when both are requested. Without these information sources, CIGAR alone cannot resolve individual mismatches in `M` blocks.

### SV batch windows

BED records provide 0-based target intervals. For VCF records, the program uses available `SVTYPE`, `END`, and `SVLEN` information, as well as REF/ALT lengths, to infer plotting intervals. The target area is shaded within a larger window. `--padding N` overrides automatic padding on both sides. Insertions and small VCF variants normally receive compact flanking windows, whereas larger SVs receive larger windows.

**Scope:** SV figures show alignment evidence; they do not validate an SV call by themselves. No dedicated cross-chromosome breakend/translocation drawing is provided. Imprecise confidence intervals such as `CIPOS`/`CIEND` are not used to draw separate uncertainty bounds; complex or inconsistent SV fields should be reviewed before batch generation.

### VCF Highlight mode

The `highlight` command adds a VCF Highlights annotation track. At targeted SNV positions, observed read bases are emphasized while non-target differences are visually muted. A phased sample requires phased `GT` (with `|`) and a valid `PS` to define phase-set-specific haplotype rows. Unphased records use REF/ALT-style rows.

Within each phase block, a read is assigned to a haplotype only if its observed informative alleles uniquely support that haplotype. Conflicting, non-informative, or missing site evidence does not produce a forced assignment. Each phase block is assessed independently: a read spanning several blocks may have different haplotype states across segments. The layout groups reads by these block-specific signatures, rather than claiming a globally phased haplotype for the whole read.

Highlight-specific controls:

- `--no-hap-sort` retains coordinate-based stacking instead of variant-signature clustering.
- `--no-hap-filter` keeps reads that do not overlap target VCF sites.
- `--hap-layout packed` (default) reuses rows among nonoverlapping reads within a cluster; `--hap-layout row` keeps one alignment per row.
- `--highlight-samples` restricts the displayed VCF samples, not the BAM records by read group.
- `--mode rna` uses the RNA rendering pipeline for Highlight mode.

The command also reports Highlight summary counts (including reads covering valid target SNVs, haplotype-assigned reads, mixed-block patterns, and unassigned reads) to help interpret what is actually represented in the figure. These are read-level display/assignment summaries, not haplotype calls.

## Command-line options

Use `python bin/bamsnap-lrs <mode> --help` for the authoritative current set. All three modes share the common flags below; the Highlight mode adds its own flags.

### Common options (`dna`, `rna`, and `highlight`)

| Flag | Purpose | Default |
|---|---|---|
| `--bam FILE [FILE ...]` | One or more BAM/CRAM alignment tracks; can also be repeated | Required |
| `--pos CHR:START-END` | Render one region (or `CHR:POS`) | Mutually exclusive with `--regions` |
| `--regions FILE` | Batch regions in BED or VCF | Mutually exclusive with `--pos` |
| `--out FILE` | Single-region output `.svg`, `.pdf`, or `.png` | Required with `--pos` |
| `--out-prefix PREFIX` | Batch output prefix and extension | Required with `--regions` |
| `--padding N` | Fixed flanks for batch regions | Automatic per region |
| `--fa FILE` | Reference FASTA for mismatch comparison/CRAM | Optional for BAM; recommended for CRAM |
| `--mapq N` | Minimum mapping quality | `0` |
| `--show-supp` / `--show-secondary` | Include supplementary / secondary records | Off |
| `--width N` | Output width in pixels | `1200` |
| `--read-height N` | Read-body height | `6` |
| `--detail low|mid|high` | Base/indel detail level | `mid` (`high` in Highlight) |
| `--overview-detail hide|show` | Simplify or retain small features at overview scale | `hide` |
| `--show-axis` | Coordinate axis | Off (`on` in Highlight) |
| `--show-coverage` / `--no-coverage` | Per-track coverage display | On |
| `--coverage-height N` | Height of coverage track | `15` |
| `--coverage-max-depth N` | Set coverage plot depth maximum | Automatic |
| `--show-composition` | Add a base-composition track | Off |
| `--comp-height N` / `--comp-max-depth N` | Composition track sizing/depth | Automatic |
| `--track-title TITLE [TITLE ...]` | Override one or more track labels | Filename stem |
| `--gff FILE` / `--bed FILE` | Annotation input | Mutually exclusive |
| `--gff-display gene|transcript` | Gene or transcript-level models for GFF/GTF | `gene` |
| `--style default|jbrowse` | Visual style | `jbrowse` |
| `--color-by type|base|strand|mapq` | Read coloring | `type` |
| `--use-md` / `--use-cs` | Request MD/cs mismatch evidence | Off |
| `--hide-insertions` | De-emphasize insertions | Off |
| `--show-insertion-labels` / `--no-insertion-labels` | Insertion length labels | On |

### Highlight-specific options

| Flag | Purpose | Default |
|---|---|---|
| `--highlight-vcf FILE` | Target variant file | Required |
| `--highlight-samples SAMPLE [SAMPLE ...]` | VCF samples to display | All |
| `--mode dna|rna` | Alignment display mode | `dna` |
| `--no-hap-sort` | Disable observed-signature clustering | Sorting enabled |
| `--no-hap-filter` | Retain reads not intersecting target sites | Filtering enabled |
| `--hap-layout packed|row` | Packed or one-alignment-per-row layout | `packed` |

> **Removed from the documented CLI:** `--max-reads` and `--downsample-strategy` appear in older README revisions but are currently commented out in `src/bamsnap_lrs/cli.py`. They should not be used unless reintroduced to the command-line parser.

## Examples and gallery

The repository includes example commands in [`run_example.sh`](run_example.sh) and assets under [`example/`](example/).

**Visual comparison with JBrowse**

JBrowse browser view:

![JBrowse comparison screenshot](example/ScreenShot_jbrowse.png)

BamSnap-LRS static view:

![BamSnap-LRS DNA example](example/chrM_1000_3000.gff.bamsnapLRS.svg)

**RNA view**

![BamSnap-LRS RNA example](example/rna/test.rna.svg)

**Structural variants**

![Insertion example](example/SV_show/insertion.onebam_chrX_26021820_26021920.svg)

![Deletion example](example/SV_show/deletion.onebam_chr1_24314927_24320831.svg)

![Inversion example](example/SV_show/inversion_chrX_5770557_5773392.svg)

![Duplication example](example/SV_show/duplication_chr11_54871000_54876965.svg)

**VCF Highlight view**

![Highlight example](example/highlight/highlight_exam/highlight.one.bam.svg)

For a broader and less selectively illustrated set of loci, see the [BamSnap-LRS visualization gallery](https://wanyi0309.github.io/Bamsnap-LRS_visualization_gallery/).

## Limitations and troubleshooting

- **Missing reads or empty region:** verify contig names, reference assembly, interval, `--mapq`, whether supplementary/secondary alignments need to be included, and whether the BAM/CRAM index is available.
- **CRAM decoding:** provide the matching reference with `--fa` if decoding fails; pysam can sometimes locate references via CRAM header or reference environment settings.
- **No displayed mismatches:** supply `--fa` or request available `--use-cs`/`--use-md` tags. For broad windows, try `--overview-detail show`.
- **Annotation inputs:** `.gz` text files are supported for BED/GFF/GTF/VCF, but annotation and regions readers scan them sequentially; indexes are neither needed nor used for these files.
- **VCF semantics:** batch visualization is optimized for within-chromosome examples; it is not a full SV normalization engine or an inter-chromosomal BND viewer. Inspect records with multiple ALT alleles, missing or inconsistent `END`/`SVLEN`, or uncertain breakpoints.
- **Static output:** this is a figure-generation tool, not an interactive genome browser.

## License

The existing project README states that BamSnap-LRS is MIT-licensed. **Repository maintainers should add a `LICENSE` file** containing the applicable license text if it is not already present, so that the distribution terms are explicit.

## Acknowledgments

- Inspired by [JBrowse](https://jbrowse.org/)
- Inspired by the original [BamSnap](https://github.com/parklab/bamsnap)
- Uses [pysam](https://github.com/pysam-developers/pysam) for BAM/CRAM access
