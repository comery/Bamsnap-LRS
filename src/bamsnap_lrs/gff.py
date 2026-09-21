import os
import gzip
from dataclasses import dataclass, field
from functools import lru_cache
from typing import List, Dict


@dataclass
class Exon:
    start: int
    end: int


@dataclass
class Transcript:
    id: str
    name: str
    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    exons: List[Exon] = field(default_factory=list)
    cds: List[Exon] = field(default_factory=list)


@dataclass
class Gene:
    id: str
    name: str
    chrom: str
    start: int
    end: int
    strand: str
    exons: List[Exon] = field(default_factory=list)
    cds: List[Exon] = field(default_factory=list)
    transcripts: List[Transcript] = field(default_factory=list)


def _open_annotation_file(path: str):
    """Open plain-text or gzip-compressed GFF/GTF annotation files."""
    if path.lower().endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_attributes(attr_str: str) -> Dict[str, str]:
    """Parse both GFF3 (key=value) and GTF (key \"value\") attributes."""
    attrs: Dict[str, str] = {}
    for part in attr_str.split(';'):
        part = part.strip()
        if not part:
            continue

        if '=' in part:
            key, val = part.split('=', 1)
        else:
            fields = part.split(' ', 1)
            if len(fields) != 2:
                continue
            key, val = fields

        attrs[key.strip()] = val.strip().strip('"')
    return attrs


def _get_or_create_gene(
    genes: Dict[str, Gene],
    gene_id: str,
    gene_name: str,
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> Gene:
    gene = genes.get(gene_id)
    if gene is None:
        gene = Gene(
            id=gene_id,
            name=gene_name or gene_id,
            chrom=chrom,
            start=start,
            end=end,
            strand=strand,
        )
        genes[gene_id] = gene
    else:
        gene.start = min(gene.start, start)
        gene.end = max(gene.end, end)
        if gene_name and (not gene.name or gene.name == gene.id):
            gene.name = gene_name
        if gene.strand not in ('+', '-') and strand in ('+', '-'):
            gene.strand = strand
    return gene


def _get_or_create_transcript(
    transcripts: Dict[str, Transcript],
    genes: Dict[str, Gene],
    transcript_id: str,
    transcript_name: str,
    gene_id: str,
    gene_name: str,
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> Transcript:
    tx = transcripts.get(transcript_id)
    if tx is None:
        tx = Transcript(
            id=transcript_id,
            name=transcript_name or transcript_id,
            gene_id=gene_id,
            gene_name=gene_name or gene_id,
            chrom=chrom,
            start=start,
            end=end,
            strand=strand,
        )
        transcripts[transcript_id] = tx
        gene = _get_or_create_gene(
            genes, gene_id, gene_name, chrom, start, end, strand
        )
        gene.transcripts.append(tx)
    else:
        tx.start = min(tx.start, start)
        tx.end = max(tx.end, end)
        if transcript_name and (not tx.name or tx.name == tx.id):
            tx.name = transcript_name
        if gene_name and (not tx.gene_name or tx.gene_name == tx.gene_id):
            tx.gene_name = gene_name
    return tx


@lru_cache(maxsize=2)
def load_gff(gff_path: str) -> Dict[str, List[Gene]]:
    """Parse a GFF/GTF file once and cache all gene/transcript annotations.

    The cached structure is grouped by chromosome. GFF/GTF coordinates are
    converted from 1-based inclusive to 0-based half-open coordinates.
    """
    if not os.path.exists(gff_path):
        return {}

    genes: Dict[str, Gene] = {}
    transcripts: Dict[str, Transcript] = {}
    transcript_to_gene: Dict[str, str] = {}

    with _open_annotation_file(gff_path) as f:
        for line in f:
            if not line or line.startswith('#'):
                continue

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue

            r_chrom = parts[0]
            r_type = parts[2].lower()
            r_strand = parts[6]

            try:
                r_start = int(parts[3]) - 1
                r_end = int(parts[4])
            except ValueError:
                continue

            attrs = parse_attributes(parts[8])

            # Common GTF fields.
            gtf_gene_id = attrs.get('gene_id')
            gtf_gene_name = attrs.get('gene_name')
            gtf_tx_id = attrs.get('transcript_id')
            gtf_tx_name = attrs.get('transcript_name')

            # GFF3 fields.
            gff_id = attrs.get('ID')
            parent_raw = attrs.get('Parent')
            parent_ids = [p.strip() for p in parent_raw.split(',')] if parent_raw else []

            if r_type == 'gene':
                gene_id = gtf_gene_id or gff_id
                if not gene_id:
                    continue
                gene_name = gtf_gene_name or attrs.get('Name') or attrs.get('gene') or gene_id
                _get_or_create_gene(
                    genes, gene_id, gene_name, r_chrom, r_start, r_end, r_strand
                )
                continue

            if r_type in ('transcript', 'mrna'):
                transcript_id = gtf_tx_id or gff_id
                if not transcript_id:
                    continue

                # GTF provides gene_id directly. In GFF3 the mRNA/transcript
                # Parent points to the gene.
                gene_id = gtf_gene_id or (parent_ids[0] if parent_ids else None)
                if not gene_id:
                    continue

                gene_name = (
                    gtf_gene_name
                    or attrs.get('gene')
                    or attrs.get('gene_name')
                    or (genes[gene_id].name if gene_id in genes else gene_id)
                )
                transcript_name = gtf_tx_name or attrs.get('Name') or transcript_id

                _get_or_create_transcript(
                    transcripts,
                    genes,
                    transcript_id,
                    transcript_name,
                    gene_id,
                    gene_name,
                    r_chrom,
                    r_start,
                    r_end,
                    r_strand,
                )
                transcript_to_gene[transcript_id] = gene_id
                continue

            if r_type not in ('exon', 'cds'):
                continue

            # A GTF feature normally gives both gene_id and transcript_id.
            # A GFF3 exon/CDS normally uses Parent=<transcript ID>. Parent can
            # contain multiple transcript IDs, so attach the feature to each.
            feature_tx_ids: List[str] = []
            if gtf_tx_id:
                feature_tx_ids = [gtf_tx_id]
            elif parent_ids:
                feature_tx_ids = parent_ids

            gene_id_from_record = gtf_gene_id
            gene_name_from_record = gtf_gene_name or attrs.get('gene')

            attached_gene_ids = set()

            for transcript_id in feature_tx_ids:
                gene_id = gene_id_from_record or transcript_to_gene.get(transcript_id)
                if not gene_id:
                    continue

                gene_name = (
                    gene_name_from_record
                    or (genes[gene_id].name if gene_id in genes else gene_id)
                )

                tx = _get_or_create_transcript(
                    transcripts,
                    genes,
                    transcript_id,
                    transcript_id,
                    gene_id,
                    gene_name,
                    r_chrom,
                    r_start,
                    r_end,
                    r_strand,
                )
                transcript_to_gene[transcript_id] = gene_id

                feature = Exon(r_start, r_end)
                if r_type == 'exon':
                    tx.exons.append(feature)
                else:
                    tx.cds.append(feature)

                attached_gene_ids.add(gene_id)

            # Keep gene-level feature lists for the existing collapsed gene
            # display. Each source feature line is added once per gene even if
            # a GFF3 Parent field lists multiple transcripts from that gene.
            if gene_id_from_record:
                attached_gene_ids.add(gene_id_from_record)
            elif not attached_gene_ids:
                for parent_id in parent_ids:
                    if parent_id in genes:
                        attached_gene_ids.add(parent_id)

            for gene_id in attached_gene_ids:
                gene_name = (
                    gene_name_from_record
                    or (genes[gene_id].name if gene_id in genes else gene_id)
                )
                gene = _get_or_create_gene(
                    genes, gene_id, gene_name, r_chrom, r_start, r_end, r_strand
                )
                feature = Exon(r_start, r_end)
                if r_type == 'exon':
                    gene.exons.append(feature)
                else:
                    gene.cds.append(feature)

    # Sort nested features and transcripts once, before caching.
    for tx in transcripts.values():
        tx.exons.sort(key=lambda x: (x.start, x.end))
        tx.cds.sort(key=lambda x: (x.start, x.end))

    genes_by_chrom: Dict[str, List[Gene]] = {}
    for gene in genes.values():
        gene.exons.sort(key=lambda x: (x.start, x.end))
        gene.cds.sort(key=lambda x: (x.start, x.end))
        gene.transcripts.sort(key=lambda t: (t.start, t.end, t.id))
        genes_by_chrom.setdefault(gene.chrom, []).append(gene)

    for chrom_genes in genes_by_chrom.values():
        chrom_genes.sort(key=lambda g: (g.start, g.end, g.id))

    return genes_by_chrom


def query_gff(
    annotation: Dict[str, List[Gene]],
    chrom: str,
    start: int,
    end: int,
) -> List[Gene]:
    """Return genes overlapping [start, end)."""
    return [
        gene
        for gene in annotation.get(chrom, [])
        if gene.end >= start and gene.start <= end
    ]


def parse_gff(gff_path: str, chrom: str, start: int, end: int) -> List[Gene]:
    """Backward-compatible GFF/GTF region query.

    The annotation file is parsed only on the first call for a given path;
    subsequent regions reuse the in-memory cache.
    """
    return query_gff(load_gff(gff_path), chrom, start, end)


def clear_gff_cache() -> None:
    """Clear the in-memory annotation cache."""
    load_gff.cache_clear()


def gff_cache_info():
    """Return cache statistics for testing/debugging."""
    return load_gff.cache_info()
