import os
from dataclasses import dataclass, field
from typing import List, Dict, Optional


@dataclass
class Exon:
    start: int
    end: int


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


def parse_gff(gff_path: str, chrom: str, start: int, end: int) -> List[Gene]:
    """Parse GFF/GTF file and extract genes within range"""
    if not os.path.exists(gff_path):
        return []

    genes: Dict[str, Gene] = {}
    transcripts: Dict[str, str] = {}  # transcript_id -> gene_id

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            r_chrom = parts[0]
            if r_chrom != chrom:
                continue
            
            r_start = int(parts[3]) - 1
            r_end = int(parts[4])
            
            # Check range overlap
            if r_end < start or r_start > end:
                continue
            
            r_type = parts[2].lower()
            r_strand = parts[6]
            attrs = parse_attributes(parts[8])
            
            gene_id = attrs.get('gene_id') or attrs.get('ID')
            gene_name = attrs.get('gene_name') or attrs.get('Name') or gene_id
            transcript_id = attrs.get('transcript_id') or attrs.get('ID')
            parent_id = attrs.get('Parent')

            if r_type in ['gene', 'transcript', 'mrna']:
                if gene_id and gene_id not in genes:
                    genes[gene_id] = Gene(
                        id=gene_id,
                        name=gene_name,
                        chrom=r_chrom,
                        start=r_start,
                        end=r_end,
                        strand=r_strand
                    )
                if transcript_id and gene_id:
                    transcripts[transcript_id] = gene_id
            
            elif r_type == 'exon':
                target_gene_id = None
                if parent_id in transcripts:
                    target_gene_id = transcripts[parent_id]
                elif parent_id in genes:
                    target_gene_id = parent_id
                elif gene_id in genes:
                    target_gene_id = gene_id
                
                if target_gene_id and target_gene_id in genes:
                    genes[target_gene_id].exons.append(Exon(r_start, r_end))
                    # Update gene boundaries if not set
                    genes[target_gene_id].start = min(genes[target_gene_id].start, r_start)
                    genes[target_gene_id].end = max(genes[target_gene_id].end, r_end)

            elif r_type == 'cds':
                target_gene_id = None
                if parent_id in transcripts:
                    target_gene_id = transcripts[parent_id]
                elif parent_id in genes:
                    target_gene_id = parent_id
                elif gene_id in genes:
                    target_gene_id = gene_id
                
                if target_gene_id and target_gene_id in genes:
                    genes[target_gene_id].cds.append(Exon(r_start, r_end))

    # Sort exons and CDS for each gene
    for gene in genes.values():
        gene.exons.sort(key=lambda x: x.start)
        gene.cds.sort(key=lambda x: x.start)

    # Return genes that overlap with the range
    result = [g for g in genes.values() if g.end >= start and g.start <= end]
    # Sort genes by start position
    result.sort(key=lambda x: x.start)
    return result


def parse_attributes(attr_str: str) -> Dict[str, str]:
    """Parse GFF/GTF attribute string"""
    attrs = {}
    # Handle both GFF3 (key=value) and GTF (key "value") formats
    for part in attr_str.split(';'):
        part = part.strip()
        if not part:
            continue
        
        if '=' in part:
            key, val = part.split('=', 1)
        else:
            # GTF format: key "value"
            parts = part.split(' ', 1)
            if len(parts) == 2:
                key, val = parts
            else:
                continue
        
        key = key.strip()
        val = val.strip().strip('"')
        attrs[key] = val
    return attrs

