#!/usr/bin/env python3
"""Convert GFF file to BED format"""

import sys
from collections import defaultdict

def parse_attributes(attr_str):
    """Parse GFF attribute string"""
    attrs = {}
    for part in attr_str.split(';'):
        part = part.strip()
        if not part:
            continue
        if '=' in part:
            key, val = part.split('=', 1)
            attrs[key.strip()] = val.strip()
    return attrs

def gff_to_bed(gff_path, bed_path):
    """Convert GFF file to BED format"""
    genes = {}  # gene_id -> gene info
    gene_exons = defaultdict(list)  # gene_id -> list of exons
    gene_cds = defaultdict(list)  # gene_id -> list of CDS
    transcript_to_gene = {}  # transcript_id -> gene_id
    
    # First pass: collect gene and transcript information
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom = parts[0]
            feature_type = parts[2].lower()
            start = int(parts[3])  # GFF is 1-based
            end = int(parts[4])
            strand = parts[6]
            attrs = parse_attributes(parts[8])
            
            # Extract gene_id (remove parentheses if present)
            gene_id = attrs.get('gene_id')
            if gene_id:
                gene_id = gene_id.split('(')[0]  # Remove parentheses part
            else:
                gene_id_from_id = attrs.get('ID', '')
                if gene_id_from_id.startswith('gene_'):
                    gene_id = gene_id_from_id[5:].split('(')[0]
                else:
                    gene_id = gene_id_from_id.split('(')[0]
            
            gene_name = attrs.get('Name') or gene_id or attrs.get('ID', '')
            
            if feature_type in ['gene', 'ncrna_gene']:
                if gene_id and gene_id not in genes:
                    genes[gene_id] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'name': gene_name
                    }
            elif feature_type in ['transcript', 'mrna', 'trna', 'rrna']:
                # Map transcript to gene
                transcript_id_full = attrs.get('ID', '')
                transcript_id = transcript_id_full.replace('transcript_', '').split('(')[0]
                parent = attrs.get('Parent', '')
                
                # Clean gene_id (remove parentheses if present)
                clean_gene_id = gene_id.split('(')[0] if gene_id else None
                
                if clean_gene_id:
                    transcript_to_gene[transcript_id] = clean_gene_id
                elif parent:
                    # Try to extract gene_id from parent
                    if parent.startswith('gene_'):
                        clean_parent = parent[5:].split('(')[0]
                        transcript_to_gene[transcript_id] = clean_parent
                    else:
                        clean_parent = parent.split('(')[0]
                        transcript_to_gene[transcript_id] = clean_parent
            elif feature_type == 'exon':
                parent = attrs.get('Parent', '')
                if parent:
                    # Extract transcript name from parent
                    transcript_name = parent.replace('transcript_', '').split('(')[0]
                    # Find gene_id for this transcript
                    target_gene_id = transcript_to_gene.get(transcript_name)
                    if not target_gene_id:
                        # Try direct match with gene names
                        for gid, gene_info in genes.items():
                            if gene_info['name'] == transcript_name or transcript_name.startswith(gene_info['name']):
                                target_gene_id = gid
                                break
                    if target_gene_id and target_gene_id in genes:
                        gene_exons[target_gene_id].append((start, end))
            elif feature_type == 'cds':
                parent = attrs.get('Parent', '')
                if parent:
                    transcript_name = parent.replace('transcript_', '').split('(')[0]
                    target_gene_id = transcript_to_gene.get(transcript_name)
                    if not target_gene_id:
                        for gid, gene_info in genes.items():
                            if gene_info['name'] == transcript_name or transcript_name.startswith(gene_info['name']):
                                target_gene_id = gid
                                break
                    if target_gene_id and target_gene_id in genes:
                        gene_cds[target_gene_id].append((start, end))
    
    # Second pass: write BED file
    with open(bed_path, 'w') as out:
        for gene_id, gene_info in sorted(genes.items(), key=lambda x: x[1]['start']):
            chrom = gene_info['chrom']
            gff_start = gene_info['start']  # GFF is 1-based
            gff_end = gene_info['end']
            strand = gene_info['strand'] if gene_info['strand'] in ['+', '-'] else '.'
            name = gene_info['name']
            
            # Convert to BED format (0-based, end exclusive)
            bed_start = gff_start - 1
            bed_end = gff_end
            
            # Get exons for this gene
            exons = gene_exons.get(gene_id, [])
            if not exons:
                # If no exons, use gene boundaries as single exon
                exons = [(gff_start, gff_end)]
            
            # Sort exons by start position and remove duplicates
            exons = sorted(set(exons), key=lambda x: x[0])
            
            # Get CDS for thickStart/thickEnd
            cds_list = gene_cds.get(gene_id, [])
            if cds_list:
                cds_list = sorted(set(cds_list), key=lambda x: x[0])
                thick_start = cds_list[0][0] - 1  # Convert to 0-based
                thick_end = cds_list[-1][1]
            else:
                # If no CDS, use gene boundaries
                thick_start = bed_start
                thick_end = bed_end
            
            # Calculate blocks (exons relative to gene start)
            block_count = len(exons)
            block_sizes = []
            block_starts = []
            
            for exon_start, exon_end in exons:
                block_start = exon_start - gff_start  # Relative to gene start
                block_size = exon_end - exon_start
                block_starts.append(str(block_start))
                block_sizes.append(str(block_size))
            
            # Default color: blue for genes, green for ncRNA
            if 'ncRNA' in gene_id or 'trn' in name.lower() or 'rrn' in name.lower():
                rgb = "0,128,0"  # Green for ncRNA
            else:
                rgb = "0,0,255"  # Blue for protein-coding
            
            # Score: use 0 as default
            score = "0"
            
            # Write BED line
            bed_line = [
                chrom,
                str(bed_start),
                str(bed_end),
                name,
                score,
                strand,
                str(thick_start),
                str(thick_end),
                rgb,
                str(block_count),
                ",".join(block_sizes) + ",",
                ",".join(block_starts) + ","
            ]
            
            out.write("\t".join(bed_line) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_gff_to_bed.py <input.gff> <output.bed>")
        sys.exit(1)
    
    gff_to_bed(sys.argv[1], sys.argv[2])
    print(f"Converted {sys.argv[1]} to {sys.argv[2]}")

