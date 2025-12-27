"""Parse BED and VCF files to extract genomic regions"""
import os
import gzip
from typing import List, Tuple, Optional


def _parse_vcf_info(info_str: str) -> dict:
    """Parse VCF INFO field into a dictionary"""
    info_dict = {}
    if not info_str or info_str == '.':
        return info_dict
    
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    
    return info_dict


def _calculate_vcf_region(chrom: str, pos: int, ref: str, alt: str, info: Optional[dict] = None, 
                         default_padding: int = 250) -> Tuple[str, int, int]:
    """
    Calculate genomic region for a VCF variant.
    
    Args:
        chrom: Chromosome name
        pos: 1-based position
        ref: Reference allele
        alt: Alternate allele(s), comma-separated
        info: INFO field dictionary
        default_padding: Default padding around variant (bp)
    
    Returns:
        (chrom, start, end) tuple
    """
    # Check for END tag in INFO (for structural variants)
    end_pos = None
    if info:
        if 'END' in info:
            try:
                end_pos = int(info['END'])
            except (ValueError, TypeError):
                pass
        
        # For structural variants, use SVLEN if available
        if 'SVTYPE' in info and 'SVLEN' in info:
            try:
                svlen = abs(int(info['SVLEN']))
                end_pos = pos + svlen
            except (ValueError, TypeError):
                pass
    
    # Calculate region boundaries
    if end_pos and end_pos > pos:
        # Structural variant with explicit END
        start = max(0, pos - default_padding)
        end = end_pos + default_padding
    else:
        # SNV or small indel
        # Calculate variant length
        ref_len = len(ref)
        
        # For ALT, take the longest allele if multiple
        alt_alleles = alt.split(',') if alt else ['']
        max_alt_len = max(len(a) for a in alt_alleles) if alt_alleles else 0
        
        # Variant spans from pos to pos + max(ref_len, alt_len) - 1
        variant_end = pos + max(ref_len, max_alt_len) - 1
        
        # Add padding
        start = max(0, pos - default_padding)
        end = variant_end + default_padding
    
    return (chrom, start, end)


def parse_regions_file(regions_path: str, vcf_padding: int = 250) -> List[Tuple[str, int, int]]:
    """
    Parse BED or VCF file to extract genomic regions.
    
    Args:
        regions_path: Path to BED or VCF file
        vcf_padding: Padding around VCF variants (bp), default 250
        
    Returns:
        List of (chrom, start, end) tuples
        
    BED format:
        - Columns: chrom, start, end, ...
        - Start is 0-based, end is exclusive
        - Returns: (chrom, start, end)
        
    VCF format:
        - Columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, ...
        - POS is 1-based
        - Supports standard VCF and VCFv4.3+ formats
        - Handles structural variants with END tag in INFO field
        - Handles SVTYPE and SVLEN for structural variants
        - Calculates region based on REF/ALT lengths
        - Returns: (chrom, start, end)
    """
    if not os.path.exists(regions_path):
        raise FileNotFoundError(f"Regions file not found: {regions_path}")
    
    regions: List[Tuple[str, int, int]] = []
    
    # Detect file format by extension
    is_vcf = regions_path.lower().endswith('.vcf') or regions_path.lower().endswith('.vcf.gz')
    is_gzipped = regions_path.lower().endswith('.gz')
    
    # Open file (handle gzip if needed)
    if is_gzipped:
        open_func = gzip.open
        mode = 'rt'  # text mode for gzip
    else:
        open_func = open
        mode = 'r'
    
    with open_func(regions_path, mode) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            # Skip VCF header lines (all lines starting with ## or #)
            if line.startswith('#'):
                # Check for VCF header line (starts with #CHROM)
                if is_vcf and line.startswith('#CHROM'):
                    # This is the column header line, continue to data
                    continue
                continue
            
            parts = line.split('\t')
            
            if is_vcf:
                # VCF format: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, ...
                if len(parts) < 5:  # At least CHROM, POS, ID, REF, ALT
                    continue
                
                chrom = parts[0]
                try:
                    pos = int(parts[1])  # 1-based position
                    ref = parts[3] if len(parts) > 3 else ''
                    alt = parts[4] if len(parts) > 4 else ''
                    
                    # Parse INFO field if present
                    info = None
                    if len(parts) > 7:
                        info = _parse_vcf_info(parts[7])
                    
                    # Calculate region
                    chrom, start, end = _calculate_vcf_region(chrom, pos, ref, alt, info, vcf_padding)
                    regions.append((chrom, start, end))
                except (ValueError, IndexError) as e:
                    # Skip malformed lines
                    continue
            else:
                # BED format: chrom, start, end, ...
                if len(parts) < 3:
                    continue
                
                chrom = parts[0]
                try:
                    start = int(parts[1])  # 0-based
                    end = int(parts[2])    # exclusive
                    if start < 0 or end <= start:
                        continue
                    regions.append((chrom, start, end))
                except ValueError:
                    continue
    
    return regions

