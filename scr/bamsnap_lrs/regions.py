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


def _first_int(value) -> Optional[int]:
    """Return the first integer from a VCF INFO value such as '123' or '-123,456'."""
    if value is None:
        return None
    try:
        return int(str(value).split(',')[0])
    except (ValueError, TypeError):
        return None


def _round_half_length_to_250(length: int, minimum: int = 250) -> int:
    """Use about half of the region/SV length, rounded to a nearby 500-bp step."""
    if length <= 0:
        return minimum
    half_len = length / 2.0
    rounded = int((half_len + 250) // 250) * 250
    return max(minimum, rounded)


def _infer_bed_padding(target_start: int, target_end: int) -> int:
    """Auto padding for BED regions when --padding is not provided."""
    region_len = max(1, target_end - target_start)
    return _round_half_length_to_250(region_len)


def _infer_vcf_padding(pos: int, ref: str, alt: str, info: Optional[dict]) -> int:
    """Infer padding for one VCF record when --padding is not provided.

    Rules:
      1. Insertions and small SNV/indel records use 150 bp.
      2. Other SVs use about half of SV length, rounded to a nearby 250-bp step.
    """
    info = info or {}
    svtype = str(info.get('SVTYPE', '')).upper()
    end_pos = _first_int(info.get('END'))
    svlen = _first_int(info.get('SVLEN'))
    abs_svlen = abs(svlen) if svlen is not None else 0

    ref_len = len(ref) if ref else 1
    alt_alleles = alt.split(',') if alt else ['']
    max_alt_len = max((len(a) for a in alt_alleles), default=0)
    is_symbolic_ins = any(a.upper() == '<INS>' for a in alt_alleles)

    # Long-read SV VCFs often encode insertions as SVTYPE=INS and END=POS.
    if svtype == 'INS' or is_symbolic_ins or (end_pos is not None and end_pos == pos and abs_svlen > 0):
        return 150

    sv_len = 0
    if abs_svlen > 0:
        sv_len = abs_svlen
    elif end_pos is not None and end_pos > pos:
        sv_len = end_pos - pos

    is_sv = svtype in {'INV', 'DEL', 'DUP', 'CNV', 'BND'} or sv_len > 50
    if is_sv and sv_len > 0:
        return _round_half_length_to_250(sv_len)

    # Small SNV / small indel: use the same compact window as insertion records.
    if max(ref_len, max_alt_len) <= 50:
        return 150

    # Fallback for long non-symbolic alleles without SVTYPE.
    return _round_half_length_to_250(max(ref_len, max_alt_len))


def _calculate_vcf_region(chrom: str, pos: int, ref: str, alt: str, info: Optional[dict] = None,
                         padding: Optional[int] = None) -> Tuple[str, int, int]:
    """
    Calculate plotting region for a VCF variant.

    If padding is provided by the user, use it directly. If padding is None,
    infer padding for each VCF record:
      - insertion / POS == END INS: 150 bp
      - small SNV / small indel: 150 bp
      - other SVs: roughly half of SV length, rounded to a nearby 500-bp step
    """
    LARGE_VARIANT_THRESHOLD = 50
    plot_padding = padding if padding is not None else _infer_vcf_padding(pos, ref, alt, info)

    end_pos = None
    if info:
        svlen = _first_int(info.get('SVLEN'))
        if info.get('SVTYPE') and svlen is not None:
            end_pos = pos + abs(svlen)

        parsed_end = _first_int(info.get('END'))
        if parsed_end is not None:
            end_pos = parsed_end

    use_sv_region = False
    if end_pos is not None:
        variant_length = abs(_first_int(info.get('SVLEN')) or 0) if info else 0
        svtype = str(info.get('SVTYPE', '')).upper() if info else ''
        is_symbolic_alt = alt.startswith('<') and alt.endswith('>')

        if (
            end_pos > pos
            or end_pos == pos
            or variant_length > LARGE_VARIANT_THRESHOLD
            or svtype in {'INV', 'DEL', 'DUP', 'CNV', 'INS', 'BND'}
            or is_symbolic_alt
        ):
            use_sv_region = True

    if use_sv_region:
        # VCF POS and END are 1-based; convert to 0-based before applying padding
        # so the plot window aligns with pysam / renderer 0-based coordinates.
        start = max(0, (pos - 1) - plot_padding)
        end = (end_pos - 1) + plot_padding
    else:
        # SNV or small indel
        ref_len = len(ref) if ref else 1
        alt_alleles = alt.split(',') if alt else ['']
        max_alt_len = max(len(a) for a in alt_alleles) if alt_alleles else 0
        variant_end = pos + max(ref_len, max_alt_len) - 1

        start = max(0, (pos - 1) - plot_padding)
        end = (variant_end - 1) + plot_padding

    return (chrom, start, end)

def _calculate_vcf_target_region(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    info: Optional[dict] = None,
) -> Tuple[str, int, int]:
    """
    Calculate the original target region of a VCF record (without padding).

    Returns:
        (chrom, target_start, target_end)
    """
    LARGE_VARIANT_THRESHOLD = 50

    end_pos = None
    if info:
        if 'SVTYPE' in info and 'SVLEN' in info:
            try:
                svlen = abs(int(str(info['SVLEN']).split(',')[0]))
                end_pos = pos + svlen
            except (ValueError, TypeError):
                pass

        if 'END' in info:
            try:
                end_pos = int(info['END'])
            except (ValueError, TypeError):
                pass

    variant_length = 0
    if info and 'SVLEN' in info:
        try:
            variant_length = abs(int(str(info['SVLEN']).split(',')[0]))
        except (ValueError, TypeError):
            variant_length = 0

    svtype = str(info.get('SVTYPE', '')).upper() if info else ''
    is_symbolic_alt = alt.startswith('<') and alt.endswith('>')

    use_sv_region = False
    if end_pos is not None:
        if (
            end_pos > pos
            or end_pos == pos
            or variant_length > LARGE_VARIANT_THRESHOLD
            or svtype in {'INV', 'DEL', 'DUP', 'CNV', 'INS', 'BND'}
            or is_symbolic_alt
        ):
            use_sv_region = True

    if use_sv_region:
        # VCF POS and END are 1-based; convert to 0-based half-open for the
        # renderer which uses pysam / 0-based coordinates throughout.
        target_start = pos - 1
        target_end = max(end_pos, pos + 1) - 1
    else:
        ref_len = len(ref) if ref else 1
        alt_alleles = alt.split(',') if alt else ['']
        max_alt_len = max(len(a) for a in alt_alleles) if alt_alleles else 1
        variant_end = pos + max(ref_len, max_alt_len) - 1

        # VCF POS is 1-based; convert to 0-based so focus_region aligns with
        # insertion / mismatch markers which are drawn at 0-based ref positions.
        target_start = pos - 1
        target_end = max(variant_end, pos + 1) - 1

    return chrom, target_start, target_end

def parse_regions_file(regions_path: str, padding: Optional[int] = None) -> List[Tuple[str, int, int, int, int]]:
    """
    Parse BED or VCF file to extract genomic regions.
    
    Args:
        regions_path: Path to BED or VCF file
        padding: Optional user-defined padding added to both sides of each BED/VCF region.
                 If None, padding is inferred from each BED/VCF record.
        
    Returns:
        List of (chrom, plot_start, plot_end, target_start, target_end) tuples
        
    BED format:
        - Columns: chrom, start, end, ...
        - Start is 0-based, end is exclusive
        - Returns: (chrom, plot_start, plot_end, target_start, target_end)
        
    VCF format:
        - Columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, ...
        - POS is 1-based
        - Supports standard VCF and VCFv4.3+ formats
        - Handles structural variants with END tag in INFO field
        - Handles SVTYPE and SVLEN for structural variants
        - Calculates region based on REF/ALT lengths
        - Returns: (chrom, plot_start, plot_end, target_start, target_end)
    """
    if not os.path.exists(regions_path):
        raise FileNotFoundError(f"Regions file not found: {regions_path}")
    
    regions: List[Tuple[str, int, int, int, int]] = []
    
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
                    chrom, plot_start, plot_end = _calculate_vcf_region(chrom, pos, ref, alt, info, padding)
                    _, target_start, target_end = _calculate_vcf_target_region(chrom, pos, ref, alt, info)
                    regions.append((chrom, plot_start, plot_end, target_start, target_end))
                except (ValueError, IndexError) as e:
                    # Skip malformed lines
                    continue
            else:
                # BED format: chrom, start, end, ...
                if len(parts) < 3:
                    continue
                
                chrom = parts[0]
                try:
                    target_start = int(parts[1])   # 0-based
                    target_end = int(parts[2])     # exclusive
                    if target_start < 0 or target_end <= target_start:
                        continue

                    plot_padding = padding if padding is not None else _infer_bed_padding(target_start, target_end)
                    plot_start = max(0, target_start - plot_padding)
                    plot_end = target_end + plot_padding

                    regions.append((chrom, plot_start, plot_end, target_start, target_end))
                except ValueError:
                    continue
    
    return regions

