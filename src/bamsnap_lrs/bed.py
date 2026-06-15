import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple


@dataclass
class BedBlock:
    """Represents a block (exon) in a BED feature"""
    start: int
    end: int


@dataclass
class BedFeature:
    """Represents a BED feature"""
    chrom: str
    start: int  # 0-based, start position
    end: int    # 0-based, end position (exclusive)
    name: Optional[str] = None
    score: Optional[int] = None  # 0-1000
    strand: Optional[str] = None  # + or -
    thick_start: Optional[int] = None  # thickStart
    thick_end: Optional[int] = None    # thickEnd
    item_rgb: Optional[Tuple[int, int, int]] = None  # RGB color
    blocks: List[BedBlock] = field(default_factory=list)  # blockCount, blockSizes, blockStarts


def parse_bed(bed_path: str, chrom: str, start: int, end: int) -> List[BedFeature]:
    """Parse BED file and extract features within range"""
    if not os.path.exists(bed_path):
        return []

    features: List[BedFeature] = []

    with open(bed_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            # Required columns: chrom, start, end
            r_chrom = parts[0]
            if r_chrom != chrom:
                continue
            
            try:
                r_start = int(parts[1])  # BED is 0-based
                r_end = int(parts[2])    # BED end is exclusive
            except ValueError:
                continue
            
            # Check range overlap
            if r_end <= start or r_start >= end:
                continue
            
            # Parse optional columns
            name = parts[3] if len(parts) > 3 else None
            score = None
            if len(parts) > 4 and parts[4] != '.':
                try:
                    score = int(float(parts[4]))  # BED score can be float
                except ValueError:
                    pass
            
            strand = parts[5] if len(parts) > 5 and parts[5] in ['+', '-', '.'] else None
            if strand == '.':
                strand = None
            
            thick_start = None
            thick_end = None
            if len(parts) > 6 and parts[6] != '.':
                try:
                    thick_start = int(parts[6])
                except ValueError:
                    pass
            if len(parts) > 7 and parts[7] != '.':
                try:
                    thick_end = int(parts[7])
                except ValueError:
                    pass
            
            # Parse RGB color
            item_rgb = None
            if len(parts) > 8 and parts[8] != '.':
                try:
                    rgb_parts = parts[8].split(',')
                    if len(rgb_parts) == 3:
                        item_rgb = (int(rgb_parts[0]), int(rgb_parts[1]), int(rgb_parts[2]))
                except (ValueError, IndexError):
                    pass
            
            # Parse blocks (exons)
            blocks: List[BedBlock] = []
            if len(parts) > 11:  # blockCount, blockSizes, blockStarts
                try:
                    block_count = int(parts[9])
                    block_sizes = [int(x) for x in parts[10].rstrip(',').split(',') if x]
                    block_starts = [int(x) for x in parts[11].rstrip(',').split(',') if x]
                    
                    if len(block_sizes) == block_count and len(block_starts) == block_count:
                        for i in range(block_count):
                            block_start = r_start + block_starts[i]
                            block_end = block_start + block_sizes[i]
                            blocks.append(BedBlock(block_start, block_end))
                except (ValueError, IndexError):
                    pass
            
            feature = BedFeature(
                chrom=r_chrom,
                start=r_start,
                end=r_end,
                name=name,
                score=score,
                strand=strand,
                thick_start=thick_start,
                thick_end=thick_end,
                item_rgb=item_rgb,
                blocks=blocks
            )
            features.append(feature)
    
    # Sort features by start position
    features.sort(key=lambda x: x.start)
    return features

