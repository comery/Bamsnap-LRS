"""
Highlight module - Parses VCF files for target-site highlighting.

This module extracts SNP sites (and small variants) from a user-provided VCF
so the renderer can:
  1. Display a dedicated 'VCF track' above the reads, showing per-sample
     REF/ALT alleles as colored bars (IGV-style).
  2. On the reads themselves, only color mismatches that overlap these
     highlighted positions (other mismatches are rendered in pale gray to
     reduce visual clutter, making it easier to trace SNP linkage on reads).

Data model:
  - HighlightSite: one genomic position (0-based) that came from the VCF.
    Stores ref allele and a per-sample list of ALT alleles (only alleles that
    differ from REF are kept). For a multi-sample VCF, the samples list
    preserves the column order from the VCF header.
"""
import os
import gzip
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class SampleCall:
    """Resolved genotype for one sample at one site.

    Four GT states (used by REF/ALT row rendering when unphased):
      * hom-REF  (0/0)   -> has_ref=True,  alt_base=None
      * het      (0/1)   -> has_ref=True,  alt_base="<ALT>"
      * hom-ALT  (1/1)   -> has_ref=False, alt_base="<ALT>"
      * missing  (./.)   -> has_ref=False, alt_base=None

    Phasing-aware fields (set when GT uses '|' AND PS is provided):
      is_phased  : True iff GT used '|' (any number of '|'-separated alleles)
                   AND PS tag is present
      hap_bases  : list of allele bases in PS order, one per '|'-separated
                   slot. e.g. "0|1|2" with ALT="T,G" -> ["A", "T", "G"]
                   None entries mean '.' (missing on that hap).
                   Length is the sample's ploidy at this site, so a diploid
                   has 2 entries, a triploid has 3, a haploid has 1.
      phase_set  : PS tag value, or None when unphased

    Legacy fields (still present for diploid back-compat / unphased row):
      hap1_base  = hap_bases[0] if available else None
      hap2_base  = hap_bases[1] if available else None
    """
    has_ref: bool = False
    alt_base: Optional[str] = None
    # Phasing-aware fields
    is_phased: bool = False
    hap_bases: List[Optional[str]] = field(default_factory=list)
    phase_set: Optional[str] = None

    # Back-compat convenience properties
    @property
    def hap1_base(self) -> Optional[str]:
        return self.hap_bases[0] if len(self.hap_bases) > 0 else None

    @property
    def hap2_base(self) -> Optional[str]:
        return self.hap_bases[1] if len(self.hap_bases) > 1 else None

    @property
    def ploidy(self) -> int:
        """Number of alleles in this sample's GT at this site (>=1)."""
        return max(1, len(self.hap_bases))


@dataclass
class HighlightSite:
    """A single highlighted genomic position from the VCF.

    Attributes
    ----------
    chrom : str
        Chromosome name.
    pos : int
        0-based genomic position (VCF POS - 1).
    ref : str
        REF allele (uppercase, single base for SNPs).
    sample_calls : Dict[str, SampleCall]
        Per-sample resolved genotype. See SampleCall for the four possible
        states. Previously this field was ``sample_alts`` mapping to a single
        allele string — that could not distinguish hom-REF from hom-ALT from
        missing, which mattered for the REF/ALT row semantics.
    is_snv : bool
        True when both REF and all ALTs are single bases (simple SNV).
        Used to decide whether to draw REF/ALT base-color cells (SNV only)
        or just a location-marker bar (indel).
    """
    chrom: str
    pos: int
    ref: str
    sample_calls: Dict[str, SampleCall] = field(default_factory=dict)
    is_snv: bool = True

    # ---- Backwards-compat shim ---------------------------------------------
    # Older code / external callers may still look at `.sample_alts` (maps to
    # the ALT allele or None). Provide it as a read-only derived property so
    # any third-party code that imported HighlightSite doesn't break.
    @property
    def sample_alts(self) -> Dict[str, Optional[str]]:
        return {k: v.alt_base for k, v in self.sample_calls.items()}


def _open_maybe_gzip(path: str):
    """Open a plain or gzipped text file."""
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _resolve_sample_call(gt_field: str, ref: str, alt_alleles: List[str],
                         phase_set: Optional[str] = None) -> SampleCall:
    """Resolve a VCF GT string into a SampleCall.

    Splits the GT on '|' (phased) or '/' (unphased). For phased GTs with a
    PS, EVERY hap allele is recovered into SampleCall.hap_bases — including
    triploid and higher (your trio data uses 0|1|0 in the child sample).

    Examples (ref="A", alts=["T","G"]):
      "0/0"           -> hap_bases=[]              (unphased; no per-hap info)
      "0|1"   PS=37   -> hap_bases=["A","T"]
      "0|1|0" PS=37   -> hap_bases=["A","T","A"]   (triploid)
      "1/2"           -> hap_bases=[]
      "1"             -> hap_bases=["T"]           (haploid; unphased single)
      "./."           -> hap_bases=[]
    """
    if not gt_field or gt_field == ".":
        return SampleCall(False, None)

    is_phased_sep = "|" in gt_field and phase_set is not None
    raw = gt_field.replace("|", "/")
    parts = raw.split("/")

    def _idx_to_base(s: str) -> Optional[str]:
        if s in ("", "."):
            return None
        try:
            i = int(s)
        except ValueError:
            return None
        if i == 0:
            return ref.upper()
        if 1 <= i <= len(alt_alleles):
            return alt_alleles[i - 1].upper()
        return None

    # has_ref + alt_base (legacy fields, used for unphased rendering)
    has_ref = False
    alt_idx_seen = None
    for p in parts:
        if p in ("", "."):
            continue
        try:
            idx = int(p)
        except ValueError:
            continue
        if idx == 0:
            has_ref = True
        elif alt_idx_seen is None:
            alt_idx_seen = idx
    alt_base = None
    if alt_idx_seen is not None:
        pos = alt_idx_seen - 1
        if 0 <= pos < len(alt_alleles):
            alt_base = alt_alleles[pos].upper()

    # Per-hap bases for ALL ploidy levels (only when phased+PS)
    hap_bases: List[Optional[str]] = []
    if is_phased_sep and len(parts) >= 1:
        hap_bases = [_idx_to_base(p) for p in parts]

    return SampleCall(
        has_ref=has_ref,
        alt_base=alt_base,
        is_phased=bool(is_phased_sep and any(h is not None for h in hap_bases)),
        hap_bases=hap_bases,
        phase_set=phase_set if is_phased_sep else None,
    )


def parse_highlight_vcf(
    vcf_path: str,
    chrom: str,
    start: int,
    end: int,
) -> Tuple[List[HighlightSite], List[str]]:
    """Parse a VCF and return highlight sites falling within a region.

    Parameters
    ----------
    vcf_path : str
        Path to the .vcf or .vcf.gz file.
    chrom, start, end : str, int, int
        Region of interest (0-based, half-open: [start, end)).

    Returns
    -------
    (sites, samples)
        sites   : list of HighlightSite sorted by position.
        samples : list of sample names, in VCF-column order. Empty list if
                  the VCF has no sample columns (sites-only VCF); in that
                  case a single synthetic sample "VCF" is added to each site
                  carrying the first ALT, so the track still renders.
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"Highlight VCF not found: {vcf_path}")

    sites: List[HighlightSite] = []
    samples: List[str] = []

    with _open_maybe_gzip(vcf_path) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            # Header line with sample names
            if line.startswith("#CHROM"):
                # Columns 0..8 are fixed VCF fields; columns 9+ are samples.
                cols = line.lstrip("#").split("\t")
                if len(cols) > 9:
                    samples = cols[9:]
                continue
            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                continue

            v_chrom = parts[0]
            if v_chrom != chrom:
                continue
            try:
                v_pos_1b = int(parts[1])  # VCF is 1-based
            except ValueError:
                continue
            v_pos_0b = v_pos_1b - 1
            if v_pos_0b < start or v_pos_0b >= end:
                continue

            ref = parts[3].upper()
            alt_field = parts[4]
            if alt_field in (".", ""):
                continue  # no variant reported
            alt_alleles = [a.upper() for a in alt_field.split(",") if a]
            if not alt_alleles:
                continue

            is_snv = len(ref) == 1 and all(len(a) == 1 for a in alt_alleles)

            # Per-sample genotype extraction
            sample_calls: Dict[str, SampleCall] = {}
            if samples and len(parts) >= 10:
                # FORMAT is parts[8]; find GT and (optionally) PS index.
                fmt_keys = parts[8].split(":") if len(parts) > 8 else []
                try:
                    gt_idx = fmt_keys.index("GT")
                except ValueError:
                    gt_idx = 0
                ps_idx = None
                if "PS" in fmt_keys:
                    ps_idx = fmt_keys.index("PS")
                for i, sname in enumerate(samples):
                    col = 9 + i
                    if col >= len(parts):
                        sample_calls[sname] = SampleCall(False, None)
                        continue
                    sample_fields = parts[col].split(":")
                    gt_val = sample_fields[gt_idx] if gt_idx < len(sample_fields) else "."
                    ps_val = None
                    if ps_idx is not None and ps_idx < len(sample_fields):
                        raw_ps = sample_fields[ps_idx]
                        if raw_ps and raw_ps != "." and raw_ps != "0":
                            ps_val = raw_ps
                    sample_calls[sname] = _resolve_sample_call(gt_val, ref, alt_alleles, ps_val)
            else:
                # Sites-only VCF: treat as unphased het carrying first ALT
                sample_calls["VCF"] = SampleCall(has_ref=True, alt_base=alt_alleles[0])

            sites.append(HighlightSite(
                chrom=v_chrom,
                pos=v_pos_0b,
                ref=ref,
                sample_calls=sample_calls,
                is_snv=is_snv,
            ))

    # If the VCF had no sample columns we still want a consistent samples list
    # for the renderer to iterate over.
    if not samples and sites:
        samples = ["VCF"]

    sites.sort(key=lambda s: s.pos)
    return sites, samples


def build_highlight_index(sites: List[HighlightSite]) -> Dict[int, HighlightSite]:
    """Build an O(1) position lookup: 0-based pos -> HighlightSite.

    Used by the renderer to quickly test whether a given reference position
    is a highlighted site while walking over read segments.
    """
    return {s.pos: s for s in sites}
