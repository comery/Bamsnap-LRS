"""
make_example.py
===============
Generate a small, self-contained Bamsnap-LRS *highlight* example.

Outputs (all written next to this script, or to --outdir):
    example_ref.fa        500 bp single-contig reference  (chr1)
    example.vcf           single phased sample, 2 phase blocks + scattered SNPs
    example_reads.sam     long reads spanning the two blocks and their junction

Design (all coordinates 0-based unless noted; VCF POS is 1-based):
    chrom            : chr1, length 500
    Phase block A    : PS=1001, SNP sites at 0-based 60,95,130,165,200
                       sample GT "1|0"  -> hap1 = ALT, hap2 = REF
    Phase block B    : PS=2001, SNP sites at 0-based 300,335,370,405,440
                       sample GT "0|1"  -> hap1 = REF, hap2 = ALT
    Scattered SNPs   : 0-based 25, 250, 470   (unphased "0/1", no PS)
    Junction         : the A/B boundary sits around 200..300; most reads are
                       built to span this so the two blocks' linkage is visible.

The reads carry deliberate haplotype combinations across the two blocks
(A-hap1 + B-hap1, A-hap2 + B-hap2, and the two "switched" combinations) plus a
few short single-block reads and a couple of noisy reads with extra non-site
mismatches / a small deletion, so the renderer changes can be checked.
"""
import os
import random
import argparse

CHROM = "scaffold1"
REF_LEN = 500
SEED = 7

# (0-based pos, REF, ALT, phase_set or None, gt) ; gt is the sample GT string.
# For phased blocks GT slot order is hap1|hap2.
SITES = [
    # scattered (unphased) -------------------------------------------------
    (25,  None, None, None, "0/1"),
    # block A (PS=60) ----------------------------------------------------
    (60,  None, None, "60", "1|0"),
    (95,  None, None, "60", "1|0"),
    (130, None, None, "60", "1|0"),
    (165, None, None, "60", "1|0"),
    (200, None, None, "60", "1|0"),
    # scattered (unphased), sits in the junction ---------------------------
    (250, None, None, None, "0/1"),
    # block B (PS=300) ----------------------------------------------------
    (300, None, None, "300", "0|1"),
    (335, None, None, "300", "0|1"),
    (370, None, None, "300", "0|1"),
    (405, None, None, "300", "0|1"),
    (440, None, None, "300", "0|1"),
    # scattered (unphased) -------------------------------------------------
    (470, None, None, None, "0/1"),
]

BASES = "ACGT"


def _rand_ref(n, rng):
    return "".join(rng.choice(BASES) for _ in range(n))


def _pick_alt(ref_base, rng):
    return rng.choice([b for b in BASES if b != ref_base])


def build_sites(ref, rng):
    """Fill in REF/ALT for every site from the reference sequence."""
    out = []
    for pos, _, _, ps, gt in SITES:
        ref_base = ref[pos]
        alt_base = _pick_alt(ref_base, rng)
        out.append({
            "pos": pos, "ref": ref_base, "alt": alt_base,
            "ps": ps, "gt": gt,
        })
    return out


def hap_base_of(site, hap_idx):
    """Return the allele base for hap slot `hap_idx` of a phased site."""
    slots = site["gt"].replace("|", "/").split("/")
    if hap_idx >= len(slots):
        return None
    a = slots[hap_idx]
    if a == "0":
        return site["ref"]
    if a == "1":
        return site["alt"]
    return None


def write_reference(path, ref):
    with open(path, "w") as f:
        f.write(f">{CHROM}\n")
        for i in range(0, len(ref), 60):
            f.write(ref[i:i + 60] + "\n")


def write_vcf(path, sites, sample="SAMPLE1"):
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={CHROM},length={REF_LEN}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                   "FILTER", "INFO", "FORMAT", sample]),
    ]
    for s in sites:
        if s["ps"]:
            fmt, val = "GT:PS", f'{s["gt"]}:{s["ps"]}'
        else:
            fmt, val = "GT", s["gt"]
        lines.append("\t".join([
            CHROM, str(s["pos"] + 1), ".", s["ref"], s["alt"],
            ".", "PASS", ".", fmt, val,
        ]))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------
# Each read spec: (qname, ref_start, ref_end, hapA, hapB, n_noise, del_spec)
#   hapA / hapB : which hap (1, 2) the read carries over block A / block B
#                 sites it covers; None means "carry REF everywhere".
#   n_noise     : number of extra random non-site mismatches (tests graying).
#   del_spec    : optional (del_start, del_len) small deletion, or None.
READ_SPECS = [
    # The TRUE linkage in this dataset across the two phase blocks is:
    #   block A hap1  <->  block B hap2
    #   block A hap2  <->  block B hap1
    # so every junction-spanning read carries one of those two combinations.
    # Spans are spread across the whole 0..500 scaffold, including reads that
    # run almost edge to edge.
    # ---- spanning, A-hap1 linked to B-hap2 ----
    ("rA1B2_1",  0,  500, 1, 2, 0, None),
    ("rA1B2_2",  5,  480, 1, 2, 1, None),
    ("rA1B2_3", 20,  470, 1, 2, 0, (220, 4)),
    ("rA1B2_4", 40,  495, 1, 2, 0, None),
    ("rA1B2_5", 55,  460, 1, 2, 2, None),
    ("rA1B2_6", 10,  490, 1, 2, 0, None),
    ("rA1B2_7", 70,  455, 1, 2, 0, None),
    ("rA1B2_8", 30,  475, 1, 2, 1, None),
    # ---- spanning, A-hap2 linked to B-hap1 ----
    ("rA2B1_1",  0,  499, 2, 1, 0, None),
    ("rA2B1_2",  8,  485, 2, 1, 1, None),
    ("rA2B1_3", 25,  470, 2, 1, 0, None),
    ("rA2B1_4", 45,  498, 2, 1, 0, (240, 6)),
    ("rA2B1_5", 50,  465, 2, 1, 2, None),
    ("rA2B1_6", 15,  492, 2, 1, 0, None),
    ("rA2B1_7", 80,  450, 2, 1, 0, None),
    ("rA2B1_8", 35,  488, 2, 1, 1, None),
    # ---- short reads covering only block A ----
    ("rA1_only_1",  0, 230, 1, None, 0, None),
    ("rA1_only_2", 20, 215, 1, None, 1, None),
    ("rA2_only_1", 10, 240, 2, None, 0, None),
    ("rA2_only_2", 35, 225, 2, None, 0, None),
    # ---- short reads covering only block B ----
    ("rB1_only_1", 270, 499, None, 1, 0, None),
    ("rB1_only_2", 290, 485, None, 1, 1, None),
    ("rB2_only_1", 275, 498, None, 2, 0, None),
    ("rB2_only_2", 295, 480, None, 2, 0, None),
    # ---- reads that match NO single hap (conflicting alleles in block A) ----
    ("rMixed_1",  20,  480, "mix", None, 1, None),
    ("rMixed_2",  60,  470, "mix", None, 0, None),
]


def build_read_seq(ref, spec, sites, rng):
    qname, rs, re, hapA, hapB, n_noise, del_spec = spec
    seq = list(ref[rs:re])  # read covers ref[rs:re)

    def set_site(site, hap, mix_counter=[0]):
        if hap == "mix":
            # Alternate hap1 / hap2 alleles across the block's sites so the
            # read is consistent with neither hap.
            slot = mix_counter[0] % 2
            mix_counter[0] += 1
            base = hap_base_of(site, slot)
            if base is None:
                base = site["ref"]
        elif hap is None:
            base = site["ref"]            # carry REF
        else:
            base = hap_base_of(site, hap - 1)  # hap is 1/2 -> slot 0/1
            if base is None:
                base = site["ref"]
        idx = site["pos"] - rs
        if 0 <= idx < len(seq):
            seq[idx] = base

    for site in sites:
        if not (rs <= site["pos"] < re):
            continue
        if site["ps"] == "1001":
            set_site(site, hapA)
        elif site["ps"] == "2001":
            set_site(site, hapB)
        else:
            # scattered unphased site: carry REF (kept gray by the renderer,
            # but still colored as a highlight cell at its position)
            set_site(site, None)

    # Extra random non-site mismatches (to test that non-target mismatches
    # render in the same gray as matches).
    site_pos = {s["pos"] for s in sites}
    placed = 0
    attempts = 0
    while placed < n_noise and attempts < 200:
        attempts += 1
        gp = rng.randint(rs + 2, re - 3)
        if gp in site_pos:
            continue
        idx = gp - rs
        seq[idx] = _pick_alt(seq[idx], rng)
        placed += 1

    seq = "".join(seq)

    # CIGAR: full-length match, optionally with a small deletion.
    if del_spec is not None:
        dstart, dlen = del_spec
        # local read-coordinate of the deletion (deletion consumes ref only)
        left = dstart - rs
        right_ref = (re - rs) - (left + dlen)
        # remove the deleted ref bases from the read sequence
        seq = seq[:left] + seq[left + dlen:]
        cigar = f"{left}M{dlen}D{right_ref}M"
    else:
        cigar = f"{re - rs}M"

    return qname, rs, re, seq, cigar


def write_sam(path, ref, sites, rng):
    header = [
        "@HD\tVN:1.6\tSO:coordinate",
        f"@SQ\tSN:{CHROM}\tLN:{REF_LEN}",
        "@PG\tID:make_example\tPN:make_example",
    ]
    rows = []
    for spec in READ_SPECS:
        qname, rs, re, seq, cigar = build_read_seq(ref, spec, sites, rng)
        flag = 0  # forward, primary
        rows.append("\t".join([
            qname, str(flag), CHROM, str(rs + 1), "60", cigar,
            "*", "0", "0", seq, "*",
        ]))
    rows.sort(key=lambda r: int(r.split("\t")[3]))
    with open(path, "w") as f:
        f.write("\n".join(header + rows) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(SEED)

    ref = _rand_ref(REF_LEN, rng)
    sites = build_sites(ref, rng)

    write_reference(os.path.join(args.outdir, "example_ref.fa"), ref)
    write_vcf(os.path.join(args.outdir, "example.vcf"), sites)
    write_sam(os.path.join(args.outdir, "example_reads.sam"), ref, sites, rng)

    print("Wrote example_ref.fa, example.vcf, example_reads.sam to", args.outdir)
    print("\nSite summary (0-based pos | REF>ALT | PS | GT | hap1/hap2):")
    for s in sites:
        h1 = hap_base_of(s, 0) if s["ps"] else "-"
        h2 = hap_base_of(s, 1) if s["ps"] else "-"
        print(f"  {s['pos']:>3}  {s['ref']}>{s['alt']}  PS={s['ps'] or '.':<5} "
              f"GT={s['gt']:<4} hap1={h1} hap2={h2}")


if __name__ == "__main__":
    main()
