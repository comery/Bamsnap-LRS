#!/usr/bin/env python3

import argparse
import gzip
import random


def main():
    parser = argparse.ArgumentParser(
        description="Randomly assign phase-set blocks to an existing VCF.gz."
    )
    parser.add_argument("-i", "--input", required=True, help="Input VCF.gz")
    parser.add_argument("-o", "--output", required=True, help="Output VCF")
    parser.add_argument(
        "--min-snps",
        type=int,
        default=3,
        help="Minimum SNP count per block",
    )
    parser.add_argument(
        "--max-snps",
        type=int,
        default=8,
        help="Maximum SNP count per block",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=208,
        help="Random seed",
    )
    args = parser.parse_args()

    if args.min_snps < 1:
        raise ValueError("--min-snps must be at least 1")

    if args.max_snps < args.min_snps:
        raise ValueError("--max-snps must be >= --min-snps")

    rng = random.Random(args.seed)

    current_chrom = None
    current_ps = None
    block_size = 0
    snps_in_block = 0
    ps_header_found = False

    with gzip.open(args.input, "rt") as fin, open(
        args.output, "w", encoding="utf-8"
    ) as fout:

        for line in fin:
            if line.startswith("##FORMAT=<ID=PS,"):
                ps_header_found = True
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                if not ps_header_found:
                    fout.write(
                        '##FORMAT=<ID=PS,Number=1,Type=Integer,'
                        'Description="Artificial phase set; position of '
                        'the first SNP in the block">\n'
                    )

                fout.write(line)
                continue

            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = fields[0]
            pos = fields[1]

            # Start a new block at a new chromosome or when the
            # current block has reached its randomly selected size.
            if chrom != current_chrom or snps_in_block >= block_size:
                current_chrom = chrom
                current_ps = pos
                block_size = rng.randint(args.min_snps, args.max_snps)
                snps_in_block = 0

            format_fields = fields[8].split(":")

            if "PS" in format_fields:
                ps_index = format_fields.index("PS")
            else:
                format_fields.append("PS")
                ps_index = len(format_fields) - 1

            fields[8] = ":".join(format_fields)

            # Add or replace PS in every sample column.
            for sample_index in range(9, len(fields)):
                sample_fields = fields[sample_index].split(":")

                while len(sample_fields) < len(format_fields):
                    sample_fields.append(".")

                sample_fields[ps_index] = current_ps
                fields[sample_index] = ":".join(sample_fields)

            fout.write("\t".join(fields) + "\n")
            snps_in_block += 1


if __name__ == "__main__":
    main()