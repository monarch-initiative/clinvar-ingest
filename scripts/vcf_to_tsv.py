#!/usr/bin/env python3
"""
Convert ClinVar VCF to TSV format for Koza processing.

Uses cyvcf2 to parse VCF and expand INFO fields into separate columns.
"""

import argparse
import sys
from pathlib import Path

from cyvcf2 import VCF


# INFO fields to extract (based on transform.yaml columns)
INFO_FIELDS = [
    "AF_ESP",
    "AF_EXAC",
    "AF_TGP",
    "ALLELEID",
    "CLNDN",
    "CLNDNINCL",
    "CLNDISDB",
    "CLNDISDBINCL",
    "CLNHGVS",
    "CLNREVSTAT",
    "CLNSIG",
    "CLNSIGCONF",
    "CLNSIGINCL",
    "CLNVC",
    "CLNVCSO",
    "CLNVI",
    "DBVARID",
    "GENEINFO",
    "MC",
    "ONCDN",
    "ONCDNINCL",
    "ONCDISDB",
    "ONCDISDBINCL",
    "ONC",
    "ONCINCL",
    "ONCREVSTAT",
    "ONCCONF",
    "ORIGIN",
    "RS",
    "SCIDN",
    "SCIDNINCL",
    "SCIDISDB",
    "SCIDISDBINCL",
    "SCIREVSTAT",
    "SCI",
    "SCIINCL",
]

# Standard VCF columns
VCF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]


def get_info_value(variant, field):
    """Extract INFO field value, returning '.' if missing."""
    try:
        value = variant.INFO.get(field)
        if value is None:
            return "."
        if isinstance(value, tuple):
            return "|".join(str(v) for v in value)
        return str(value)
    except KeyError:
        return "."


def convert_vcf_to_tsv(vcf_path: Path, output_path: Path):
    """Convert VCF file to TSV with expanded INFO fields."""
    vcf = VCF(str(vcf_path))

    # Build header
    header = VCF_COLUMNS + INFO_FIELDS

    with open(output_path, "w") as out:
        # Write header
        out.write("\t".join(header) + "\n")

        # Process each variant
        for variant in vcf:
            row = [
                variant.CHROM,
                str(variant.POS),
                variant.ID or ".",
                variant.REF,
                ",".join(variant.ALT) if variant.ALT else ".",
                str(variant.QUAL) if variant.QUAL else ".",
                variant.FILTER or "PASS",
            ]

            # Add INFO fields
            for field in INFO_FIELDS:
                row.append(get_info_value(variant, field))

            out.write("\t".join(row) + "\n")

    vcf.close()


def main():
    parser = argparse.ArgumentParser(
        description="Convert ClinVar VCF to TSV format for Koza processing"
    )
    parser.add_argument("vcf", type=Path, help="Input VCF file (can be gzipped)")
    parser.add_argument("output", type=Path, help="Output TSV file")

    args = parser.parse_args()

    if not args.vcf.exists():
        print(f"Error: VCF file not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)

    convert_vcf_to_tsv(args.vcf, args.output)
    print(f"Converted {args.vcf} -> {args.output}")


if __name__ == "__main__":
    main()
