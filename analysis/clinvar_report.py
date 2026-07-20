#!/usr/bin/env python3
"""Exploratory report on ClinVar review-star cutoffs and classification mix.

Combines two independent analyses into one HTML report:

1. Review-star cutoff impact (per-submission evidence, from
   submission_summary.txt): `clinvar_helpers.variant_records_to_disease()`
   drops any submission record whose review status maps to fewer than
   `star_min` stars (production currently hardcodes `var2disease_star_min =
   3`). Section 1 re-runs that same disease mapping at star_min in
   {0, 1, 2, 3, 4} and reports, for each threshold, distinct variants and
   distinct (gene, disease) pairs. These counts intentionally ignore the
   separate HPO-overlap requirement in process_row() (which gates whether
   the production transform emits ANY edge at all) so the star cutoff's
   effect can be isolated.

   Note: star_min=2 and star_min=3 always produce identical results here
   because submission_summary.txt only ever uses six literal ReviewStatus
   values (practice guideline / reviewed by expert panel / criteria
   provided, single submitter / no assertion criteria provided / no
   classification provided / flagged submission) -- "2 stars" (criteria
   provided, multiple submitters, no conflicts) is an aggregate,
   cross-submitter status that only appears in ClinVar's variant-level
   files (see section 3, which uses CLNREVSTAT from the VCF directly and
   *does* see this tier), never per-submission-record here.

2. Multi-submitter concordance rescue: gene-disease pairs whose ONLY
   supporting evidence is star 0/1 (never reach star_min=2) are additionally
   checked for >=2 distinct Submitters independently reporting the *same*
   MONDO disease with the *exact same* ClinicalSignificance (Pathogenic/
   Likely pathogenic family only), regardless of each submitter's own
   review status. This reconstructs a proxy for the missing per-record
   "2 star" tier from raw submitter agreement. Each such pair is also
   annotated with whether any contributing record used
   CollectionMethod == "literature only" (a publication proxy) --
   informational only, not used as a filter.

3. Clinical significance & review-status crossfilter (variant-level, from
   clinvar.vcf.gz directly): CLNREVSTAT here IS ClinVar's aggregate
   per-variant review status (so all 5 star levels are real, including
   "2 stars"). CLNSIG is a free-form, sometimes pipe/slash-combined string
   (100+ distinct values observed) bucketed via classify_clnsig() into
   P / LP / P-LP / B / LB / B-LB / VUS / Conflicting / Other / Not
   classified. An interactive crossfilter widget lets you toggle
   significance / star / variant-type / size / literature / concordance /
   STRchive filters and see live variant counts. The STRchive dimension flags
   whether a variant's genomic footprint overlaps one of STRchive's ~80
   curated pathogenic short-tandem-repeat loci (github.com/dashnowlab/
   STRchive, fetched and cached on first run), independent of ClinVar's own
   CLNVC "Microsatellite" label.

4. Structural variants / CNVs -- what's NOT in the VCF: clinvar.vcf.gz (the
   only input sections 1-3 use) requires a fixed genomic position + REF/ALT,
   which structurally excludes copy-number gain/loss, translocations, and
   other large rearrangements. ClinVar publishes those separately, in
   variant_summary.txt.gz (same tab_delimited/ directory as
   submission_summary.txt, Start/Stop coordinate ranges instead of a fixed
   REF/ALT). This section downloads that file (~440MB, auto-fetched on first
   run) and reports real counts, size distributions, and disease-linkage
   rates for those SV/CNV types -- including the finding that only ~6-10%
   of copy-number gain/loss entries resolve to a MONDO/MedGen/OMIM disease
   id at all (the rest are ClinGen-style dosage-sensitivity regions with no
   single named condition), and that GeneSymbol degrades to unparseable
   prose ("covers 42 genes, none of which curated...") for large CNVs,
   unlike the clean delimited gene lists small variants get.

Requires data downloaded/preprocessed per the top-level pipeline (run from
the repo root):
    just download
    just preprocess

Usage (run from analysis/, matching analysis/justfile's working directory):
    cd analysis && PYTHONPATH=../src uv run --project .. python clinvar_report.py
    # or: just --justfile analysis/justfile report
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import requests

from clinvar_helpers import (
    format_id_to_map,
    make_genes_from_row,
    make_medgen_to_mondo_map,
    make_mondo_map,
    make_variant_record_map,
    predicate_map,
    review_star_map,
    variant_records_to_disease,
)

STAR_LEVELS = [0, 1, 2, 3, 4]
PRODUCTION_STAR_MIN = 3
EXPERT_STAR_MIN = 2  # threshold above which a pair no longer needs "rescuing"
# star levels enumerated in the gene-disease-pair detail table (the full
# breakdown for 0/1 is tens of thousands of rows -- not useful to browse)
ENUMERATED_LEVELS = [2, 4]
MIN_CONCORDANT_SUBMITTERS = 2
MAX_VARIANT_SAMPLE = 8

CLNSIG_BUCKETS = ["P", "LP", "P/LP", "VUS", "LB", "B", "B/LB", "Conflicting", "Other", "Not classified"]
NOT_CLASSIFIED_VALUES = {
    ".",
    "not_provided",
    "no_classification_provided",
    "no_classification_for_the_single_variant",
    "no_classifications_from_unflagged_records",
}

LOC_GENE_RE = re.compile(r"^LOC\d+$", re.IGNORECASE)


def drop_loc_genes(gene_ids: list, gene_symbols: list) -> tuple[list, list]:
    """Drop NCBI placeholder "LOC" gene symbols (LOC + Entrez GeneID, e.g.
    LOC109461479) assigned to loci with no curated official gene symbol.
    ClinVar's GENEINFO field lists every locus overlapping a variant's
    position, so a repeat-expansion variant in a real gene like HTT also
    tags neighboring LOC-symbol predicted loci -- not useful signal, always
    excluded here rather than left as a toggle."""
    kept = [(gid, sym) for gid, sym in zip(gene_ids, gene_symbols) if not LOC_GENE_RE.match(sym)]
    if not kept:
        return [], []
    ids, syms = zip(*kept)
    return list(ids), list(syms)


# ---------------------------------------------------------------------------
# Section 0: illustrative examples (variant types; one gene, many diseases)
# ---------------------------------------------------------------------------

EXAMPLE_VARIANT_IDS = {
    "snv": "55630",  # BRCA1 NC_000017.11:g.43045711G>C, Pathogenic, 4-star
    "indel_del": "53158",  # CFTR NC_000007.14:g.117480086_117480108del, Pathogenic, 4-star
    "indel_ins": "973964",  # RPE65 NC_000001.11:g.68431505_68431506insCAGC, Pathogenic, 3-star
    "str": "183387",  # FMR1 CGG repeat expansion NC_000023.11:g.147912051CGG[201], Fragile X syndrome
}
# Structural variants -- not in clinvar.vcf.gz, only in variant_summary.txt.gz (see section 4)
EXAMPLE_SV_IDS = {
    "cnv_del_large": "2579266",  # 22q11.21(chr22:18985739-21081116)x1, ~2.1Mb, DiGeorge syndrome
    "cnv_dup_large": "4846746",  # 22q11.21(chr22:18932238-20324016)x3, ~1.4Mb, Congenital heart disease
    "cnv_del_gene": "2579216",  # DMD Xp21.1(chrX:31601428-31727780)x0, ~126kb, Duchenne muscular dystrophy
    "cnv_dup_gene": "58614",  # DMD Xp21.1(chrX:32662366-32758964)x2, ~97kb, phenotype not resolved in this record
    "inversion": "90611",  # MSH2 inversion, ~38kb, Lynch syndrome, 3-star (expert panel)
}

GENE_MODEL_EXAMPLES = [
    {
        "gene": "LMNA",
        "heading": "One gene, many diseases: LMNA",
        "blurb": (
            "<strong>LMNA</strong> encodes a nuclear envelope protein and is one of the clearest examples "
            "of <em>locus heterogeneity going the other direction</em> &mdash; instead of many genes causing "
            "one disease, different variants in this single gene cause strikingly different conditions "
            '("laminopathies"): muscular dystrophy, dilated cardiomyopathy, peripheral neuropathy '
            "(Charcot-Marie-Tooth), lipodystrophy, and the premature-aging disorder Hutchinson-Gilford "
            "progeria syndrome, among others. Below, every Pathogenic/Likely pathogenic SNV in LMNA "
            "({count:,} variants) is plotted at its real genomic position against the gene's actual exon "
            "structure (GRCh38, via Ensembl), colored by the disease it's associated with. Hover a point for "
            "details, click to open its ClinVar record."
        ),
    },
    {
        "gene": "SMCHD1",
        "diseases": {"MONDO:0008031", "MONDO:0011323"},  # FSHD2, BAMS only
        "heading": "Same gene, opposite mutation patterns: SMCHD1",
        "blurb": (
            "<strong>SMCHD1</strong> illustrates a different kind of pattern: two diseases from the same gene, "
            "with opposite mutation <em>distributions</em>. Facioscapulohumeral muscular dystrophy type 2 "
            "(FSHD2) is caused by loss-of-function variants scattered across almost the entire gene &mdash; "
            "disease arises from haploinsufficiency, so nearly any variant that breaks the protein works, "
            "wherever it falls. Bosma arhinia microphthalmia syndrome (BAMS), in contrast, is caused by "
            "missense variants tightly clustered in the N-terminal ATPase (GHKL) domain, encoded by a small "
            "run of early exons &mdash; only a change to that specific region produces the effect behind BAMS. "
            "Below, {count:,} Pathogenic/Likely pathogenic SNVs are plotted the same way as above: real "
            "genomic position, real exon structure, colored by disease. In the real downloaded data, FSHD2 "
            "variants span the region's full width while BAMS variants stay confined to a narrow band of "
            "early exons, nested well inside the FSHD2 range."
        ),
    },
]
ENSEMBL_REST_URL = "https://rest.ensembl.org"
# top N diseases get their own legend color; everything else buckets into "Other"
GENE_MODEL_TOP_DISEASES = 8
GENE_MODEL_COLORS = [
    "#4e79a7",
    "#f28e2b",
    "#e15759",
    "#76b7b2",
    "#59a14f",
    "#edc948",
    "#b07aa1",
    "#ff9da7",
]
GENE_MODEL_OTHER_COLOR = "#9ca3af"


def load_example_variants(clinvar_tsv: Path) -> dict:
    """Pull the full row (HGVS, review status, disease names, etc.) for each
    hand-picked EXAMPLE_VARIANT_IDS entry straight out of clinvar.tsv, keyed by
    label ("snv"/"indel_del"/...), so the mutation-type gallery always reflects
    the actual downloaded release rather than a hardcoded description that
    could drift from the real record."""
    wanted = {v: k for k, v in EXAMPLE_VARIANT_IDS.items()}
    found = {}
    with open(clinvar_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["ID"] in wanted:
                found[wanted[row["ID"]]] = row
            if len(found) == len(wanted):
                break
    return found


# ---------------------------------------------------------------------------
# Mutation-type cartoon gallery -- simplified schematic diagrams (not real
# genomic-coordinate plots, unlike the gene models below) illustrating the
# mechanism of each broad class of variation ClinVar covers, each paired with
# a real hand-picked example (EXAMPLE_VARIANT_IDS / EXAMPLE_SV_IDS above).
# ---------------------------------------------------------------------------

CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE = "#e2e8f0", "#94a3b8"
CARTOON_SUB_FILL, CARTOON_SUB_STROKE = "#fef3c7", "#d97706"
CARTOON_INS_FILL, CARTOON_INS_STROKE = "#dcfce7", "#16a34a"
CARTOON_DEL_FILL, CARTOON_DEL_STROKE = "#fee2e2", "#dc2626"
CARTOON_DUP_FILL, CARTOON_DUP_STROKE = "#dbeafe", "#2563eb"
CARTOON_INV_FILL, CARTOON_INV_STROKE = "#f3e8ff", "#9333ea"
CARTOON_STR_FILL, CARTOON_STR_STROKE = "#cffafe", "#0891b2"


def _cartoon_svg(width=260, height=130):
    return (
        f"<svg viewBox='0 0 {width} {height}' width='{width}' height='{height}' "
        f"style='width:100%; height:auto; max-width:{width}px; display:block;' "
        f"xmlns='http://www.w3.org/2000/svg' font-family='ui-sans-serif, system-ui'>"
    )


def _cartoon_box(x, y, w, h, label, fill, stroke, text_color="#0f172a", font_size=11, rx=4, dashed=False):
    dash = " stroke-dasharray='3,2'" if dashed else ""
    return (
        f"<rect x='{x:.1f}' y='{y:.1f}' width='{w:.1f}' height='{h:.1f}' fill='{fill}' "
        f"stroke='{stroke}' stroke-width='1.5' rx='{rx}'{dash}/>"
        f"<text x='{x + w / 2:.1f}' y='{y + h / 2 + font_size / 3:.1f}' font-size='{font_size}' "
        f"text-anchor='middle' fill='{text_color}'>{label}</text>"
    )


def _cartoon_label(x, y, text, size=11, color="#64748b", anchor="start", weight="normal"):
    return (
        f"<text x='{x}' y='{y}' font-size='{size}' fill='{color}' text-anchor='{anchor}' "
        f"font-weight='{weight}'>{text}</text>"
    )


def _cartoon_caption(text, width=260, y=122):
    return _cartoon_label(width / 2, y, text, size=10.5, color="#64748b", anchor="middle")


def cartoon_snv():
    parts = [_cartoon_svg()]
    seq = ["A", "T", "G", "C", "G", "A", "T"]
    alt_base = "T"
    alt_idx = 3
    box_w, box_h, gap = 24, 26, 3
    start_x = 70
    y_ref, y_alt = 18, 64
    parts.append(_cartoon_label(10, y_ref + box_h / 2 + 4, "Ref"))
    parts.append(_cartoon_label(10, y_alt + box_h / 2 + 4, "Alt"))
    for i, base in enumerate(seq):
        x = start_x + i * (box_w + gap)
        fill, stroke = (
            (CARTOON_SUB_FILL, CARTOON_SUB_STROKE) if i == alt_idx else (CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE)
        )
        parts.append(_cartoon_box(x, y_ref, box_w, box_h, base, fill, stroke))
        alt_label = alt_base if i == alt_idx else base
        parts.append(_cartoon_box(x, y_alt, box_w, box_h, alt_label, fill, stroke))
    cx = start_x + alt_idx * (box_w + gap) + box_w / 2
    parts.append(
        f"<line x1='{cx}' y1='{y_ref + box_h + 2}' x2='{cx}' y2='{y_alt - 4}' stroke='{CARTOON_SUB_STROKE}' "
        f"stroke-width='1.5' marker-end='url(#arrow-sub)'/>"
        f"<defs><marker id='arrow-sub' markerWidth='6' markerHeight='6' refX='3' refY='3' orient='auto'>"
        f"<path d='M0,0 L6,3 L0,6 Z' fill='{CARTOON_SUB_STROKE}'/></marker></defs>"
    )
    parts.append(_cartoon_caption("single-base substitution (G&gt;C)"))
    parts.append("</svg>")
    return "".join(parts)


def cartoon_indel(kind):
    parts = [_cartoon_svg()]
    box_w, box_h, gap = 20, 26, 2.5
    start_x = 55
    y_ref, y_alt = 18, 64
    if kind == "deletion":
        seq = ["A", "T", "G", "C", "G", "A", "T", "C"]
        del_range = (3, 5)  # indices 3,4 removed
        parts.append(_cartoon_label(10, y_ref + box_h / 2 + 4, "Ref"))
        parts.append(_cartoon_label(10, y_alt + box_h / 2 + 4, "Alt"))
        for i, base in enumerate(seq):
            x = start_x + i * (box_w + gap)
            highlighted = del_range[0] <= i < del_range[1]
            fill, stroke = (
                (CARTOON_DEL_FILL, CARTOON_DEL_STROKE)
                if highlighted
                else (CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE)
            )
            parts.append(_cartoon_box(x, y_ref, box_w, box_h, base, fill, stroke, dashed=highlighted))
        vi = 0
        for i, base in enumerate(seq):
            if del_range[0] <= i < del_range[1]:
                continue
            x = start_x + vi * (box_w + gap)
            parts.append(_cartoon_box(x, y_alt, box_w, box_h, base, CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
            vi += 1
        gap_x = start_x + del_range[0] * (box_w + gap) - gap / 2
        parts.append(
            f"<line x1='{gap_x:.1f}' y1='{y_alt - 4}' x2='{gap_x:.1f}' y2='{y_alt + box_h + 4}' "
            f"stroke='{CARTOON_DEL_STROKE}' stroke-width='1.5' stroke-dasharray='2,2'/>"
        )
        parts.append(_cartoon_caption("2 bp deleted"))
    else:
        seq = ["A", "T", "G", "C", "G", "A"]
        ins_after = 3
        inserted = ["T", "T"]
        parts.append(_cartoon_label(10, y_ref + box_h / 2 + 4, "Ref"))
        parts.append(_cartoon_label(10, y_alt + box_h / 2 + 4, "Alt"))
        for i, base in enumerate(seq):
            x = start_x + i * (box_w + gap)
            parts.append(_cartoon_box(x, y_ref, box_w, box_h, base, CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
        vi = 0
        for i, base in enumerate(seq):
            x = start_x + vi * (box_w + gap)
            parts.append(_cartoon_box(x, y_alt, box_w, box_h, base, CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
            vi += 1
            if i == ins_after - 1:
                for ib in inserted:
                    x = start_x + vi * (box_w + gap)
                    parts.append(_cartoon_box(x, y_alt, box_w, box_h, ib, CARTOON_INS_FILL, CARTOON_INS_STROKE))
                    vi += 1
        parts.append(_cartoon_caption("2 bp inserted"))
    parts.append("</svg>")
    return "".join(parts)


def cartoon_str():
    parts = [_cartoon_svg()]
    box_w, box_h, gap = 30, 24, 3
    start_x = 70
    y_ref, y_alt = 18, 64
    parts.append(_cartoon_label(10, y_ref + box_h / 2 + 4, "Ref"))
    parts.append(_cartoon_label(10, y_alt + box_h / 2 + 4, "Exp"))
    for i in range(5):
        x = start_x + i * (box_w + gap)
        parts.append(_cartoon_box(x, y_ref, box_w, box_h, "CGG", CARTOON_STR_FILL, CARTOON_STR_STROKE, font_size=9))
    for i in range(3):
        x = start_x + i * (box_w + gap)
        parts.append(_cartoon_box(x, y_alt, box_w, box_h, "CGG", CARTOON_STR_FILL, CARTOON_STR_STROKE, font_size=9))
    dots_x = start_x + 3 * (box_w + gap)
    parts.append(_cartoon_label(dots_x + box_w / 2, y_alt + box_h / 2 + 4, "&#8226;&#8226;&#8226;", anchor="middle"))
    parts.append(_cartoon_caption("~5 repeats &rarr; 200+ repeats"))
    parts.append("</svg>")
    return "".join(parts)


def cartoon_large_cnv(kind):
    parts = [_cartoon_svg()]
    box_w, box_h, gap = 20, 22, 3
    lcr_w = 7
    y_ref, y_alt = 22, 68
    x0 = 45
    genes = ["g1", "g2", "g3", "g4", "g5"]

    def lcr(x, y):
        return f"<path d='M{x},{y} L{x + lcr_w},{y + box_h / 2} L{x},{y + box_h} Z' fill='#94a3b8'/>"

    x = x0
    parts.append(_cartoon_label(10, y_ref + box_h / 2 + 4, "Ref"))
    parts.append(_cartoon_box(x, y_ref, box_w, box_h, genes[0], CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
    x += box_w + gap
    lcr1_x = x
    parts.append(lcr(x, y_ref))
    x += lcr_w + gap
    for g in genes[1:4]:
        parts.append(_cartoon_box(x, y_ref, box_w, box_h, g, CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
        x += box_w + gap
    lcr2_x = x
    parts.append(lcr(x, y_ref))
    x += lcr_w + gap
    parts.append(_cartoon_box(x, y_ref, box_w, box_h, genes[4], CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
    parts.append(
        _cartoon_label(
            (lcr1_x + lcr2_x) / 2 + lcr_w, y_ref - 6, "flanking low-copy repeats", size=8.5, anchor="middle"
        )
    )

    if kind == "del":
        x = x0
        parts.append(_cartoon_label(10, y_alt + box_h / 2 + 4, "Alt"))
        parts.append(_cartoon_box(x, y_alt, box_w, box_h, genes[0], CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
        x += box_w + gap
        parts.append(lcr(x, y_alt))
        x += lcr_w + gap
        parts.append(_cartoon_box(x, y_alt, box_w, box_h, genes[4], CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
        parts.append(_cartoon_caption("g2–g4 deleted (NAHR at flanking repeats)"))
    else:
        color, stroke = CARTOON_DUP_FILL, CARTOON_DUP_STROKE
        x = x0
        parts.append(_cartoon_label(10, y_alt + box_h / 2 + 4, "Alt"))
        parts.append(_cartoon_box(x, y_alt, box_w, box_h, genes[0], CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
        x += box_w + gap
        parts.append(lcr(x, y_alt))
        x += lcr_w + gap
        for g in genes[1:4]:
            parts.append(_cartoon_box(x, y_alt, box_w, box_h, g, color, stroke))
            x += box_w + gap
        for g in genes[1:4]:
            parts.append(_cartoon_box(x, y_alt, box_w, box_h, g, color, stroke))
            x += box_w + gap
        parts.append(lcr(x, y_alt))
        x += lcr_w + gap
        parts.append(_cartoon_box(x, y_alt, box_w, box_h, genes[4], CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
        parts.append(_cartoon_caption("g2–g4 duplicated (extra copy)"))
    parts.append("</svg>")
    return "".join(parts)


def cartoon_exon_cnv(kind):
    parts = [_cartoon_svg()]
    box_w, box_h, gap = 15, 20, 2
    y_ref, y_alt = 18, 68
    x0 = 45
    n = 8
    affected = (3, 6)  # exons 4,5,6 (0-idx 3,4,5)

    def track(y, exons, colors):
        out = [
            f"<line x1='{x0 - 6}' y1='{y + box_h / 2}' x2='{x0 + len(exons) * (box_w + gap)}' "
            f"y2='{y + box_h / 2}' stroke='#94a3b8' stroke-width='1.5'/>"
        ]
        x = x0
        for i, label in exons:
            fill, stroke = colors.get(i, (CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
            out.append(_cartoon_box(x, y, box_w, box_h, label, fill, stroke, font_size=7.5))
            x += box_w + gap
        return "".join(out)

    parts.append(_cartoon_label(10, y_ref + box_h / 2 + 4, "Ref"))
    ref_exons = [(i, str(i + 1)) for i in range(n)]
    ref_colors = {
        i: (CARTOON_DEL_FILL, CARTOON_DEL_STROKE)
        if kind == "del" and affected[0] <= i < affected[1]
        else (CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE)
        for i in range(n)
    }
    parts.append(track(y_ref, ref_exons, ref_colors))

    parts.append(_cartoon_label(10, y_alt + box_h / 2 + 4, "Alt"))
    if kind == "del":
        alt_exons = [(i, str(i + 1)) for i in range(n) if not (affected[0] <= i < affected[1])]
        alt_colors = {i: (CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE) for i, _ in alt_exons}
        parts.append(track(y_alt, alt_exons, alt_colors))
        parts.append(_cartoon_caption("exons 4–6 deleted"))
    else:
        alt_exons = [(i, str(i + 1)) for i in range(n)]
        for i in range(*affected):
            alt_exons.append((i, str(i + 1) + "'"))
        alt_colors = {}
        alt_exons_indexed = []
        for idx, (_orig_i, label) in enumerate(alt_exons):
            alt_exons_indexed.append((idx, label))
            alt_colors[idx] = (
                (CARTOON_DUP_FILL, CARTOON_DUP_STROKE) if "'" in label else (CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE)
            )
        parts.append(track(y_alt, alt_exons_indexed, alt_colors))
        parts.append(_cartoon_caption("exons 4–6 duplicated"))
    parts.append("</svg>")
    return "".join(parts)


def cartoon_inversion():
    parts = [_cartoon_svg()]
    box_w, box_h, gap = 32, 24, 4
    y_ref, y_alt = 18, 68
    x0 = 40
    labels = ["A", "B", "C", "D"]

    def chevron(x, y, color, flip=False):
        d = (
            f"M{x},{y} L{x + 6},{y + box_h / 2} L{x},{y + box_h}"
            if not flip
            else f"M{x + 6},{y} L{x},{y + box_h / 2} L{x + 6},{y + box_h}"
        )
        return f"<path d='{d}' fill='none' stroke='{color}' stroke-width='2'/>"

    parts.append(_cartoon_label(10, y_ref + box_h / 2 + 4, "Ref"))
    x = x0
    for lab in labels:
        parts.append(_cartoon_box(x, y_ref, box_w, box_h, lab, CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE))
        parts.append(chevron(x + box_w + 1, y_ref, "#94a3b8"))
        x += box_w + gap + 8

    parts.append(_cartoon_label(10, y_alt + box_h / 2 + 4, "Alt"))
    order = ["A", "C", "B", "D"]
    x = x0
    for lab in order:
        inverted = lab in ("B", "C")
        fill, stroke = (CARTOON_INV_FILL, CARTOON_INV_STROKE) if inverted else (CARTOON_NEUTRAL_FILL, CARTOON_NEUTRAL_STROKE)
        parts.append(_cartoon_box(x, y_alt, box_w, box_h, lab, fill, stroke))
        parts.append(chevron(x + box_w + 1, y_alt, CARTOON_INV_STROKE if inverted else "#94a3b8", flip=inverted))
        x += box_w + gap + 8
    parts.append(
        f"<path d='M120,10 Q145,-2 170,10' fill='none' stroke='{CARTOON_INV_STROKE}' stroke-width='1.5' "
        f"marker-end='url(#arrow-inv)'/>"
        f"<defs><marker id='arrow-inv' markerWidth='6' markerHeight='6' refX='3' refY='3' orient='auto'>"
        f"<path d='M0,0 L6,3 L0,6 Z' fill='{CARTOON_INV_STROKE}'/></marker></defs>"
    )
    parts.append(_cartoon_caption("B–C segment inverted"))
    parts.append("</svg>")
    return "".join(parts)


MUTATION_GALLERY = [
    {"key": "snv", "kind": "tsv", "title": "SNV", "cartoon": cartoon_snv},
    {"key": "indel_del", "kind": "tsv", "title": "Indel — deletion", "cartoon": lambda: cartoon_indel("deletion")},
    {"key": "indel_ins", "kind": "tsv", "title": "Indel — insertion", "cartoon": lambda: cartoon_indel("insertion")},
    {"key": "str", "kind": "tsv", "title": "Microsatellite / STR expansion", "cartoon": cartoon_str},
    {
        "key": "cnv_del_large",
        "kind": "sv",
        "title": "CNV — large recurrent deletion",
        "cartoon": lambda: cartoon_large_cnv("del"),
        "gene_label": "22q11.21",
    },
    {
        "key": "cnv_dup_large",
        "kind": "sv",
        "title": "CNV — large recurrent duplication",
        "cartoon": lambda: cartoon_large_cnv("dup"),
        "gene_label": "22q11.21",
    },
    {
        "key": "cnv_del_gene",
        "kind": "sv",
        "title": "CNV — single-gene deletion (few exons)",
        "cartoon": lambda: cartoon_exon_cnv("del"),
        "gene_label": "DMD",
    },
    {
        "key": "cnv_dup_gene",
        "kind": "sv",
        "title": "CNV — single-gene duplication (few exons)",
        "cartoon": lambda: cartoon_exon_cnv("dup"),
        "gene_label": "DMD",
    },
    {"key": "inversion", "kind": "sv", "title": "Inversion", "cartoon": cartoon_inversion, "gene_label": "MSH2"},
]


def ensure_ensembl_gene_model(data_dir: Path, gene_symbol: str) -> Path:
    """Fetch and cache this gene's canonical-transcript exon structure (GRCh38,
    matching clinvar.tsv's own coordinate system) from Ensembl's REST API --
    used to draw the gene model in the "one gene, many diseases" illustration."""
    path = data_dir / f"ensembl_{gene_symbol}.json"
    if path.exists():
        return path
    url = f"{ENSEMBL_REST_URL}/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1;content-type=application/json"
    print(f"Downloading {url} (one-time)...")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    path.write_bytes(resp.content)
    return path


def load_gene_model(data_dir: Path, gene_symbol: str) -> dict:
    """Canonical-transcript exon coordinates (GRCh38) for gene_symbol, as a
    {chrom, strand, tx_start, tx_end, exons: [(start, end), ...]} dict."""
    path = ensure_ensembl_gene_model(data_dir, gene_symbol)
    gene = json.loads(path.read_text())
    canonical_id = gene["canonical_transcript"].split(".")[0]
    transcript = next(t for t in gene["Transcript"] if t["id"] == canonical_id)
    exons = sorted((e["start"], e["end"]) for e in transcript["Exon"])
    return {
        "chrom": gene["seq_region_name"],
        "strand": gene["strand"],
        "tx_start": transcript["start"],
        "tx_end": transcript["end"],
        "exons": exons,
    }


def build_gene_model_variants(
    clinvar_tsv: Path,
    var_records: dict,
    map_to_mondo: dict,
    mondo_labels: dict,
    gene_symbol: str,
    disease_filter: set | None = None,
) -> list[dict]:
    """Every Pathogenic/Likely-pathogenic SNV in gene_symbol with >=1 resolved
    MONDO disease (star_min=0 -- any submission evidence, matching section 3's
    own "any P/LP evidence" pair-building pass), each tagged with one primary
    disease (the lexicographically-first MONDO id among those mapped, or the
    first one that's also in disease_filter if given) purely so a single-gene,
    single-color-per-point plot has something deterministic to color by; a
    handful of variants map to >1 disease and this ignores the rest for
    plotting purposes only. If disease_filter is given (a set of MONDO ids),
    variants whose resolved diseases don't intersect it are dropped entirely
    -- used to zoom a multi-disease gene down to just the diseases being
    contrasted (e.g. SMCHD1's FSHD2 vs. BAMS), ignoring rarer/near-duplicate
    disease terms that would otherwise clutter the legend."""
    variants = []
    with open(clinvar_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["CLNVC"] != "single_nucleotide_variant":
                continue
            _, gene_syms = make_genes_from_row(row["GENEINFO"])
            if gene_symbol not in gene_syms:
                continue
            records = var_records.get(row["ID"])
            if not records:
                continue
            dis, _, _ = variant_records_to_disease(records, map_to_mondo, star_min=0)
            if not dis:
                continue
            if disease_filter is not None:
                dis = {d for d in dis if d in disease_filter}
                if not dis:
                    continue
            mondo_id = sorted(dis)[0]
            variants.append(
                {
                    "id": row["ID"],
                    "pos": int(row["POS"]),
                    "ref": row["REF"],
                    "alt": row["ALT"],
                    "mondo": mondo_id,
                    "disease_name": mondo_labels.get(mondo_id, mondo_id),
                }
            )
    return variants


def assign_gene_model_colors(variants: list[dict]) -> tuple[dict, list[dict]]:
    """Color the GENE_MODEL_TOP_DISEASES most-frequent diseases individually;
    bucket the long tail into one "Other" color. Returns (mondo_id -> color,
    legend rows sorted by count) -- the legend row for "Other" (if any) is
    always last regardless of its count, since it isn't one disease."""
    counts = Counter(v["mondo"] for v in variants)
    top = counts.most_common(GENE_MODEL_TOP_DISEASES)
    color_map = {mondo_id: GENE_MODEL_COLORS[i] for i, (mondo_id, _) in enumerate(top)}
    names = {v["mondo"]: v["disease_name"] for v in variants}
    legend = [{"mondo": mondo_id, "name": names[mondo_id], "color": color_map[mondo_id], "n": n} for mondo_id, n in top]
    other_n = sum(n for mondo_id, n in counts.items() if mondo_id not in color_map)
    if other_n:
        legend.append({"mondo": None, "name": "Other", "color": GENE_MODEL_OTHER_COLOR, "n": other_n})
    return color_map, legend


def render_gene_model_svg(gene_symbol: str, gene_model: dict, variants: list[dict], color_map: dict) -> str:
    """SVG gene model (exon boxes on a genomic coordinate line) with one
    color-coded, clickable "lollipop" per SNV at its real genomic position,
    linking out to the variant's ClinVar page. Variants that round to the same
    pixel column (recurrent-hotspot exons can have dozens within a few bp)
    are stacked vertically rather than overplotted, so each stays individually
    visible and clickable -- stack height doubles as a crude density read."""
    width, left_pad, right_pad = 920, 70, 30
    exon_height = 16
    dot_r = 3
    dot_spacing = 7
    plot_w = width - left_pad - right_pad

    tx_start, tx_end = gene_model["tx_start"], gene_model["tx_end"]
    span = tx_end - tx_start

    def x_of(pos: int) -> float:
        return left_pad + (pos - tx_start) / span * plot_w

    # bucket by rounded pixel column so near-coincident genomic positions stack
    # instead of overlapping; stable order (position, id) keeps repeat renders identical
    buckets: dict[int, list[dict]] = {}
    for v in sorted(variants, key=lambda v: (v["pos"], v["id"])):
        buckets.setdefault(round(x_of(v["pos"])), []).append(v)
    max_stack = max((len(b) for b in buckets.values()), default=1)

    stack_area = max_stack * dot_spacing + 20
    track_y = stack_area + 20
    height = track_y + exon_height / 2 + 40

    parts = [
        f"<svg viewBox='0 0 {width} {height:.0f}' width='{width}' height='{height:.0f}' "
        f"style='width:100%; height:auto; max-width:{width}px; display:block;' "
        f"xmlns='http://www.w3.org/2000/svg' font-family='ui-sans-serif, system-ui'>"
    ]

    # backbone (introns) + exon boxes
    parts.append(f"<line x1='{left_pad}' y1='{track_y}' x2='{width - right_pad}' y2='{track_y}' stroke='#94a3b8' stroke-width='2'/>")
    for e_start, e_end in gene_model["exons"]:
        x0, x1 = x_of(max(e_start, tx_start)), x_of(min(e_end, tx_end))
        parts.append(
            f"<rect x='{x0:.1f}' y='{track_y - exon_height / 2:.1f}' width='{max(x1 - x0, 1):.1f}' "
            f"height='{exon_height}' fill='#475569' rx='2'/>"
        )

    strand_label = "5' → 3'" if gene_model["strand"] == 1 else "3' ← 5'"
    parts.append(
        f"<text x='{left_pad}' y='{track_y + 34}' font-size='12' fill='#334155'>"
        f"{gene_symbol} &#8226; chr{gene_model['chrom']}:{tx_start:,}-{tx_end:,} &#8226; {strand_label}</text>"
    )

    # lollipops, stacked within each pixel bucket
    stem_base = track_y - exon_height / 2
    for x, bucket in buckets.items():
        for i, v in enumerate(bucket):
            y = stem_base - 4 - i * dot_spacing
            color = color_map.get(v["mondo"], GENE_MODEL_OTHER_COLOR)
            url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{v['id']}/"
            title = f"{gene_symbol}:g.{v['pos']}{v['ref']}&gt;{v['alt']} &#8226; {v['disease_name']} &#8226; VariationID {v['id']}"
            stem_y1 = stem_base if i == 0 else stem_base - (i - 1) * dot_spacing - 4
            parts.append(
                f"<a href='{url}' target='_blank' rel='noopener'><title>{title}</title>"
                f"<line x1='{x}' y1='{stem_y1:.1f}' x2='{x}' y2='{y:.1f}' "
                f"stroke='{color}' stroke-width='1' opacity='0.5'/>"
                f"<circle cx='{x}' cy='{y:.1f}' r='{dot_r}' fill='{color}' fill-opacity='0.85' "
                f"stroke='{color}' stroke-width='1'/></a>"
            )

    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# Section 1 & 2: review-star cutoff + multi-submitter concordance rescue
# ---------------------------------------------------------------------------


def load_maps(data_dir: Path):
    var_records = make_variant_record_map(str(data_dir / "submission_summary.txt.gz"))
    map_to_mondo = make_mondo_map(str(data_dir / "mondo.sssom.tsv"))
    medgen_to_mondo = make_medgen_to_mondo_map(str(data_dir / "MedGenIDMappings.txt.gz"))
    map_to_mondo.update(medgen_to_mondo)
    return var_records, map_to_mondo


def load_mondo_labels(data_dir: Path) -> dict:
    """mondo.sssom.tsv's subject_id/subject_label columns give a human-readable
    name for each MONDO id -- the same subject_id values used as disease ids
    throughout this report (see make_mondo_map())."""
    labels: dict = {}
    with open(data_dir / "mondo.sssom.tsv") as f:
        header = None
        for line in f:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if header is None:
                header = {k: i for i, k in enumerate(cols)}
                continue
            subject_id = cols[header["subject_id"]]
            subject_label = cols[header["subject_label"]]
            if subject_id.startswith("MONDO:") and subject_label:
                labels.setdefault(subject_id, subject_label)
    return labels


def concordance_groups(record_list: list, map_to_mondo: dict) -> dict:
    """Group Pathogenic/Likely-pathogenic records (any review status) by
    (mondo_disease_id, exact ClinicalSignificance) and collect the distinct
    Submitters that agree, plus whether any contributing record used a
    literature-only collection method. Mirrors variant_records_to_disease()'s
    per-record ReportedPhenotypeInfo -> SubmittedPhenotypeInfo fallback, but
    ignores review stars entirely (concordance is a star-independent signal).
    """
    groups: dict = {}
    for rec in record_list:
        clinsig = rec["ClinicalSignificance"]
        if clinsig not in predicate_map:
            continue

        mondo_ids = set()
        for mg_mapping in rec["ReportedPhenotypeInfo"].split(";"):
            mg_map = "MedGen:{}".format(mg_mapping.split(":")[0])
            if mg_map in map_to_mondo:
                mondo_ids.update(map_to_mondo[mg_map].keys())

        if not mondo_ids:
            for dis_id in rec["SubmittedPhenotypeInfo"].split(";"):
                dis_id = format_id_to_map(dis_id)
                if dis_id in map_to_mondo:
                    mondo_ids.update(map_to_mondo[dis_id].keys())
                elif dis_id is not None and "MONDO:" in dis_id:
                    mondo_ids.add(dis_id)

        is_literature = rec["CollectionMethod"] == "literature only"
        for mondo_id in mondo_ids:
            key = (mondo_id, clinsig)
            entry = groups.setdefault(key, {"submitters": set(), "has_literature": False})
            entry["submitters"].add(rec["Submitter"])
            entry["has_literature"] = entry["has_literature"] or is_literature

    return groups


def compute_star_data(clinvar_tsv: Path, var_records: dict, map_to_mondo: dict, mondo_labels: dict):
    """Single pass over clinvar.tsv computing summary counts (all star
    levels), per-pair variant sets (ENUMERATED_LEVELS), and the 0/1-star
    multi-submitter-concordance rescue analysis."""
    variant_sets = {s: set() for s in STAR_LEVELS}
    pair_sets = {s: set() for s in STAR_LEVELS}
    pair_variant_ids = {s: {} for s in ENUMERATED_LEVELS}

    # (gene_sym, gene_id, mondo_id) -> {"variants": set, "submitter_count": int, "has_literature": bool}
    sub_expert_pairs: dict = {}

    with open(clinvar_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            varid = row["ID"]
            if varid not in var_records:
                continue

            gene_ids, gene_symbols = drop_loc_genes(*make_genes_from_row(row["GENEINFO"]))
            records = var_records[varid]

            dis_by_star = {}
            for star_min in STAR_LEVELS:
                dis, _, _ = variant_records_to_disease(records, map_to_mondo, star_min=star_min)
                dis_by_star[star_min] = dis
                if not dis:
                    continue
                variant_sets[star_min].add(varid)
                for gene_id, gene_sym in zip(gene_ids, gene_symbols):
                    for disease_id in dis:
                        pair_sets[star_min].add((gene_id, disease_id))
                        if star_min in pair_variant_ids:
                            key = (gene_sym, gene_id, disease_id)
                            pair_variant_ids[star_min].setdefault(key, set()).add(varid)

            # diseases only reachable via star 0/1 (never reach the expert tier)
            sub_expert_diseases = set(dis_by_star[0]) - set(dis_by_star[EXPERT_STAR_MIN])
            if not sub_expert_diseases:
                continue

            groups = concordance_groups(records, map_to_mondo)
            disease_max_submitters: dict = {}
            disease_has_literature: dict = {}
            for (mondo_id, _clinsig), entry in groups.items():
                n = len(entry["submitters"])
                if n > disease_max_submitters.get(mondo_id, 0):
                    disease_max_submitters[mondo_id] = n
                disease_has_literature[mondo_id] = (
                    disease_has_literature.get(mondo_id, False) or entry["has_literature"]
                )

            for gene_id, gene_sym in zip(gene_ids, gene_symbols):
                for mondo_id in sub_expert_diseases:
                    key = (gene_sym, gene_id, mondo_id)
                    entry = sub_expert_pairs.setdefault(
                        key,
                        {
                            "variants": set(),
                            "submitter_count": 0,
                            "has_literature": False,
                        },
                    )
                    entry["variants"].add(varid)
                    entry["submitter_count"] = max(entry["submitter_count"], disease_max_submitters.get(mondo_id, 0))
                    entry["has_literature"] = entry["has_literature"] or disease_has_literature.get(mondo_id, False)

    counts = {s: {"variants": len(variant_sets[s]), "gene_disease_pairs": len(pair_sets[s])} for s in STAR_LEVELS}

    # Merge into one row per unique pair (at the lowest enumerated level),
    # with a variant count at each enumerated star level (None if absent).
    base_level = ENUMERATED_LEVELS[0]
    pairs = []
    for key in sorted(pair_variant_ids[base_level]):
        gene_sym, gene_id, mondo_id = key
        entry = {"gene": gene_sym, "gene_id": gene_id, "mondo": mondo_id}
        for level in ENUMERATED_LEVELS:
            varids = pair_variant_ids[level].get(key)
            entry[f"v{level}"] = len(varids) if varids else None
        pairs.append(entry)

    rescue_rows = []
    for (gene_sym, gene_id, mondo_id), entry in sub_expert_pairs.items():
        rescue_rows.append(
            {
                "gene": gene_sym,
                "gene_id": gene_id,
                "mondo": mondo_id,
                "disease_name": mondo_labels.get(mondo_id, ""),
                "n_variants": len(entry["variants"]),
                "submitter_count": entry["submitter_count"],
                "has_literature": entry["has_literature"],
                "rescued": entry["submitter_count"] >= MIN_CONCORDANT_SUBMITTERS,
            }
        )
    rescue_rows.sort(key=lambda r: (-r["rescued"], -r["submitter_count"], r["gene"]))

    rescued_rows = [r for r in rescue_rows if r["rescued"]]
    rescued_variant_ids = {
        vid for r in rescued_rows for vid in sub_expert_pairs[(r["gene"], r["gene_id"], r["mondo"])]["variants"]
    }
    rescue_summary = {
        "sub_expert_pairs_total": len(rescue_rows),
        "sub_expert_pairs_rescued": len(rescued_rows),
        "sub_expert_variants_total": len(variant_sets[0] - variant_sets[EXPERT_STAR_MIN]),
        "sub_expert_variants_rescued": len(rescued_variant_ids),
    }

    return counts, pairs, rescue_summary, rescue_rows


# ---------------------------------------------------------------------------
# Section 3: CLNSIG / CLNREVSTAT / CLNVC crossfilter
# ---------------------------------------------------------------------------


def classify_clnsig(value: str) -> str:
    if value in NOT_CLASSIFIED_VALUES:
        return "Not classified"

    tokens = {t for part in value.split("|") for t in part.split("/")}

    if any(t.startswith("Conflicting") for t in tokens):
        return "Conflicting"

    has_p = any(t == "Pathogenic" or t.startswith("Pathogenic,") for t in tokens)
    has_lp = any(t == "Likely_pathogenic" or t.startswith("Likely_pathogenic,") for t in tokens)
    has_b = any(t == "Benign" or t.startswith("Benign,") for t in tokens)
    has_lb = any(t == "Likely_benign" or t.startswith("Likely_benign,") for t in tokens)
    has_vus = any(
        t == "Uncertain_significance" or t.startswith("VUS") or t in ("Uncertain_risk_allele", "Likely_risk_allele")
        for t in tokens
    )

    if has_p and has_lp:
        return "P/LP"
    if has_p:
        return "P"
    if has_lp:
        return "LP"
    if has_b and has_lb:
        return "B/LB"
    if has_b:
        return "B"
    if has_lb:
        return "LB"
    if has_vus:
        return "VUS"
    return "Other"


SIZE_BUCKETS = ["1 bp (SNV)", "2-10 bp", "11-100 bp", "101-1,000 bp", "1,001-10,000 bp", "10,000+ bp"]


def classify_size(ref: str, alt: str) -> str:
    """Bucket by max(len(REF), len(ALT)) -- the physical allele footprint,
    independent of ClinVar's own Type/CLNVC label (which classifies by
    reference-sequence *context*, e.g. a 2bp deletion inside a tandem repeat
    gets called "Microsatellite" rather than "Deletion" -- see conversation).
    This lets you filter by actual size regardless of that label.
    """
    n = max(len(ref), len(alt))
    if n <= 1:
        return "1 bp (SNV)"
    if n <= 10:
        return "2-10 bp"
    if n <= 100:
        return "11-100 bp"
    if n <= 1000:
        return "101-1,000 bp"
    if n <= 10000:
        return "1,001-10,000 bp"
    return "10,000+ bp"


STRCHIVE_URL = "https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/STRchive-loci.json"


def ensure_strchive_downloaded(data_dir: Path) -> Path:
    """STRchive (github.com/dashnowlab/STRchive) curates ~80 short tandem repeat
    loci with established Mendelian-disease associations -- fetch its loci table
    once and cache it alongside the other downloaded inputs."""
    path = data_dir / "STRchive-loci.json"
    if path.exists():
        return path
    print(f"Downloading {STRCHIVE_URL} (one-time)...")
    resp = requests.get(STRCHIVE_URL, timeout=60)
    resp.raise_for_status()
    path.write_bytes(resp.content)
    return path


def load_strchive_intervals(data_dir: Path) -> dict[str, list[tuple[int, int]]]:
    """hg38 (start_hg38/stop_hg38, 1-based) intervals for each STRchive locus,
    grouped by bare chromosome name -- matching clinvar.tsv's own CHROM column,
    which (like the source VCF) carries no "chr" prefix."""
    path = ensure_strchive_downloaded(data_dir)
    loci = json.loads(path.read_text())
    intervals: dict[str, list[tuple[int, int]]] = {}
    for locus in loci:
        chrom = locus.get("chrom")
        start, stop = locus.get("start_hg38"), locus.get("stop_hg38")
        if not chrom or start is None or stop is None:
            continue
        intervals.setdefault(chrom.removeprefix("chr"), []).append((start, stop))
    return intervals


def in_strchive(chrom: str, pos: int, ref_len: int, intervals: dict) -> bool:
    """Does this variant's genomic footprint (POS .. POS+len(REF)-1) overlap any
    STRchive-curated pathogenic repeat locus on the same chromosome? Only ~80
    loci total, so a linear scan per variant is cheap."""
    var_start, var_end = pos, pos + ref_len - 1
    for locus_start, locus_end in intervals.get(chrom, ()):
        if var_start <= locus_end and var_end >= locus_start:
            return True
    return False


def build_clnsig_cube(
    clinvar_tsv: Path, var_records: dict, map_to_mondo: dict, mondo_labels: dict, strchive_intervals: dict
) -> tuple[list[dict], list[str], list[dict]]:
    """Full pass over clinvar.tsv (every row, not just those with a
    submission_summary match) tallying (star, clnsig_bucket, clnvc,
    has_literature, has_concordance), plus -- for rows with a
    submission_summary match -- the set of (gene, MONDO disease) pairs
    implied by that variant (any Pathogenic/Likely-pathogenic submission
    record, any star; same disease-mapping logic as section 1/2's star_min=0
    pass, since P/LP-family classification is a prerequisite for any disease
    mapping regardless of the variant's own aggregate CLNSIG/star shown
    here). Non-Pathogenic/Likely-pathogenic cells (VUS, B, LB, Conflicting,
    Other, Not classified) will always have an empty pair set -- production
    never creates disease edges for those, so that's the expected (not a
    bug) result of toggling those buckets on in the crossfilter.

    has_literature: does ANY submission record for this variant (any
    ClinicalSignificance) use CollectionMethod == "literature only"?
    has_concordance: does this variant have >=MIN_CONCORDANT_SUBMITTERS
    distinct Submitters agreeing on the same disease + exact same
    ClinicalSignificance (Pathogenic/Likely-pathogenic family only) -- the
    same signal used for section 2's rescue analysis. Variants absent from
    submission_summary.txt get False for both (undetermined).
    in_strchive: does this variant's genomic footprint overlap one of
    STRchive's ~80 curated pathogenic short-tandem-repeat loci (hg38
    coordinates)? See in_strchive() -- informational only, not used to gate
    disease-pair inclusion.

    Each cell also tracks, per pair, the exact contributing-variant count and
    a capped sample of ClinVar VariationIDs (MAX_VARIANT_SAMPLE) -- some
    pairs have thousands of contributing variants (BRCA1/BRCA2 top out
    around 3-5k), so embedding every id would bloat the report; the sample
    is still enough to click through to a few real examples, and the count
    is exact (not itself capped) so "+N more" in the UI is never wrong.
    """
    counts: Counter = Counter()
    clnvc_types: set = set()
    cell_pairs: dict = {}
    pair_index: dict = {}
    pair_gene_symbol: dict = {}

    with open(clinvar_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            star = review_star_map.get(row["CLNREVSTAT"])
            if star is None:
                star = 0
            clnsig = classify_clnsig(row["CLNSIG"])
            clnvc = row["CLNVC"] if row["CLNVC"] != "." else "unspecified"
            clnvc_types.add(clnvc)
            size = classify_size(row["REF"], row["ALT"])
            in_strchive_flag = in_strchive(row["CHROM"], int(row["POS"]), len(row["REF"]), strchive_intervals)

            records = var_records.get(row["ID"])
            has_literature = False
            has_concordance = False
            dis: dict = {}
            if records is not None:
                has_literature = any(rec["CollectionMethod"] == "literature only" for rec in records)
                groups = concordance_groups(records, map_to_mondo)
                has_concordance = any(
                    len(entry["submitters"]) >= MIN_CONCORDANT_SUBMITTERS for entry in groups.values()
                )
                dis, _, _ = variant_records_to_disease(records, map_to_mondo, star_min=0)

            cell_key = (star, clnsig, clnvc, size, has_literature, has_concordance, in_strchive_flag)
            counts[cell_key] += 1

            if not dis:
                continue
            gene_ids, gene_symbols = drop_loc_genes(*make_genes_from_row(row["GENEINFO"]))
            cell_map = cell_pairs.setdefault(cell_key, {})
            for gene_id, gene_sym in zip(gene_ids, gene_symbols):
                pair_gene_symbol.setdefault(gene_id, gene_sym)
                for mondo_id in dis:
                    # keyed by (gene_id, disease) only, matching section 1/2's canonical pair
                    # identity -- NCBIGene id is what production actually emits; a handful of
                    # ids (135 observed) carry >1 symbol spelling across records (gene renames),
                    # which would otherwise inflate this count vs section 1's headline number.
                    pkey = (gene_id, mondo_id)
                    idx = pair_index.setdefault(pkey, len(pair_index))
                    entry = cell_map.setdefault(idx, {"n": 0, "sample": []})
                    entry["n"] += 1
                    if len(entry["sample"]) < MAX_VARIANT_SAMPLE:
                        entry["sample"].append(row["ID"])

    cube = [
        {
            "star": star,
            "clnsig": clnsig,
            "clnvc": clnvc,
            "size": size,
            "literature": literature,
            "concordant": concordant,
            "strchive": strchive,
            "count": n,
            "pairs": [
                {"i": idx, "n": e["n"], "v": e["sample"]}
                for idx, e in sorted(
                    cell_pairs.get((star, clnsig, clnvc, size, literature, concordant, strchive), {}).items()
                )
            ],
        }
        for (star, clnsig, clnvc, size, literature, concordant, strchive), n in sorted(counts.items())
    ]

    pair_list = [None] * len(pair_index)
    for (gene_id, mondo_id), idx in pair_index.items():
        pair_list[idx] = {
            "gene": pair_gene_symbol.get(gene_id, ""),
            "gene_id": gene_id,
            "mondo": mondo_id,
            "disease_name": mondo_labels.get(mondo_id, ""),
        }

    return cube, sorted(clnvc_types), pair_list


# ---------------------------------------------------------------------------
# Section 4: structural variants / CNVs -- what's not in the VCF
# ---------------------------------------------------------------------------

VARIANT_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
# Type values from variant_summary.txt that never appear in the VCF's CLNVC field --
# the VCF requires a fixed genomic position + REF/ALT, which structurally excludes these.
SV_TYPES = ["copy number gain", "copy number loss", "Translocation", "Complex", "fusion", "Tandem duplication", "Inversion"]


def ensure_variant_summary_downloaded(data_dir: Path) -> Path:
    """variant_summary.txt.gz isn't part of the production download.yaml (it's not
    consumed by the ingest itself, only by this exploratory section) -- fetch it once,
    ~440MB, if not already present."""
    path = data_dir / "variant_summary.txt.gz"
    if path.exists():
        return path
    print(f"Downloading {VARIANT_SUMMARY_URL} (~440MB, one-time)...")
    with requests.get(VARIANT_SUMMARY_URL, stream=True, timeout=300) as resp:
        resp.raise_for_status()
        with open(path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)
    return path


def build_sv_summary(data_dir: Path) -> dict:
    """Full pass over variant_summary.txt.gz (GRCh38 rows only) tallying the Type
    distribution and, for the SV_TYPES specifically: count, how many resolve to a
    disease id (PhenotypeIDS not '-'/empty) vs not, genomic span (Stop - Start), and
    every row (resolved AND unresolved, any ClinicalSignificance) for the full
    paginated table -- not just a curated top-N. Each row carries a `resolved`
    boolean so the report can offer a "SV resolved to disease" filter checkbox
    rather than silently only ever showing the resolved subset.
    """
    path = ensure_variant_summary_downloaded(data_dir)

    type_counts: Counter = Counter()
    sv_stats = {t: {"count": 0, "resolved": 0, "spans": []} for t in SV_TYPES}
    sv_rows = []

    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        hcols = {k: i for i, k in enumerate(header)}
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if cols[hcols["Assembly"]] != "GRCh38":
                continue

            vtype = cols[hcols["Type"]]
            type_counts[vtype] += 1
            if vtype not in SV_TYPES:
                continue

            stat = sv_stats[vtype]
            stat["count"] += 1
            pheno_ids = cols[hcols["PhenotypeIDS"]]
            resolved = pheno_ids not in ("-", "")
            if resolved:
                stat["resolved"] += 1

            try:
                span = int(cols[hcols["Stop"]]) - int(cols[hcols["Start"]])
            except ValueError:
                span = None
            if span is not None:
                stat["spans"].append(span)

            if span is not None:
                sv_rows.append(
                    {
                        "name": cols[hcols["Name"]],
                        "type": vtype,
                        "allele_id": cols[hcols["AlleleID"]],
                        "variation_id": cols[hcols["VariationID"]],
                        "clinsig": cols[hcols["ClinicalSignificance"]],
                        "phenotype": cols[hcols["PhenotypeList"]],
                        "phenotype_ids": pheno_ids,
                        "resolved": resolved,
                        "span": span,
                        "review_status": cols[hcols["ReviewStatus"]],
                        "num_submitters": int(cols[hcols["NumberSubmitters"]]),
                    }
                )

    sv_rows.sort(key=lambda e: -e["span"])

    sv_type_rows = []
    for t in SV_TYPES:
        s = sv_stats[t]
        spans = sorted(s["spans"])
        median_span = spans[len(spans) // 2] if spans else None
        sv_type_rows.append(
            {
                "type": t,
                "count": s["count"],
                "resolved": s["resolved"],
                "resolved_pct": (100 * s["resolved"] / s["count"]) if s["count"] else 0,
                "span_min": min(spans) if spans else None,
                "span_median": median_span,
                "span_max": max(spans) if spans else None,
            }
        )

    return {
        "type_counts": dict(type_counts.most_common()),
        "sv_type_rows": sv_type_rows,
        "sv_rows": sv_rows,
        "total_variants": sum(type_counts.values()),
    }


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------


def render_html(
    results: dict,
    pairs: list,
    rescue_summary: dict,
    rescue_rows: list,
    cube: list,
    clnvc_types: list,
    pair_list: list,
    sv_summary: dict,
    example_variants: dict,
    gene_model_results: list[dict],
) -> str:
    max_variants = max(r["variants"] for r in results.values()) or 1
    max_pairs = max(r["gene_disease_pairs"] for r in results.values()) or 1

    chart_width, bar_height, gap, left_pad = 640, 32, 16, 90

    def bars(metric_key: str, max_val: int) -> str:
        rows = []
        for i, star in enumerate(STAR_LEVELS):
            value = results[star][metric_key]
            bar_w = (value / max_val) * (chart_width - left_pad - 60) if max_val else 0
            y = i * (bar_height + gap)
            highlight = " fill='#f59e0b'" if star == PRODUCTION_STAR_MIN else " fill='#2563eb'"
            rows.append(
                f"<text x='{left_pad - 10}' y='{y + bar_height / 2 + 5}' text-anchor='end' "
                f"font-size='13' fill='#334155'>{star} star{'s' if star != 1 else ''}</text>"
                f"<rect x='{left_pad}' y='{y}' width='{bar_w:.1f}' height='{bar_height}'{highlight} rx='3'/>"
                f"<text x='{left_pad + bar_w + 8}' y='{y + bar_height / 2 + 5}' font-size='13' "
                f"fill='#0f172a'>{value:,}</text>"
            )
        svg_height = len(STAR_LEVELS) * (bar_height + gap)
        return (
            f"<svg viewBox='0 0 {chart_width} {svg_height}' xmlns='http://www.w3.org/2000/svg'>"
            + "".join(rows)
            + "</svg>"
        )

    row_class = ' class="prod"'
    table_rows = "".join(
        f"<tr{row_class if star == PRODUCTION_STAR_MIN else ''}>"
        f"<td>{star}{' (production default)' if star == PRODUCTION_STAR_MIN else ''}</td>"
        f"<td>{results[star]['variants']:,}</td>"
        f"<td>{results[star]['gene_disease_pairs']:,}</td>"
        f"</tr>"
        for star in STAR_LEVELS
    )

    UNRESOLVED_PHENOTYPES = {"See cases", "not provided", "-", ""}

    def mutation_gallery_card(entry: dict) -> str:
        row = example_variants[entry["key"]]
        if entry["kind"] == "tsv":
            gene = make_genes_from_row(row["GENEINFO"])[1][0]
            notation = row["CLNHGVS"]
            disease = row["CLNDN"].split("|")[0].replace("_", " ")
            stars = review_star_map.get(row["CLNREVSTAT"], 0)
            variant_id = row["ID"]
            note_html = ""
        else:
            gene = entry.get("gene_label", "")
            notation = row["name"]
            stars = review_star_map.get(row["review_status"].replace(" ", "_"), 0)
            variant_id = row["variation_id"]
            if row["phenotype"] in UNRESOLVED_PHENOTYPES:
                disease = (
                    f"<span style='color:#94a3b8;'>phenotype not resolved in this record "
                    f"(ClinVar placeholder &ldquo;{row['phenotype'] or 'not provided'}&rdquo;)</span>"
                )
            else:
                disease = row["phenotype"]
            note_html = (
                f"<div style='color:#64748b; font-size:11px; margin-top:2px;'>"
                f"~{row['span']:,} bp &mdash; not in clinvar.vcf.gz, see section 4</div>"
            )

        star_html = "&#9733;" * stars + "&#9734;" * (4 - stars)
        cartoon_html = entry["cartoon"]()
        return f"""<div class="xf-panel">
  <div style="font-size:11px; font-weight:700; color:#2563eb; text-transform:uppercase; letter-spacing:0.03em;">{entry['title']}</div>
  <div style="margin-top:8px;">{cartoon_html}</div>
  <div style="font-size:14px; font-weight:600; margin-top:8px;">{gene}</div>
  <div style="font-family:ui-monospace,monospace; font-size:11px; color:#475569; margin-top:4px; word-break:break-all;">{notation}</div>
  <div style="font-size:13px; margin-top:6px;">{disease}</div>
  <div style="font-size:12px; color:#b45309; margin-top:4px;">{star_html}</div>
  {note_html}
  <div style="margin-top:8px;">
    <a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{variant_id}/" target="_blank" rel="noopener">ClinVar VariationID {variant_id} &#8599;</a>
  </div>
</div>"""

    example_cards_html = "".join(mutation_gallery_card(e) for e in MUTATION_GALLERY)

    def gene_model_section(result: dict) -> str:
        legend_html = "".join(
            f"<div style='display:flex; align-items:center; gap:6px; font-size:12px; color:#334155; padding:2px 0;'>"
            f"<span style='width:10px; height:10px; border-radius:2px; background:{row['color']}; flex:none;'></span>"
            f"<span>{row['name']}</span><span style='color:#94a3b8;'>({row['n']})</span></div>"
            for row in result["legend"]
        )
        blurb = result["blurb"].format(count=result["count"])
        return f"""<h3 style="margin-top:2rem;">{result['heading']}</h3>
<p class="subtitle">{blurb}</p>
<div class="table-wrap" style="overflow-x:auto; padding:12px;">
{result['svg']}
</div>
<div style="display:flex; flex-wrap:wrap; gap:4px 24px; margin-top:8px;">
{legend_html}
</div>"""

    gene_model_sections_html = "".join(gene_model_section(r) for r in gene_model_results)

    pairs_json = json.dumps(pairs)
    rescue_json = json.dumps(rescue_rows)
    cube_json = json.dumps(cube)
    clnsig_buckets_json = json.dumps(CLNSIG_BUCKETS)
    clnvc_types_json = json.dumps(clnvc_types)
    size_buckets_json = json.dumps(SIZE_BUCKETS)
    pair_list_json = json.dumps(pair_list)
    sv_rows_json = json.dumps(sv_summary["sv_rows"])
    total_variants = sum(r["count"] for r in cube)

    rescue_pct = (
        100 * rescue_summary["sub_expert_pairs_rescued"] / rescue_summary["sub_expert_pairs_total"]
        if rescue_summary["sub_expert_pairs_total"]
        else 0
    )

    vcf_only_types = [
        "single nucleotide variant",
        "Deletion",
        "Duplication",
        "Microsatellite",
        "Indel",
        "Insertion",
        "Inversion",
    ]
    sv_type_table_rows = "".join(
        f"<tr>"
        f"<td>{r['type']}</td>"
        f"<td class='num'>{r['count']:,}</td>"
        f"<td class='num'>{r['resolved']:,} ({r['resolved_pct']:.1f}%)</td>"
        f"<td class='num'>{r['span_min']:,}</td>"
        f"<td class='num'>{r['span_median']:,}</td>"
        f"<td class='num'>{r['span_max']:,}</td>"
        f"</tr>"
        if r["count"]
        else f"<tr><td>{r['type']}</td><td class='num'>0</td><td class='num'>&mdash;</td>"
        f"<td class='num'>&mdash;</td><td class='num'>&mdash;</td><td class='num'>&mdash;</td></tr>"
        for r in sv_summary["sv_type_rows"]
    )
    vcf_type_table_rows = "".join(
        f"<tr><td>{t}</td><td class='num'>{sv_summary['type_counts'].get(t, 0):,}</td></tr>" for t in vcf_only_types
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>ClinVar exploration report</title>
<style>
  body {{ font-family: -apple-system, sans-serif; max-width: 1000px; margin: 40px auto; color: #0f172a; }}
  body {{ padding: 0 16px; }}
  h1 {{ font-size: 1.4rem; }}
  h2 {{ font-size: 1.1rem; margin-top: 2.5rem; }}
  h3 {{ font-size: 1rem; margin: 0 0 8px; color: #334155; }}
  p.subtitle {{ color: #64748b; }}
  table {{ border-collapse: collapse; width: 100%; margin-top: 1rem; font-size: 13px; }}
  th, td {{ text-align: left; padding: 6px 12px; border-bottom: 1px solid #e2e8f0; }}
  tr.prod td {{ font-weight: 600; background: #fffbeb; }}
  .chart {{ margin-top: 0.5rem; }}
  .controls {{ display: flex; gap: 8px; align-items: center; margin-top: 1rem; flex-wrap: wrap; }}
  .search-box {{ flex: 1; min-width: 200px; padding: 6px 10px; border-radius: 6px; }}
  .search-box {{ border: 1px solid #d1d5db; font-size: 13px; }}
  .count-label {{ font-size: 12px; color: #64748b; white-space: nowrap; }}
  .table-wrap {{ max-height: 480px; overflow: auto; border: 1px solid #e2e8f0; border-radius: 8px; margin-top: 8px; }}
  .table-wrap table {{ margin-top: 0; }}
  .table-wrap thead th {{ position: sticky; top: 0; background: #f8fafc; cursor: pointer; white-space: nowrap; }}
  .table-wrap td.mono {{ font-family: ui-monospace, monospace; font-size: 12px; color: #475569; }}
  .table-wrap td.num {{ text-align: right; }}
  .yes {{ color: #b45309; font-weight: 600; }}
  .no {{ color: #cbd5e1; }}
  .summary-box {{ background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 8px; margin-top: 1rem; }}
  .summary-box {{ padding: 12px 16px; }}
  .summary-box strong {{ color: #0f172a; }}
  .lit-badge {{ display: inline-block; padding: 1px 6px; border-radius: 4px; font-size: 11px; font-weight: 600; }}
  .lit-yes {{ background: #dbeafe; color: #1e40af; }}
  .lit-no {{ background: #f1f5f9; color: #94a3b8; }}
  .total-boxes {{ display: flex; gap: 16px; margin: 1.5rem 0; flex-wrap: wrap; }}
  .total-box {{ background: #0f172a; color: #fff; border-radius: 10px; padding: 20px 24px; flex: 1; min-width: 220px; }}
  .total-box .n {{ font-size: 2.2rem; font-weight: 700; }}
  .total-box .label {{ font-size: 13px; color: #94a3b8; }}
  .xf-panels {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 16px; margin-top: 1rem; }}
  @media (max-width: 800px) {{ .xf-panels {{ grid-template-columns: 1fr; }} }}
  .xf-panel {{ border: 1px solid #e2e8f0; border-radius: 10px; padding: 14px; }}
  .xf-row {{ display: flex; align-items: center; gap: 6px; padding: 3px 0; cursor: pointer; user-select: none; }}
  .xf-row input {{ cursor: pointer; }}
  .xf-row .swatch {{ width: 10px; height: 10px; border-radius: 2px; flex: none; }}
  .xf-bar-track {{ flex: 1; height: 16px; background: #f1f5f9; border-radius: 4px; overflow: hidden; }}
  .xf-bar-fill {{ height: 100%; background: #2563eb; }}
  .xf-bar-fill.off {{ background: #cbd5e1; }}
  .xf-row .count {{ font-size: 11px; color: #64748b; width: 64px; text-align: right; flex: none; }}
  .xf-actions {{ display: flex; gap: 8px; margin-top: 8px; }}
  .xf-actions button {{ font-size: 11px; padding: 3px 8px; border-radius: 5px; cursor: pointer; }}
  .xf-actions button {{ border: 1px solid #cbd5e1; background: #fff; }}
  .pagination {{ display: flex; align-items: center; gap: 12px; margin-top: 8px; }}
  .pagination button {{ font-size: 12px; padding: 4px 10px; border-radius: 5px; cursor: pointer; }}
  .pagination button {{ border: 1px solid #cbd5e1; background: #fff; }}
  .pagination button:disabled {{ opacity: 0.4; cursor: not-allowed; }}
  .pagination .page-info {{ font-size: 12px; color: #64748b; }}
  .btn-primary {{ font-size: 13px; padding: 7px 14px; border-radius: 6px; cursor: pointer; }}
  .btn-primary {{ border: 1px solid #1d4ed8; background: #2563eb; color: #fff; font-weight: 600; }}
  .btn-primary:hover {{ background: #1d4ed8; }}
</style>
</head>
<body>
<h1>ClinVar exploration report</h1>
<div class="summary-box">
  <strong><a href="https://www.ncbi.nlm.nih.gov/clinvar/" target="_blank" rel="noopener">ClinVar</a></strong>
  is a free, public archive run by NCBI (part of NIH) that aggregates reported relationships between human
  genetic variants and health conditions. Clinical testing laboratories, research groups, and expert panels
  submit their own interpretation of a variant &mdash; e.g. <em>Pathogenic</em>, <em>Likely pathogenic</em>,
  <em>Uncertain significance</em>, <em>Likely benign</em>, or <em>Benign</em> &mdash; along with the disease
  it's associated with and the evidence behind that call. Because submitters vary widely in rigor, every
  submission carries a <strong>review status</strong> ("star rating", 0&ndash;4) reflecting how it was
  assessed: from a single submitter with no assessed criteria (0&#9733;) up to a formal expert panel or
  practice guideline (4&#9733;). Multiple independent submitters can (and often do) submit interpretations
  for the same variant, sometimes agreeing and sometimes conflicting. This ingest, and the analyses below,
  are about deciding which of that submitted evidence is trustworthy enough to promote into the knowledge
  graph as a variant-disease association.
</div>

<h2>Illustrative examples</h2>
<p class="subtitle">
  ClinVar covers several fundamentally different kinds of variation, from a single altered base up to
  megabase-scale rearrangements. Each card below is a simplified schematic (not to scale, and not a real
  genomic-coordinate plot -- see the gene models further down for that) of the underlying mechanism, paired
  with a real example pulled from the actual downloaded release and linked to its ClinVar record. The two
  "large recurrent" CNVs both come from the same 22q11.21 region (the classic DiGeorge/velocardiofacial
  locus) to show a deletion and its reciprocal duplication side by side; the two "single-gene" CNVs both
  come from DMD, whose exon-level deletions and duplications are the textbook example of this category.
  Structural types (CNVs, inversion) aren't in <code>clinvar.vcf.gz</code> at all -- see section 4.
</p>
<div class="xf-panels" style="grid-template-columns: repeat(3, 1fr);">
{example_cards_html}
</div>

{gene_model_sections_html}

<p class="subtitle">
  Two analyses over the downloaded ClinVar release: (1&ndash;2) review-star cutoff impact on the
  production disease-mapping filter, and (3) a crossfilter over ClinVar's own variant-level
  classification and review-status fields.
</p>

<h2>1. Review-star cutoff impact on variant &amp; gene-disease-pair counts</h2>
<p class="subtitle">
  Effect of <code>var2disease_star_min</code> on variant &amp; gene-disease-pair counts, using
  per-submission evidence from <code>submission_summary.txt</code>. Production currently uses
  <strong>{PRODUCTION_STAR_MIN} stars</strong> (highlighted below).
</p>

<h3>Distinct variants with &ge;1 disease mapping</h3>
<div class="chart">{bars("variants", max_variants)}</div>

<h3>Distinct gene-disease pairs</h3>
<div class="chart">{bars("gene_disease_pairs", max_pairs)}</div>

<table>
<tr><th>star_min</th><th>Variants</th><th>Gene-disease pairs</th></tr>
{table_rows}
</table>

<h3>Gene &rarr; MONDO disease pairs (star_min &ge;2 / &ge;3, and the &ge;4 subset)</h3>
<p class="subtitle">
  star_min=2 and star_min=3 always produce identical pairs here (see script docstring), so they're
  shown as one column. The &ge;4 column highlights which of those pairs also survive the strictest cutoff.
  ClinVar's <code>GENEINFO</code> field lists every locus overlapping a variant's position, which would
  otherwise surface NCBI's placeholder "LOC" gene symbols (<code>LOC</code> + Entrez GeneID, e.g.
  <code>LOC109461479</code>, assigned when a locus has no curated official gene symbol) as extra rows
  alongside the real gene for the same variant/disease -- these are dropped everywhere in this report
  (see <code>drop_loc_genes()</code>).
</p>
<div class="controls">
  <input id="pairs-search" class="search-box" type="text" placeholder="Filter by gene or MONDO id...">
  <label style="font-size:13px; color:#475569; display:flex; align-items:center; gap:4px; white-space:nowrap;">
    <input id="pairs-only4" type="checkbox"> only &ge;4&#9733; survivors
  </label>
  <span id="pairs-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th data-k="gene">Gene &#8597;</th>
  <th data-k="gene_id">Gene ID &#8597;</th>
  <th data-k="mondo">MONDO disease &#8597;</th>
  <th data-k="v2">Variants (&ge;2&#9733;/&ge;3&#9733;) &#8597;</th>
  <th data-k="v4">Variants (&ge;4&#9733;) &#8597;</th>
</tr>
</thead>
<tbody id="pairs-rows"></tbody>
</table>
</div>

<h2>2. 0/1&#9733;-only pairs: multi-submitter concordance rescue</h2>
<p class="subtitle">
  Gene-disease pairs whose <strong>only</strong> supporting evidence is star 0/1
  (never reach star_min={EXPERT_STAR_MIN}) are additionally checked for
  &ge;{MIN_CONCORDANT_SUBMITTERS} distinct Submitters independently reporting the
  <strong>same MONDO disease with the exact same ClinicalSignificance</strong> (Pathogenic/Likely pathogenic family
  only), regardless of any individual submitter's own review status. Each pair is also annotated with whether any
  contributing record used <code>CollectionMethod == "literature only"</code> (a publication proxy) &mdash;
  informational only, not used as a filter here.
</p>
<div class="summary-box">
  <strong>{rescue_summary["sub_expert_pairs_rescued"]:,}</strong> of
  <strong>{rescue_summary["sub_expert_pairs_total"]:,}</strong> 0/1&#9733;-only gene-disease pairs
  ({rescue_pct:.1f}%) have &ge;{MIN_CONCORDANT_SUBMITTERS} concordant submitters.<br>
  <strong>{rescue_summary["sub_expert_variants_rescued"]:,}</strong> of
  <strong>{rescue_summary["sub_expert_variants_total"]:,}</strong> variants with
  no expert-tier (&ge;{EXPERT_STAR_MIN}&#9733;) evidence at all gain at least one rescued pair.
</div>
<div class="controls">
  <input id="rescue-search" class="search-box" type="text" placeholder="Filter by gene, MONDO id, or disease name...">
  <label style="font-size:13px; color:#475569; display:flex; align-items:center; gap:4px; white-space:nowrap;">
    <input id="rescue-only-rescued" type="checkbox" checked> only rescued
  </label>
  <span id="rescue-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th data-k="gene">Gene &#8597;</th>
  <th data-k="gene_id">Gene ID &#8597;</th>
  <th data-k="mondo">MONDO disease &#8597;</th>
  <th data-k="disease_name">Disease name &#8597;</th>
  <th data-k="n_variants"># variants &#8597;</th>
  <th data-k="submitter_count">Concordant submitters &#8597;</th>
  <th data-k="has_literature">Has literature submission &#8597;</th>
</tr>
</thead>
<tbody id="rescue-rows"></tbody>
</table>
</div>

<h2>3. Clinical significance &amp; review-status crossfilter</h2>
<p class="subtitle">
  Built from all {total_variants:,} variants in <code>clinvar.vcf.gz</code>, using ClinVar's own
  aggregate, variant-level fields (CLNSIG / CLNREVSTAT / CLNVC) &mdash; not the per-submission
  evidence used in sections 1&ndash;2. <code>CLNREVSTAT</code> here is ClinVar's real aggregate
  review status, so it genuinely includes the "2 star" tier that's absent per-record in section 1.
  Toggle checkboxes in any panel; the other panels re-filter to match, and the totals update live.
  Each panel's own bars stay visible when toggled off (dimmed) so you can see what you're excluding.
  Gene-disease pairs are only ever produced from Pathogenic/Likely-pathogenic evidence (same as
  sections 1&ndash;2) &mdash; filtering to VUS/Benign/Conflicting/Other/Not-classified alone will
  always show 0 pairs, which is expected: production never creates disease edges for those.
  "Has literature submission" and "Multiple concordant submitters" are per-variant flags computed
  from <code>submission_summary.txt</code> (same source as sections 1&ndash;2, not the VCF fields
  above) &mdash; variants absent from that file show as "No" for both. "Overlaps STRchive locus"
  checks the variant's genomic footprint (POS through POS+len(REF)-1) against the hg38 coordinates
  of <a href="https://strchive.org/loci/" target="_blank" rel="noopener">STRchive</a>'s ~80 curated,
  disease-associated short-tandem-repeat loci &mdash; informational only, independent of ClinVar's
  own CLNVC "Microsatellite" label (see below). "Size" buckets by
  <code>max(len(REF), len(ALT))</code> &mdash; the actual physical allele footprint, independent of
  ClinVar's own Type/CLNVC label. Those can disagree: a 2bp deletion inside a tandem repeat gets
  labeled "Microsatellite" rather than "Deletion" by ClinVar's own convention (see e.g.
  <a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/1323792/" target="_blank"
  rel="noopener">VariationID 1323792</a>), so filtering by CLNVC alone can hide or include variants
  you wouldn't expect by size &mdash; this panel lets you filter on actual size regardless of label.
</p>

<div class="total-boxes">
  <div class="total-box">
    <div class="n" id="xf-total-count">-</div>
    <div class="label">variants matching current filters</div>
  </div>
  <div class="total-box">
    <div class="n" id="xf-pairs-count">-</div>
    <div class="label">gene-disease pairs matching current filters</div>
  </div>
</div>

<div class="xf-panels">
  <div class="xf-panel">
    <h3>Clinical significance</h3>
    <div id="xf-clnsig-rows"></div>
    <div class="xf-actions">
      <button data-dim="clnsig" data-action="all">All</button>
      <button data-dim="clnsig" data-action="none">None</button>
    </div>
  </div>
  <div class="xf-panel">
    <h3>Review status (stars)</h3>
    <div id="xf-star-rows"></div>
    <div class="xf-actions">
      <button data-dim="star" data-action="all">All</button>
      <button data-dim="star" data-action="none">None</button>
    </div>
  </div>
  <div class="xf-panel">
    <h3>Variant type</h3>
    <div id="xf-clnvc-rows"></div>
    <div class="xf-actions">
      <button data-dim="clnvc" data-action="all">All</button>
      <button data-dim="clnvc" data-action="none">None</button>
    </div>
  </div>
  <div class="xf-panel">
    <h3>Size (max of REF/ALT length)</h3>
    <div id="xf-size-rows"></div>
    <div class="xf-actions">
      <button data-dim="size" data-action="all">All</button>
      <button data-dim="size" data-action="none">None</button>
    </div>
  </div>
  <div class="xf-panel">
    <h3>Has literature submission</h3>
    <div id="xf-literature-rows"></div>
    <div class="xf-actions">
      <button data-dim="literature" data-action="all">All</button>
      <button data-dim="literature" data-action="none">None</button>
    </div>
  </div>
  <div class="xf-panel">
    <h3>Multiple concordant submitters (&ge;{MIN_CONCORDANT_SUBMITTERS})</h3>
    <div id="xf-concordant-rows"></div>
    <div class="xf-actions">
      <button data-dim="concordant" data-action="all">All</button>
      <button data-dim="concordant" data-action="none">None</button>
    </div>
  </div>
  <div class="xf-panel">
    <h3>Overlaps STRchive locus</h3>
    <div id="xf-strchive-rows"></div>
    <div class="xf-actions">
      <button data-dim="strchive" data-action="all">All</button>
      <button data-dim="strchive" data-action="none">None</button>
    </div>
  </div>
</div>

<h3>Gene-disease pairs matching current filters</h3>
<p class="subtitle">
  Adjusting checkboxes above updates the live totals immediately, but the pair list below only
  recomputes when you click <strong>Apply filters</strong> &mdash; listing distinct (gene, disease)
  pairs can be a larger computation than the aggregate counts, so it's on-demand rather than live.
  Each pair also lists a sample of the contributing ClinVar variant IDs (hyperlinked to the actual
  ClinVar record), capped at 5 per row since some pairs -- BRCA1/BRCA2 top out around 3-5k
  contributing variants -- have far too many to list in full; "# variants" is always the exact count.
  NCBI's placeholder "LOC" gene symbols (<code>LOC</code> + Entrez GeneID, e.g. <code>LOC109461479</code>,
  assigned when a locus has no curated official gene symbol -- ClinVar's <code>GENEINFO</code> would
  otherwise surface these as extra rows alongside the real gene for the same variant/disease, e.g. the
  HTT repeat-expansion variant also tagging LOC109461479 and LOC129929027) are dropped everywhere in
  this report (see <code>drop_loc_genes()</code>).
</p>
<div class="controls">
  <button id="xf-apply" class="btn-primary">Apply filters</button>
  <input id="xf-pairs-table-search" class="search-box" type="text"
    placeholder="Filter by gene, MONDO id, or disease name...">
  <span id="xf-pairs-table-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th data-sort-for="xf-pairs-table-rows" data-sort-key="gene">Gene &#8597;</th>
  <th data-sort-for="xf-pairs-table-rows" data-sort-key="gene_id">Gene ID &#8597;</th>
  <th data-sort-for="xf-pairs-table-rows" data-sort-key="mondo">MONDO disease &#8597;</th>
  <th data-sort-for="xf-pairs-table-rows" data-sort-key="disease_name">Disease name &#8597;</th>
  <th data-sort-for="xf-pairs-table-rows" data-sort-key="variant_count"># variants &#8597;</th>
  <th>ClinVar IDs</th>
</tr>
</thead>
<tbody id="xf-pairs-table-rows"></tbody>
</table>
</div>
<div class="pagination">
  <button id="xf-pairs-table-prev">&larr; Prev</button>
  <span id="xf-pairs-table-page-info" class="page-info"></span>
  <button id="xf-pairs-table-next">Next &rarr;</button>
</div>

<h2>4. Structural variants &amp; CNVs &mdash; what's not in the VCF</h2>
<p class="subtitle">
  Sections 1&ndash;3 above are built entirely from <code>clinvar.vcf.gz</code> and
  <code>submission_summary.txt</code> &mdash; both of which require a fixed genomic
  position plus REF/ALT sequence. That structurally excludes copy-number gain/loss,
  translocations, and other large rearrangements. ClinVar publishes those separately in
  <code>variant_summary.txt.gz</code> (Start/Stop coordinate ranges instead of a fixed
  REF/ALT), from the same <code>tab_delimited/</code> directory as
  <code>submission_summary.txt</code>. This section downloads that file once (~440MB)
  and reports real counts from the GRCh38 rows.
</p>

<h3>Types absent from the VCF's CLNVC field</h3>
<div class="table-wrap">
<table>
<thead>
<tr><th>Type</th><th># variants</th><th># disease-resolved (has a MONDO/MedGen/OMIM id)</th>
<th>Span min (bp)</th><th>Span median (bp)</th><th>Span max (bp)</th></tr>
</thead>
<tbody>{sv_type_table_rows}</tbody>
</table>
</div>
<p class="subtitle">
  Only <strong>{sv_summary["sv_type_rows"][0]["resolved_pct"]:.1f}%</strong> of copy-number-gain
  and <strong>{sv_summary["sv_type_rows"][1]["resolved_pct"]:.1f}%</strong> of copy-number-loss
  entries resolve to any disease id at all &mdash; the rest are ClinGen-style regional
  dosage-sensitivity assertions (ClinVar reports their <code>PhenotypeList</code> as the literal
  placeholder <code>"See cases"</code>, with no MONDO/MedGen/OMIM id), not tied to one named
  condition. Spans reach into the tens of millions of bp &mdash; three to four orders of
  magnitude past the ~10kb ceiling the VCF's fixed-REF/ALT format can represent at all.
</p>

<h3>For comparison: types the VCF already carries (same file, same assembly)</h3>
<div class="table-wrap">
<table>
<thead><tr><th>Type</th><th># variants</th></tr></thead>
<tbody>{vcf_type_table_rows}</tbody>
</table>
</div>

<h3>All SV/CNV variants ({len(sv_summary["sv_rows"]):,} rows)</h3>
<p class="subtitle">
  Every GRCh38 row across copy number gain/loss, Translocation, and Complex, any
  ClinicalSignificance. Only ~6-10% of copy-number gain/loss rows resolve to a disease id at all
  (see the table above) -- uncheck "SV resolved to disease" to see the unresolved majority, whose
  <code>PhenotypeList</code> is the literal ClinVar placeholder <code>"See cases"</code>. Note the
  gene list this ingest already builds (<code>make_genes_from_row()</code>, splitting
  <code>GENEINFO</code> on <code>|</code>) has no equivalent here &mdash;
  <code>variant_summary.txt</code>'s own <code>GeneSymbol</code> column degrades to unparseable
  prose for large CNVs (e.g. <em>"covers 42 genes, none of which curated to show dosage
  sensitivity"</em>), so pulling an actual gene list for these would require a separate
  genomic-coordinate lookup, not just parsing a delimited column.
</p>
<div class="controls">
  <input id="sv-table-search" class="search-box" type="text"
    placeholder="Filter by region, disease, or disease id...">
  <label style="font-size:13px; color:#475569; display:flex; align-items:center; gap:4px; white-space:nowrap;">
    <input id="sv-resolved-only" type="checkbox" checked> SV resolved to disease
  </label>
  <span id="sv-table-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th data-sort-for="sv-table-rows" data-sort-key="variation_id">ClinVar ID &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="allele_id">AlleleID &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="name">Name / nomenclature &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="type">Type &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="clinsig">Clinical significance &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="resolved">Resolved &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="phenotype">Disease &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="phenotype_ids">Disease ID(s) &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="span">Span (bp) &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="num_submitters"># submitters &#8597;</th>
  <th data-sort-for="sv-table-rows" data-sort-key="review_status">Review status &#8597;</th>
</tr>
</thead>
<tbody id="sv-table-rows"></tbody>
</table>
</div>
<div class="pagination">
  <button id="sv-table-prev">&larr; Prev</button>
  <span id="sv-table-page-info" class="page-info"></span>
  <button id="sv-table-next">Next &rarr;</button>
</div>

<script>
function createPaginatedTable(config) {{
  // config: {{ tbodyId, countElId, searchElId, prevBtnId, nextBtnId, pageInfoElId,
  //   pageSize, columns: [{{key, label, className, format(v, row)}}], searchFields: [...],
  //   defaultSortKey, defaultSortDir, checkboxFilterElId, checkboxFilterPredicate(row) }}
  var data = [];
  var query = "";
  var sortKey = config.defaultSortKey;
  var sortDir = config.defaultSortDir || 1;
  var page = 0;
  var pageSize = config.pageSize || 25;

  var tbody = document.getElementById(config.tbodyId);
  var countEl = config.countElId ? document.getElementById(config.countElId) : null;
  var searchEl = config.searchElId ? document.getElementById(config.searchElId) : null;
  var pageInfoEl = document.getElementById(config.pageInfoElId);
  var prevBtn = document.getElementById(config.prevBtnId);
  var nextBtn = document.getElementById(config.nextBtnId);
  var checkboxEl = config.checkboxFilterElId ? document.getElementById(config.checkboxFilterElId) : null;

  function filteredSorted() {{
    var rows = data;
    if (checkboxEl && checkboxEl.checked && config.checkboxFilterPredicate) {{
      rows = rows.filter(config.checkboxFilterPredicate);
    }}
    if (query) {{
      var q = query.toLowerCase();
      rows = rows.filter(function(r) {{
        return config.searchFields.some(function(f) {{
          return String(r[f] == null ? "" : r[f]).toLowerCase().indexOf(q) !== -1;
        }});
      }});
    }}
    rows = rows.slice().sort(function(a, b) {{
      var av = a[sortKey], bv = b[sortKey];
      if (typeof av === "string") return av.localeCompare(bv) * sortDir;
      return ((av || 0) - (bv || 0)) * sortDir;
    }});
    return rows;
  }}

  function render() {{
    var rows = filteredSorted();
    var totalPages = Math.max(1, Math.ceil(rows.length / pageSize));
    if (page >= totalPages) page = totalPages - 1;
    if (page < 0) page = 0;
    var pageRows = rows.slice(page * pageSize, (page + 1) * pageSize);

    tbody.innerHTML = pageRows.map(function(r) {{
      return "<tr>" + config.columns.map(function(c) {{
        var v = c.format ? c.format(r[c.key], r) : r[c.key];
        return "<td class='" + (c.className || "") + "'>" + v + "</td>";
      }}).join("") + "</tr>";
    }}).join("");

    if (countEl) countEl.textContent = rows.length.toLocaleString() + " rows";
    if (pageInfoEl) pageInfoEl.textContent = rows.length ? ("Page " + (page + 1) + " of " + totalPages) : "";
    if (prevBtn) prevBtn.disabled = page <= 0;
    if (nextBtn) nextBtn.disabled = page >= totalPages - 1;
  }}

  if (searchEl) {{
    searchEl.addEventListener("input", function(e) {{
      query = e.target.value;
      page = 0;
      render();
    }});
  }}
  if (checkboxEl) {{
    checkboxEl.addEventListener("change", function() {{
      page = 0;
      render();
    }});
  }}
  if (prevBtn) prevBtn.addEventListener("click", function() {{ page--; render(); }});
  if (nextBtn) nextBtn.addEventListener("click", function() {{ page++; render(); }});

  document.querySelectorAll("[data-sort-for='" + config.tbodyId + "']").forEach(function(th) {{
    th.addEventListener("click", function() {{
      var k = th.getAttribute("data-sort-key");
      if (sortKey === k) {{
        sortDir *= -1;
      }} else {{
        sortKey = k;
        sortDir = 1;
      }}
      page = 0;
      render();
    }});
  }});

  return {{
    setData: function(newData) {{
      data = newData;
      page = 0;
      render();
    }},
  }};
}}

(function() {{
  var data = {pairs_json};
  var state = {{ sortKey: "v2", sortDir: -1, query: "" }};
  var tbody = document.getElementById("pairs-rows");
  var countEl = document.getElementById("pairs-count");
  var only4El = document.getElementById("pairs-only4");

  function render() {{
    var rows = data.filter(function(r) {{
      if (only4El.checked && r.v4 === null) return false;
      if (!state.query) return true;
      var q = state.query.toLowerCase();
      return r.gene.toLowerCase().indexOf(q) !== -1 ||
             r.gene_id.toLowerCase().indexOf(q) !== -1 ||
             r.mondo.toLowerCase().indexOf(q) !== -1;
    }});

    rows.sort(function(a, b) {{
      var av = a[state.sortKey], bv = b[state.sortKey];
      if (av === null) av = -1;
      if (bv === null) bv = -1;
      if (typeof av === "string") return av.localeCompare(bv) * state.sortDir;
      return (av - bv) * state.sortDir;
    }});

    countEl.textContent = rows.length + " of " + data.length + " pairs";

    tbody.innerHTML = rows.map(function(r) {{
      var v4cell = r.v4 === null
        ? "<span class='no'>&mdash;</span>"
        : "<span class='yes'>" + r.v4 + "</span>";
      return "<tr>" +
        "<td>" + r.gene + "</td>" +
        "<td class='mono'>" + r.gene_id + "</td>" +
        "<td class='mono'>" + r.mondo + "</td>" +
        "<td class='num'>" + r.v2 + "</td>" +
        "<td class='num'>" + v4cell + "</td>" +
      "</tr>";
    }}).join("");
  }}

  document.getElementById("pairs-search").addEventListener("input", function(e) {{
    state.query = e.target.value;
    render();
  }});
  only4El.addEventListener("change", render);
  document.querySelectorAll(".table-wrap th[data-k]").forEach(function(th) {{
    if (th.closest("table").querySelector("#pairs-rows") === null) return;
    th.addEventListener("click", function() {{
      var k = th.getAttribute("data-k");
      if (state.sortKey === k) {{
        state.sortDir *= -1;
      }} else {{
        state.sortKey = k;
        state.sortDir = (k === "gene" || k === "gene_id" || k === "mondo") ? 1 : -1;
      }}
      render();
    }});
  }});

  render();
}})();

(function() {{
  var data = {rescue_json};
  var state = {{ sortKey: "submitter_count", sortDir: -1, query: "" }};
  var tbody = document.getElementById("rescue-rows");
  var countEl = document.getElementById("rescue-count");
  var onlyRescuedEl = document.getElementById("rescue-only-rescued");

  function render() {{
    var rows = data.filter(function(r) {{
      if (onlyRescuedEl.checked && !r.rescued) return false;
      if (!state.query) return true;
      var q = state.query.toLowerCase();
      return r.gene.toLowerCase().indexOf(q) !== -1 ||
             r.gene_id.toLowerCase().indexOf(q) !== -1 ||
             r.mondo.toLowerCase().indexOf(q) !== -1 ||
             r.disease_name.toLowerCase().indexOf(q) !== -1;
    }});

    rows.sort(function(a, b) {{
      var av = a[state.sortKey], bv = b[state.sortKey];
      if (typeof av === "boolean") {{ av = av ? 1 : 0; bv = bv ? 1 : 0; }}
      if (typeof av === "string") return av.localeCompare(bv) * state.sortDir;
      return (av - bv) * state.sortDir;
    }});

    countEl.textContent = rows.length + " of " + data.length + " pairs";

    tbody.innerHTML = rows.map(function(r) {{
      var litBadge = r.has_literature
        ? "<span class='lit-badge lit-yes'>literature</span>"
        : "<span class='lit-badge lit-no'>&mdash;</span>";
      var submitterCell = r.rescued
        ? "<span class='yes'>" + r.submitter_count + "</span>"
        : "<span class='no'>" + r.submitter_count + "</span>";
      return "<tr>" +
        "<td>" + r.gene + "</td>" +
        "<td class='mono'>" + r.gene_id + "</td>" +
        "<td class='mono'>" + r.mondo + "</td>" +
        "<td>" + (r.disease_name || "<span class='no'>&mdash;</span>") + "</td>" +
        "<td class='num'>" + r.n_variants + "</td>" +
        "<td class='num'>" + submitterCell + "</td>" +
        "<td class='num'>" + litBadge + "</td>" +
      "</tr>";
    }}).join("");
  }}

  document.getElementById("rescue-search").addEventListener("input", function(e) {{
    state.query = e.target.value;
    render();
  }});
  onlyRescuedEl.addEventListener("change", render);
  document.querySelectorAll(".table-wrap th[data-k]").forEach(function(th) {{
    if (th.closest("table").querySelector("#rescue-rows") === null) return;
    th.addEventListener("click", function() {{
      var k = th.getAttribute("data-k");
      if (state.sortKey === k) {{
        state.sortDir *= -1;
      }} else {{
        state.sortKey = k;
        state.sortDir = (k === "gene" || k === "gene_id" || k === "mondo") ? 1 : -1;
      }}
      render();
    }});
  }});

  render();
}})();

(function() {{
  var cube = {cube_json};
  var pairList = {pair_list_json};
  var clnsigBuckets = {clnsig_buckets_json};
  var clnvcTypes = {clnvc_types_json};
  var sizeBuckets = {size_buckets_json};
  var starLevels = [0, 1, 2, 3, 4];

  var clnsigColors = {{
    "P": "#b91c1c", "LP": "#ea580c", "P/LP": "#c2410c",
    "VUS": "#a16207", "LB": "#16a34a", "B": "#15803d", "B/LB": "#166534",
    "Conflicting": "#7c3aed", "Other": "#64748b", "Not classified": "#94a3b8"
  }};

  var boolValues = [true, false];
  var ALL_DIMS = ["clnsig", "star", "clnvc", "size", "literature", "concordant", "strchive"];

  var state = {{
    clnsig: new Set(clnsigBuckets),
    star: new Set(starLevels),
    clnvc: new Set(clnvcTypes),
    size: new Set(sizeBuckets),
    literature: new Set(boolValues),
    concordant: new Set(boolValues),
    strchive: new Set(boolValues)
  }};

  function otherDims(dim) {{
    return ALL_DIMS.filter(function(d) {{ return d !== dim; }});
  }}

  function matches(row, dims) {{
    for (var i = 0; i < dims.length; i++) {{
      var d = dims[i];
      if (!state[d].has(row[d])) return false;
    }}
    return true;
  }}

  function total() {{
    var sum = 0;
    cube.forEach(function(row) {{ if (matches(row, ALL_DIMS)) sum += row.count; }});
    return sum;
  }}

  function pairsTotal() {{
    var seen = new Set();
    cube.forEach(function(row) {{
      if (!matches(row, ALL_DIMS)) return;
      row.pairs.forEach(function(p) {{ seen.add(p.i); }});
    }});
    return seen.size;
  }}

  function breakdown(dim, otherDims) {{
    var totals = {{}};
    cube.forEach(function(row) {{
      if (!matches(row, otherDims)) return;
      totals[row[dim]] = (totals[row[dim]] || 0) + row.count;
    }});
    return totals;
  }}

  function renderPanel(dim, values, otherDims, containerId, labelFn, colorFn) {{
    var totals = breakdown(dim, otherDims);
    var maxVal = Math.max.apply(null, values.map(function(v) {{ return totals[v] || 0; }})) || 1;
    var container = document.getElementById(containerId);
    container.innerHTML = values.map(function(v) {{
      var n = totals[v] || 0;
      var on = state[dim].has(v);
      var pct = (n / maxVal) * 100;
      var color = colorFn ? colorFn(v) : "#2563eb";
      return "<div class='xf-row' data-dim='" + dim + "' data-value='" + v + "'>" +
        "<input type='checkbox'" + (on ? " checked" : "") + " data-dim='" + dim + "' data-value='" + v + "'>" +
        "<span class='swatch' style='background:" + color + "'></span>" +
        "<span style='width:170px; font-size:12px; color:#334155;'>" + labelFn(v) + "</span>" +
        "<span class='xf-bar-track'><span class='xf-bar-fill" + (on ? "" : " off") +
          "' style='width:" + pct.toFixed(1) + "%; background:" + (on ? color : "#cbd5e1") + ";'></span></span>" +
        "<span class='count'>" + n.toLocaleString() + "</span>" +
      "</div>";
    }}).join("");
  }}

  function renderAll() {{
    document.getElementById("xf-total-count").textContent = total().toLocaleString();
    document.getElementById("xf-pairs-count").textContent = pairsTotal().toLocaleString();
    renderPanel("clnsig", clnsigBuckets, otherDims("clnsig"), "xf-clnsig-rows",
      function(v) {{ return v; }}, function(v) {{ return clnsigColors[v] || "#2563eb"; }});
    renderPanel("star", starLevels, otherDims("star"), "xf-star-rows",
      function(v) {{ return v + (v === 1 ? " star" : " stars"); }}, function() {{ return "#2563eb"; }});
    renderPanel("clnvc", clnvcTypes, otherDims("clnvc"), "xf-clnvc-rows",
      function(v) {{ return v.replace(/_/g, " "); }}, function() {{ return "#2563eb"; }});
    renderPanel("size", sizeBuckets, otherDims("size"), "xf-size-rows",
      function(v) {{ return v; }}, function() {{ return "#2563eb"; }});
    renderPanel("literature", boolValues, otherDims("literature"), "xf-literature-rows",
      function(v) {{ return v ? "Yes" : "No"; }}, function() {{ return "#2563eb"; }});
    renderPanel("concordant", boolValues, otherDims("concordant"), "xf-concordant-rows",
      function(v) {{ return v ? "Yes" : "No"; }}, function() {{ return "#2563eb"; }});
    renderPanel("strchive", boolValues, otherDims("strchive"), "xf-strchive-rows",
      function(v) {{ return v ? "Yes" : "No"; }}, function() {{ return "#2563eb"; }});
  }}

  document.body.addEventListener("change", function(e) {{
    if (e.target.tagName !== "INPUT" || !e.target.hasAttribute("data-dim")) return;
    var dim = e.target.getAttribute("data-dim");
    if (!(dim in state)) return;
    var raw = e.target.getAttribute("data-value");
    var value;
    if (dim === "star") {{
      value = parseInt(raw, 10);
    }} else if (dim === "literature" || dim === "concordant" || dim === "strchive") {{
      value = raw === "true";
    }} else {{
      value = raw;
    }}
    if (e.target.checked) {{
      state[dim].add(value);
    }} else {{
      state[dim].delete(value);
    }}
    renderAll();
  }});

  document.body.addEventListener("click", function(e) {{
    var btn = e.target.closest(".xf-actions button[data-dim]");
    if (!btn) return;
    var dim = btn.getAttribute("data-dim");
    var action = btn.getAttribute("data-action");
    var values = dim === "clnsig" ? clnsigBuckets
      : dim === "star" ? starLevels
      : dim === "clnvc" ? clnvcTypes
      : dim === "size" ? sizeBuckets
      : boolValues;
    state[dim] = action === "all" ? new Set(values) : new Set();
    renderAll();
  }});

  var PAIRS_TABLE_VARIANT_DISPLAY_CAP = 5;

  function clinvarLink(id) {{
    var url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/" + id + "/";
    return "<a href='" + url + "' target='_blank' rel='noopener'>" + id + "</a>";
  }}

  var pairsTable = createPaginatedTable({{
    tbodyId: "xf-pairs-table-rows",
    countElId: "xf-pairs-table-count",
    searchElId: "xf-pairs-table-search",
    prevBtnId: "xf-pairs-table-prev",
    nextBtnId: "xf-pairs-table-next",
    pageInfoElId: "xf-pairs-table-page-info",
    pageSize: 25,
    defaultSortKey: "gene",
    defaultSortDir: 1,
    searchFields: ["gene", "gene_id", "mondo", "disease_name"],
    columns: [
      {{ key: "gene" }},
      {{ key: "gene_id", className: "mono" }},
      {{ key: "mondo", className: "mono" }},
      {{ key: "disease_name", format: function(v) {{ return v || "<span class='no'>&mdash;</span>"; }} }},
      {{ key: "variant_count", className: "num" }},
      {{ key: "variant_sample", format: function(v, row) {{
        var shown = v.slice(0, PAIRS_TABLE_VARIANT_DISPLAY_CAP);
        var html = shown.map(clinvarLink).join(", ");
        var remaining = row.variant_count - shown.length;
        if (remaining > 0) html += " <span class='count-label'>+" + remaining + " more</span>";
        return html;
      }} }},
    ],
  }});

  document.getElementById("xf-apply").addEventListener("click", function() {{
    var merged = {{}};
    cube.forEach(function(row) {{
      if (!matches(row, ALL_DIMS)) return;
      row.pairs.forEach(function(p) {{
        var m = merged[p.i];
        if (!m) {{ m = merged[p.i] = {{ n: 0, v: [] }}; }}
        m.n += p.n;
        m.v = m.v.concat(p.v);
      }});
    }});
    var rows = Object.keys(merged).map(function(k) {{
      var idx = parseInt(k, 10);
      var base = pairList[idx];
      var m = merged[idx];
      return {{
        gene: base.gene,
        gene_id: base.gene_id,
        mondo: base.mondo,
        disease_name: base.disease_name,
        variant_count: m.n,
        variant_sample: m.v,
      }};
    }});
    pairsTable.setData(rows);
  }});

  renderAll();
}})();

(function() {{
  var svTable = createPaginatedTable({{
    tbodyId: "sv-table-rows",
    countElId: "sv-table-count",
    searchElId: "sv-table-search",
    prevBtnId: "sv-table-prev",
    nextBtnId: "sv-table-next",
    pageInfoElId: "sv-table-page-info",
    pageSize: 25,
    defaultSortKey: "span",
    defaultSortDir: -1,
    searchFields: ["name", "phenotype", "phenotype_ids", "variation_id", "allele_id"],
    checkboxFilterElId: "sv-resolved-only",
    checkboxFilterPredicate: function(row) {{ return row.resolved; }},
    columns: [
      {{ key: "variation_id", className: "mono", format: function(v) {{
        var url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/" + v + "/";
        return "<a href='" + url + "' target='_blank' rel='noopener'>" + v + "</a>";
      }} }},
      {{ key: "allele_id", className: "mono" }},
      {{ key: "name" }},
      {{ key: "type" }},
      {{ key: "clinsig" }},
      {{ key: "resolved", format: function(v) {{
        return v ? "<span class='yes'>Yes</span>" : "<span class='no'>No</span>";
      }} }},
      {{ key: "phenotype" }},
      {{ key: "phenotype_ids", className: "mono" }},
      {{ key: "span", className: "num", format: function(v) {{ return v.toLocaleString(); }} }},
      {{ key: "num_submitters", className: "num" }},
      {{ key: "review_status" }},
    ],
  }});
  svTable.setData({sv_rows_json});
}})();
</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("../data"))
    parser.add_argument("--output", type=Path, default=Path("clinvar_report.html"))
    args = parser.parse_args()

    var_records, map_to_mondo = load_maps(args.data_dir)
    mondo_labels = load_mondo_labels(args.data_dir)
    clinvar_tsv = args.data_dir / "clinvar.tsv"
    results, pairs, rescue_summary, rescue_rows = compute_star_data(
        clinvar_tsv, var_records, map_to_mondo, mondo_labels
    )

    print(f"{'star_min':>8} | {'variants':>10} | {'gene_disease_pairs':>19}")
    for star in STAR_LEVELS:
        r = results[star]
        print(f"{star:>8} | {r['variants']:>10,} | {r['gene_disease_pairs']:>19,}")

    print(
        f"\n0/1-star-only pairs: {rescue_summary['sub_expert_pairs_total']:,} total, "
        f"{rescue_summary['sub_expert_pairs_rescued']:,} rescued by "
        f">={MIN_CONCORDANT_SUBMITTERS} concordant submitters"
    )
    print(
        f"0/1-star-only variants (no >={EXPERT_STAR_MIN}-star evidence at all): "
        f"{rescue_summary['sub_expert_variants_total']:,} total, "
        f"{rescue_summary['sub_expert_variants_rescued']:,} gain a rescued pair"
    )

    strchive_intervals = load_strchive_intervals(args.data_dir)
    cube, clnvc_types, pair_list = build_clnsig_cube(
        clinvar_tsv, var_records, map_to_mondo, mondo_labels, strchive_intervals
    )
    total = sum(r["count"] for r in cube)
    by_clnsig: Counter = Counter()
    for r in cube:
        by_clnsig[r["clnsig"]] += r["count"]
    print(f"\nCLNSIG distribution (all {total:,} variants):")
    for bucket in CLNSIG_BUCKETS:
        n = by_clnsig[bucket]
        print(f"{bucket:>15} | {n:>10,} | {100 * n / total:>5.1f}%")
    all_pairs = {p["i"] for r in cube for p in r["pairs"]}
    print(f"\nGene-disease pairs across all filters: {len(all_pairs):,}")
    n_strchive_loci = sum(len(v) for v in strchive_intervals.values())
    n_strchive_variants = sum(r["count"] for r in cube if r["strchive"])
    print(f"\nSTRchive loci loaded: {n_strchive_loci:,}; variants overlapping a locus: {n_strchive_variants:,}")

    sv_summary = build_sv_summary(args.data_dir)
    print(f"\nStructural variant types (variant_summary.txt, GRCh38, {sv_summary['total_variants']:,} rows):")
    for r in sv_summary["sv_type_rows"]:
        if r["count"]:
            print(f"{r['type']:>20} | {r['count']:>7,} | resolved={r['resolved']:>6,} ({r['resolved_pct']:.1f}%)")
        else:
            print(f"{r['type']:>20} | {r['count']:>7,}")
    resolved_count = sum(1 for r in sv_summary["sv_rows"] if r["resolved"])
    print(f"Total SV/CNV rows: {len(sv_summary['sv_rows']):,} ({resolved_count:,} disease-resolved)")

    example_variants = load_example_variants(clinvar_tsv)
    sv_rows_by_id = {r["variation_id"]: r for r in sv_summary["sv_rows"]}
    for key, variation_id in EXAMPLE_SV_IDS.items():
        example_variants[key] = sv_rows_by_id[variation_id]

    gene_model_results = []
    for example in GENE_MODEL_EXAMPLES:
        gene = example["gene"]
        gene_model = load_gene_model(args.data_dir, gene)
        gene_model_variants = build_gene_model_variants(
            clinvar_tsv, var_records, map_to_mondo, mondo_labels, gene, example.get("diseases")
        )
        gene_model_colors, gene_model_legend = assign_gene_model_colors(gene_model_variants)
        gene_model_svg = render_gene_model_svg(gene, gene_model, gene_model_variants, gene_model_colors)
        print(f"\n{gene} illustrative gene model: {len(gene_model_variants):,} SNVs, {len(gene_model_legend):,} disease colors")
        gene_model_results.append(
            {
                "heading": example["heading"],
                "blurb": example["blurb"],
                "svg": gene_model_svg,
                "legend": gene_model_legend,
                "count": len(gene_model_variants),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        render_html(
            results,
            pairs,
            rescue_summary,
            rescue_rows,
            cube,
            clnvc_types,
            pair_list,
            sv_summary,
            example_variants,
            gene_model_results,
        )
    )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
