#!/usr/bin/env python3
"""Exploratory report on ClinVar review-star cutoffs and classification mix.

Section numbers below match the rendered report's headings (1 Purpose,
2 Illustrative examples, 3 How ClinVar curation works, then:)

4. Review-star cutoff impact (per-submission evidence, from
   submission_summary.txt): `clinvar_helpers.variant_records_to_disease()`
   drops any submission record whose review status maps to fewer than
   `star_min` stars (production currently hardcodes `var2disease_star_min =
   3`). Section 6 re-runs that same disease mapping at star_min in
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
   files (see section {S.crossfilter}, which uses CLNREVSTAT from the VCF directly and
   *does* see this tier), never per-submission-record here.

4b. Multi-submitter concordance rescue (same section {S.star_cutoff} block): gene-disease pairs whose ONLY
   supporting evidence is star 0/1 (never reach star_min=2) are additionally
   checked for >=2 distinct Submitters independently reporting the *same*
   MONDO disease with the *exact same* ClinicalSignificance (Pathogenic/
   Likely pathogenic family only), regardless of each submitter's own
   review status. This reconstructs a proxy for the missing per-record
   "2 star" tier from raw submitter agreement. Each such pair is also
   annotated with whether any contributing record used
   CollectionMethod == "literature only" (a publication proxy) --
   informational only, not used as a filter.

5. Monarch KG gene-disease associations vs ClinVar: reconciles the pairs
   section {S.star_cutoff} derives against the curated gene-disease edges in the Monarch
   KG (monarch-kg.tar.gz, ~330MB, auto-fetched and reduced to small derived
   extracts on first run). Monarch's edges come from OMIM / Orphanet /
   ClinGen and are assertions about the GENE; ClinVar's pairs are DERIVED by
   crossing (variant->gene) with (variant->disease), so this measures how
   much submitted variant evidence stands behind each curated association
   and what ClinVar implies that no curator has asserted. Because MONDO
   carries several co-existing terms over one clinical area, exact-id
   mismatches are additionally classified by whether Monarch links the same
   gene to an ancestor or descendant of ClinVar's term. That classification
   is a *candidate* list, not a duplicate count -- see the SCN1A worked
   example, where nearby terms split the evidence but MONDO holds them
   apart deliberately.

6. Clinical significance & review-status crossfilter (variant-level, from
   clinvar.vcf.gz directly): CLNREVSTAT here IS ClinVar's aggregate
   per-variant review status (so all 5 star levels are real, including
   "2 stars"). CLNSIG is a free-form, sometimes pipe/slash-combined string
   (100+ distinct values observed) bucketed via classify_clnsig() into
   P / LP / P-LP / B / LB / B-LB / VUS / Conflicting / Other / Not
   classified. An interactive crossfilter widget lets you toggle
   significance / star / variant-type / size / literature / concordance /
   STRchive / production-filter dimensions and see live variant counts. That
   last dimension re-runs section {S.star_cutoff}'s exact inclusion criterion
   (>=PRODUCTION_STAR_MIN stars on some individual submission record, OR
   >=MIN_CONCORDANT_SUBMITTERS concordant submitters) per variant, so the
   crossfilter can select the population production actually ingests --
   necessary because the star dimension here is CLNREVSTAT (variant-level)
   while production filters on per-record ReviewStatus, so ticking ">=3 stars"
   in the star panel does NOT reproduce the production filter.
   The STRchive dimension flags
   whether a variant's genomic footprint overlaps one of STRchive's ~80
   curated pathogenic short-tandem-repeat loci (github.com/dashnowlab/
   STRchive, fetched and cached on first run), independent of ClinVar's own
   CLNVC "Microsatellite" label.

7. New ClinVar ingest recommendation: consolidates sections {S.phenotype_terms}-7 into four
   concrete, prioritised changes to src/clinvar_helpers.py -- decide
   inclusion per (gene, disease) pair rather than per variant; test
   concordance over the MONDO hierarchy instead of exact ids; attribute a
   variant to its causal gene rather than every locus in GENEINFO; and drop
   deprecated MONDO terms. None are implemented; the section reports the
   projected pair counts for each threshold so the tradeoff is explicit.

8. Structural variants / CNVs -- what's NOT in the VCF: clinvar.vcf.gz (the
   only input sections {S.phenotype_terms}-7 use) requires a fixed genomic position + REF/ALT,
   which structurally excludes copy-number gain/loss, translocations, and
   other large rearrangements. ClinVar publishes those separately, in
   variant_summary.txt.gz (same tab_delimited/ directory as
   submission_summary.txt, Start/Stop coordinate ranges instead of a fixed
   REF/ALT). This section downloads that file (~440MB, auto-fetched on first
   run) and reports real counts, size distributions, and disease-linkage
   rates for those SV/CNV types. Note the build split matters enormously
   here: of ~49.8k distinct CNV variants only ~14.4k have a GRCh38 row (and
   just 8% of those resolve to a disease id -- they are largely ClinGen-style
   dosage-sensitivity regions with no single named condition), while the
   ~35.4k that ClinVar has only ever placed on GRCh37 are ~70% resolved.
   Filtering to GRCh38 alone would therefore discard most of ClinVar's
   curated CNV-disease evidence; see ASSEMBLY_PREFERENCE. Also that
   GeneSymbol degrades to unparseable
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
    ASSEMBLY_PREFERENCE,
    map_CLNDISDB_to_mondo,
    map_mondo_to_hp,
    parse_CLNDISDB,
    format_id_to_map,
    make_genes_from_row,
    make_variant_gene_map,
    make_medgen_to_mondo_map,
    make_mondo_map,
    make_variant_record_map,
    predicate_map,
    aggregate_star_min,
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

def variant_genes_for(variant_genes: dict, varid: str) -> tuple[list, list]:
    """The gene ClinVar attributes this variant to, as parallel id/symbol lists.

    Mirrors what the production transform does: gene attribution comes from
    variant_summary.txt.gz (clinvar_helpers.make_variant_gene_map), NOT from the
    VCF's GENEINFO field. GENEINFO is populated positionally -- every locus
    overlapping the variant -- so building pairs from it invented associations for
    antisense transcripts (CFTR-AS1), readthrough fusions (NPHP3-ACAD11), locus
    control regions (GH-LCR) and NCBI LOC placeholders, each inheriting the real
    gene's whole disease and submitter roster. That was a defect in the ingest and
    is fixed; this report reads the same corrected source so its pair counts match
    what the transform actually emits.

    Returns empty lists where ClinVar declines to attribute a gene (GeneID == -1) --
    no GENEINFO fallback, matching process_row().
    """
    entry = variant_genes.get(varid)
    if entry is None:
        return [], []
    return [entry[0]], [entry[1]]


# ---------------------------------------------------------------------------
# Section 2: illustrative examples (variant types; one gene, many diseases)
# ---------------------------------------------------------------------------

EXAMPLE_VARIANT_IDS = {
    "snv": "55630",  # BRCA1 NC_000017.11:g.43045711G>C, Pathogenic, 4-star
    "indel_del": "53158",  # CFTR NC_000007.14:g.117480086_117480108del, Pathogenic, 4-star
    "indel_ins": "973964",  # RPE65 NC_000001.11:g.68431505_68431506insCAGC, Pathogenic, 3-star
    "str": "183387",  # FMR1 CGG repeat expansion NC_000023.11:g.147912051CGG[201], Fragile X syndrome
}
# Structural variants -- not in clinvar.vcf.gz, only in variant_summary.txt.gz (see section the structural-variants section)
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
    {
        "gene": "STAT3",
        "diseases": {"MONDO:0007818", "MONDO:0014414"},  # HIES (LOF), STAT3 GOF disease
        "heading": "Gain vs. loss of function, same gene: STAT3",
        "blurb": (
            "<strong>STAT3</strong> shows the LMNA/SMCHD1 pattern from a different angle: instead of two "
            "diseases telling apart by <em>where</em> a variant falls, these two are opposite <em>directions</em> "
            "of the same molecular event. Dominant-negative <strong>loss-of-function</strong> variants disrupt "
            "STAT3's DNA-binding ability and cause autosomal dominant hyper-IgE (Job) syndrome &mdash; recurrent "
            "infections, eczema, and markedly elevated IgE from impaired immune signaling. "
            "<strong>Gain-of-function</strong> variants push STAT3 transcriptional activity <em>up</em> instead, "
            "and cause an almost opposite clinical picture: early-onset multisystem autoimmunity and "
            "lymphoproliferation, from an immune system that won't stop signaling. Below, {count:,} "
            "Pathogenic/Likely pathogenic SNVs are plotted the same way as above. In UniProt's own curated "
            "variant calls, the loss-of-function/HIES variants cluster tightly in the DNA-binding domain "
            "(&asymp;aa 380-465), while the gain-of-function variants scatter far more broadly, from the "
            "coiled-coil domain through the SH2 and transactivation domains &mdash; the reverse clustering "
            "direction from SMCHD1 above, a reminder that there's no universal rule for which mechanism "
            "clusters and which spreads."
        ),
    },
    {
        "gene": "CARD11",
        "diseases": {"MONDO:0014081", "MONDO:0014645"},  # SCID due to CARD11 deficiency (LOF), BENTA disease (GOF)
        "heading": "Gain vs. loss of function, same gene: CARD11",
        "blurb": (
            "<strong>CARD11</strong> is part of the CBM complex (CARD11-BCL10-MALT1) that activates NF-&kappa;B "
            "signaling in lymphocytes. <strong>Loss-of-function</strong> variants (biallelic, or dominant-negative "
            "in milder cases) prevent the complex from assembling properly and cause severe combined "
            "immunodeficiency &mdash; the signal never turns on. <strong>Gain-of-function</strong> variants "
            "instead lock the complex into constitutive, ligand-independent activation, causing BENTA disease "
            "(B cell Expansion with NF-&kappa;B and T cell Anergy) &mdash; the signal never turns off. Below, "
            "{count:,} Pathogenic/Likely pathogenic SNVs are plotted the same way as above. Real BENTA-associated "
            "positions cluster right at the CARD domain/coiled-coil boundary (&asymp;aa 123-134); "
            "immunodeficiency-associated variants are more spread, including large C-terminal deletions "
            "removing the guanylate-kinase-like domain entirely."
        ),
    },
    {
        "gene": "KCNQ2",
        "heading": "Same clinical label, opposite mechanisms: KCNQ2",
        "blurb": (
            "<strong>KCNQ2</strong> is the sharpest example here of why gene-level (or even position-level) "
            "interpretation isn't enough. Unlike the examples above, ClinVar's own clinical diagnosis labels "
            "for KCNQ2 &mdash; benign familial neonatal seizures, developmental and epileptic encephalopathy "
            "&mdash; do <em>not</em> separate by mechanism: both <strong>loss-of-function</strong> (the classic, "
            "common mechanism, from mild self-limited neonatal seizures to severe dominant-negative "
            "encephalopathy) and rare <strong>gain-of-function</strong> variants can produce encephalopathy, and "
            "their variant positions are extensively intermixed below, not separated by domain the way the other "
            "three examples are. UniProt directly curates one confirmed gain-of-function variant, p.Arg201Cys, "
            "in the S4 voltage-sensor segment &mdash; \"results in loss of voltage-dependent channel gating and "
            "highly increased potassium currents,\" the mechanistic opposite of the surrounding loss-of-function "
            "variants (highlighted in purple below; everything else shown is the real, unfiltered ClinVar "
            "population, colored by ClinVar's own diagnosis label). "
            "<strong>This is also the one with a real prescribing consequence:</strong> in published cohorts, "
            "loss-of-function KCNQ2-DEE has responded well to sodium-channel blockers (carbamazepine, phenytoin, "
            "oxcarbazepine), while gain-of-function variants have been proposed to need the "
            "<em>opposite</em> approach &mdash; a Kv7-channel <em>blocker</em> (amitriptyline) rather than a "
            "channel opener. A drug chosen for the wrong mechanism could plausibly make either patient worse, "
            "which is precisely the argument for representing molecular mechanism, not just gene identity, "
            "in the knowledge graph."
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



# ---------------------------------------------------------------------------
# SCN1A: MONDO's real disease hierarchy, the gene, and variant sharing
# ---------------------------------------------------------------------------

SCN1A_SUBGRAPH_FILE = "scn1a_mondo_subgraph.json"

# Row assignment follows MONDO's actual subclass depth (fetched from OLS4); only the
# horizontal placement is chosen for legibility.
SCN1A_LAYOUT = {
    "MONDO:0005027": (350, 30, 240),
    "MONDO:0015650": (20, 120, 260),
    "MONDO:0100620": (620, 120, 320),
    "MONDO:0100022": (20, 220, 260),
    "MONDO:0100062": (620, 235, 320),
    "MONDO:0800490": (20, 330, 300),
    "MONDO:0800491": (150, 460, 280),
    "MONDO:0100135": (470, 460, 220),
    "MONDO:0100079": (720, 460, 220),
    "MONDO:0030268": (720, 555, 220),
}
SCN1A_GENE_BOX = (150, 780, 280)

# average glyph width as a fraction of font-size for the sans stack used here; used to
# wrap label text so it stays inside its box rather than spilling past the border
_CHAR_W = 0.53


def _fit_lines(text: str, box_w: int, font_size: float, max_lines: int = 3) -> list:
    """Greedy word-wrap to fit box_w, with an ellipsis if it still will not fit."""
    per_line = max(8, int((box_w - 16) / (font_size * _CHAR_W)))
    words, lines, cur = text.split(), [], ""
    for word in words:
        trial = f"{cur} {word}".strip()
        if len(trial) <= per_line:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = word
            if len(lines) == max_lines:
                break
    if cur and len(lines) < max_lines:
        lines.append(cur)
    if not lines:
        return [text[:per_line]]
    consumed = len(" ".join(lines))
    if consumed < len(text) - 1:
        last = lines[-1]
        lines[-1] = (last[: per_line - 1] + "\u2026") if len(last) + 1 > per_line else last + "\u2026"
    return lines


def load_scn1a_subgraph(data_dir: Path) -> dict:
    """MONDO's real hierarchy for SCN1A's main disease terms, cached from OLS4, with
    per-term variant/submitter counts and the variant overlap between terms."""
    path = data_dir / SCN1A_SUBGRAPH_FILE
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def _short_term(label: str) -> str:
    """Compact a MONDO label for one-line use, keeping what distinguishes it. Several of
    SCN1A's terms share the prefix "developmental and epileptic encephalopathy", so a plain
    truncation renders different diseases identically."""
    short = label.replace("developmental and epileptic encephalopathy", "DEE")
    short = short.replace("neonatal/infantile-onset epilepsy syndrome with", "neonatal/infantile")
    return short if len(short) <= 44 else short[:43] + "\u2026"


def diagram_scn1a_hierarchy(sub: dict) -> str:
    """SCN1A's disease terms as MONDO actually relates them, with the gene attached and
    the variant overlap between terms drawn on top.

    The structure is MONDO's own (OLS4 `hierarchicalParents`), not a simplification. The
    point it makes is the opposite of what a quick look suggests: these are *distinct*
    syndromes MONDO deliberately separates -- MONDO:0800491 is Ohtahara syndrome / early
    infantile epileptic encephalopathy (onset <=3 months, suppression-burst EEG), and the
    Dravet entry carries an explicit note that it is a distinct class from DEE 6A. Variants
    shared between them are therefore evidence of SCN1A's phenotypic spectrum, not of
    duplicate terminology.
    """
    if not sub:
        return "<p class='subtitle'>SCN1A subgraph not cached -- see load_scn1a_subgraph().</p>"

    labels, stats, parents = sub["labels"], sub["stats"], sub["parents"]
    shares = sub.get("shares", {})
    boxes, edges, notes, share_labels = [], [], [], []

    def centre(m):
        x, y, w = SCN1A_LAYOUT[m]
        return x + w / 2, y

    heights = {}
    for _m, (_x, _y, _w) in SCN1A_LAYOUT.items():
        _st = stats.get(_m, {})
        _n = len(_fit_lines(labels.get(_m, _m), _w, 11))
        _h = 8 + 13 * _n + 13 + (15 if _st.get("variants") else 14)
        if _st.get("submitters") == 1 and _st.get("lead"):
            _h += 12
        heights[_m] = _h + 6

    for child, ps in parents.items():
        if child not in SCN1A_LAYOUT:
            continue
        cx, cy = centre(child)
        for i, parent in enumerate(ps):
            if parent not in SCN1A_LAYOUT:
                continue
            px, py, pw = SCN1A_LAYOUT[parent]
            x2, y2 = px + pw / 2, py + heights.get(parent, 62)
            dash = " stroke-dasharray='5,4'" if i else ""
            edges.append(
                f"<path d='M{cx:.0f},{cy:.0f} C{cx:.0f},{(cy + y2) / 2:.0f} "
                f"{x2:.0f},{(cy + y2) / 2:.0f} {x2:.0f},{y2:.0f}' fill='none' "
                f"stroke='#94a3b8' stroke-width='1.6'{dash} marker-end='url(#arr-scn)'/>"
            )

    for mondo, (x, y, w) in SCN1A_LAYOUT.items():
        st = stats.get(mondo, {})
        n_var, n_sub = st.get("variants", 0), st.get("submitters", 0)
        lead = st.get("lead", "")
        if n_var == 0:
            fill, stroke, txt = "#f1f5f9", "#cbd5e1", "#475569"
        elif n_sub == 1:
            fill, stroke, txt = "#fee2e2", "#b91c1c", "#7f1d1d"
        elif n_sub >= 20:
            fill, stroke, txt = "#dcfce7", "#15803d", "#14532d"
        else:
            fill, stroke, txt = "#fef9c3", "#a16207", "#713f12"

        h = heights[mondo]
        boxes.append(
            f"<rect x='{x}' y='{y}' width='{w}' height='{h}' rx='6' fill='{fill}' "
            f"stroke='{stroke}' stroke-width='1.5'/>"
        )
        cy = y + 16
        for line in _fit_lines(labels.get(mondo, mondo), w, 11):
            boxes.append(
                f"<text x='{x + 8}' y='{cy}' font-size='11' font-weight='600' fill='{txt}'>"
                f"{line}</text>"
            )
            cy += 13
        boxes.append(
            f"<text x='{x + 8}' y='{cy}' font-size='9' font-family='ui-monospace,monospace' "
            f"fill='#64748b'>{mondo}</text>"
        )
        cy += 15
        if n_var:
            boxes.append(
                f"<text x='{x + 8}' y='{cy}' font-size='10' fill='{txt}'>"
                f"<tspan font-weight='700'>{n_var:,}</tspan> variants &#183; "
                f"<tspan font-weight='700'>{n_sub}</tspan> lab{'' if n_sub == 1 else 's'}</text>"
            )
            if n_sub == 1 and lead:
                cy += 12
                for line in _fit_lines("only: " + lead, w, 8.5, max_lines=1):
                    boxes.append(
                        f"<text x='{x + 8}' y='{cy}' font-size='8.5' fill='#b91c1c'>{line}</text>"
                    )
        else:
            boxes.append(
                f"<text x='{x + 8}' y='{cy}' font-size='9.5' fill='#94a3b8'>"
                f"no SCN1A variants here</text>"
            )

    # OMIM / Orphanet coverage, drawn above each box in its own layer so it can be
    # toggled off. Which source vocabulary a term carries explains which labs land on it.
    xref_layer = []
    xrefs = sub.get("xrefs", {})
    for mondo, (x, y, w) in SCN1A_LAYOUT.items():
        xr = xrefs.get(mondo, {})
        bits = []
        if xr.get("OMIM"):
            bits.append(("#7c3aed", ", ".join(xr["OMIM"])))
        if xr.get("Orphanet"):
            bits.append(("#0891b2", ", ".join(xr["Orphanet"])))
        if not bits:
            bits.append(("#94a3b8", "no OMIM / Orphanet id"))
        bx = x
        for colour, text in bits:
            pw = len(text) * 5.2 + 12
            xref_layer.append(
                f"<rect x='{bx:.0f}' y='{y - 17}' width='{pw:.0f}' height='14' rx='7' "
                f"fill='#ffffff' stroke='{colour}' stroke-width='1'/>"
                f"<text x='{bx + 6:.0f}' y='{y - 6.5}' font-size='8.5' fill='{colour}'>{text}</text>"
            )
            bx += pw + 4

    # the gene, wired to every term that carries its variants
    gx, gy, gw = SCN1A_GENE_BOX
    boxes.append(
        f"<rect x='{gx}' y='{gy}' width='{gw}' height='46' rx='6' fill='#1e293b' stroke='#0f172a'/>"
        f"<text x='{gx + 12}' y='{gy + 20}' font-size='13' font-weight='700' fill='#f8fafc'>SCN1A</text>"
        f"<text x='{gx + 12}' y='{gy + 36}' font-size='9.5' fill='#cbd5e1'>"
        f"NCBIGene:6323 &#183; one gene, many syndromes</text>"
    )
    for mondo in SCN1A_LAYOUT:
        if not stats.get(mondo, {}).get("variants"):
            continue
        tx, ty, tw = SCN1A_LAYOUT[mondo]
        edges.append(
            f"<path d='M{gx + gw / 2:.0f},{gy} C{gx + gw / 2:.0f},{gy - 40} "
            f"{tx + tw / 2:.0f},{ty + 120} {tx + tw / 2:.0f},{ty + heights.get(mondo, 62)}' fill='none' "
            f"stroke='#334155' stroke-width='1' stroke-dasharray='2,3' opacity='0.55'/>"
        )

    # Variant sharing drawn as connectors, laid out so none of them cross.
    #
    # Two rules make that work. Anchors along a box's bottom edge are ordered by how far
    # away the partner is -- the most distant partner takes the outermost anchor -- so
    # spans nest rather than partially overlap. Lanes are then assigned widest-span-first,
    # putting the widest connector furthest from the boxes. A vertical therefore never
    # drops through another connector's horizontal.
    SHARE_COLOR = "#ea580c"
    top = [(pair, n) for pair, n in sorted(shares.items(), key=lambda kv: -kv[1])
           if all(t in SCN1A_LAYOUT for t in pair.split("|"))][:3]
    max_n = max((n for _p, n in top), default=1)

    def box_centre(t):
        x, _y, w = SCN1A_LAYOUT[t]
        return x + w / 2

    # per box: partners, ordered most-distant first, then anchored outside-in
    anchors = {}
    for term in {t for pair, _n in top for t in pair.split("|")}:
        partners = []
        for pair, _n in top:
            a, b = pair.split("|")
            if term in (a, b):
                partners.append((pair, b if term == a else a))
        cx = box_centre(term)
        # a partner to the right wants a left-hand anchor when it is the far one, so the
        # wider span encloses the narrower
        partners.sort(key=lambda pp: -abs(box_centre(pp[1]) - cx))
        x, _y, w = SCN1A_LAYOUT[term]
        n_p = len(partners)
        for i, (pair, partner) in enumerate(partners):
            if box_centre(partner) > cx:
                frac = (i + 1) / (n_p + 1)          # further partner -> further left
            else:
                frac = 1 - (i + 1) / (n_p + 1)      # mirrored for left-hand partners
            anchors[(term, pair)] = x + w * frac

    # widest span gets the lowest lane
    spans = {}
    for pair, _n in top:
        a, b = pair.split("|")
        xa, xb = anchors[(a, pair)], anchors[(b, pair)]
        spans[pair] = abs(xa - xb)
    ordered = sorted(top, key=lambda pn: -spans[pn[0]])
    lane_of = {}
    lane_y = 660 + 40 * (len(ordered) - 1)
    for pair, _n in ordered:
        lane_of[pair] = lane_y
        lane_y -= 40

    for pair, n in top:
        a, b = pair.split("|")
        xa, xb = anchors[(a, pair)], anchors[(b, pair)]
        abot = SCN1A_LAYOUT[a][1] + heights.get(a, 62)
        bbot = SCN1A_LAYOUT[b][1] + heights.get(b, 62)
        ly = lane_of[pair]
        width = 3 + 6 * (n / max_n)
        notes.append(
            f"<path d='M{xa:.0f},{abot:.0f} L{xa:.0f},{ly} L{xb:.0f},{ly} "
            f"L{xb:.0f},{bbot:.0f}' fill='none' stroke='{SHARE_COLOR}' "
            f"stroke-width='{width:.1f}' stroke-linecap='round' stroke-linejoin='round'/>"
        )
        mid = (xa + xb) / 2
        label = f"{n:,} variants reported to both"
        box_w = len(label) * 5.6 + 16
        share_labels.append(
            f"<rect x='{mid - box_w / 2:.0f}' y='{ly - 10}' width='{box_w:.0f}' "
            f"height='20' rx='10' fill='#ffffff' stroke='{SHARE_COLOR}' stroke-width='1.4'/>"
            f"<text x='{mid:.0f}' y='{ly + 4}' font-size='10.5' font-weight='700' "
            f"text-anchor='middle' fill='#9a3412'>{label}</text>"
            f"<text x='{mid:.0f}' y='{ly + 24}' font-size='9' text-anchor='middle' "
            f"fill='#9a3412'>{_short_term(labels.get(a, a))} &#8596; "
            f"{_short_term(labels.get(b, b))}</text>"
        )

    legend = (
        "<rect x='20' y='10' width='920' height='58' rx='6' fill='#f8fafc' stroke='#e2e8f0'/>"
        "<text x='32' y='27' font-size='10.5' font-weight='700' fill='#334155'>Legend</text>"
        # edge kinds
        "<line x1='90' y1='38' x2='125' y2='38' stroke='#94a3b8' stroke-width='1.6'/>"
        "<text x='131' y='41' font-size='9.5' fill='#475569'>subclass_of (first parent)</text>"
        "<line x1='285' y1='38' x2='320' y2='38' stroke='#94a3b8' stroke-width='1.6' "
        "stroke-dasharray='5,4'/>"
        "<text x='326' y='41' font-size='9.5' fill='#475569'>second parent (MONDO is a DAG)</text>"
        "<line x1='545' y1='38' x2='580' y2='38' stroke='#334155' stroke-width='1' "
        "stroke-dasharray='2,3'/>"
        "<text x='586' y='41' font-size='9.5' fill='#475569'>SCN1A variants reported here</text>"
        "<line x1='90' y1='57' x2='125' y2='57' stroke='#ea580c' stroke-width='7' "
        "stroke-linecap='round'/>"
        "<text x='131' y='60' font-size='9.5' fill='#9a3412'>variants reported to BOTH terms "
        "&#8212; SCN1A's phenotypic spectrum, not duplicate terminology</text>"
        # colour key
        "<rect x='545' y='51' width='12' height='11' rx='2' fill='#fee2e2' stroke='#b91c1c'/>"
        "<text x='561' y='60' font-size='9.5' fill='#475569'>1 lab</text>"
        "<rect x='600' y='51' width='12' height='11' rx='2' fill='#fef9c3' stroke='#a16207'/>"
        "<text x='616' y='60' font-size='9.5' fill='#475569'>&lt;20 labs</text>"
        "<rect x='670' y='51' width='12' height='11' rx='2' fill='#dcfce7' stroke='#15803d'/>"
        "<text x='686' y='60' font-size='9.5' fill='#475569'>20+ labs</text>"
        "<rect x='745' y='51' width='12' height='11' rx='2' fill='#f1f5f9' stroke='#cbd5e1'/>"
        "<text x='761' y='60' font-size='9.5' fill='#475569'>no SCN1A variants</text>"
    )

    # z-order: connectors and hierarchy edges first, boxes on top so nothing is drawn
    # over them, then the connector labels which must stay readable above everything
    return (
        "<svg viewBox='0 0 960 960' xmlns='http://www.w3.org/2000/svg' "
        "style='width:100%; height:auto;'>"
        "<defs><marker id='arr-scn' viewBox='0 0 10 10' refX='9' refY='5' markerWidth='5' "
        "markerHeight='5' orient='auto-start-reverse'>"
        "<path d='M 0 0 L 10 5 L 0 10 z' fill='#94a3b8'/></marker></defs>"
        + legend
        + "<g transform='translate(0,90)'>"
        + "".join(edges)
        + "<g id='scn1a-shares'>" + "".join(notes) + "</g>"
        + "".join(boxes)
        + "<g id='scn1a-shares-labels'>" + "".join(share_labels) + "</g>"
        + "<g id='scn1a-xrefs' style='display:none'>" + "".join(xref_layer) + "</g>"
        + "</g>"
        + "<text x='30' y='950' font-size='9.5' fill='#64748b'>"
        "hierarchy from MONDO via OLS4 &#183; counts from this ClinVar release</text>"
        "</svg>"
    )


def diagram_bams_fshd2_venn():
    """Venn diagram of hallmark clinical features for BAMS vs. FSHD2 -- the
    two SMCHD1-linked diseases plotted above. Colors match that plot's own
    legend (FSHD2 = GENE_MODEL_COLORS[0], BAMS = GENE_MODEL_COLORS[1]) for
    visual continuity. The two circles barely overlap on purpose: despite
    coming from the same gene, these are clinically distinct, largely
    non-overlapping phenotypes (craniofacial/endocrine developmental disorder
    vs. progressive muscular dystrophy) -- the intersection is left labeled
    "none known" rather than populated with anything, since that's the real,
    literature-consistent answer, not an artifact of leaving it blank."""
    width, height = 620, 400
    fshd2_color, bams_color = "#4e79a7", "#f28e2b"
    lcx, lcy, rcx, rcy, r = 210, 190, 410, 190, 150
    line_h = 17

    parts = [_cartoon_svg(width, height)]
    parts.append(f"<circle cx='{lcx}' cy='{lcy}' r='{r}' fill='{fshd2_color}' fill-opacity='0.18' stroke='{fshd2_color}' stroke-width='2'/>")
    parts.append(f"<circle cx='{rcx}' cy='{rcy}' r='{r}' fill='{bams_color}' fill-opacity='0.18' stroke='{bams_color}' stroke-width='2'/>")

    # Center every line on its circle's own cx, and vertically center the whole block on
    # cy -- rather than left-anchoring at a fixed x, which let long lines run past the
    # circle's edge wherever it narrows away from the vertical middle.
    parts.append(_cartoon_label(lcx, 55, "FSHD2", size=16, color=fshd2_color, anchor="middle", weight="700"))
    fshd2_features = [
        "Progressive facial,",
        "scapular & humeral weakness",
        "Asymmetric muscle weakness",
        "Elevated creatine kinase",
        "Sensorineural hearing loss",
        "(subset of patients)",
    ]
    fshd2_start_y = lcy - (len(fshd2_features) - 1) * line_h / 2
    for i, line in enumerate(fshd2_features):
        parts.append(_cartoon_label(lcx, fshd2_start_y + i * line_h, line, size=10.5, color="#1e3a5f", anchor="middle"))

    parts.append(_cartoon_label(rcx, 55, "BAMS", size=16, color="#b45309", anchor="middle", weight="700"))
    bams_features = [
        "Arhinia (absent nose)",
        "Microphthalmia",
        "Choanal atresia",
        "Hypogonadotropic",
        "hypogonadism",
        "Anosmia",
    ]
    bams_start_y = rcy - (len(bams_features) - 1) * line_h / 2
    for i, line in enumerate(bams_features):
        parts.append(_cartoon_label(rcx, bams_start_y + i * line_h, line, size=10.5, color="#7c2d12", anchor="middle"))

    parts.append(_cartoon_label((lcx + rcx) / 2, lcy - 15, "no shared", size=10.5, color="#475569", anchor="middle", weight="700"))
    parts.append(_cartoon_label((lcx + rcx) / 2, lcy, "hallmark features", size=9.5, color="#475569", anchor="middle"))
    parts.append(_cartoon_label((lcx + rcx) / 2, lcy + 15, "(none known)", size=9, color="#94a3b8", anchor="middle"))

    parts.append(_cartoon_caption("Same gene, clinically non-overlapping phenotypes", width=width, y=height - 24))
    parts.append(_cartoon_caption("muscular dystrophy vs. craniofacial/endocrine syndrome", width=width, y=height - 10))
    parts.append("</svg>")
    return "".join(parts)


def diagram_smchd1_domains():
    """Linear protein domain map for human SMCHD1 (UniProt A6NHR9, 2005 aa).
    Domain boundaries (ATPase activity domain 111-702, SMC hinge 1720-1847)
    and the variant positions plotted below are real values pulled from
    UniProt's own feature/disease-variant annotations, not estimated --
    see https://rest.uniprot.org/uniprotkb/A6NHR9.json. This is protein-space,
    which is what actually explains BAMS's clustering: those variants land
    inside the ATPase domain specifically, while FSHD2 not only scatters
    point mutations across the whole protein but also includes large
    truncating deletions that remove everything downstream of their start
    point -- including the SMC hinge domain entirely."""
    width, height = 900, 300
    left_pad, right_pad = 60, 40
    bar_w = width - left_pad - right_pad
    protein_len = 2005

    def x_of(pos):
        return left_pad + (pos - 1) / (protein_len - 1) * bar_w

    bar_y, bar_h = 110, 18
    atpase_color, hinge_color = "#6366f1", "#a855f7"
    bams_color, fshd2_color = "#f28e2b", "#4e79a7"

    parts = [_cartoon_svg(width, height)]

    parts.append(
        f"<rect x='{left_pad}' y='{bar_y}' width='{bar_w}' height='{bar_h}' rx='4' "
        f"fill='{CARTOON_NEUTRAL_FILL}' stroke='{CARTOON_NEUTRAL_STROKE}' stroke-width='1.5'/>"
    )

    ax0, ax1 = x_of(111), x_of(702)
    parts.append(f"<rect x='{ax0:.1f}' y='{bar_y}' width='{ax1 - ax0:.1f}' height='{bar_h}' rx='4' fill='{atpase_color}' fill-opacity='0.8'/>")
    parts.append(_cartoon_label((ax0 + ax1) / 2, bar_y - 10, "ATPase (GHKL) domain", size=11, color=atpase_color, anchor="middle", weight="700"))
    parts.append(_cartoon_label((ax0 + ax1) / 2, bar_y + bar_h + 16, "111–702", size=9.5, color="#64748b", anchor="middle"))

    hx0, hx1 = x_of(1720), x_of(1847)
    parts.append(f"<rect x='{hx0:.1f}' y='{bar_y}' width='{max(hx1 - hx0, 4):.1f}' height='{bar_h}' rx='4' fill='{hinge_color}' fill-opacity='0.8'/>")
    parts.append(_cartoon_label((hx0 + hx1) / 2, bar_y - 10, "SMC hinge", size=11, color=hinge_color, anchor="middle", weight="700"))
    parts.append(_cartoon_label((hx0 + hx1) / 2, bar_y + bar_h + 16, "1720–1847", size=9.5, color="#64748b", anchor="middle"))

    parts.append(_cartoon_label(left_pad, bar_y + bar_h + 34, "aa 1", size=9.5, color="#94a3b8", anchor="start"))
    parts.append(_cartoon_label(left_pad + bar_w, bar_y + bar_h + 34, "aa 2005", size=9.5, color="#94a3b8", anchor="end"))

    bams_positions = [107, 129, 134, 135, 136, 137, 139, 141, 171, 242, 342, 345, 348, 400, 420, 473, 518, 523, 524, 552]
    bams_y = bar_y - 34
    for pos in bams_positions:
        x = x_of(pos)
        parts.append(f"<circle cx='{x:.1f}' cy='{bams_y}' r='3.5' fill='{bams_color}' fill-opacity='0.65' stroke='{bams_color}' stroke-width='1'/>")
    parts.append(_cartoon_label(left_pad - 8, bams_y + 3, "BAMS", size=10.5, color=bams_color, anchor="end", weight="700"))

    fshd2_positions = [110, 137, 194, 263, 353, 478, 527, 615, 690, 716]
    fshd2_y = bar_y + bar_h + 55
    for pos in fshd2_positions:
        x = x_of(pos)
        parts.append(f"<circle cx='{x:.1f}' cy='{fshd2_y}' r='3.5' fill='{fshd2_color}' fill-opacity='0.65' stroke='{fshd2_color}' stroke-width='1'/>")
    parts.append(_cartoon_label(left_pad - 8, fshd2_y + 3, "FSHD2", size=10.5, color=fshd2_color, anchor="end", weight="700"))

    deletions = [(138, 2005), (195, 2005), (344, 2005), (434, 2005)]
    for i, (start, end) in enumerate(deletions):
        y = fshd2_y + 16 + i * 12
        dx0, dx1 = x_of(start), x_of(end)
        parts.append(
            f"<line x1='{dx0:.1f}' y1='{y}' x2='{dx1:.1f}' y2='{y}' stroke='{fshd2_color}' "
            f"stroke-width='2.5' stroke-linecap='round' opacity='0.5'/>"
        )
    parts.append(_cartoon_label(left_pad - 8, fshd2_y + 16 + 1.5 * 12, "(+truncating", size=8.5, color="#94a3b8", anchor="end"))
    parts.append(_cartoon_label(left_pad - 8, fshd2_y + 16 + 1.5 * 12 + 11, "deletions)", size=8.5, color="#94a3b8", anchor="end"))

    parts.append(_cartoon_caption("SMCHD1 (UniProt A6NHR9, 2005 aa) — real domain boundaries + variant positions", width=width, y=height - 24))
    parts.append(_cartoon_caption("BAMS clusters inside the ATPase domain; FSHD2 spreads across it and truncates the hinge", width=width, y=height - 10))
    parts.append("</svg>")
    return "".join(parts)


def diagram_protein_domains(gene_symbol, uniprot_id, protein_len, domain_boxes, variant_series, caption_lines):
    """Generic linear protein domain map, generalized from diagram_smchd1_domains
    (which predates this and keeps its own bespoke version) so STAT3/CARD11/KCNQ2/
    MSH2 don't each need a hand-written copy. All positions are real, UniProt-sourced
    (formal DOMAIN/TRANSMEM features and/or curated Natural variant calls) -- see each
    gene's GENE_MODEL_EXAMPLES blurb for the exact UniProt accession and citations.

    domain_boxes: list of (start, end, label, color) drawn directly on the backbone.
    variant_series: list of dicts, each a labeled row of dots:
      {"label", "color", "positions": [int,...], "side": "above"|"below" (default
      "above"), "highlight_positions": optional set of positions to draw bigger,
      outlined, and more opaque -- for calling out one or two specific variants
      (e.g. a single confirmed gain-of-function residue) against a muted background
      series, rather than a uniform population.}
    caption_lines: 1-2 strings, rendered centered at the bottom.
    """
    width = 900
    left_pad, right_pad = 175, 40
    bar_w = width - left_pad - right_pad
    bar_h = 18
    row_h = 22

    def x_of(pos):
        return left_pad + (pos - 1) / (protein_len - 1) * bar_w

    above = [s for s in variant_series if s.get("side", "above") == "above"]
    below = [s for s in variant_series if s.get("side") == "below"]
    # 60px reserves room for the two staggered domain-label/range-text tiers
    # around the backbone, regardless of how many variant-series rows stack above/below it.
    above_h = len(above) * row_h + 60
    below_h = len(below) * row_h + 60
    track_y = above_h
    height = track_y + bar_h + below_h + 24 * len(caption_lines) + 20

    parts = [_cartoon_svg(width, height)]
    parts.append(
        f"<rect x='{left_pad}' y='{track_y}' width='{bar_w}' height='{bar_h}' rx='4' "
        f"fill='{CARTOON_NEUTRAL_FILL}' stroke='{CARTOON_NEUTRAL_STROKE}' stroke-width='1.5'/>"
    )
    for i, (start, end, label, color) in enumerate(domain_boxes):
        x0, x1 = x_of(start), x_of(end)
        tier = i % 2
        parts.append(
            f"<rect x='{x0:.1f}' y='{track_y}' width='{max(x1 - x0, 3):.1f}' height='{bar_h}' rx='3' "
            f"fill='{color}' fill-opacity='0.85'/>"
        )
        parts.append(
            _cartoon_label((x0 + x1) / 2, track_y - 8 - tier * 16, label, size=8.5, color=color, anchor="middle", weight="700")
        )
        parts.append(
            _cartoon_label((x0 + x1) / 2, track_y + bar_h + 14 + tier * 13, f"{start}-{end}", size=7.5, color="#64748b", anchor="middle")
        )

    aa_line_y = track_y + bar_h + 44
    parts.append(_cartoon_label(left_pad, aa_line_y, "aa 1", size=9, color="#94a3b8", anchor="start"))
    parts.append(_cartoon_label(left_pad + bar_w, aa_line_y, f"aa {protein_len}", size=9, color="#94a3b8", anchor="end"))
    parts.append(
        _cartoon_label(
            left_pad + bar_w / 2, height - 24 * len(caption_lines) - 4,
            f"{gene_symbol} (UniProt {uniprot_id}, {protein_len} aa)", size=8.5, color="#94a3b8", anchor="middle",
        )
    )

    for i, s in enumerate(above):
        y = track_y - 44 - i * row_h
        parts.append(_cartoon_label(left_pad - 10, y + 3, s["label"], size=9.5, color=s["color"], anchor="end", weight="700"))
        hl = s.get("highlight_positions", set())
        for pos in s["positions"]:
            x = x_of(pos)
            is_hl = pos in hl
            r = 5.5 if is_hl else 3.2
            stroke = "#0f172a" if is_hl else s["color"]
            parts.append(
                f"<circle cx='{x:.1f}' cy='{y}' r='{r}' fill='{s['color']}' "
                f"fill-opacity='{0.9 if is_hl else 0.55}' stroke='{stroke}' stroke-width='{2 if is_hl else 1}'/>"
            )

    for i, s in enumerate(below):
        y = track_y + bar_h + 60 + i * row_h
        parts.append(_cartoon_label(left_pad - 10, y + 3, s["label"], size=9.5, color=s["color"], anchor="end", weight="700"))
        hl = s.get("highlight_positions", set())
        for pos in s["positions"]:
            x = x_of(pos)
            is_hl = pos in hl
            r = 5.5 if is_hl else 3.2
            stroke = "#0f172a" if is_hl else s["color"]
            parts.append(
                f"<circle cx='{x:.1f}' cy='{y}' r='{r}' fill='{s['color']}' "
                f"fill-opacity='{0.9 if is_hl else 0.55}' stroke='{stroke}' stroke-width='{2 if is_hl else 1}'/>"
            )

    for i, line in enumerate(caption_lines):
        parts.append(_cartoon_caption(line, width=width, y=height - 10 - (len(caption_lines) - 1 - i) * 14))

    parts.append("</svg>")
    return "".join(parts)


# Per-gene domain-diagram registry consumed by main()'s GENE_MODEL_EXAMPLES loop.
# Domain boundaries and variant positions below are pulled directly from each
# gene's UniProt entry (rest.uniprot.org/uniprotkb/<id>.json): "Domain"/"Coiled
# coil"/"Transmembrane" features for domain_boxes, "Natural variant" features
# (filtered by disease tag / functional-effect text) for variant_series. Where
# UniProt doesn't formally tag a domain (STAT3's coiled-coil/DNA-binding/
# transactivation domains; MSH2's five-domain MutS architecture), boundaries
# follow the structural-biology literature cited in each caption instead.
GENE_DOMAIN_DIAGRAMS = {
    "SMCHD1": {
        "svg_fn": diagram_smchd1_domains,
        "uniprot_id": "A6NHR9",
        "caption": (
            "In protein space, not just genomic space, the same story holds: BAMS variants sit inside the "
            "ATPase domain specifically, while FSHD2 includes large truncating deletions that remove the SMC "
            "hinge domain entirely."
        ),
    },
    "STAT3": {
        "svg_fn": lambda: diagram_protein_domains(
            "STAT3",
            "P40763",
            770,
            domain_boxes=[
                (138, 320, "Coiled-coil", "#6366f1"),
                (321, 465, "DNA-binding", "#f28e2b"),
                (580, 670, "SH2", "#16a34a"),
                (705, 770, "Transactivation", "#a855f7"),
            ],
            variant_series=[
                {
                    "label": "HIES1 (LOF)",
                    "color": "#2563eb",
                    "side": "above",
                    "positions": [382, 384, 389, 395, 423, 425, 437, 463, 611, 621, 622, 637, 644, 657],
                },
                {
                    "label": "ADMIO1 (GOF)",
                    "color": "#dc2626",
                    "side": "below",
                    "positions": [
                        70, 103, 107, 152, 166, 174, 218, 260, 278, 302, 325, 330, 331, 344, 353, 357,
                        361, 387, 389, 392, 394, 408, 415, 419, 420, 421, 422, 443, 447, 448, 506, 546,
                        640, 643, 646, 658, 663, 664, 703, 715, 716,
                    ],
                },
            ],
            caption_lines=[
                "STAT3 (UniProt P40763, 770 aa) — SH2 domain boundary is UniProt-annotated; CCD/DBD/TAD follow published structural boundaries",
                "HIES1 (loss-of-function) clusters in the DNA-binding domain; ADMIO1 (gain-of-function) scatters from the coiled-coil through the SH2/TAD",
            ],
        ),
        "uniprot_id": "P40763",
        "caption": (
            "In protein space, HIES1 (hyper-IgE syndrome, loss-of-function) variants cluster tightly around the "
            "DNA-binding domain, while ADMIO1 (autoimmune disease, gain-of-function) variants scatter across the "
            "coiled-coil, linker, SH2, and transactivation domains &mdash; the reverse clustering direction from "
            "SMCHD1 above."
        ),
    },
    "CARD11": {
        "svg_fn": lambda: diagram_protein_domains(
            "CARD11",
            "Q9BXL7",
            1154,
            domain_boxes=[
                (18, 110, "CARD", "#6366f1"),
                (130, 449, "Coiled-coil", "#f28e2b"),
                (667, 755, "PDZ", "#16a34a"),
                (973, 1140, "GUK", "#a855f7"),
            ],
            variant_series=[
                {
                    "label": "BENTA (GOF)",
                    "color": "#dc2626",
                    "side": "above",
                    "positions": [123, 134],
                    "highlight_positions": {123, 134},
                },
                {
                    "label": "IMD11 (LOF)",
                    "color": "#2563eb",
                    "side": "below",
                    "positions": [57, 194, 945, 975],
                },
            ],
            caption_lines=[
                "CARD11 (UniProt Q9BXL7, 1154 aa) — real domain boundaries + disease-variant positions",
                "BENTA (gain-of-function) sits right at the CARD/coiled-coil boundary; IMD11 (loss-of-function) spans into the GUK domain",
            ],
        ),
        "uniprot_id": "Q9BXL7",
        "caption": (
            "In protein space, the two confirmed BENTA (gain-of-function) positions sit right at the CARD "
            "domain/coiled-coil boundary, while IMD11 (loss-of-function) variants are more spread &mdash; "
            "including one inside the C-terminal guanylate-kinase-like domain."
        ),
    },
    "KCNQ2": {
        "svg_fn": lambda: diagram_protein_domains(
            "KCNQ2",
            "O43526",
            872,
            domain_boxes=[
                (91, 113, "S1", "#94a3b8"),
                (124, 145, "S2", "#94a3b8"),
                (164, 183, "S3", "#94a3b8"),
                (197, 215, "S4", "#6366f1"),
                (228, 253, "S5", "#94a3b8"),
                (264, 276, "H5", "#f28e2b"),
                (288, 314, "S6", "#94a3b8"),
            ],
            variant_series=[
                {
                    "label": "BFNS1 (LOF)",
                    "color": "#2563eb",
                    "side": "above",
                    "positions": [114, 154, 159, 196, 204, 208, 214, 217, 228, 243, 284, 333, 353, 358, 448, 547, 588, 637],
                },
                {
                    "label": "DEE7 (LOF)",
                    "color": "#f28e2b",
                    "side": "below",
                    "positions": [207, 210, 213, 234, 247, 266, 268, 276, 291, 294, 301, 306, 554, 561, 578, 581],
                },
                {
                    "label": "R201C (GOF)",
                    "color": "#a855f7",
                    "side": "above",
                    "positions": [201],
                    "highlight_positions": {201},
                },
            ],
            caption_lines=[
                "KCNQ2 (UniProt O43526, 872 aa) — real S1-S6 transmembrane topology + disease-variant positions",
                "BFNS1/DEE7 (loss-of-function) variants intermix across every segment; the single confirmed GOF variant (R201C) sits in S4",
            ],
        ),
        "uniprot_id": "O43526",
        "caption": (
            "In protein space, benign familial neonatal seizures (BFNS1) and developmental and epileptic "
            "encephalopathy (DEE7) &mdash; both predominantly loss-of-function &mdash; intermix across every "
            "transmembrane segment, unlike the clean domain separation in the examples above. The one "
            "UniProt-confirmed gain-of-function variant, p.Arg201Cys, sits in the S4 voltage-sensor segment "
            "(highlighted in purple), the mechanistic opposite of its loss-of-function neighbors."
        ),
    },
}


# ---------------------------------------------------------------------------
# Closing proposal: current vs. proposed Biolink model for variants
# ---------------------------------------------------------------------------

DIAGRAM_NODE_FILL = "#1e293b"
DIAGRAM_NODE_TEXT = "#f8fafc"
DIAGRAM_EDGE_COLOR = "#64748b"


def _diagram_node(x, y, w, h, lines, fill=DIAGRAM_NODE_FILL, text_color=DIAGRAM_NODE_TEXT, font_size=12.5):
    parts = [f"<rect x='{x}' y='{y}' width='{w}' height='{h}' rx='8' fill='{fill}'/>"]
    n = len(lines)
    line_h = font_size + 4.5
    start_y = y + h / 2 - (n - 1) * line_h / 2 + font_size / 3
    for i, line in enumerate(lines):
        weight = "font-weight='700'" if i == 0 else ""
        parts.append(
            f"<text x='{x + w / 2}' y='{start_y + i * line_h}' font-size='{font_size}' {weight} "
            f"text-anchor='middle' fill='{text_color}'>{line}</text>"
        )
    return "".join(parts)


def _diagram_arrow(x1, y1, x2, y2, color=DIAGRAM_EDGE_COLOR, marker="url(#arrow-diag)"):
    return f"<line x1='{x1}' y1='{y1}' x2='{x2}' y2='{y2}' stroke='{color}' stroke-width='1.5' marker-end='{marker}'/>"


def _diagram_edge_label(x, y, lines, color="#334155", anchor="start", font_size=10.5):
    parts = []
    for i, line in enumerate(lines):
        fill = color if i == 0 else "#64748b"
        parts.append(f"<text x='{x}' y='{y + i * (font_size + 4)}' font-size='{font_size}' fill='{fill}' text-anchor='{anchor}'>{line}</text>")
    return "".join(parts)


ARROW_MARKER_DEFS = (
    "<defs><marker id='arrow-diag' markerWidth='8' markerHeight='8' refX='7' refY='4' orient='auto'>"
    "<path d='M0,0 L8,4 L0,8 Z' fill='#64748b'/></marker>"
    "<marker id='arrow-diag-red' markerWidth='8' markerHeight='8' refX='7' refY='4' orient='auto'>"
    "<path d='M0,0 L8,4 L0,8 Z' fill='#dc2626'/></marker></defs>"
)


def diagram_sequencevariant_hierarchy():
    """Biolink Model class hierarchy for SequenceVariant, per
    https://biolink.github.io/biolink-model/SequenceVariant/. Shown here because
    production emits every ClinVar node -- SNV, indel, or STR alike -- as the
    parent SequenceVariant category, never the more specific Snv subclass. That's
    not just true of the ingest code (src/clinvar_helpers.py); it's confirmed live
    in the production Monarch graph too, e.g. a real query against
    https://api-v3.monarchinitiative.org/v3/api/entity/CLINVAR:55630 (a BRCA1 SNV)
    returns "category":"biolink:SequenceVariant", not biolink:Snv."""
    width, height = 660, 280
    parts = [_cartoon_svg(width, height), ARROW_MARKER_DEFS]

    parts.append(_diagram_node(20, 12, 210, 42, ["NamedThing"], font_size=12, fill="#475569"))
    parts.append(
        _diagram_node(
            20, 90, 210, 55,
            ["BiologicalEntity", "(+ GenomicEntity, PhysicalEssence,", "OntologyClass mixins)"],
            font_size=10, fill="#475569",
        )
    )
    parts.append(
        _diagram_node(
            20, 185, 260, 55,
            ["SequenceVariant", "← what production emits,", "for every ClinVar record"],
            font_size=11, fill="#2563eb",
        )
    )
    parts.append(
        _diagram_node(
            380, 185, 260, 55,
            ["Snv (Single Nucleotide Variant)", "direct subclass —", "never used by this ingest"],
            font_size=11, fill="#94a3b8",
        )
    )

    parts.append(_diagram_arrow(125, 54, 125, 90))
    parts.append(_diagram_arrow(125, 145, 125, 185))
    parts.append(_diagram_arrow(280, 212, 380, 212))
    parts.append(_diagram_edge_label(288, 205, ["is_a"], font_size=9.5))

    parts.append(
        _cartoon_caption(
            "Biolink class hierarchy (biolink-model docs) — node category stays SequenceVariant either way",
            width=width, y=height - 8,
        )
    )
    parts.append("</svg>")
    return "".join(parts)


def diagram_current_model():
    width, height = 720, 340
    parts = [_cartoon_svg(width, height), ARROW_MARKER_DEFS]

    parts.append(_diagram_node(30, 130, 190, 80, ["SequenceVariant", "type: SO consequence terms", "(e.g. missense_variant)"]))
    parts.append(_diagram_node(470, 25, 210, 55, ["Gene", "(NCBIGene) &mdash; ×N, identical", "predicate every time"], font_size=11.5))
    parts.append(_diagram_node(470, 145, 210, 50, ["Disease", "(MONDO)"]))
    parts.append(_diagram_node(470, 260, 210, 50, ["PhenotypicFeature", "(HPO)"]))

    parts.append(_diagram_arrow(220, 150, 470, 55))
    parts.append(_diagram_edge_label(230, 100, ["is_sequence_variant_of"]))
    parts.append(_diagram_arrow(220, 165, 470, 168))
    parts.append(_diagram_edge_label(230, 195, ["causes /", "associated_with_increased_likelihood_of"]))
    parts.append(_diagram_arrow(220, 190, 470, 280))
    parts.append(_diagram_edge_label(230, 230, ["contributes_to"]))

    parts.append(
        f"<rect x='30' y='250' width='190' height='55' rx='6' fill='#fee2e2' stroke='#dc2626' stroke-width='1.5' stroke-dasharray='4,3'/>"
        f"<text x='125' y='270' font-size='11' font-weight='700' text-anchor='middle' fill='#b91c1c'>Not representable:</text>"
        f"<text x='125' y='284' font-size='10' text-anchor='middle' fill='#b91c1c'>CNV, translocation,</text>"
        f"<text x='125' y='297' font-size='10' text-anchor='middle' fill='#b91c1c'>inversion, complex SV</text>"
    )
    parts.append(
        _cartoon_caption(
            "Today: one node type, one undifferentiated gene edge, SNV/indel/microsatellite only", width=width, y=height - 8
        )
    )
    parts.append("</svg>")
    return "".join(parts)


def diagram_proposed_model():
    width, height = 780, 500
    parts = [_cartoon_svg(width, height), ARROW_MARKER_DEFS]

    parts.append(
        _diagram_node(
            20,
            190,
            180,
            90,
            ["SequenceVariant", "type: VARIANT CLASS (SO)", "SNV / indel / STR /", "CNV-gain / CNV-loss /", "inversion / complex"],
            font_size=11,
        )
    )

    gene_defs = [
        {
            "y": 30,
            "label_y": 28,
            "label": ["Gene A", "(fully contained)"],
            "edge": ["impact: whole_gene", "consequence: frameshift_variant", "morph: amorph (activity: decreased)"],
            "color": CARTOON_DEL_STROKE,
        },
        {
            "y": 175,
            "label_y": 108,
            "label": ["Gene B", "(breakpoint inside)"],
            "edge": ["impact: partial_gene (exons 4-6)", "consequence: exon_loss_variant", "morph: unknown / TBD"],
            "color": CARTOON_SUB_STROKE,
        },
        {
            "y": 320,
            "label_y": 270,
            "label": ["Gene C", "(no exonic overlap)"],
            "edge": ["impact: regulatory_only", "consequence: n/a (no CDS overlap)", "morph: n/a"],
            "color": "#64748b",
        },
    ]
    for g in gene_defs:
        parts.append(_diagram_node(430, g["y"], 190, 65, g["label"], font_size=12))
        parts.append(_diagram_arrow(200, 235, 430, g["y"] + 32, color=g["color"], marker=f"url(#arrow-diag)"))
        parts.append(_diagram_edge_label(210, g["label_y"], g["edge"], color=g["color"]))

    parts.append(_diagram_node(20, 340, 150, 55, ["Disease", "(MONDO)"], font_size=12))
    parts.append(_diagram_node(195, 340, 150, 55, ["PhenotypicFeature", "(HPO)"], font_size=12))
    parts.append(_diagram_arrow(70, 280, 85, 340, color=DIAGRAM_EDGE_COLOR))
    parts.append(_diagram_arrow(150, 280, 250, 340, color=DIAGRAM_EDGE_COLOR))
    parts.append(_diagram_edge_label(20, 305, ["causes / associated_with"], anchor="start", font_size=9.5))
    parts.append(_diagram_edge_label(230, 305, ["contributes_to"], anchor="start", font_size=9.5))

    parts.append(
        _cartoon_caption(
            "Proposed: variant class stays on the node; consequence + morph + impact move onto each per-gene edge",
            width=width,
            y=height - 8,
        )
    )
    parts.append("</svg>")
    return "".join(parts)


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
    variant_genes: dict,
    gene_symbol: str,
    disease_filter: set | None = None,
) -> list[dict]:
    """Every Pathogenic/Likely-pathogenic SNV in gene_symbol with >=1 resolved
    MONDO disease (star_min=0 -- any submission evidence, matching section {S.crossfilter}'s
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
            _, gene_syms = variant_genes_for(variant_genes, row["ID"])
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
# Section 6: review-star cutoff + multi-submitter concordance rescue
# ---------------------------------------------------------------------------


def load_maps(data_dir: Path):
    var_records = make_variant_record_map(str(data_dir / "submission_summary.txt.gz"))
    map_to_mondo = make_mondo_map(str(data_dir / "mondo.sssom.tsv"))
    medgen_to_mondo = make_medgen_to_mondo_map(str(data_dir / "MedGenIDMappings.txt.gz"))
    map_to_mondo.update(medgen_to_mondo)
    return var_records, map_to_mondo


def sssom_label_collisions(data_dir: Path) -> dict:
    """Source terms that more than one MONDO class claims.

    If MONDO:A and MONDO:B both assert skos:exactMatch to source terms carrying the same
    name, transitivity says they should be one class. This is a far better duplicate
    *detector* than the subclass graph -- two terms sharing a parent are usually just
    clinically related, whereas two terms claiming the same external identity are a
    genuine inconsistency somewhere.

    It is a detector and nothing more. Every mapping in the file carries
    semapv:UnspecifiedMatching, so nothing records how any of them was made and there is
    no provenance to adjudicate a collision with. See recommendation 2: report these for
    curation, never merge on them. The SCN1A case is exactly why -- "dravet syndrome" is
    claimed by three MONDO classes, and reading mondo#745 shows the OMIM:607208 ->
    MONDO:0100079 mapping is a deliberate curatorial decision, not the bug it looks like.
    A collision means "a human should look at this", nothing stronger.
    """
    path = data_dir / "mondo.sssom.tsv"
    predicates: Counter = Counter()
    justifications: Counter = Counter()
    by_label: dict = {}
    total = 0
    with open(path) as f:
        header = None
        for line in f:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if header is None:
                header = {k: i for i, k in enumerate(cols)}
                continue
            if len(cols) <= max(header.values()):
                continue
            total += 1
            predicates[cols[header["predicate_id"]]] += 1
            justifications[cols[header["mapping_justification"]]] += 1
            label = cols[header["object_label"]].strip().lower()
            if label:
                by_label.setdefault(label, set()).add(cols[header["subject_id"]])

    collisions = {lab: sorted(t) for lab, t in by_label.items() if len(t) > 1}
    exact = predicates.get("skos:exactMatch", 0)
    unspec = justifications.get("semapv:UnspecifiedMatching", 0)
    return {
        "total": total,
        "n_labels": len(by_label),
        "n_collisions": len(collisions),
        "exact_pct": 100 * exact / (total or 1),
        "unspecified_pct": 100 * unspec / (total or 1),
        "dravet": collisions.get("dravet syndrome", []),
    }


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


def summarize_review_status(var_records: dict) -> list[dict]:
    """Real distribution of submission_summary.txt's per-record ReviewStatus
    field, across every record for every variant (var_records is the full,
    unfiltered file contents -- see make_variant_record_map()). Includes every
    status review_star_map knows how to score, even ones with zero observed
    records, so it's visible at a glance which of the ten defined statuses
    (notably "criteria provided, multiple submitters, no conflicts", the 2-star
    tier) never actually appear per-record in this file -- see 'How ClinVar
    curation works' above and section {S.star_cutoff}'s own docstring."""
    counts: Counter = Counter()
    for records in var_records.values():
        for rec in records:
            counts[rec["ReviewStatus"]] += 1

    known_statuses = {k.replace("_", " "): v for k, v in review_star_map.items() if k != "."}
    all_statuses = known_statuses.keys() | counts.keys()
    rows = [
        {
            "status": status,
            "count": counts.get(status, 0),
            "stars": known_statuses.get(status),
        }
        for status in all_statuses
    ]
    rows.sort(key=lambda r: (-(r["stars"] if r["stars"] is not None else -1), -r["count"]))
    return rows


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


def compute_star_data(clinvar_tsv: Path, var_records: dict, map_to_mondo: dict, mondo_labels: dict, variant_genes: dict):
    """Single pass over clinvar.tsv computing summary counts (all star levels,
    including "2-star, computed") and a three-way mutually-exclusive partition
    of every gene-disease pair by evidence tier (pairs_ge3 / pairs_2star_concordance
    / pairs_remaining, section {S.star_cutoff}). Also builds pair_submitters (every pair's full
    submitter roster, any star level), used by the ClinGen-coverage analysis."""
    variant_sets = {s: set() for s in STAR_LEVELS}
    pair_sets = {s: set() for s in STAR_LEVELS}
    pair_variant_ids = {s: {} for s in [0] + ENUMERATED_LEVELS}
    # "2-star, computed": since no individual submission record ever actually carries the aggregate
    # 2-star ReviewStatus (see 'How ClinVar curation works' above), a raw star_min=2 filter over
    # per-record data is identical to star_min=3 -- this instead reconstructs what a true >=2-star
    # population would be by rescuing pairs via >=MIN_CONCORDANT_SUBMITTERS concordant submitters,
    # the same proxy mechanism section the star-cutoff section uses. Not a raw filter; computed.
    variant_set_computed2: set = set()
    pair_set_computed2: set = set()
    pair_variant_ids_computed2: dict = {}

    # (gene_sym, gene_id, mondo_id) -> {"submitters": set, "variants": set} -- every P/LP-family gene-disease
    # pair's full submitter roster, any star level, regardless of whether it needs concordance rescue. Used
    # for the ClinGen-coverage analysis: which pairs' only submitters are ClinGen-affiliated.
    pair_submitters: dict = {}

    with open(clinvar_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            varid = row["ID"]
            if varid not in var_records:
                continue

            gene_ids, gene_symbols = variant_genes_for(variant_genes, varid)
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

            dis_computed2, _, _ = variant_records_to_disease(
                records, map_to_mondo, star_min=EXPERT_STAR_MIN, rescue_min_submitters=MIN_CONCORDANT_SUBMITTERS
            )
            if dis_computed2:
                variant_set_computed2.add(varid)
                for gene_id, gene_sym in zip(gene_ids, gene_symbols):
                    for disease_id in dis_computed2:
                        pair_set_computed2.add((gene_id, disease_id))
                        key = (gene_sym, gene_id, disease_id)
                        pair_variant_ids_computed2.setdefault(key, set()).add(varid)

            groups = concordance_groups(records, map_to_mondo) if gene_ids else {}
            for (mondo_id, _clinsig), entry in groups.items():
                for gene_id, gene_sym in zip(gene_ids, gene_symbols):
                    key = (gene_sym, gene_id, mondo_id)
                    pe = pair_submitters.setdefault(key, {"submitters": set(), "variants": set()})
                    pe["submitters"].update(entry["submitters"])
                    pe["variants"].add(varid)

    counts = {s: {"variants": len(variant_sets[s]), "gene_disease_pairs": len(pair_sets[s])} for s in STAR_LEVELS}
    counts["2c"] = {"variants": len(variant_set_computed2), "gene_disease_pairs": len(pair_set_computed2)}

    # Three-way, mutually-exclusive partition of every gene-disease pair with >=1
    # Pathogenic/Likely-pathogenic submission record (pair_sets[0]), by evidence tier.
    #
    # Group A: raw >=3-star evidence exists for this pair (2 and 3 are identical raw
    # thresholds, see docstring above) -- the strongest, unambiguous tier.
    # Pooled submitter roster per canonical (gene_id, mondo_id) pair -- the union across ALL of
    # that pair's variants, collapsing the handful of gene ids carrying >1 symbol spelling. This
    # is pair-level support: how many independent submitters have said anything P/LP about this
    # gene-disease relationship at all. The tier assignment above is a per-VARIANT test (does
    # some single variant clear >=3 stars or draw >=MIN_CONCORDANT_SUBMITTERS submitters by
    # itself), so the two can diverge sharply -- see the tier-3 blurb for a worked example.
    pooled_submitters: dict = {}
    for (_gene_sym, gene_id, mondo_id), entry in pair_submitters.items():
        pooled_submitters.setdefault((gene_id, mondo_id), set()).update(entry["submitters"])

    def n_submitters(gene_id: str, mondo_id: str) -> int:
        return len(pooled_submitters.get((gene_id, mondo_id), ()))

    # "# variants" means something different in each tier list below (tier 1 counts only its
    # >=3-star variants, tier 2 only its concordance-rescued ones, tier 3 all of them), which
    # makes those columns incomparable across tiers. This is the tier-independent measure --
    # every variant with >=1 P/LP record mapping to the pair -- used wherever pairs from
    # different tiers appear side by side (section the monarch-kg section).
    all_pair_variants: dict = {}
    for (_gene_sym, gene_id, mondo_id), varids in pair_variant_ids[0].items():
        all_pair_variants.setdefault((gene_id, mondo_id), set()).update(varids)

    def n_variants_all(gene_id: str, mondo_id: str) -> int:
        return len(all_pair_variants.get((gene_id, mondo_id), ()))

    base_level = ENUMERATED_LEVELS[0]
    pairs_ge3 = []
    for key in sorted(pair_variant_ids[base_level]):
        gene_sym, gene_id, mondo_id = key
        entry = {
            "gene": gene_sym,
            "gene_id": gene_id,
            "mondo": mondo_id,
            "disease_name": mondo_labels.get(mondo_id, ""),
            "n_submitters": n_submitters(gene_id, mondo_id),
            "n_variants_all": n_variants_all(gene_id, mondo_id),
        }
        for level in ENUMERATED_LEVELS:
            varids = pair_variant_ids[level].get(key)
            entry[f"v{level}"] = len(varids) if varids else None
        pairs_ge3.append(entry)

    # Group B: no raw >=3-star evidence anywhere for this pair, but reachable via
    # >=MIN_CONCORDANT_SUBMITTERS concordant submitters ("2-star, computed"). Mutually
    # exclusive of Group A by construction: pair_set_computed2's raw-star check is the
    # same star_min=3 threshold as pair_sets[3], so excluding pair_sets[3] here removes
    # exactly (and only) the pairs Group A already covers.
    pairs_2star_concordance = []
    for (gene_sym, gene_id, mondo_id), varids in pair_variant_ids_computed2.items():
        if (gene_id, mondo_id) in pair_sets[3]:
            continue
        pairs_2star_concordance.append(
            {
                "gene": gene_sym,
                "gene_id": gene_id,
                "mondo": mondo_id,
                "disease_name": mondo_labels.get(mondo_id, ""),
                "n_variants": len(varids),
                "n_submitters": n_submitters(gene_id, mondo_id),
                "n_variants_all": n_variants_all(gene_id, mondo_id),
            }
        )

    # Group C: everything else -- >=1 P/LP submission record, but neither raw >=3-star
    # evidence nor concordance rescue. Mutually exclusive of A and B by construction:
    # excluding pair_set_computed2 (which is exactly A union B) leaves only the remainder.
    pairs_remaining = []
    for (gene_sym, gene_id, mondo_id), varids in pair_variant_ids[0].items():
        if (gene_id, mondo_id) in pair_set_computed2:
            continue
        pairs_remaining.append(
            {
                "gene": gene_sym,
                "gene_id": gene_id,
                "mondo": mondo_id,
                "disease_name": mondo_labels.get(mondo_id, ""),
                "n_variants": len(varids),
                "n_submitters": n_submitters(gene_id, mondo_id),
                "n_variants_all": n_variants_all(gene_id, mondo_id),
            }
        )

    # Per-term variant sets for the SCN1A hierarchy diagram: which variants each of its
    # disease terms claims, so the overlap between them can be drawn rather than asserted.
    scn1a_sets: dict = {}
    for (_gene_sym, gene_id, mondo_id), varids in pair_variant_ids[0].items():
        if gene_id == "NCBIGene:6323":
            scn1a_sets.setdefault(mondo_id, set()).update(varids)

    return (
        counts,
        pairs_ge3,
        pairs_2star_concordance,
        pairs_remaining,
        pair_submitters,
        pair_sets[0],
        scn1a_sets,
    )


CLINGEN_SUBMITTER_RE = re.compile(r"clingen", re.IGNORECASE)


def summarize_clingen_coverage(pair_submitters: dict, mondo_labels: dict) -> tuple[dict, list]:
    """Classify every gene-disease pair's full submitter roster (pair_submitters
    from compute_star_data -- any Pathogenic/Likely-pathogenic evidence, any
    star) by whether it includes a ClinGen-affiliated submitter (a Variant
    Curation Expert Panel, the ACMG Interpretation Working Group, etc.,
    detected as a "clingen" substring in the Submitter field -- the only
    practical proxy available, since submission_summary.txt carries no
    separate submitter-type/affiliation column) and whether any OTHER,
    non-ClinGen submitter also independently reports the same pair.

    single_submitter_clingen_pairs answers "of pairs with exactly one
    submitter, how many of those lone submitters are ClinGen?" -- these can
    never be concordance-rescued (rescue needs >=2 submitters), so today they
    only ever enter production if that single submission itself reaches
    3+ stars, i.e. IS a ClinGen expert-panel review.

    clingen_only_rows is the broader question: pairs where ClinGen is present
    and ALL of that pair's submitters are ClinGen-affiliated (covers e.g. two
    different VCEPs both weighing in on the same pair, not just literal
    single-submitter cases) -- these would have zero P/LP evidence at all,
    from anyone, if ClinGen's own submissions were removed.
    """
    single_submitter_pairs = 0
    single_submitter_clingen_pairs = 0
    clingen_only_rows = []

    for (gene_sym, gene_id, mondo_id), entry in pair_submitters.items():
        submitters = entry["submitters"]
        clingen_submitters = {s for s in submitters if CLINGEN_SUBMITTER_RE.search(s)}
        other_submitters = submitters - clingen_submitters

        if len(submitters) == 1:
            single_submitter_pairs += 1
            if clingen_submitters:
                single_submitter_clingen_pairs += 1

        if clingen_submitters and not other_submitters:
            clingen_only_rows.append(
                {
                    "gene": gene_sym,
                    "gene_id": gene_id,
                    "mondo": mondo_id,
                    "disease_name": mondo_labels.get(mondo_id, ""),
                    "clingen_submitters": sorted(clingen_submitters),
                    "n_clingen_submitters": len(clingen_submitters),
                    "n_variants": len(entry["variants"]),
                    "variant_sample": sorted(entry["variants"], key=int)[:MAX_VARIANT_SAMPLE],
                }
            )

    clingen_only_rows.sort(key=lambda r: (-r["n_variants"], r["gene"]))
    summary = {
        "total_pairs": len(pair_submitters),
        "single_submitter_pairs": single_submitter_pairs,
        "single_submitter_clingen_pairs": single_submitter_clingen_pairs,
        "clingen_only_pairs": len(clingen_only_rows),
    }
    return summary, clingen_only_rows


# ---------------------------------------------------------------------------
# Section 9: CLNSIG / CLNREVSTAT / CLNVC crossfilter
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


# ---------------------------------------------------------------------------
# Section 8: Monarch KG gene-disease associations vs ClinVar
# ---------------------------------------------------------------------------

MONARCH_KG_URL = "https://data.monarchinitiative.org/monarch-kg/latest/monarch-kg.tar.gz"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
# Monarch KG gene nodes are HGNC-keyed; ClinVar's GENEINFO is Entrez. HGNC's own
# complete set is the authoritative crosswalk (the KG's gene node xrefs carry
# ENSEMBL/OMIM but not NCBIGene).
MONARCH_DERIVED = ("mk_gene_disease.tsv", "mk_mondo_subclass.tsv", "mk_mondo_terms.tsv")


def ensure_monarch_kg(data_dir: Path) -> None:
    """Fetch the Monarch KG release and reduce it to three small derived files.

    monarch-kg.tar.gz is ~330MB and expands to a 5.3GB edges TSV, so everything
    downstream reads the derived extracts instead:
      mk_gene_disease.tsv   subject, object, primary_knowledge_source, original_subject
                            for every *GeneToDisease* edge (~16k rows)
      mk_mondo_subclass.tsv child, parent for every MONDO subclass_of edge (~46k)
      mk_mondo_terms.tsv    id, name, deprecated for every MONDO node
    """
    if all((data_dir / f).exists() for f in MONARCH_DERIVED):
        return

    tar_path = data_dir / "monarch-kg.tar.gz"
    if not tar_path.exists():
        print(f"Downloading {MONARCH_KG_URL} (~330MB, one-time)...")
        with requests.get(MONARCH_KG_URL, stream=True, timeout=1800) as resp:
            resp.raise_for_status()
            with open(tar_path, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)

    edges = data_dir / "monarch-kg_edges.tsv"
    nodes = data_dir / "monarch-kg_nodes.tsv"
    if not (edges.exists() and nodes.exists()):
        import tarfile

        print("Extracting monarch-kg_{nodes,edges}.tsv (~5.7GB)...")
        with tarfile.open(tar_path, "r:gz") as tf:
            tf.extractall(data_dir, filter="data")

    print("Reducing Monarch KG to derived extracts...")
    with (
        open(edges, newline="") as fh,
        open(data_dir / "mk_gene_disease.tsv", "w") as gd,
        open(data_dir / "mk_mondo_subclass.tsv", "w") as sc,
    ):
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if "GeneToDisease" in row["category"]:
                gd.write(f"{row['subject']}\t{row['object']}\t{row['primary_knowledge_source']}\t{row['original_subject']}\n")
            elif (
                row["predicate"] == "biolink:subclass_of"
                and row["subject"].startswith("MONDO:")
                and row["object"].startswith("MONDO:")
            ):
                sc.write(f"{row['subject']}\t{row['object']}\n")

    with open(nodes, newline="") as fh, open(data_dir / "mk_mondo_terms.tsv", "w") as out:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["id"].startswith("MONDO:"):
                out.write(f"{row['id']}\t{row['name']}\t{row['deprecated']}\n")


def ensure_hgnc(data_dir: Path) -> Path:
    path = data_dir / "hgnc_complete_set.txt"
    if path.exists():
        return path
    print(f"Downloading {HGNC_URL} (one-time)...")
    resp = requests.get(HGNC_URL, timeout=300)
    resp.raise_for_status()
    path.write_bytes(resp.content)
    return path


def load_monarch_gene_disease(data_dir: Path) -> dict:
    """Monarch KG's curated gene-disease associations, keyed the same way this
    report keys ClinVar pairs: (NCBIGene id, MONDO id).

    Returns pairs, per-source pair sets, the MONDO parent map (for the
    ancestor-aware reconciliation below), MONDO labels and the deprecated set.
    """
    ensure_monarch_kg(data_dir)
    hgnc_path = ensure_hgnc(data_dir)

    hgnc_to_entrez = {}
    entrez_to_symbol = {}
    with open(hgnc_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("entrez_id"):
                gene_id = "NCBIGene:{}".format(row["entrez_id"])
                hgnc_to_entrez[row["hgnc_id"]] = gene_id
                entrez_to_symbol[gene_id] = row["symbol"]

    pairs: set = set()
    by_source: dict = {}
    pair_sources: dict = {}
    with open(data_dir / "mk_gene_disease.tsv") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            subject, obj, source = parts[0], parts[1], parts[2]
            original = parts[3] if len(parts) > 3 else ""
            if not obj.startswith("MONDO:"):
                continue
            # OMIM/ClinGen edges carry an HGNC subject; Orphanet edges additionally
            # preserve the Entrez id in original_subject, which is used as a fallback.
            gene_id = hgnc_to_entrez.get(subject)
            if gene_id is None and original.startswith("NCBIGene:"):
                gene_id = original
            if gene_id is None:
                continue
            pairs.add((gene_id, obj))
            by_source.setdefault(source, set()).add((gene_id, obj))
            pair_sources.setdefault((gene_id, obj), set()).add(source.replace("infores:", ""))

    parents: dict = {}
    with open(data_dir / "mk_mondo_subclass.tsv") as fh:
        for line in fh:
            child, parent = line.rstrip("\n").split("\t")
            parents.setdefault(child, set()).add(parent)

    labels: dict = {}
    deprecated: set = set()
    with open(data_dir / "mk_mondo_terms.tsv") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            labels[parts[0]] = parts[1] if len(parts) > 1 else ""
            if len(parts) > 2 and parts[2].strip().lower() in ("true", "1"):
                deprecated.add(parts[0])

    return {
        "pairs": pairs,
        "by_source": by_source,
        "pair_sources": pair_sources,
        "entrez_to_symbol": entrez_to_symbol,
        "parents": parents,
        "labels": labels,
        "deprecated": deprecated,
        "known_terms": set(labels),
    }


SUPPORT_BUCKETS = ["1", "2", "3-4", "5-9", "10+"]

# ClinVar's GENEINFO lists every locus overlapping a variant's position, not just the
# causal gene, so a heavily-submitted CFTR variant also tags CFTR-AS1 and CFTR-AS2 and
# the ingest minted a gene-disease pair for each (now fixed -- see section the ingest-recommendation section). These two
# LOC<GeneID> placeholders; these two patterns catch the other systematic offenders --
# antisense/divergent transcripts overlapping a real disease gene, and uncharacterized
# open reading frames. Flagged rather than dropped here, so the report can show how much
# of the "novel" ClinVar-only signal is positional overlap rather than new biology.
ANTISENSE_RE = re.compile(r"-(AS\d*|DT|IT\d*|OT\d*)$")
UNCHAR_ORF_RE = re.compile(r"^C\d+orf\d+$", re.IGNORECASE)


def classify_gene_symbol(symbol: str) -> str:
    if ANTISENSE_RE.search(symbol):
        return "antisense/divergent transcript"
    if UNCHAR_ORF_RE.match(symbol):
        return "uncharacterized ORF"
    return "protein-coding"


def support_bucket(n: int) -> str:
    if n >= 10:
        return "10+"
    if n >= 5:
        return "5-9"
    if n >= 3:
        return "3-4"
    if n == 2:
        return "2"
    return "1"


def build_monarch_comparison(
    pairs_ge3: list, pairs_2star: list, pairs_remaining: list, monarch: dict,
    pair_submitters: dict | None = None,
) -> dict:
    """Reconcile ClinVar's implied gene-disease pairs against the Monarch KG's
    curated ones.

    Note the two are NOT the same kind of statement. Monarch's gene-disease
    edges come from OMIM / Orphanet / ClinGen curation of the *gene*; ClinVar's
    pairs here are derived (variant -> gene) x (variant -> disease) from
    individual variant submissions. ClinVar is not itself a source of
    gene-disease edges in the KG, so this measures how much independent
    submitted variant evidence stands behind each curated association, and what
    ClinVar implies that no curator has asserted.

    Matching is on exact (NCBIGene id, MONDO id). Because MONDO carries several
    co-existing terms over one clinical area (see the SCN1A worked example),
    ClinVar-only pairs are additionally checked for an *ancestor-aware* match:
    does Monarch associate the same gene with an ancestor or descendant of
    ClinVar's disease term?

    That check produces a *candidate* list, not a duplicate count, and the
    difference matters. Hierarchical proximity means two terms are clinically
    related; it does not make them redundant. In the SCN1A example both nearby
    pairs turn out to be distinctions MONDO maintains deliberately -- one at
    ClinGen's explicit request -- so this column must never be used to collapse
    rows automatically.
    """
    clinvar: dict = {}
    for tier, rows in ((1, pairs_ge3), (2, pairs_2star), (3, pairs_remaining)):
        for r in rows:
            key = (r["gene_id"], r["mondo"])
            n_var = r.get("n_variants_all", 0)
            prev = clinvar.get(key)
            # a gene id carrying >1 symbol spelling can yield duplicate rows; keep the
            # strongest tier and the larger variant count
            if prev is None or tier < prev["tier"]:
                clinvar[key] = {
                    "gene": r["gene"],
                    "gene_id": r["gene_id"],
                    "mondo": r["mondo"],
                    "disease_name": r["disease_name"],
                    "n_variants": n_var,
                    "n_submitters": r["n_submitters"],
                    "tier": tier,
                }
            elif n_var > prev["n_variants"]:
                prev["n_variants"] = n_var

    mk_pairs = monarch["pairs"]
    parents = monarch["parents"]
    labels = monarch["labels"]

    clinvar_keys = set(clinvar)
    both = clinvar_keys & mk_pairs
    clinvar_only = clinvar_keys - mk_pairs
    monarch_only = mk_pairs - clinvar_keys

    # Monarch diseases per gene, for the ancestor-aware reconciliation
    mk_by_gene: dict = {}
    for gene_id, mondo_id in mk_pairs:
        mk_by_gene.setdefault(gene_id, set()).add(mondo_id)

    # How well is each Monarch association supported by ClinVar submissions?
    supported_rows = []
    support_xtab: Counter = Counter()
    for key in both:
        c = clinvar[key]
        support_xtab[(c["tier"], support_bucket(c["n_submitters"]))] += 1
        supported_rows.append(
            {
                "gene": c["gene"],
                "gene_id": c["gene_id"],
                "mondo": c["mondo"],
                "disease_name": c["disease_name"] or labels.get(c["mondo"], ""),
                "n_variants": c["n_variants"],
                "n_submitters": c["n_submitters"],
                "tier": c["tier"],
            }
        )

    # In ClinVar, absent from Monarch -- split by whether an ancestor/descendant of the
    # ClinVar term is already associated with that same gene in Monarch (i.e. the same
    # relationship recorded at a different level of the hierarchy) or whether the gene is
    # absent from Monarch's gene-disease set entirely.
    anc_cache: dict = {}

    def ancestors(term: str) -> set:
        if term not in anc_cache:
            anc_cache[term] = mondo_ancestors(term, parents)
        return anc_cache[term]

    only_rows = []
    only_xtab: Counter = Counter()
    kinship_counts: Counter = Counter()
    for key in clinvar_only:
        c = clinvar[key]
        gene_id, mondo_id = key
        mk_diseases = mk_by_gene.get(gene_id, set())
        if not mk_diseases:
            kinship = "gene absent from Monarch"
        else:
            if ancestors(mondo_id) & mk_diseases:
                kinship = "ancestor in Monarch"
            elif any(mondo_id in ancestors(d) for d in mk_diseases):
                kinship = "descendant in Monarch"
            else:
                kinship = "gene in Monarch, unrelated term"
        kinship_counts[kinship] += 1
        only_xtab[(c["tier"], support_bucket(c["n_submitters"]))] += 1
        only_rows.append(
            {
                "gene": c["gene"],
                "gene_id": c["gene_id"],
                "mondo": c["mondo"],
                "disease_name": c["disease_name"] or labels.get(c["mondo"], ""),
                "n_variants": c["n_variants"],
                "n_submitters": c["n_submitters"],
                "tier": c["tier"],
                "kinship": kinship,
            }
        )

    # Monarch associations with no ClinVar P/LP evidence at all.
    #
    # mondo_status answers "is this term actually in MONDO?" -- three outcomes, because the
    # KG can carry an association whose disease term MONDO has since retired (Orphanet still
    # asserts SCN1A -> MONDO:0011794 "obsolete Dravet syndrome") or one whose id has no MONDO
    # node in this release at all. Exposed as a filter so the genuinely-current associations
    # can be separated from ones the ontology has moved on from.
    pair_sources = monarch["pair_sources"]
    symbols = monarch["entrez_to_symbol"]
    known_terms = monarch["known_terms"]
    monarch_only_rows = []
    mondo_status_counts: Counter = Counter()
    for gene_id, mondo_id in monarch_only:
        if mondo_id in monarch["deprecated"]:
            status = "deprecated in MONDO"
        elif mondo_id in known_terms:
            status = "current in MONDO"
        else:
            status = "not in MONDO"
        mondo_status_counts[status] += 1
        monarch_only_rows.append(
            {
                "gene": symbols.get(gene_id, ""),
                "gene_id": gene_id,
                "mondo": mondo_id,
                "disease_name": labels.get(mondo_id, ""),
                "sources": ", ".join(sorted(pair_sources.get((gene_id, mondo_id), ()))),
                "mondo_status": status,
            }
        )

    deprecated_clinvar = sorted(
        (clinvar[k] for k in clinvar_keys if k[1] in monarch["deprecated"]),
        key=lambda r: -r["n_variants"],
    )

    # SCN1A worked example -- computed live rather than hardcoded, so the numbers in the
    # narrative can never drift from the data the rest of the section is built from.
    SCN1A_GENE = "NCBIGene:6323"
    # the single lab behind a term, where there is only one -- used to label the diagram
    scn1a_leads = {}
    for (gene_sym, gene_id, mondo_id), entry in pair_submitters.items():
        if gene_id == SCN1A_GENE and len(entry["submitters"]) == 1:
            scn1a_leads[mondo_id] = next(iter(entry["submitters"]))

    scn1a_rows = sorted(
        (
            {
                "mondo": mondo_id,
                "disease_name": c["disease_name"] or labels.get(mondo_id, ""),
                "n_variants": c["n_variants"],
                "n_submitters": c["n_submitters"],
                "tier": c["tier"],
                "in_monarch": (gene_id, mondo_id) in mk_pairs,
                "sources": ", ".join(sorted(pair_sources.get((gene_id, mondo_id), ()))),
            }
            for (gene_id, mondo_id), c in clinvar.items()
            if gene_id == SCN1A_GENE
        ),
        key=lambda r: -r["n_variants"],
    )


    # Multi-variant support per pair -- the fourth kind of "2 star" in the two-star
    # section: several independent variants tying a gene to a disease, which is evidence
    # ClinVar's own aggregate never considers because it is computed per variant.
    multi_variant_pairs = 0
    multi_variant_single_submitter = 0
    for c in clinvar.values():
        if c["n_variants"] >= 2:
            multi_variant_pairs += 1
            if c["n_submitters"] <= 1:
                multi_variant_single_submitter += 1
    # Pair-level support profile over EVERY ClinVar pair (not just the Monarch-matched
    # ones) -- the evidence behind section the ingest-recommendation section's threshold recommendation. tier 1/2 are what
    # production ingests today; the pooled submitter count is the pair-level alternative.
    all_xtab: Counter = Counter()
    for c in clinvar.values():
        all_xtab[(c["tier"], support_bucket(c["n_submitters"]))] += 1

    ingested_today = sum(n for (tier, _b), n in all_xtab.items() if tier in (1, 2))
    thin_ingested = all_xtab[(1, "1")] + all_xtab[(2, "1")]
    dropped_ge3 = sum(all_xtab[(3, b)] for b in ("3-4", "5-9", "10+"))
    dropped_ge2 = dropped_ge3 + all_xtab[(3, "2")]
    # proposed rule: any raw >=3-star evidence (tier 1) OR >=N pooled distinct submitters
    projected = {}
    for n_min in (2, 3):
        projected[n_min] = sum(
            1
            for c in clinvar.values()
            if c["tier"] == 1 or c["n_submitters"] >= n_min
        )

    return {
        "multi_variant_pairs": multi_variant_pairs,
        "multi_variant_single_submitter": multi_variant_single_submitter,
        "all_xtab": all_xtab,
        "ingested_today": ingested_today,
        "thin_ingested": thin_ingested,
        "dropped_ge3": dropped_ge3,
        "dropped_ge2": dropped_ge2,
        "projected": projected,
        "n_clinvar": len(clinvar_keys),
        "n_monarch": len(mk_pairs),
        "n_both": len(both),
        "n_clinvar_only": len(clinvar_only),
        "n_monarch_only": len(monarch_only),
        "by_source": {k: len(v) for k, v in sorted(monarch["by_source"].items())},
        "support_xtab": support_xtab,
        "only_xtab": only_xtab,
        "kinship_counts": kinship_counts,
        "mondo_status_counts": mondo_status_counts,
        "scn1a_rows": scn1a_rows,
        "scn1a_leads": scn1a_leads,
        "supported_rows": sorted(supported_rows, key=lambda r: -r["n_submitters"]),
        "only_rows": sorted(only_rows, key=lambda r: -r["n_submitters"]),
        "monarch_only_rows": sorted(monarch_only_rows, key=lambda r: (r["gene"], r["mondo"])),
        "deprecated_rows": deprecated_clinvar,
    }


def summarize_submission_file(var_records: dict) -> dict:
    """Field-level profile of submission_summary.txt.gz, computed from the already
    in-memory var_records so it costs no extra I/O. Feeds section {S.input_files}'s input-file
    breakdown -- in particular how much of ReportedPhenotypeInfo is the
    "C3661900:not provided" placeholder, which maps to no MONDO term and so makes
    those records invisible to every disease-mapping and concordance test."""
    n_records = 0
    placeholder = 0
    submitters: set = set()
    methods: Counter = Counter()
    for records in var_records.values():
        for rec in records:
            n_records += 1
            reported = rec["ReportedPhenotypeInfo"]
            if "C3661900" in reported or "not provided" in reported.lower():
                placeholder += 1
            submitters.add(rec["Submitter"])
            methods[rec["CollectionMethod"]] += 1
    return {
        "n_records": n_records,
        "n_variants": len(var_records),
        "n_submitters": len(submitters),
        "placeholder_phenotype": placeholder,
        "methods": methods,
    }


def compare_gene_attribution(clinvar_tsv: Path, data_dir: Path) -> dict:
    """GENEINFO's every-overlapping-locus list vs variant_summary's curated GeneID.

    clinvar.vcf.gz's GENEINFO is documented as "Gene(s) for the variant" and is
    populated positionally -- every gene whose span covers the variant, so
    antisense transcripts (CFTR-AS1), divergent transcripts (MFF-DT) and
    uncharacterized ORFs (C11orf65) ride along with the real gene and inherit
    its whole disease/submitter roster via make_genes_from_row(). ClinVar's own
    per-variant attribution lives in variant_summary.txt.gz (GeneID /
    GeneSymbol / HGNC_ID). This quantifies the swap: how many attributions each
    source makes, what share of each is a nomenclature-identifiable artifact,
    how often a multi-gene GENEINFO collapses to a single curated gene, and how
    much recall would be lost to GeneID=-1 or a missing row.
    """
    ensure_variant_summary_downloaded(data_dir)

    vs_symbols: dict = {}
    vs_sym_text: dict = {}
    vs_field: Counter = Counter()
    vcf_field: Counter = Counter()
    # one row per variant per build -- keep only the most-preferred build present for each
    # VariationID, so nothing is counted twice and GRCh37-only variants are not dropped
    vs_rank = {name: i for i, name in enumerate(ASSEMBLY_PREFERENCE)}
    vs_best: dict = {}
    with gzip.open(data_dir / "variant_summary.txt.gz", "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rank = vs_rank.get(row["Assembly"])
            if rank is None:
                continue
            vid = row["VariationID"]
            seen = vs_best.get(vid)
            if seen is not None and seen <= rank:
                continue
            if seen is not None:
                # a better build supersedes what a lower-preference row contributed
                vs_symbols.pop(vid, None)
                vs_sym_text.pop(vid, None)
            vs_best[vid] = rank
            vs_symbols.setdefault(vid, set()).add(row["GeneID"])
            vs_sym_text.setdefault(vid, set()).add(row["GeneSymbol"])


    # All four counted from the final per-variant selection rather than incrementally:
    # a variant present on both builds passes the rank test twice on its way to GRCh38,
    # so incrementing inside the loop double-counts it.
    vs_field["rows_kept"] = len(vs_best)
    for vid, rank in vs_best.items():
        vs_field[f"build_{ASSEMBLY_PREFERENCE[rank]}"] += 1
        if vs_symbols.get(vid) == {"-1"}:
            vs_field["geneid_minus1"] += 1
        # large CNVs degrade GeneSymbol to prose ("covers 42 genes, none of which...")
        if any(" " in sym or len(sym) > 40 for sym in vs_sym_text.get(vid, ())):
            vs_field["symbol_prose"] += 1

    stats = {
        "geneinfo_attributions": 0,
        "geneinfo_artifact": 0,
        "vs_attributions": 0,
        "vs_artifact": 0,
        "multi_gene_variants": 0,
        "multi_resolved_to_one": 0,
        "no_vs_row": 0,
        "vs_minus1": 0,
        "geneinfo_hist": Counter(),
        "vs_hist": Counter(),
        "vcf_field": vcf_field,
        "vs_field": vs_field,
    }

    with open(clinvar_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            vcf_field["rows"] += 1
            for col in ("CLNDISDB", "GENEINFO", "MC", "RS"):
                if row[col] == ".":
                    vcf_field[f"empty_{col}"] += 1
            gene_ids, gene_syms = make_genes_from_row(row["GENEINFO"])
            if not gene_ids:
                continue
            stats["geneinfo_hist"][min(len(gene_ids), 5)] += 1
            stats["geneinfo_attributions"] += len(gene_ids)
            stats["geneinfo_artifact"] += sum(1 for s in gene_syms if classify_gene_symbol(s) != "protein-coding")

            gids = vs_symbols.get(row["ID"])
            if gids is None:
                stats["no_vs_row"] += 1
                continue
            if gids == {"-1"}:
                stats["vs_minus1"] += 1
            syms = vs_sym_text.get(row["ID"], set())
            stats["vs_hist"][min(len(gids), 5)] += 1
            stats["vs_attributions"] += len(gids)
            stats["vs_artifact"] += sum(1 for s in syms if classify_gene_symbol(s) != "protein-coding")
            if len(gene_ids) > 1:
                stats["multi_gene_variants"] += 1
                if len(gids) == 1:
                    stats["multi_resolved_to_one"] += 1

    return stats


def mondo_ancestors(term: str, parents: dict, max_depth: int = 12) -> set:
    """Transitive subclass_of closure of a MONDO term (excluding the term itself)."""
    seen: set = set()
    frontier = {term}
    for _ in range(max_depth):
        nxt: set = set()
        for node in frontier:
            for parent in parents.get(node, ()):
                if parent not in seen:
                    seen.add(parent)
                    nxt.add(parent)
        if not nxt:
            break
        frontier = nxt
    return seen


# Section order for the rendered report. Numbers, the nav bar and every "section N"
# cross-reference are generated from this list -- renumbering used to mean editing ~50
# hardcoded places and was repeatedly got wrong. To reorder, move a line.
SECTION_ORDER = [
    ("purpose", "Purpose"),
    ("illustrative-examples", "Illustrative examples"),
    ("clinvar-curation", "How ClinVar curation works"),
    ("input-files", "Input files: what the ingest reads"),
    ("scope-decisions", "Scope: what this ingest keeps, decided up front"),
    ("phenotype-terms", "Phenotype (HPO) terms on ClinVar variants"),
    ("two-star", "The four kinds of 2-star evidence"),
    ("pairing", "Variant:disease vs gene:disease pairing"),
    ("star-cutoff", "Review-star cutoff impact on variant & gene-disease-pair counts"),
    ("evidence-tiers", "Evidence tiers and a predicate split by evidence"),
    ("monarch-kg", "Monarch KG gene-disease associations vs ClinVar"),
    ("ingest-recommendation", "New ClinVar ingest recommendation"),
    ("ingest-compare", "Previous vs new ingest"),
    ("crossfilter", "Multi-class variant clinical significance"),
    ("structural-variants", "Structural variants & CNVs \u2014 what's not in the VCF"),
    ("biolink-proposal", "Beyond SNVs: reframing the Biolink model for variants"),
]

# Short nav labels where the full title is too long for the bar
NAV_LABELS = {
    "illustrative-examples": "Illustrative examples",
    "clinvar-curation": "ClinVar curation",
    "input-files": "Input files",
    "scope-decisions": "Scope decisions",
    "phenotype-terms": "Phenotype terms",
    "two-star": "2-star evidence",
    "pairing": "Pairing",
    "star-cutoff": "Star cutoff & pair tiers",
    "evidence-tiers": "Evidence tiers",
    "monarch-kg": "Monarch KG vs ClinVar",
    "ingest-recommendation": "Ingest recommendation",
    "ingest-compare": "Previous vs new ingest",
    "crossfilter": "Multi-class significance",
    "structural-variants": "Structural variants",
    "biolink-proposal": "Biolink proposal",
}

SECTION_NUMBER = {anchor: i for i, (anchor, _t) in enumerate(SECTION_ORDER, 1)}


class _SectionRefs:
    """`S.monarch_kg` -> the number of the monarch-kg section. Underscores map to
    hyphens so it reads naturally inside an f-string."""

    def __getattr__(self, name):
        return SECTION_NUMBER[name.replace("_", "-")]


S = _SectionRefs()


def section_heading(anchor: str) -> str:
    title = dict(SECTION_ORDER)[anchor]
    return f'<h2 id="{anchor}">{SECTION_NUMBER[anchor]}. {title}</h2>'


def section_nav() -> str:
    links = "".join(
        f'  <a href="#{a}">{SECTION_NUMBER[a]}. {NAV_LABELS.get(a, t)}</a>\n'
        for a, t in SECTION_ORDER
    )
    return '<nav class="section-nav">\n' + links + '</nav>'



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


# ---------------------------------------------------------------------------
# Section 12: previous vs new ingest, with selectable filter parameters
# ---------------------------------------------------------------------------

CONCORDANCE_LEVELS = [1, 2, 3, 4, 5]  # 5 means ">=5"




def load_emitted_summary(output_dir: Path) -> dict:
    """Counts read straight out of the emitted KGX files.

    Every headline figure in this report comes from here rather than being written into
    prose, so the report cannot disagree with the artifacts it describes. Returns zeros
    if the transform has not been run yet.
    """
    nodes_path = output_dir / "clinvar_variant_nodes.tsv"
    edges_path = output_dir / "clinvar_variant_edges.tsv"
    summary = {
        "nodes": 0,
        "edges": 0,
        "by_predicate": Counter(),
        "genes": 0,
        "diseases": 0,
        "classes": Counter(),
        "available": False,
    }
    if not (nodes_path.exists() and edges_path.exists()):
        return summary

    summary["available"] = True
    with open(nodes_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            summary["nodes"] += 1
            if row.get("type"):
                summary["classes"][row["type"]] += 1

    genes: set = set()
    diseases: set = set()
    with open(edges_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            summary["edges"] += 1
            summary["by_predicate"][row["predicate"]] += 1
            if row["predicate"].endswith("is_sequence_variant_of"):
                genes.add(row["object"])
            else:
                diseases.add(row["object"])
    summary["genes"] = len(genes)
    summary["diseases"] = len(diseases)
    return summary


def load_hp_labels(data_dir: Path) -> dict:
    """HP term labels from the Monarch KG node dump, cached to a small extract.
    Used to name phenotype terms in section {S.star_cutoff} and to spot obsoleted ones (HPO keeps
    retired terms with an "obsolete " label prefix)."""
    cache = data_dir / "mk_hp_terms.tsv"
    if not cache.exists():
        nodes = data_dir / "monarch-kg_nodes.tsv"
        if not nodes.exists():
            return {}
        with open(nodes, newline="") as fh, open(cache, "w") as out:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["id"].startswith("HP:"):
                    out.write(f"{row['id']}\t{row['name']}\n")
    labels = {}
    with open(cache) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 1:
                labels[parts[0]] = parts[1]
    return labels


PHENO_HIST_CAP = 8


def build_phenotype_profile(
    clinvar_tsv: Path, var_records: dict, map_to_mondo: dict, mondo_labels: dict, hp_labels: dict
) -> dict:
    """How ClinVar attaches HPO terms to variants, and whether those terms vary between
    variants of the same disease.

    ClinVar does not record observed patient phenotypes per variant. HPO ids appear inside
    CLNDISDB groups alongside MONDO / MedGen / OMIM / Orphanet ids for the *same* condition
    -- they are cross-references naming the condition in HPO's vocabulary. This function
    measures the consequence: the per-variant term-count distribution (including variants
    with none), and, per disease, how much the HPO set actually differs across that
    disease's variants.

    consistency = (variants carrying the disease's most common HPO set) / (variants with
    any HPO set for it). A value at or near 1.0 means the terms are a property of the
    disease, not of the variant.
    """
    hist_all: Counter = Counter()
    hist_plp: Counter = Counter()
    hp_usage: Counter = Counter()
    disease_sets: dict = {}

    with open(clinvar_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            groups = parse_CLNDISDB(row["CLNDISDB"])
            hp_terms = {hp for g in groups for hp in g["HP"]}
            hist_all[min(len(hp_terms), PHENO_HIST_CAP)] += 1
            for hp in hp_terms:
                hp_usage[hp] += 1

            records = var_records.get(row["ID"])
            if records is None:
                continue
            dis, _, _ = variant_records_to_disease(records, map_to_mondo, star_min=0)
            if not dis:
                continue
            hist_plp[min(len(hp_terms), PHENO_HIST_CAP)] += 1

            mapped, _ = map_CLNDISDB_to_mondo(groups, map_to_mondo)
            for mondo_id, hps in map_mondo_to_hp(mapped, dis).items():
                disease_sets.setdefault(mondo_id, Counter())[frozenset(hps)] += 1

    rows = []
    for mondo_id, sets in disease_sets.items():
        non_empty = {k: v for k, v in sets.items() if k}
        n_with = sum(non_empty.values())
        if n_with < 2:
            continue
        modal_set, modal_n = max(non_empty.items(), key=lambda kv: kv[1])
        rows.append(
            {
                "mondo": mondo_id,
                "disease_name": mondo_labels.get(mondo_id, ""),
                "n_variants": sum(sets.values()),
                "n_with_hpo": n_with,
                "n_distinct_sets": len(non_empty),
                "consistency": round(modal_n / n_with, 4),
                "modal_terms": ", ".join(
                    f"{h} {hp_labels.get(h, '')}".strip() for h in sorted(modal_set)
                )[:180],
            }
        )
    rows.sort(key=lambda r: -r["n_variants"])

    consistency_hist: Counter = Counter()
    for r in rows:
        consistency_hist[min(int(r["consistency"] * 10), 10)] += 1

    obsolete = [(h, n) for h, n in hp_usage.most_common() if hp_labels.get(h, "").startswith("obsolete")]
    return {
        "hist_all": hist_all,
        "hist_plp": hist_plp,
        "rows": rows,
        "consistency_hist": consistency_hist,
        "n_distinct_hp": len(hp_usage),
        "top_hp": [(h, hp_labels.get(h, ""), n) for h, n in hp_usage.most_common(10)],
        "obsolete": obsolete[:10],
        "n_obsolete": len(obsolete),
        "obsolete_variants": sum(n for _h, n in obsolete),
    }


# ---------------------------------------------------------------------------
# Section 6: evidence tiers and the proposed predicate split
# ---------------------------------------------------------------------------

VAR_CITATIONS_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt"


def ensure_var_citations(data_dir: Path) -> Path:
    """ClinVar's per-variant literature links (VariationID -> PubMed / PMC / BookShelf)."""
    path = data_dir / "var_citations.txt"
    if path.exists():
        return path
    print(f"Downloading {VAR_CITATIONS_URL} (~220MB, one-time)...")
    with requests.get(VAR_CITATIONS_URL, stream=True, timeout=1800) as resp:
        resp.raise_for_status()
        with open(path, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
    return path


def load_pubmed_variants(data_dir: Path) -> set:
    """VariationIDs with at least one PubMed citation. PubMedCentral and NCBIBookShelf
    rows are excluded -- BookShelf entries are GeneReviews-style overviews rather than
    primary reports, and PMC duplicates PubMed."""
    path = ensure_var_citations(data_dir)
    cited: set = set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["citation_source"] == "PubMed":
                cited.add(row["VariationID"])
    return cited


# Proposed predicate assignment, by evidence tier rather than by ClinicalSignificance:
#   >=2 star (ClinVar aggregate) or inferred 2 star (concordance)  ->  causes
#   1 star with a supporting publication                           ->  associated_with
TIER_CAUSES = "causes"
TIER_ASSOCIATED = "associated_with"


def build_evidence_tiers(
    clinvar_tsv: Path,
    var_records: dict,
    map_to_mondo: dict,
    variant_genes: dict,
    pubmed_variants: set,
) -> dict:
    """Staged view of what each evidence tier would contribute, over every small variant
    in clinvar.vcf.gz (SNVs and indels -- structural variants are not in this file at all).

    Stages are cumulative:
      1. every variant in the VCF
      2. ...with >=1 Pathogenic / Likely pathogenic / Pathogenic-Likely pathogenic record
      3. ...whose ClinVar aggregate CLNREVSTAT reaches >=2 stars          -> causes
      4. ...plus those rescued by >=MIN_CONCORDANT_SUBMITTERS concordant
         submitters, i.e. an inferred 2 star                              -> causes
      5. ...plus 1-star records on a variant with a PubMed citation       -> associated_with

    Stage 3 subsumes the per-record >=3-star path: a variant with an expert-panel record
    aggregates to >=3 stars, so CLNREVSTAT already reflects it.

    A (gene, disease) pair is counted at `causes` if any of its variants qualifies via
    stage 3 or 4, and at `associated_with` only if it qualifies solely via stage 5 --
    the stronger tier wins where both apply.
    """
    stage_variants = {k: set() for k in range(1, 6)}
    stage_genes = {k: set() for k in range(1, 6)}
    stage_diseases = {k: set() for k in range(1, 6)}
    pairs_causes: set = set()
    pairs_assoc: set = set()
    lit_only_variants: set = set()
    clnvc_counts: Counter = Counter()

    with open(clinvar_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            varid = row["ID"]
            stage_variants[1].add(varid)
            clnvc_counts[row["CLNVC"]] += 1

            records = var_records.get(varid)
            if records is None:
                continue
            gene_entry = variant_genes.get(varid)
            gene_id = gene_entry[0] if gene_entry else None

            dis_any, _, _ = variant_records_to_disease(records, map_to_mondo, star_min=0)
            if not dis_any:
                continue
            stage_variants[2].add(varid)
            if gene_id:
                stage_genes[2].add(gene_id)
            stage_diseases[2].update(dis_any)

            agg_star = review_star_map.get(row["CLNREVSTAT"], 0)
            groups = concordance_groups(records, map_to_mondo)
            conc_diseases = {
                mondo_id
                for (mondo_id, _cs), entry in groups.items()
                if len(entry["submitters"]) >= MIN_CONCORDANT_SUBMITTERS
            }
            has_pub = varid in pubmed_variants
            if any(r["CollectionMethod"] == "literature only" for r in records):
                lit_only_variants.add(varid)

            # stage 3: ClinVar's own aggregate
            causes_dis = set(dis_any) if agg_star >= aggregate_star_min else set()
            if causes_dis:
                stage_variants[3].add(varid)
                if gene_id:
                    stage_genes[3].add(gene_id)
                stage_diseases[3].update(causes_dis)

            # stage 4: inferred 2 star via concordance
            causes_dis |= conc_diseases
            if causes_dis:
                stage_variants[4].add(varid)
                if gene_id:
                    stage_genes[4].add(gene_id)
                stage_diseases[4].update(causes_dis)

            # stage 5: 1-star with a publication, for whatever is left
            assoc_dis = (set(dis_any) - causes_dis) if (has_pub and agg_star >= 1) else set()
            if causes_dis or assoc_dis:
                stage_variants[5].add(varid)
                if gene_id:
                    stage_genes[5].add(gene_id)
                stage_diseases[5].update(causes_dis | assoc_dis)

            if gene_id:
                for d in causes_dis:
                    pairs_causes.add((gene_id, d))
                for d in assoc_dis:
                    pairs_assoc.add((gene_id, d))

    # the stronger tier wins where a pair qualifies both ways
    pairs_assoc -= pairs_causes

    return {
        "stage_variants": {k: len(v) for k, v in stage_variants.items()},
        "stage_genes": {k: len(v) for k, v in stage_genes.items()},
        "stage_diseases": {k: len(v) for k, v in stage_diseases.items()},
        "pairs_causes": len(pairs_causes),
        "pairs_assoc": len(pairs_assoc),
        "n_pubmed": len(pubmed_variants),
        "n_lit_only": len(lit_only_variants),
        "clnvc": clnvc_counts,
    }


def build_filter_cube(clinvar_tsv: Path, var_records: dict, map_to_mondo: dict, variant_genes: dict) -> dict:
    """Cubes supporting the interactive filter comparison in section {S.structural_variants}.

    The production inclusion rule is a disjunction of three monotone conditions on a
    (variant, disease):

        per-record stars >= star_min  OR  concordant submitters >= N  OR  aggregate stars >= agg_min

    Because each condition alone suffices, both counts below are exact under any
    parameter combination rather than approximations:

      variant_cube  keyed (max per-record star over all this variant's diseases,
                    max concordant submitters over all its diseases, aggregate star,
                    has a ClinVar-asserted gene) -> variant count.
                    A variant is included iff any one of its diseases is, and
                    max_star >= star_min implies some disease reaches star_min.
      edge_cube     keyed (per-disease max star, per-disease concordance, aggregate
                    star) -> (variant, disease) count, i.e. disease edges.

    NOT modelled: the CLNDISDB-echo gate (stage 5 in FILTERING.md). It is all-or-nothing
    per variant and depends on which diseases survive, so it cannot be folded into a
    static cube. It is very nearly inert -- under the previous filter it dropped 2 of
    a couple of variants in a hundred thousand -- and section {S.ingest_compare} reports the cube's
    prediction against the emitted counts so the size of that approximation stays visible.
    emitted counts so the size of that approximation stays visible.
    """
    variant_cube: Counter = Counter()
    edge_cube: Counter = Counter()

    def conc_bucket(n: int) -> int:
        return min(n, CONCORDANCE_LEVELS[-1])

    with open(clinvar_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            varid = row["ID"]
            records = var_records.get(varid)
            if records is None:
                continue
            agg_star = review_star_map.get(row["CLNREVSTAT"], 0)

            # per-disease max review-star among supporting P/LP records
            star_by_disease: dict = {}
            for rec in records:
                clinsig = rec["ClinicalSignificance"]
                if clinsig not in predicate_map:
                    continue
                stars = review_star_map[rec["ReviewStatus"].replace(" ", "_")]
                mondo_ids = set()
                for mg in rec["ReportedPhenotypeInfo"].split(";"):
                    key = "MedGen:{}".format(mg.split(":")[0])
                    if key in map_to_mondo:
                        mondo_ids.update(map_to_mondo[key])
                if not mondo_ids:
                    for dis_id in rec["SubmittedPhenotypeInfo"].split(";"):
                        dis_id = format_id_to_map(dis_id)
                        if dis_id in map_to_mondo:
                            mondo_ids.update(map_to_mondo[dis_id])
                        elif dis_id is not None and "MONDO:" in dis_id:
                            mondo_ids.add(dis_id)
                for mondo_id in mondo_ids:
                    if stars > star_by_disease.get(mondo_id, -1):
                        star_by_disease[mondo_id] = stars

            if not star_by_disease:
                continue

            # per-disease max concordant-submitter count, across ClinicalSignificance values
            groups = concordance_groups(records, map_to_mondo)
            conc_by_disease: dict = {}
            for (mondo_id, _clinsig), entry in groups.items():
                n = len(entry["submitters"])
                if n > conc_by_disease.get(mondo_id, 0):
                    conc_by_disease[mondo_id] = n

            for mondo_id, stars in star_by_disease.items():
                edge_cube[(stars, conc_bucket(conc_by_disease.get(mondo_id, 1)), agg_star)] += 1

            variant_cube[
                (
                    max(star_by_disease.values()),
                    conc_bucket(max(conc_by_disease.values()) if conc_by_disease else 1),
                    agg_star,
                    variant_genes.get(varid) is not None,
                )
            ] += 1

    return {
        "variants": [{"s": k[0], "c": k[1], "a": k[2], "g": k[3], "n": v} for k, v in sorted(variant_cube.items())],
        "edges": [{"s": k[0], "c": k[1], "a": k[2], "n": v} for k, v in sorted(edge_cube.items())],
    }


def build_clnsig_cube(
    clinvar_tsv: Path, var_records: dict, map_to_mondo: dict, mondo_labels: dict, strchive_intervals: dict,
    variant_genes: dict,
) -> tuple[list[dict], list[str], list[dict]]:
    """Full pass over clinvar.tsv (every row, not just those with a
    submission_summary match) tallying (star, clnsig_bucket, clnvc,
    has_literature, has_concordance), plus -- for rows with a
    submission_summary match -- the set of (gene, MONDO disease) pairs
    implied by that variant (any Pathogenic/Likely-pathogenic submission
    record, any star; same disease-mapping logic as section {S.star_cutoff}'s star_min=0
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
    same signal used for section {S.star_cutoff}'s rescue analysis. Variants absent from
    submission_summary.txt get False for both (undetermined).
    in_strchive: does this variant's genomic footprint overlap one of
    STRchive's ~80 curated pathogenic short-tandem-repeat loci (hg38
    coordinates)? See in_strchive() -- informational only, not used to gate
    disease-pair inclusion.
    in_production: does this variant survive production's actual inclusion
    filter -- i.e. does variant_records_to_disease() at
    star_min=PRODUCTION_STAR_MIN with the
    >=MIN_CONCORDANT_SUBMITTERS concordance rescue still yield >=1 MONDO
    disease? This is exactly the criterion behind section {S.star_cutoff}'s headline
    counts (and so, like them, it deliberately ignores process_row()'s
    separate HPO-overlap requirement). It is computed from per-submission
    ReviewStatus, NOT from the variant-level CLNREVSTAT star shown in the
    "star" dimension -- the two star systems are different (see module
    docstring), which is why this is its own dimension rather than something
    you can reconstruct by ticking star boxes.

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
            in_production = False
            dis: dict = {}
            if records is not None:
                has_literature = any(rec["CollectionMethod"] == "literature only" for rec in records)
                groups = concordance_groups(records, map_to_mondo)
                has_concordance = any(
                    len(entry["submitters"]) >= MIN_CONCORDANT_SUBMITTERS for entry in groups.values()
                )
                dis, _, _ = variant_records_to_disease(records, map_to_mondo, star_min=0)
                prod_dis, _, _ = variant_records_to_disease(
                    records,
                    map_to_mondo,
                    star_min=PRODUCTION_STAR_MIN,
                    rescue_min_submitters=MIN_CONCORDANT_SUBMITTERS,
                )
                in_production = bool(prod_dis)

            cell_key = (star, clnsig, clnvc, size, has_literature, has_concordance, in_strchive_flag, in_production)
            counts[cell_key] += 1

            if not dis:
                continue
            gene_ids, gene_symbols = variant_genes_for(variant_genes, row["ID"])
            cell_map = cell_pairs.setdefault(cell_key, {})
            for gene_id, gene_sym in zip(gene_ids, gene_symbols):
                pair_gene_symbol.setdefault(gene_id, gene_sym)
                for mondo_id in dis:
                    # keyed by (gene_id, disease) only, matching section the star-cutoff section's canonical pair
                    # identity -- NCBIGene id is what production actually emits; a handful of
                    # ids (135 observed) carry >1 symbol spelling across records (gene renames),
                    # which would otherwise inflate this count vs section the star-cutoff section's headline number.
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
            "production": production,
            "count": n,
            "pairs": [
                {"i": idx, "n": e["n"], "v": e["sample"]}
                for idx, e in sorted(
                    cell_pairs.get(
                        (star, clnsig, clnvc, size, literature, concordant, strchive, production), {}
                    ).items()
                )
            ],
        }
        for (star, clnsig, clnvc, size, literature, concordant, strchive, production), n in sorted(counts.items())
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
# Section 12: structural variants / CNVs -- what's not in the VCF
# ---------------------------------------------------------------------------

VARIANT_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
# Type values from variant_summary.txt that never appear in the VCF's CLNVC field --
# the VCF requires a fixed genomic position + REF/ALT, which structurally excludes these.
SV_TYPES = ["copy number gain", "copy number loss", "Translocation", "Complex", "fusion", "Tandem duplication", "Inversion"]


def resolve_phenotype_ids_to_mondo(phenotype_ids: str, map_to_mondo: dict) -> set:
    """Parse variant_summary.txt's PhenotypeIDS field ("|"-separated phenotypes,
    each a comma-separated list of cross-ref ids) into MONDO ids. Prefers a
    direct "MONDO:MONDO:0000000"-style tag -- ClinVar's own double-prefixed
    format, present on most disease-resolved CNV/SV rows -- and falls back to
    the same MedGen/OMIM -> MONDO lookup (map_to_mondo) used everywhere else in
    this report for rows that only carry MedGen/OMIM cross-refs."""
    mondo_ids = set()
    for group in phenotype_ids.split("|"):
        for token in group.split(","):
            token = token.strip()
            if token.startswith("MONDO:MONDO:"):
                mondo_ids.add("MONDO:" + token.rsplit(":", 1)[-1])
            elif token in map_to_mondo:
                mondo_ids.update(map_to_mondo[token].keys())
    return mondo_ids


def is_single_clean_gene_symbol(gene_symbol: str) -> bool:
    """variant_summary.txt's GeneSymbol is ';'-delimited for a handful of genes
    or degrades to prose ("subset of 2646 genes: ..." / "covers 42 genes...")
    for large regional CNVs -- and its GeneID column is the sentinel "-1" for
    any of those multi-gene rows, so there's no reliable per-gene id to pair
    with a disease. Restrict to rows unambiguously about exactly one gene."""
    return gene_symbol not in ("-", "") and ";" not in gene_symbol and not gene_symbol.startswith(("subset of", "covers "))


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


def build_sv_summary(data_dir: Path, map_to_mondo: dict) -> dict:
    """Two passes over variant_summary.txt.gz tallying the Type
    distribution and, for the SV_TYPES specifically: count, how many resolve to a
    disease id (PhenotypeIDS not '-'/empty) vs not, genomic span (Stop - Start), and
    every row (resolved AND unresolved, any ClinicalSignificance) for the full
    paginated table -- not just a curated top-N. Each row carries a `resolved`
    boolean so the report can offer a "SV resolved to disease" filter checkbox
    rather than silently only ever showing the resolved subset.

    Also builds gene_disease_pairs: for Pathogenic/Likely-pathogenic SV_TYPES
    rows unambiguously about exactly one gene (is_single_clean_gene_symbol),
    the (gene_sym, gene_id, mondo_id) pairs implied by that gene + its resolved
    disease(s) -- used to find pairs whose only ClinVar evidence is structural,
    invisible to production's SNV/indel-only pipeline (clinvar.vcf.gz has no
    way to represent a CNV/inversion at all).

    The file carries one row per variant per genome build. The first pass records
    which build should represent each VariationID (ASSEMBLY_PREFERENCE: GRCh38,
    else GRCh37) and the second tallies exactly one row per variant, so a variant
    present on both builds is never counted twice and one ClinVar has only placed
    on GRCh37 is not silently dropped.
    """
    path = ensure_variant_summary_downloaded(data_dir)

    type_counts: Counter = Counter()
    sv_stats = {t: {"count": 0, "resolved": 0, "spans": []} for t in SV_TYPES}
    sv_rows = []
    gene_disease_pairs: dict = {}

    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        hcols = {k: i for i, k in enumerate(header)}
        # First pass: which build should represent each variant. Buffering the rows
        # themselves would mean holding millions of split lines in memory, so this
        # re-reads the file instead and keeps only a rank per VariationID.
        sv_rank = {name: i for i, name in enumerate(ASSEMBLY_PREFERENCE)}
        sv_best: dict = {}
        sv_unplaced: set = set()
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            varid = cols[hcols["VariationID"]]
            rank = sv_rank.get(cols[hcols["Assembly"]])
            if rank is None:
                # NCBI36 / "na" -- unusable. Track SV-typed ones so the count of variants
                # this section cannot place is reported rather than silently dropped.
                if cols[hcols["Type"]] in SV_TYPES:
                    sv_unplaced.add(varid)
                continue
            seen = sv_best.get(varid)
            if seen is None or rank < seen:
                sv_best[varid] = rank

    sv_unplaced -= set(sv_best)
    sv_build_stats: dict = {name: {"count": 0, "resolved": 0} for name in ASSEMBLY_PREFERENCE}

    sv_done: set = set()
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            rank = sv_rank.get(cols[hcols["Assembly"]])
            varid = cols[hcols["VariationID"]]
            # exactly one row per variant: the preferred build, first occurrence only
            if rank is None or rank != sv_best.get(varid) or varid in sv_done:
                continue
            sv_done.add(varid)

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
            build = cols[hcols["Assembly"]]
            if build in sv_build_stats:
                sv_build_stats[build]["count"] += 1
                if resolved:
                    sv_build_stats[build]["resolved"] += 1

            try:
                span = int(cols[hcols["Stop"]]) - int(cols[hcols["Start"]])
            except ValueError:
                span = None
            if span is not None:
                stat["spans"].append(span)

            if span is not None:
                variation_id = cols[hcols["VariationID"]]
                clinsig = cols[hcols["ClinicalSignificance"]]
                sv_rows.append(
                    {
                        "name": cols[hcols["Name"]],
                        "type": vtype,
                        "allele_id": cols[hcols["AlleleID"]],
                        "variation_id": variation_id,
                        "clinsig": clinsig,
                        "phenotype": cols[hcols["PhenotypeList"]],
                        "phenotype_ids": pheno_ids,
                        "resolved": resolved,
                        "span": span,
                        "review_status": cols[hcols["ReviewStatus"]],
                        "num_submitters": int(cols[hcols["NumberSubmitters"]]),
                    }
                )

                gene_symbol = cols[hcols["GeneSymbol"]]
                if resolved and clinsig in predicate_map and is_single_clean_gene_symbol(gene_symbol):
                    gene_id = cols[hcols["GeneID"]]
                    for mondo_id in resolve_phenotype_ids_to_mondo(pheno_ids, map_to_mondo):
                        key = (gene_symbol, gene_id, mondo_id)
                        entry = gene_disease_pairs.setdefault(key, {"types": set(), "variants": set()})
                        entry["types"].add(vtype)
                        entry["variants"].add(variation_id)

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
        "gene_disease_pairs": gene_disease_pairs,
        "build_stats": sv_build_stats,
        "unplaced": len(sv_unplaced),
    }


def sv_only_gene_disease_pairs(sv_pairs: dict, snv_pair_set: set, mondo_labels: dict) -> list[dict]:
    """Gene-disease pairs with CNV/SV evidence (build_sv_summary's
    gene_disease_pairs) that have NO corresponding SNV/indel evidence anywhere
    in clinvar.vcf.gz (snv_pair_set -- section {S.star_cutoff}'s own star_min=0 pair_sets[0],
    any Pathogenic/Likely-pathogenic evidence, any review status). These pairs
    are structurally invisible to production: it only ever reads
    clinvar.vcf.gz, which has no way to represent a CNV or inversion at all,
    so a gene-disease association resting solely on structural evidence never
    produces a VariantToGeneAssociation/VariantToDiseaseAssociation edge today,
    no matter how well-reviewed that structural evidence is."""
    rows = []
    for (gene_sym, gene_id, mondo_id), entry in sv_pairs.items():
        if (f"NCBIGene:{gene_id}", mondo_id) in snv_pair_set:
            continue
        rows.append(
            {
                "gene": gene_sym,
                "gene_id": gene_id,
                "mondo": mondo_id,
                "disease_name": mondo_labels.get(mondo_id, ""),
                "types": sorted(entry["types"]),
                "n_variants": len(entry["variants"]),
                "variant_sample": sorted(entry["variants"], key=int)[:MAX_VARIANT_SAMPLE],
            }
        )
    rows.sort(key=lambda r: (-r["n_variants"], r["gene"]))
    return rows


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------


def render_html(
    results: dict,
    pairs_ge3: list,
    pairs_2star_concordance: list,
    pairs_remaining: list,
    cube: list,
    clnvc_types: list,
    pair_list: list,
    sv_summary: dict,
    example_variants: dict,
    gene_model_results: list[dict],
    clingen_summary: dict,
    clingen_only_rows: list,
    sv_only_pairs: list,
    review_status_rows: list,
    monarch_comparison: dict,
    gene_attribution: dict,
    variant_genes: dict,
    submission_profile: dict,
    filter_cube: dict,
    phenotype_profile: dict,
    evidence_tiers: dict,
    emitted: dict,
    scn1a_subgraph: dict,
) -> str:
    max_variants = max(r["variants"] for r in results.values()) or 1
    max_pairs = max(r["gene_disease_pairs"] for r in results.values()) or 1

    chart_width, bar_height, gap, left_pad = 640, 32, 16, 90
    # Display order for section the star-cutoff section's chart/table: "2c" (2-star, computed) is inserted right after the
    # raw 2-star level, since it's the reconstructed version of that same tier -- see compute_star_data().
    CHART_LEVELS = [0, 1, 2, "2c", 3, 4]

    def level_label(level) -> str:
        if level == "2c":
            return "2 stars *"
        return f"{level} star{'s' if level != 1 else ''}"

    def bars(metric_key: str, max_val: int) -> str:
        rows = []
        for i, level in enumerate(CHART_LEVELS):
            value = results[level][metric_key]
            bar_w = (value / max_val) * (chart_width - left_pad - 60) if max_val else 0
            y = i * (bar_height + gap)
            if level == "2c":
                fill = " fill='#8b5cf6' stroke='#6d28d9' stroke-width='1.5' stroke-dasharray='4,3'"
            elif level == PRODUCTION_STAR_MIN:
                fill = " fill='#f59e0b'"
            else:
                fill = " fill='#2563eb'"
            rows.append(
                f"<text x='{left_pad - 10}' y='{y + bar_height / 2 + 5}' text-anchor='end' "
                f"font-size='13' fill='#334155'>{level_label(level)}</text>"
                f"<rect x='{left_pad}' y='{y}' width='{bar_w:.1f}' height='{bar_height}'{fill} rx='3'/>"
                f"<text x='{left_pad + bar_w + 8}' y='{y + bar_height / 2 + 5}' font-size='13' "
                f"fill='#0f172a'>{value:,}</text>"
            )
        svg_height = len(CHART_LEVELS) * (bar_height + gap)
        return (
            f"<svg viewBox='0 0 {chart_width} {svg_height}' xmlns='http://www.w3.org/2000/svg'>"
            + "".join(rows)
            + "</svg>"
        )

    row_class = ' class="prod"'
    table_rows = "".join(
        (
            f"<tr style='background:#f5f3ff;'>"
            f"<td>2 <span style='color:#6d28d9; font-weight:700;'>*</span> (computed)</td>"
            f"<td>{results['2c']['variants']:,}</td>"
            f"<td>{results['2c']['gene_disease_pairs']:,}</td>"
            f"</tr>"
        )
        if level == "2c"
        else (
            f"<tr{row_class if level == PRODUCTION_STAR_MIN else ''}>"
            f"<td>{level}{' (production default)' if level == PRODUCTION_STAR_MIN else ''}</td>"
            f"<td>{results[level]['variants']:,}</td>"
            f"<td>{results[level]['gene_disease_pairs']:,}</td>"
            f"</tr>"
        )
        for level in CHART_LEVELS
    )

    UNRESOLVED_PHENOTYPES = {"See cases", "not provided", "-", ""}

    def mutation_gallery_card(entry: dict) -> str:
        row = example_variants[entry["key"]]
        if entry["kind"] == "tsv":
            gene = (variant_genes_for(variant_genes, row["ID"])[1] or ["?"])[0]
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
                f"~{row['span']:,} bp &mdash; not in clinvar.vcf.gz, see section {S.structural_variants}</div>"
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
        domain_html = ""
        if result.get("domain_svg"):
            domain_html = f"""<div class="xf-panel" style="margin-top:1rem;">
  <div style="font-size:11px; font-weight:700; color:#2563eb; text-transform:uppercase; letter-spacing:0.03em;">Protein domain view</div>
  <div style="margin-top:8px;">{result['domain_svg']}</div>
  <div style="font-size:12px; color:#475569; margin-top:6px;">
    Real domain boundaries and disease-variant positions from
    <a href="https://www.uniprot.org/uniprotkb/{result['domain_uniprot_id']}" target="_blank" rel="noopener">UniProt {result['domain_uniprot_id']}</a>.
    {result['domain_caption']}
  </div>
</div>"""
        phenotype_html = ""
        if result.get("phenotype_svg"):
            phenotype_html = f"""<div class="xf-panel" style="max-width:520px; margin-top:1rem;">
  <div style="font-size:11px; font-weight:700; color:#2563eb; text-transform:uppercase; letter-spacing:0.03em;">BAMS vs. FSHD2 phenotype overlap</div>
  <div style="margin-top:8px;">{result['phenotype_svg']}</div>
  <div style="font-size:12px; color:#475569; margin-top:6px;">
    Despite arising from the same gene, BAMS and FSHD2 have essentially no overlapping hallmark clinical
    features -- consistent with the opposite mutation-distribution story above (domain-specific missense vs.
    loss-of-function-anywhere), these read as two clinically distinct diseases that happen to share a locus.
  </div>
</div>"""
        return f"""<h3>{result['heading']}</h3>
<p class="subtitle">{blurb}</p>
<div class="table-wrap" style="overflow-x:auto; padding:12px;">
{result['svg']}
</div>
<div style="display:flex; flex-wrap:wrap; gap:4px 24px; margin-top:8px;">
{legend_html}
</div>
{domain_html}
{phenotype_html}"""

    gene_model_sections_html = "".join(gene_model_section(r) for r in gene_model_results)

    pairs_ge3_json = json.dumps(pairs_ge3)
    pairs_2star_json = json.dumps(pairs_2star_concordance)
    pairs_remaining_json = json.dumps(pairs_remaining)
    cube_json = json.dumps(cube)
    clnsig_buckets_json = json.dumps(CLNSIG_BUCKETS)
    clnvc_types_json = json.dumps(clnvc_types)
    size_buckets_json = json.dumps(SIZE_BUCKETS)
    pair_list_json = json.dumps(pair_list)
    sv_rows_json = json.dumps(sv_summary["sv_rows"])
    total_variants = sum(r["count"] for r in cube)

    clingen_single_pct = (
        100 * clingen_summary["single_submitter_clingen_pairs"] / clingen_summary["single_submitter_pairs"]
        if clingen_summary["single_submitter_pairs"]
        else 0
    )
    clingen_only_pct = (
        100 * clingen_summary["clingen_only_pairs"] / clingen_summary["total_pairs"] if clingen_summary["total_pairs"] else 0
    )
    clingen_only_json = json.dumps(clingen_only_rows)

    sv_b38 = sv_summary["build_stats"]["GRCh38"]
    sv_b37 = sv_summary["build_stats"]["GRCh37"]
    sv_gene_disease_pairs_count = len(sv_summary["gene_disease_pairs"])
    sv_only_pairs_count = len(sv_only_pairs)
    sv_only_pct = 100 * sv_only_pairs_count / sv_gene_disease_pairs_count if sv_gene_disease_pairs_count else 0
    sv_only_json = json.dumps(sv_only_pairs)

    total_review_status_records = sum(r["count"] for r in review_status_rows)
    review_status_table_rows = "".join(
        (
            f"<tr class='no'><td>{r['stars'] if r['stars'] is not None else '?'}</td>"
            f"<td>{r['status']} <em>(never observed per-record)</em></td>"
            f"<td class='num'>0</td><td class='num'>0.0%</td></tr>"
        )
        if r["count"] == 0
        else (
            f"<tr><td>{r['stars'] if r['stars'] is not None else '?'}</td>"
            f"<td>{r['status']}</td>"
            f"<td class='num'>{r['count']:,}</td>"
            f"<td class='num'>{100 * r['count'] / total_review_status_records:.1f}%</td></tr>"
        )
        for r in review_status_rows
    )

    mc = monarch_comparison
    scn1a_diagram = diagram_scn1a_hierarchy(scn1a_subgraph)
    _sh = scn1a_subgraph.get("shares", {}) if scn1a_subgraph else {}
    shared_headline = f"{max(_sh.values()):,}" if _sh else "the"

    # Headline figures come from the emitted artifacts, never from prose
    em = emitted
    em_gene = em["by_predicate"].get("biolink:is_sequence_variant_of", 0)
    em_causes = em["by_predicate"].get("biolink:causes", 0)
    em_assoc = em["by_predicate"].get("biolink:associated_with_increased_likelihood_of", 0)
    em_classes = ", ".join(f"{c} ({n:,})" for c, n in em["classes"].most_common()) or "n/a"

    # --- input files ---------------------------------------------------------
    ga = gene_attribution
    sf = submission_profile
    vcf_f = ga["vcf_field"]
    vs_f = ga["vs_field"]

    def pct(n, d):
        return f"{100 * n / d:.1f}%" if d else "-"

    vcf_field_rows = "".join(
        f"<tr><td><code>{col}</code></td><td>{desc}</td>"
        f"<td class='num'>{vcf_f['empty_' + col]:,}</td>"
        f"<td class='num'>{pct(vcf_f['empty_' + col], vcf_f['rows'])}</td></tr>"
        for col, desc in (
            ("CLNDISDB", "Disease cross-references (MONDO / MedGen / OMIM / Orphanet / HPO)"),
            ("GENEINFO", "Every gene overlapping the position &mdash; <strong>not</strong> the causal gene"),
            ("MC", "Molecular consequence, SO terms"),
            ("RS", "dbSNP rsID"),
        )
    )
    method_rows = "".join(
        f"<tr><td>{m}</td><td class='num'>{n:,}</td><td class='num'>{pct(n, sf['n_records'])}</td></tr>"
        for m, n in sf["methods"].most_common(6)
    )
    placeholder_pct = pct(sf["placeholder_phenotype"], sf["n_records"])

    input_files_html = f"""{section_heading("input-files")}
<p class="subtitle">
  Everything downstream is built from three ClinVar files plus two mapping files. They describe the
  same variants at different granularities, and the choice of which field to trust in which file is
  where most of this report's findings live. Row and field counts below are computed from the actual
  downloaded release, not quoted from documentation.
</p>

<div class="table-wrap">
<table>
<thead><tr><th>File</th><th>Grain (one row =)</th><th>Rows</th><th>What the ingest takes from it</th></tr></thead>
<tbody>
<tr>
  <td><code>clinvar.vcf.gz</code><br><span style="color:#64748b; font-size:11.5px;">&rarr; <code>clinvar.tsv</code></span></td>
  <td>one <strong>variant</strong> (ClinVar VariationID)</td>
  <td class="num">{vcf_f['rows']:,}</td>
  <td>The variant itself, its aggregate classification, and its disease cross-references</td>
</tr>
<tr>
  <td><code>submission_summary.txt.gz</code></td>
  <td>one <strong>submission</strong> (variant &times; submitting lab)</td>
  <td class="num">{sf['n_records']:,}</td>
  <td>Per-lab classification, review status and reported phenotype &mdash; the evidence the star and concordance filters run over</td>
</tr>
<tr>
  <td><code>variant_summary.txt.gz</code></td>
  <td>one variant &times; genome build</td>
  <td class="num">{vs_f['rows_kept']:,} <span style="color:#64748b;">(deduped)</span></td>
  <td>ClinVar's curated per-variant <strong>gene attribution</strong>, and the structural variants absent from the VCF (section {S.structural_variants})</td>
</tr>
<tr>
  <td><code>mondo.sssom.tsv</code></td>
  <td>one cross-reference mapping</td>
  <td class="num">&mdash;</td>
  <td>OMIM / Orphanet / MeSH &rarr; MONDO</td>
</tr>
<tr>
  <td><code>MedGenIDMappings.txt.gz</code></td>
  <td>one MedGen concept mapping</td>
  <td class="num">&mdash;</td>
  <td>MedGen CUI &rarr; MONDO, the primary disease route for submission records</td>
</tr>
</tbody>
</table>
</div>

<h3>1. <code>clinvar.vcf.gz</code> &mdash; the variant-level view</h3>
<p class="subtitle">
  {vcf_f['rows']:,} variants, 43 INFO columns after <code>scripts/vcf_to_tsv.py</code> flattens it.
  <code>ID</code> is the ClinVar VariationID and the join key to both other files. The fields that
  matter:
</p>
<div class="table-wrap">
<table>
<thead><tr><th>Field</th><th>Used for</th></tr></thead>
<tbody>
<tr><td><code>ID</code></td><td><code>SequenceVariant.id</code> (<code>CLINVAR:&lt;id&gt;</code>) and the join key to the other two files</td></tr>
<tr><td><code>CLNHGVS</code></td><td><code>SequenceVariant.name</code> &mdash; genomic HGVS</td></tr>
<tr><td><code>CLNDISDB</code></td><td>Disease cross-refs, parsed by <code>parse_CLNDISDB()</code>. This is where HPO terms come from, and the HPO overlap gates whether the transform emits anything at all</td></tr>
<tr><td><code>MC</code></td><td><code>SequenceVariant.type</code> &mdash; SO molecular-consequence terms</td></tr>
<tr><td><code>CLNREVSTAT</code></td><td>Aggregate review status. Recorded as an association qualifier and drives section {S.crossfilter}'s star panel. <strong>Not</strong> the star value the production filter uses &mdash; see below</td></tr>
<tr><td><code>RS</code></td><td><code>SequenceVariant.xref</code> &mdash; dbSNP</td></tr>
<tr><td><code>CLNSIG</code>, <code>CLNVC</code></td><td>Aggregate significance and variant type. Analysis only (section {S.crossfilter}'s crossfilter); the transform classifies from per-submission records instead</td></tr>
<tr class="no"><td><code>ONC*</code>, <code>SCI*</code></td><td>Oncogenicity and somatic-impact classifications. Unused &mdash; this ingest is germline only</td></tr>
</tbody>
</table>
</div>
<p class="subtitle">How often those fields are simply absent:</p>
<div class="table-wrap">
<table>
<thead><tr><th>Field</th><th>Meaning</th><th>Rows with <code>.</code></th><th>%</th></tr></thead>
<tbody>{vcf_field_rows}</tbody>
</table>
</div>

<div class="summary-box" style="border-left:4px solid #b91c1c;">
  <strong>&#9888; <code>GENEINFO</code> &mdash; the field that caused the antisense problem.</strong>
  Format is <code>SYMBOL:GeneID|SYMBOL:GeneID|...</code>, and it lists <em>every locus whose span covers
  the variant position</em>. That is what the field is documented to do &mdash; it is a positional
  annotation, not a causal claim. The ingest read it as a causal claim.
  <br><br>
  For F508del, <code>GENEINFO</code> is
  <code>CFTR:1080|CFTR-AS1:111082987</code>. Building a <code>VariantToGeneAssociation</code> from each
  entry gave <code>CFTR-AS1</code> its own edge, and from there the antisense transcript inherited
  CFTR's entire disease list and submitter roster &mdash; appearing as a gene-disease pair with 76
  independent submitters behind it. The same happened for <code>TTN-AS1</code> (74),
  <code>SYNGAP1-AS1</code> (71), <code>SNHG14</code>&rarr;Angelman, <code>RIF1</code>&rarr;nemaline
  myopathy 2, <code>GH-LCR</code>, <code>NPHP3-ACAD11</code>, and NCBI's <code>LOC</code>&lt;GeneID&gt;
  placeholders.
  <br><br>
  <strong>Fixed:</strong> gene attribution now comes from <code>variant_summary.txt.gz</code> instead.
  {ga['multi_gene_variants']:,} variants had &gt;1 gene in <code>GENEINFO</code>; all of them resolve to
  exactly one curated gene. See section {S.ingest_recommendation} for the full before/after.
</div>

<h3>2. <code>submission_summary.txt.gz</code> &mdash; the per-submission evidence</h3>
<p class="subtitle">
  {sf['n_records']:,} submission records across {sf['n_variants']:,} variants from
  <strong>{sf['n_submitters']:,} distinct submitters</strong> &mdash; roughly
  {sf['n_records'] // max(sf['n_variants'], 1)} records per variant on average, though the distribution is
  heavily skewed. This is the file the star cutoff and concordance rescue actually run over
  (<code>make_variant_record_map()</code>).
</p>
<div class="table-wrap">
<table>
<thead><tr><th>Field</th><th>Used for</th></tr></thead>
<tbody>
<tr><td><code>VariationID</code></td><td>Join key to the VCF</td></tr>
<tr><td><code>ClinicalSignificance</code></td><td>Gated through <code>predicate_map</code> &mdash; only the Pathogenic / Likely-pathogenic family produces disease edges. Also the concordance key</td></tr>
<tr><td><code>ReviewStatus</code></td><td>Scored by <code>review_star_map</code> &mdash; <strong>this</strong> is the star value the production filter uses</td></tr>
<tr><td><code>ReportedPhenotypeInfo</code></td><td>Primary disease route: <code>MedGen:&lt;CUI&gt;</code> &rarr; MONDO</td></tr>
<tr><td><code>SubmittedPhenotypeInfo</code></td><td>Fallback disease route when the MedGen lookup yields nothing</td></tr>
<tr><td><code>Submitter</code></td><td>Identity behind concordance and every pooled-submitter count in this report</td></tr>
<tr><td><code>CollectionMethod</code></td><td>Analysis only &mdash; flags <code>literature only</code> records</td></tr>
</tbody>
</table>
</div>

<div class="summary-box" style="border-left:4px solid #b91c1c;">
  <strong>&#9888; Two fields in this file quietly discard evidence.</strong>
  <br><br>
  <strong><code>ReportedPhenotypeInfo</code> is a placeholder {placeholder_pct} of the time.</strong>
  {sf['placeholder_phenotype']:,} of {sf['n_records']:,} records carry
  <code>C3661900:not provided</code> (or an equivalent). That MedGen concept maps to no MONDO term, so
  those records contribute <em>nothing</em> to disease mapping or concordance &mdash; no matter how
  authoritative the submitter. Some of the highest-volume clinical labs submit this way as a matter of
  course: GeneDx contributes 535 such Pathogenic/Likely-pathogenic records for SCN1A alone, CeGaT 148.
  Pooled submitter counts throughout this report therefore <em>understate</em> support, systematically
  and most severely for the busiest labs.
  <br><br>
  <strong><code>ReviewStatus</code> only ever takes 6 of the 10 documented values here</strong>, and
  "criteria provided, multiple submitters, no conflicts" (2&#9733;) is never one of them &mdash; it is an
  aggregate ClinVar computes <em>across</em> submitters and publishes only at the variant level. A
  per-record <code>star_min=2</code> filter is therefore identical to <code>star_min=3</code>. See
  section {S.clinvar_curation}'s table for the observed distribution.
</div>
<p class="subtitle">
  <code>CollectionMethod</code> distribution, for context on what kind of evidence these records represent:
</p>
<div class="table-wrap">
<table>
<thead><tr><th>CollectionMethod</th><th>Records</th><th>%</th></tr></thead>
<tbody>{method_rows}</tbody>
</table>
</div>

<h3>3. <code>variant_summary.txt.gz</code> &mdash; curated gene attribution</h3>
<p class="subtitle">
  Added to the ingest to fix the <code>GENEINFO</code> defect. One row per variant per genome build, so
  it must be filtered to <code>Assembly == GRCh38</code> &mdash; {vs_f['rows_grch38']:,} rows &mdash; or
  every variant is counted twice. Also the only place ClinVar publishes structural variants, which the
  VCF structurally cannot represent (section {S.structural_variants}).
</p>
<div class="table-wrap">
<table>
<thead><tr><th>Field</th><th>Used for</th></tr></thead>
<tbody>
<tr class="prod"><td><code>GeneID</code> / <code>GeneSymbol</code></td><td><strong>The gene ClinVar actually attributes the variant to</strong> &mdash; now the sole source of variant-gene edges (<code>make_variant_gene_map()</code>)</td></tr>
<tr><td><code>HGNC_ID</code></td><td>Not yet used. Worth noting: this is the id space the Monarch KG keys genes on, so emitting it would remove the HGNC&harr;Entrez crosswalk section {S.monarch_kg} has to build</td></tr>
<tr><td><code>Assembly</code></td><td>Build selection &mdash; GRCh38 preferred, GRCh37 only where no GRCh38 row exists, never both</td></tr>
<tr><td><code>Type</code>, <code>Start</code>, <code>Stop</code></td><td>Structural-variant analysis in section {S.structural_variants}</td></tr>
<tr><td><code>Name</code></td><td>Transcript-anchored HGVS, e.g. <code>NM_000492.3(CFTR):c.1521_1523del</code>. Independent confirmation of the attributed gene</td></tr>
<tr><td><code>NumberSubmitters</code></td><td>Not used &mdash; this report counts distinct submitters from <code>submission_summary</code> instead, since it needs them per (gene, disease) rather than per variant</td></tr>
</tbody>
</table>
</div>
<div class="summary-box">
  <strong>Known rough edges in this file.</strong>
  <code>GeneID = -1</code> on {vs_f['geneid_minus1']:,} of the selected rows means ClinVar declines to attribute a
  gene; the ingest emits no gene edge for those and deliberately does <strong>not</strong> fall back to
  <code>GENEINFO</code> ({ga['vs_minus1']:,} of them are variants present in the VCF).
  <code>GeneSymbol</code> degrades to unparseable prose on {vs_f['symbol_prose']:,} rows &mdash; large
  CNVs get text like "covers 42 genes, none of which curated" instead of a symbol, which is one reason
  section {S.structural_variants}'s CNV analysis can resolve so few of them to a named gene.
</div>
"""

    filter_cube_json = json.dumps(filter_cube)
    ingest_compare_html = f"""{section_heading("ingest-compare")}
<p class="subtitle">
  What the ingest emits, read directly from the artifacts in <code>output/</code> rather than
  written into this page, so the two cannot disagree. Structural variants and multi-gene events
  are out of scope here &mdash; see section {S.structural_variants}.
</p>

<div class="total-boxes">
  <div class="total-box"><div class="n">{em['nodes']:,}</div><div class="label"><code>SequenceVariant</code> nodes</div></div>
  <div class="total-box"><div class="n">{em['edges']:,}</div><div class="label">edges</div></div>
  <div class="total-box"><div class="n">{em['genes']:,}</div><div class="label">distinct genes reached</div></div>
  <div class="total-box"><div class="n">{em['diseases']:,}</div><div class="label">distinct diseases reached</div></div>
</div>
<div class="table-wrap">
<table>
<thead><tr><th>Predicate</th><th>Edges</th><th>Meaning</th></tr></thead>
<tbody>
<tr><td><code>biolink:is_sequence_variant_of</code></td>
    <td class="num">{em_gene:,}</td>
    <td>variant &rarr; the single gene ClinVar attributes it to</td></tr>
<tr><td><code>biolink:causes</code></td>
    <td class="num">{em_causes:,}</td>
    <td>variant &rarr; disease, corroborated at &ge;2&#9733; (ClinVar's aggregate, concordant
        submitters, or an expert panel)</td></tr>
<tr><td><code>biolink:associated_with_increased_likelihood_of</code></td>
    <td class="num">{em_assoc:,}</td>
    <td>variant &rarr; disease, &le;1&#9733; but with published support</td></tr>
</tbody>
</table>
</div>
<p class="subtitle">
  Variant classes emitted: {em_classes}. Exactly one predicate per (variant, disease) and one
  gene per variant, so nothing in the graph asserts two things about the same pair.
</p>

<h3>Filter parameters</h3>
<p class="subtitle">
  The three inclusion paths behind <code>causes</code> are a disjunction &mdash; a
  (variant, disease) is kept if <strong>any</strong> threshold is met. Set each below to see the
  population that results. Both counts are exact under any combination, not interpolated.
  <strong>Off</strong> disables that path entirely. The second metric is distinct
  <strong>(variant, disease) pairs</strong>, not emitted edges.
</p>
<div class="xf-panels" style="grid-template-columns: repeat(3, 1fr);">
  <div class="xf-panel">
    <h3>Per-record review status</h3>
    <div id="fp-star-rows"></div>
    <p style="font-size:11.5px; color:#64748b; margin:0.6rem 0 0;">
      Stars on an individual submission record. 2&#9733; never appears per-record, so &ge;2 and
      &ge;3 are identical here &mdash; see section {S.two_star}.
    </p>
  </div>
  <div class="xf-panel">
    <h3>Concordant submitters</h3>
    <div id="fp-conc-rows"></div>
    <p style="font-size:11.5px; color:#64748b; margin:0.6rem 0 0;">
      Distinct submitters agreeing on the same disease with the <em>exact same</em>
      <code>ClinicalSignificance</code>, within one variant.
    </p>
  </div>
  <div class="xf-panel">
    <h3>Aggregate review status</h3>
    <div id="fp-agg-rows"></div>
    <p style="font-size:11.5px; color:#64748b; margin:0.6rem 0 0;">
      ClinVar's variant-level <code>CLNREVSTAT</code>, computed across submitters.
    </p>
  </div>
</div>
<div class="controls" style="margin-top:1rem;">
  <label style="font-size:12.5px; color:#334155; display:inline-flex; align-items:center; gap:6px;">
    <input id="fp-gene-only" type="checkbox"> count only variants with a ClinVar-asserted gene
  </label>
</div>
<div class="total-boxes">
  <div class="total-box"><div class="n" id="fp-variants">-</div><div class="label">variants selected</div></div>
  <div class="total-box"><div class="n" id="fp-edges">-</div><div class="label">distinct (variant, disease) pairs</div></div>
</div>
<div class="table-wrap">
<table>
<thead><tr><th>Combination</th><th>Per-record</th><th>Concordance</th><th>Aggregate</th><th>Variants</th><th>(variant, disease) pairs</th></tr></thead>
<tbody id="fp-preset-rows"></tbody>
</table>
</div>
<p class="subtitle" style="font-size:12.5px;">
  This panel models the star / concordance / aggregate disjunction exactly, but not the
  CLNDISDB-echo gate (stage 5 in <code>FILTERING.md</code>), which is all-or-nothing per variant
  and cannot be folded into a static cube. That gate is very nearly inert, so the counts run
  marginally high: the combination matching the shipping filter predicts
  <span id="fp-selfcheck">-</span> against {em['nodes']:,} nodes actually emitted &mdash; though
  the emitted figure is additionally narrowed by the variant-class prune and widened by the
  &le;1&#9733; publication tier, neither of which this panel models.
</p>

"""

    # --- phenotype terms -----------------------------------------------------
    pp = phenotype_profile
    pp_hist_all = pp["hist_all"]
    pp_total_all = sum(pp_hist_all.values())
    pp_with_hpo = pp_total_all - pp_hist_all[0]

    def _pheno_hist_row(k: int) -> str:
        cls = " class='prod'" if k == 0 else ""
        label = f"{k}+" if k == PHENO_HIST_CAP else str(k)
        pct = 100 * pp_hist_all[k] / pp_total_all
        return (
            f"<tr{cls}><td>{label}</td><td class='num'>{pp_hist_all[k]:,}</td>"
            f"<td class='num'>{pct:.2f}%</td>"
            f"<td class='num'>{pp['hist_plp'].get(k, 0):,}</td></tr>"
        )

    pp_hist_rows = "".join(_pheno_hist_row(k) for k in sorted(pp_hist_all))
    pp_consist = [r["consistency"] for r in pp["rows"]]
    pp_n_dis = len(pp_consist)
    pp_mean = sum(pp_consist) / pp_n_dis if pp_n_dis else 0
    pp_perfect = sum(1 for c in pp_consist if c == 1.0)
    pp_perfect_pct = 100 * pp_perfect / pp_n_dis if pp_n_dis else 0
    pp_none_pct = 100 * pp_hist_all[0] / pp_total_all
    pp_with_pct = 100 * pp_with_hpo / pp_total_all
    pp_top_rows = "".join(
        f"<tr><td><code>{h}</code></td><td>{lbl}</td><td class='num'>{n:,}</td></tr>"
        for h, lbl, n in pp["top_hp"]
    )
    pp_obsolete_rows = "".join(
        f"<tr><td><code>{h}</code></td><td class='num'>{n:,}</td></tr>" for h, n in pp["obsolete"]
    )
    pp_rows_json = json.dumps(pp["rows"][:3000])

    phenotype_terms_html = f"""{section_heading("phenotype-terms")}
<div class="summary-box" style="border-left:4px solid #b91c1c;">
  <strong>&#9888; ClinVar does not record observed phenotypes per variant.</strong> HPO ids appear only
  inside <code>CLNDISDB</code> groups, alongside the MONDO / MedGen / OMIM / Orphanet ids for the
  <em>same</em> condition &mdash; cross-references naming that condition in HPO's vocabulary. All
  214,880 HPO-bearing <code>CLNDISDB</code> groups also carry a disease id; not one is HPO-only.
  <br><br>
  <strong>They are not supplied by submitters.</strong> Only 12 of {sf['n_records']:,} submission records put an
  HPO id in <code>ReportedPhenotypeInfo</code>, and 12,423 (0.195%) in
  <code>SubmittedPhenotypeInfo</code>. The chain is: a submitter names a condition &rarr; ClinVar
  normalises it to a MedGen concept &rarr; MedGen's cross-reference table (19,561 HPO entries in
  <code>MedGenIDMappings.txt.gz</code>) supplies the HPO id &rarr; the VCF prints it in
  <code>CLNDISDB</code>. <strong>The HPO term is a function of the disease concept, not of the variant
  or the patient</strong>, so it cannot carry variant-specific information.
</div>

<h3>The numbers</h3>
<div class="total-boxes">
  <div class="total-box"><div class="n">{pp_none_pct:.1f}%</div><div class="label">of VCF variants carry no HPO term at all ({pp_hist_all[0]:,} of {pp_total_all:,})</div></div>
  <div class="total-box"><div class="n">{pp_mean:.4f}</div><div class="label">mean HPO-set consistency across {pp_n_dis:,} diseases with &ge;2 HPO-bearing variants</div></div>
  <div class="total-box"><div class="n">{pp_perfect:,} / {pp_n_dis:,}</div><div class="label">diseases where <strong>every</strong> variant carries an identical HPO set</div></div>
</div>
<p class="subtitle">
  Only {pp_with_hpo:,} variants ({pp_with_pct:.1f}%) carry any HPO term, drawn from just
  {pp['n_distinct_hp']:,} distinct terms. <strong>{pp_perfect_pct:.1f}% of diseases show zero
  variation between their variants</strong> &mdash; inherited retinal dystrophy resolves to one
  HPO set across 4,052 variants, primary ciliary dyskinesia across 2,881, and in each case that
  set is the disease's own name in HPO (<code>HP:0012265 Ciliary dyskinesia</code>,
  <code>HP:0001639 Hypertrophic cardiomyopathy</code>). The terms are disease-determined, so a
  <code>VariantToPhenotypicFeatureAssociation</code> built from them restated the variant's
  disease edge in a second vocabulary rather than asserting a distinct phenotype. They are no
  longer emitted.
</p>

<h3>Most-used terms, and the obsolete ones</h3>
<div class="xf-panels" style="grid-template-columns: 1fr 1fr;">
  <div class="xf-panel">
    <div class="table-wrap">
    <table><thead><tr><th>HPO term</th><th>Label</th><th>Variants</th></tr></thead>
    <tbody>{pp_top_rows}</tbody></table>
    </div>
  </div>
  <div class="xf-panel">
    <p class="subtitle" style="margin-top:0;">
      {pp['n_obsolete']:,} terms HPO has retired are still referenced, across
      {pp['obsolete_variants']:,} variant mentions. They arrive through MedGen cross-references that
      were never updated &mdash; the same class of problem as the deprecated MONDO terms in section {S.crossfilter}.
    </p>
    <div class="table-wrap">
    <table><thead><tr><th>Obsolete HPO term</th><th>Variants</th></tr></thead>
    <tbody>{pp_obsolete_rows}</tbody></table>
    </div>
  </div>
</div>
"""

    scope_decisions_html = f"""{section_heading("scope-decisions")}
<p class="subtitle">
  The decisions below are applied before any of the analysis that follows, so every count in
  later sections is already inside this scope. They are stated here rather than discovered
  later because each one silently removes a large part of ClinVar.
</p>
<div class="table-wrap">
<table>
<thead><tr><th>Decision</th><th>Kept</th><th>Excluded</th><th>Where</th></tr></thead>
<tbody>
<tr>
  <td><strong>Clinical significance</strong></td>
  <td>Pathogenic, Pathogenic/Likely pathogenic, Likely pathogenic (and the low-penetrance
      form of each)</td>
  <td class="no">VUS, Benign, Likely benign, Conflicting, drug response, and everything else
      &mdash; 92% of ClinVar by variant count</td>
  <td><code>predicate_map</code></td>
</tr>
<tr>
  <td><strong>Variant class</strong></td>
  <td>SNV, Deletion, Duplication, Indel, Insertion</td>
  <td class="no">Microsatellite, Inversion and the catch-all "Variation" class. Repeat
      expansions and inversions are not well described by a fixed REF/ALT, and "Variation"
      carries no class at all</td>
  <td><code>KEPT_VARIANT_CLASSES</code></td>
</tr>
<tr>
  <td><strong>Structural variants</strong></td>
  <td>&mdash;</td>
  <td class="no">All CNVs, translocations and large rearrangements. Not a decision so much
      as a property of the input: <code>clinvar.vcf.gz</code> cannot represent them
      (section {S.structural_variants})</td>
  <td>input file</td>
</tr>
<tr>
  <td><strong>Genome build</strong></td>
  <td>GRCh38 where present, GRCh37 where it is not</td>
  <td class="no">Neither counted twice, nor GRCh37-only variants discarded</td>
  <td><code>ASSEMBLY_PREFERENCE</code></td>
</tr>
<tr>
  <td><strong>Gene attribution</strong></td>
  <td>The single gene ClinVar asserts for the variant</td>
  <td class="no">Every other locus overlapping the position &mdash; antisense transcripts,
      readthrough fusions, <code>LOC</code> placeholders</td>
  <td><code>make_variant_gene_map()</code></td>
</tr>
<tr>
  <td><strong>Phenotype edges</strong></td>
  <td>&mdash;</td>
  <td class="no">Not emitted. ClinVar's HPO ids name the same condition as the disease id
      beside them (section {S.phenotype_terms})</td>
  <td><code>process_row()</code></td>
</tr>
</tbody>
</table>
</div>
<p class="subtitle">
  Two of these are worth restating because they bound everything downstream. The
  significance filter means <strong>this ingest sees only the pathogenic tail of ClinVar</strong>
  &mdash; the VUS population, which is the majority of the archive and the part most in need
  of interpretation, is out of scope by construction. And the variant-class filter means the
  KG's ClinVar contribution is <strong>small variants only</strong>; a gene whose disease
  evidence is entirely structural is invisible here regardless of how strong that evidence is.
</p>
"""

    two_star_html = f"""{section_heading("two-star")}
<p class="subtitle">
  "2 star" is used to mean four different things in this report, and they are not
  interchangeable. Three of them exist in ClinVar; the fourth is one this ingest could compute
  but does not. Each answers a different question about corroboration.
</p>
<div class="table-wrap">
<table>
<thead><tr>
  <th>Kind</th><th>Where it lives</th><th>What it asserts</th><th>Variants / pairs</th>
</tr></thead>
<tbody>
<tr class="prod">
  <td><strong>1. ClinVar's aggregate</strong><br><span style="color:#64748b; font-size:11.5px;">used by this ingest</span></td>
  <td>the VCF's <code>CLNREVSTAT</code></td>
  <td>&ge;2 submitters applied assertion criteria and do not conflict &mdash; computed by
      ClinVar across a variant's submissions</td>
  <td class="num">662,545 variants</td>
</tr>
<tr>
  <td><strong>2. The same value, restated</strong></td>
  <td><code>variant_summary.txt</code>'s <code>ReviewStatus</code></td>
  <td>identical aggregate, published in a second file. Not independent evidence &mdash; the
      same computation</td>
  <td class="num">662,884 variants</td>
</tr>
<tr>
  <td><strong>3. Absent per-record</strong></td>
  <td><code>submission_summary.txt</code>'s <code>ReviewStatus</code></td>
  <td><em>nothing</em>. A single submission cannot be "multiple submitters", so the tier can
      never appear here. Only 6 of the 10 documented values occur, and this is not one</td>
  <td class="num no">0 records</td>
</tr>
<tr class="prod">
  <td><strong>4. Inferred per variant</strong><br><span style="color:#64748b; font-size:11.5px;">used by this ingest</span></td>
  <td>computed &mdash; <code>concordant_disease_pairs()</code></td>
  <td>&ge;{MIN_CONCORDANT_SUBMITTERS} distinct submitters independently assert the same
      disease with the same classification, <em>on one variant</em></td>
  <td class="num">51,185 variants</td>
</tr>
<tr>
  <td><strong>5. Inferred across variants</strong><br><span style="color:#b91c1c; font-size:11.5px;">not implemented</span></td>
  <td>computable &mdash; not currently computed</td>
  <td>several independent <em>variants</em> tie the same gene to the same disease. This holds
      even if every variant came from one submitter, and it is evidence none of the above
      can see, because all four are per-variant</td>
  <td class="num">{mc['multi_variant_pairs']:,} pairs</td>
</tr>
</tbody>
</table>
</div>
<div class="summary-box">
  <strong>Why kind 5 is the interesting one.</strong> Kinds 1&ndash;4 all ask "how many labs
  looked at <em>this variant</em>". Kind 5 asks "how many independent variants implicate this
  gene in this disease" &mdash; a different axis entirely, and the one that matches how a
  gene&ndash;disease relationship actually accumulates evidence.
  <br><br>
  {mc['multi_variant_pairs']:,} gene&ndash;disease pairs are carried by more than one variant.
  Of those, <strong>{mc['multi_variant_single_submitter']:,} have only a single submitter
  behind them</strong> &mdash; invisible to every multi-submitter test, yet supported by
  repeated independent observations of different variants in the same gene. SCN1A &rarr;
  early-infantile DEE (section {S.monarch_kg}) is exactly this shape: many variants, one lab.
  <br><br>
  Whether repeated variants from one lab should count as corroboration is a judgement call
  &mdash; they share that lab's methods and biases. But it is a different question from the one
  ClinVar's stars answer, and currently nothing in the pipeline asks it.
</div>
"""

    pairing_html = f"""{section_heading("pairing")}
<p class="subtitle">
  This ingest emits <strong>variant</strong>:disease edges. The Monarch KG mostly reasons about
  <strong>gene</strong>:disease relationships. The two are related by the variant&rarr;gene edge,
  but they are not the same statement, and most of the confusion in this report traces back to
  conflating them.
</p>
<div class="table-wrap">
<table>
<thead><tr><th></th><th>variant:disease</th><th>gene:disease</th></tr></thead>
<tbody>
<tr><td><strong>Emitted by this ingest?</strong></td>
    <td>yes &mdash; <code>VariantToDiseaseAssociation</code></td>
    <td class="no">no &mdash; derived downstream by joining through the gene edge</td></tr>
<tr><td><strong>What ClinVar asserts</strong></td>
    <td>directly: a submitter classified this variant for this condition</td>
    <td>indirectly: only via which gene the variant sits in</td></tr>
<tr><td><strong>Evidence unit</strong></td>
    <td>submission records on one variant</td>
    <td>every variant in the gene, pooled</td></tr>
<tr><td><strong>How this ingest decides inclusion</strong></td>
    <td class="prod">per variant &mdash; stars, concordance or aggregate</td>
    <td class="no">not decided at all; a pair exists if any one of its variants was admitted</td></tr>
<tr><td><strong>Count in this release</strong></td>
    <td class="num">193,568 (variant, disease)</td>
    <td class="num">{mc['n_clinvar']:,} (gene, disease)</td></tr>
</tbody>
</table>
</div>
<div class="summary-box" style="border-left:4px solid #b91c1c;">
  <strong>&#9888; The consequence: gene:disease pairs inherit a decision made about a single
  variant.</strong> A pair enters the KG the moment <em>one</em> of its variants clears the bar.
  Nothing ever asks whether the pair as a whole is well supported, so:
  <ul style="margin:0.5rem 0 0 1.1rem; padding:0;">
    <li>A pair backed by one expert-panel variant is admitted on the strength of that variant,
        even if its other 400 variants are single-submitter noise.</li>
    <li>A pair backed by over a thousand variants from one lab is rejected, because no single
        variant among them clears the bar (SCN1A &rarr; early-infantile DEE, section
        {S.monarch_kg}).</li>
    <li>Evidence spread thinly across many variants &mdash; the normal shape of a
        well-studied gene &mdash; never accumulates into anything.</li>
  </ul>
  <br>
  Section {S.ingest_recommendation} proposes deciding inclusion at the pair level instead, using
  the pooled submitter count that <code>compute_star_data()</code> already computes and currently
  discards.
</div>
<p class="subtitle">
  One further asymmetry worth naming: the variant&rarr;gene edge is <strong>one gene per
  variant</strong> (ClinVar's own attribution), so a variant cannot contribute to two genes'
  pairs. Before that fix it could, which is how <code>CFTR-AS1</code> acquired cystic fibrosis
  &mdash; the antisense transcript inherited every disease of the variant that overlapped it.
  Gene:disease pairing is only as sound as the gene attribution beneath it.
</p>
"""

    # --- evidence tiers ------------------------------------------------------
    et = evidence_tiers
    ET_STAGES = [
        (1, "1. All variants in <code>clinvar.vcf.gz</code> (SNVs + indels)", ""),
        (2, "2. &hellip;with &ge;1 Pathogenic / Likely pathogenic / P-LP record", ""),
        (3, "3. &hellip;ClinVar aggregate <code>CLNREVSTAT</code> &ge;2&#9733;", "causes"),
        (4, "4. &hellip;+ inferred 2&#9733; (&ge;2 concordant submitters)", "causes"),
        (5, "5. &hellip;+ 1&#9733; with a PubMed citation", "associated_with"),
    ]

    def _et_cell(v: int) -> str:
        return f"{v:,}" if v else "<span class='no'>&mdash;</span>"

    def _et_row(k: int, label: str, pred: str) -> str:
        cls = " class='prod'" if k == 4 else ""
        pred_cell = f"<code>{pred}</code>" if pred else ""
        return (
            f"<tr{cls}><td>{label}</td>"
            f"<td class='num'>{et['stage_variants'][k]:,}</td>"
            f"<td class='num'>{_et_cell(et['stage_genes'][k])}</td>"
            f"<td class='num'>{_et_cell(et['stage_diseases'][k])}</td>"
            f"<td>{pred_cell}</td></tr>"
        )

    et_rows = "".join(_et_row(k, label, pred) for k, label, pred in ET_STAGES)
    et_s4, et_s5 = et["stage_variants"][4], et["stage_variants"][5]
    et_s2 = et["stage_variants"][2]
    et_pub_pct = 100 * et_s5 / et_s2 if et_s2 else 0
    evidence_tiers_html = f"""{section_heading("evidence-tiers")}
<p class="subtitle">
  An exploration, not the current ingest: what each tier of evidence would let us associate, and what
  happens if the predicate is chosen by <em>evidence strength</em> rather than by
  <code>ClinicalSignificance</code>. The proposal tested here is
  <strong>&ge;2&#9733; (real or inferred) &rarr; <code>causes</code></strong> for Pathogenic,
  Pathogenic/Likely pathogenic and Likely pathogenic alike, and
  <strong>1&#9733; with a supporting publication &rarr;
  <code>associated_with_increased_likelihood_of</code></strong>.
</p>

<h3>What each tier reaches</h3>
<p class="subtitle">
  Cumulative, over every small variant in the VCF. Structural variants are absent from this file
  entirely (section {S.biolink_proposal}), so this is the whole SNV/indel database. The highlighted row is what the
  ingest emits today.
</p>
<div class="table-wrap">
<table>
<thead><tr><th>Stage</th><th>Variants</th><th>Distinct genes</th><th>Distinct diseases</th><th>Predicate</th></tr></thead>
<tbody>{et_rows}</tbody>
</table>
</div>

<h3>Gene&ndash;disease pairs by proposed predicate</h3>
<div class="total-boxes">
  <div class="total-box"><div class="n">{et['pairs_causes']:,}</div><div class="label">pairs at <code>causes</code> &mdash; &ge;2&#9733; real or inferred</div></div>
  <div class="total-box"><div class="n">{et['pairs_assoc']:,}</div><div class="label">pairs at <code>associated_with</code> &mdash; 1&#9733; + publication only</div></div>
  <div class="total-box"><div class="n">{et['pairs_causes'] + et['pairs_assoc']:,}</div><div class="label">total pairs under this scheme</div></div>
</div>
<p class="subtitle">
  The two columns are <strong>disjoint</strong>: a pair reachable at &ge;2&#9733; counts only at
  <code>causes</code>, and <code>associated_with</code> holds only pairs with no route to the stronger
  tier. Counting both ways would double-count the same relationship.
</p>

<div class="summary-box" style="border-left:4px solid #b91c1c;">
  <strong>&#9888; The publication signal is far weaker than it sounds.</strong>
  {et['n_pubmed']:,} variants &mdash; 55% of the whole VCF &mdash; carry at least one PubMed citation
  in <code>var_citations.txt</code>. Large cohort papers cite thousands of variants at once, and a
  citation asserts nothing about pathogenicity; the paper may well be reporting the variant as benign.
  So stage 5 admits <strong>{et_s5:,} variants, {et_pub_pct:.0f}% of every variant with any P/LP
  record at all</strong> &mdash; it is close to "accept the P/LP population" rather than a targeted
  rescue.
  <br><br>
  A much tighter alternative sits in <code>submission_summary.txt</code>:
  <code>CollectionMethod == "literature only"</code>, a submission actually <em>derived from</em>
  published evidence and tied to the pathogenicity assertion itself. Only
  <strong>{et['n_lit_only']:,} variants</strong> have one &mdash; two orders of magnitude fewer, and
  the signal points at the claim rather than merely at the variant's existence in the literature. If
  this tier is adopted, that is the better discriminator.
</div>

"""

    # --- Monarch KG reconciliation -----------------------------------------
    mk_source_rows = "".join(
        f"<tr><td><code>{src}</code></td><td class='num'>{n:,}</td></tr>"
        for src, n in sorted(mc["by_source"].items(), key=lambda kv: -kv[1])
    )

    def support_table(xtab: Counter, tier_labels: list) -> str:
        head = "".join(f"<th>{b}</th>" for b in SUPPORT_BUCKETS)
        body = []
        col_totals = Counter()
        for tier, label in tier_labels:
            cells = []
            row_total = 0
            for b in SUPPORT_BUCKETS:
                n = xtab[(tier, b)]
                row_total += n
                col_totals[b] += n
                cells.append(f"<td class='num'>{n:,}</td>")
            body.append(f"<tr><td>{label}</td>{''.join(cells)}<td class='num'><strong>{row_total:,}</strong></td></tr>")
        totals = "".join(f"<td class='num'><strong>{col_totals[b]:,}</strong></td>" for b in SUPPORT_BUCKETS)
        grand = sum(col_totals.values())
        body.append(f"<tr><td><strong>total</strong></td>{totals}<td class='num'><strong>{grand:,}</strong></td></tr>")
        return (
            f"<div class='table-wrap'><table><thead><tr><th>ClinVar evidence tier</th>{head}"
            f"<th>total</th></tr></thead><tbody>{''.join(body)}</tbody></table></div>"
        )

    TIER_LABELS = [
        (1, "Tier 1 &mdash; raw &ge;3&#9733;"),
        (2, "Tier 2 &mdash; concordance rescue"),
        (3, "Tier 3 &mdash; below threshold (not ingested)"),
    ]
    monarch_support_table = support_table(mc["support_xtab"], TIER_LABELS)
    monarch_only_table = support_table(mc["only_xtab"], TIER_LABELS)

    KINSHIP_ORDER = [
        ("ancestor in Monarch", "Monarch links this gene to a <em>broader</em> term (ClinVar is more specific)"),
        ("descendant in Monarch", "Monarch links this gene to a <em>narrower</em> term (ClinVar is more general)"),
        ("gene in Monarch, unrelated term", "Gene is in Monarch, but with no hierarchically related disease"),
        ("gene absent from Monarch", "Gene has no gene-disease association in Monarch at all"),
    ]
    kinship_total = sum(mc["kinship_counts"].values()) or 1
    kinship_rows = "".join(
        f"<tr><td>{desc}</td><td class='num'>{mc['kinship_counts'][k]:,}</td>"
        f"<td class='num'>{100 * mc['kinship_counts'][k] / kinship_total:.1f}%</td></tr>"
        for k, desc in KINSHIP_ORDER
    )
    reconcilable = mc["kinship_counts"]["ancestor in Monarch"] + mc["kinship_counts"]["descendant in Monarch"]
    reconcilable_pct = 100 * reconcilable / kinship_total

    # the collision recommendation 2 leads with, rendered from the live scan
    sssom_dravet_html = ", ".join(
        f"<code>{m}</code>" for m in mc["sssom"]["dravet"]
    ) or "no collision found in this release"

    deprecated_rows_html = "".join(
        f"<tr><td>{r['gene']}</td><td class='mono'>{r['mondo']}</td><td>{r['disease_name']}</td>"
        f"<td class='num'>{r['n_variants']:,}</td><td class='num'>{r['n_submitters']:,}</td></tr>"
        for r in mc["deprecated_rows"][:25]
    )
    n_deprecated = len(mc["deprecated_rows"])

    monarch_supported_json = json.dumps(mc["supported_rows"])
    monarch_only_json = json.dumps(mc["only_rows"])
    monarch_kg_only_json = json.dumps(mc["monarch_only_rows"])
    msc = mc["mondo_status_counts"]
    mondo_status_rows = "".join(
        f"<tr><td>{name}</td><td class='num'>{msc[name]:,}</td></tr>"
        for name in ("current in MONDO", "deprecated in MONDO", "not in MONDO")
    )
    n_not_current = msc["deprecated in MONDO"] + msc["not in MONDO"]

    def monarch_cell(r: dict) -> str:
        return f"yes &mdash; {r['sources']}" if r["in_monarch"] else "<span class='no'>no</span>"

    scn1a_rows_html = (
        "<table><thead><tr><th>MONDO term</th><th>ClinVar variants</th><th>ClinVar submitters</th>"
        "<th>ClinVar tier</th><th>In Monarch?</th></tr></thead><tbody>"
        + "".join(
            f"<tr><td><code>{r['mondo']}</code> {r['disease_name']}</td>"
            f"<td class='num'>{r['n_variants']:,}</td>"
            f"<td class='num'>{r['n_submitters']:,}</td>"
            f"<td>{r['tier']}{' (dropped)' if r['tier'] == 3 else ''}</td>"
            f"<td>{monarch_cell(r)}</td></tr>"
            for r in mc["scn1a_rows"]
        )
        + "</tbody></table>"
    )

    monarch_section_html = f"""{section_heading("monarch-kg")}
<p class="subtitle">
  The sections above derive gene-disease pairs from ClinVar submissions. This one asks how those
  line up with what the <a href="https://monarchinitiative.org/" target="_blank" rel="noopener">Monarch
  KG</a> already asserts. <strong>These are not the same kind of statement.</strong> Monarch's
  {mc['n_monarch']:,} gene&ndash;disease edges are curated assertions <em>about the gene</em>, sourced from
  OMIM, Orphanet and ClinGen; ClinVar's {mc['n_clinvar']:,} pairs here are <em>derived</em> by crossing
  (variant&rarr;gene) with (variant&rarr;disease) over individual variant submissions. ClinVar is not
  itself a source of gene-disease edges in the KG &mdash; this ingest emits variant-level edges &mdash;
  so the comparison measures how much independent submitted variant evidence stands behind each curated
  association, and what ClinVar implies that no curator has asserted.
  Matching is on exact <code>(NCBIGene id, MONDO id)</code>; gene ids are crosswalked from the KG's
  HGNC-keyed gene nodes via HGNC's own complete set. Every ClinVar gene here is the single gene ClinVar
  asserts for the variant, not the VCF's positional <code>GENEINFO</code> list &mdash; see section {S.ingest_recommendation}.
</p>

<div class="total-boxes">
  <div class="total-box"><div class="n">{mc['n_both']:,}</div><div class="label">in both (Monarch association with ClinVar variant evidence)</div></div>
  <div class="total-box"><div class="n">{mc['n_clinvar_only']:,}</div><div class="label">in ClinVar, not in Monarch</div></div>
  <div class="total-box"><div class="n">{mc['n_monarch_only']:,}</div><div class="label">in Monarch, no ClinVar P/LP evidence</div></div>
</div>

<p class="subtitle">
  Monarch's gene-disease edges by primary knowledge source:
</p>
<div class="table-wrap">
<table><thead><tr><th>Primary knowledge source</th><th># edges</th></tr></thead>
<tbody>{mk_source_rows}</tbody></table>
</div>

<h3>How well does ClinVar support Monarch's associations?</h3>
<p class="subtitle">
  Both axes describe the <em>ClinVar</em> side of the {mc['n_both']:,} pairs the two sources agree
  on &mdash; the question is how much ClinVar evidence sits behind an association Monarch has
  already curated.
</p>
<div class="table-wrap" style="max-width:760px;">
<table>
<thead><tr><th>Axis</th><th>Means</th></tr></thead>
<tbody>
<tr><td><strong>Rows &mdash; evidence tier</strong></td>
    <td>How the pair qualified. <em>Tier 1</em>: some variant has a raw &ge;3&#9733;
        (expert-panel) submission. <em>Tier 2</em>: no such variant, but one has
        &ge;{MIN_CONCORDANT_SUBMITTERS} submitters agreeing. <em>Tier 3</em>: neither &mdash;
        Monarch asserts it, this ingest contributes no variant evidence for it.</td></tr>
<tr><td><strong>Columns &mdash; pooled submitters</strong></td>
    <td>How many <em>distinct submitting labs</em> have said anything Pathogenic/Likely-pathogenic
        about that gene&ndash;disease pair, counted across all of its variants. "1" means a single
        lab; "10+" means ten or more independent labs.</td></tr>
<tr><td><strong>Each cell</strong></td>
    <td>The number of pairs with that tier and that submitter count. Reading across a row shows
        how broadly supported that tier's pairs are; reading down a column shows how those pairs
        qualified.</td></tr>
</tbody>
</table>
</div>
<p class="subtitle">
  The bottom-left corner is the one to watch: tier 3 with many submitters means Monarch curates
  the association and many labs agree, yet no <em>single variant</em> clears this ingest's bar.
</p>
{monarch_support_table}

<h3>What's in ClinVar but not in Monarch</h3>
<p class="subtitle">
  The {mc['n_clinvar_only']:,} pairs ClinVar implies that have no matching Monarch gene-disease edge,
  by tier and pooled submitter count. Most are thin single-submitter observations, as expected &mdash;
  but the high-submitter cells are candidate gene-disease relationships with substantial independent
  clinical support and no curated equivalent.
</p>
{monarch_only_table}

<h3>Are these really novel, or does Monarch already link the gene to a nearby term?</h3>
<p class="subtitle">
  MONDO carries multiple co-existing terms over one clinical area, so an exact-id mismatch does not
  mean Monarch is silent about that gene. Each ClinVar-only pair is classified below by whether Monarch
  associates the <em>same gene</em> with a term that is an ancestor or a descendant of ClinVar's term.
  <strong>{reconcilable:,} ({reconcilable_pct:.1f}%) sit at a different level of the MONDO hierarchy
  from a term Monarch already links to that gene.</strong>
</p>
<div class="summary-box" style="border-left:4px solid #b91c1c;">
  <strong>&#9888; Read this as a candidate list, not a duplicate count.</strong> Hierarchical proximity
  means two terms are <em>related</em>; it does not mean they are the same disease. The SCN1A example
  below works through two pairs that look identical in this table &mdash; both nearby in the hierarchy,
  both splitting one gene's evidence &mdash; and MONDO holds <em>both</em> apart on purpose, one of them
  at the explicit request of ClinGen's Epilepsy Gene Curation Expert Panel. Collapsing rows on the
  strength of this column alone would destroy deliberate expert distinctions.
</div>
<div class="table-wrap">
<table><thead><tr><th>Relationship to Monarch's terms for the same gene</th><th># pairs</th><th>%</th></tr></thead>
<tbody>{kinship_rows}</tbody></table>
</div>

<h2 style="font-size:1.15rem; margin-top:2rem;">The three categories, as browsable tables</h2>
<p class="subtitle">
  The same three-way split as the boxes at the top of this section, one searchable and paginated table
  each. Every row is loaded client-side &mdash; the search box filters across all listed columns and
  the column headers sort.
</p>

<h3>Table A &mdash; in both ({mc['n_both']:,} pairs)</h3>
<p class="subtitle">
  Monarch associations that ClinVar corroborates, ranked by how many independent submitters stand
  behind them.
</p>
<div class="controls">
  <input id="monarch-supported-search" class="search-box" type="text" placeholder="Filter by gene, MONDO id or disease...">
  <span id="monarch-supported-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead><tr>
  <th data-sort-for="monarch-supported-rows" data-sort-key="gene">Gene &#8597;</th>
  <th data-sort-for="monarch-supported-rows" data-sort-key="mondo">MONDO disease &#8597;</th>
  <th data-sort-for="monarch-supported-rows" data-sort-key="disease_name">Disease name &#8597;</th>
  <th data-sort-for="monarch-supported-rows" data-sort-key="n_variants"># variants &#8597;</th>
  <th data-sort-for="monarch-supported-rows" data-sort-key="n_submitters"># submitters &#8597;</th>
  <th data-sort-for="monarch-supported-rows" data-sort-key="tier">Tier &#8597;</th>
</tr></thead>
<tbody id="monarch-supported-rows"></tbody>
</table>
</div>
<div class="pagination">
  <button id="monarch-supported-prev">&larr; Prev</button>
  <span id="monarch-supported-page-info" class="page-info"></span>
  <button id="monarch-supported-next">Next &rarr;</button>
</div>

<h3>Table B &mdash; in ClinVar only ({mc['n_clinvar_only']:,} pairs)</h3>
<p class="subtitle">
  Pairs ClinVar implies with no matching Monarch gene-disease edge, sorted by pooled submitter count.
  "Kinship" is the hierarchy classification above &mdash; rows reading "ancestor"/"descendant" name a
  gene Monarch already associates with a hierarchically related term, which makes them worth a curator's
  eye rather than automatically redundant; "gene absent from Monarch" rows are the uncurated
  ones. All {mc['n_clinvar_only']:,} rows are loaded &mdash; search matches gene, MONDO id, disease
  name and kinship. Every gene here is one ClinVar actually asserts for the contributing variants;
  antisense transcripts, readthrough fusions and overlapping loci are gone (see section {S.ingest_recommendation}).
</p>
<div class="controls">
  <input id="monarch-only-search" class="search-box" type="text" placeholder="Filter by gene, MONDO id or disease...">
  <span id="monarch-only-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead><tr>
  <th data-sort-for="monarch-only-rows" data-sort-key="gene">Gene &#8597;</th>
  <th data-sort-for="monarch-only-rows" data-sort-key="mondo">MONDO disease &#8597;</th>
  <th data-sort-for="monarch-only-rows" data-sort-key="disease_name">Disease name &#8597;</th>
  <th data-sort-for="monarch-only-rows" data-sort-key="n_variants"># variants &#8597;</th>
  <th data-sort-for="monarch-only-rows" data-sort-key="n_submitters"># submitters &#8597;</th>
  <th data-sort-for="monarch-only-rows" data-sort-key="tier">Tier &#8597;</th>
  <th data-sort-for="monarch-only-rows" data-sort-key="kinship">Kinship &#8597;</th>
</tr></thead>
<tbody id="monarch-only-rows"></tbody>
</table>
</div>
<div class="pagination">
  <button id="monarch-only-prev">&larr; Prev</button>
  <span id="monarch-only-page-info" class="page-info"></span>
  <button id="monarch-only-next">Next &rarr;</button>
</div>

<h3>Table C &mdash; in Monarch only ({mc['n_monarch_only']:,} pairs)</h3>
<p class="subtitle">
  Curated Monarch associations with <strong>no</strong> Pathogenic/Likely-pathogenic ClinVar evidence
  at all &mdash; no variant in <code>clinvar.vcf.gz</code> links that gene to that disease. Mostly
  Orphanet/OMIM assertions for conditions whose molecular evidence predates ClinVar submission, or
  where the causal variants are structural (see section {S.structural_variants}).
  <br><br>
  <strong>"In MONDO" filter.</strong> An association's disease id is not automatically a live MONDO
  term: the KG can carry an edge whose term MONDO has since retired &mdash; Orphanet still asserts
  <code>SCN1A &rarr; MONDO:0011794</code> ("obsolete Dravet syndrome") &mdash; or one whose id has no
  MONDO node in this release at all. <strong>{n_not_current:,} of these
  {mc['n_monarch_only']:,} pairs are not current MONDO terms.</strong> Tick the box to show only those.
</p>
<div class="table-wrap" style="max-width:520px;">
<table><thead><tr><th>MONDO term status</th><th># pairs</th></tr></thead>
<tbody>{mondo_status_rows}</tbody></table>
</div>
<div class="controls">
  <input id="monarch-kg-only-search" class="search-box" type="text" placeholder="Filter by gene, MONDO id, disease or source...">
  <label style="font-size:12px; color:#334155; display:inline-flex; align-items:center; gap:5px;">
    <input id="monarch-kg-only-notcurrent" type="checkbox"> only terms <strong>not current</strong> in MONDO
  </label>
  <span id="monarch-kg-only-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead><tr>
  <th data-sort-for="monarch-kg-only-rows" data-sort-key="gene">Gene &#8597;</th>
  <th data-sort-for="monarch-kg-only-rows" data-sort-key="gene_id">Gene ID &#8597;</th>
  <th data-sort-for="monarch-kg-only-rows" data-sort-key="mondo">MONDO disease &#8597;</th>
  <th data-sort-for="monarch-kg-only-rows" data-sort-key="disease_name">Disease name &#8597;</th>
  <th data-sort-for="monarch-kg-only-rows" data-sort-key="sources">Knowledge source &#8597;</th>
  <th data-sort-for="monarch-kg-only-rows" data-sort-key="mondo_status">In MONDO? &#8597;</th>
</tr></thead>
<tbody id="monarch-kg-only-rows"></tbody>
</table>
</div>
<div class="pagination">
  <button id="monarch-kg-only-prev">&larr; Prev</button>
  <span id="monarch-kg-only-page-info" class="page-info"></span>
  <button id="monarch-kg-only-next">Next &rarr;</button>
</div>

<h3>Worked example: SCN1A, and why the terms don't line up</h3>
<div class="summary-box">
  <p style="margin-top:0;">
    SCN1A is the clearest illustration of why the count above is a candidate list rather than a
    duplicate count. ClinVar's variant evidence and Monarch's curated associations both concern the
    same gene and the same clinical territory &mdash; developmental and epileptic encephalopathy
    &mdash; but they land on <strong>different, deliberately distinct MONDO terms</strong>:
  </p>
  <div class="table-wrap">
  {scn1a_rows_html}
  </div>
  <div class="controls" style="margin-top:0.75rem;">
    <label style="font-size:12.5px; color:#334155; display:inline-flex; align-items:center; gap:6px;">
      <input id="scn1a-toggle-shares" type="checkbox" checked> shared-variant connectors
    </label>
    <label style="font-size:12.5px; color:#334155; display:inline-flex; align-items:center; gap:6px;">
      <input id="scn1a-toggle-xrefs" type="checkbox"> OMIM / Orphanet ids
    </label>
  </div>
  <div style="margin:0.5rem 0 1rem;">{scn1a_diagram}</div>
  <p class="subtitle" style="font-size:12.5px;">
    <strong>Which source vocabulary a term carries decides which node a lab lands on.</strong>
    Three separate MONDO terms each hold a <code>skos:exactMatch</code> to a source term named
    "Dravet syndrome": <code>OMIM:607208</code> resolves to <code>MONDO:0100079</code> (DEE 6A),
    while <code>DOID:0080422</code>, <code>ICD10CM:G40.83</code>, <code>UMLS:C0751122</code> and
    <code>MEDGEN:148243</code> all resolve to <code>MONDO:0100135</code> (Dravet syndrome), and
    <code>Orphanet:33069</code> resolves to <code>MONDO:0011794</code>, which is obsolete.
    Transitivity of <code>exactMatch</code> would say those should be one node. The counts show what
    that costs: <strong>6 variants sit on the OMIM-linked term and 590 on the MedGen-linked one</strong>.
    <br><br>
    <strong>None of that is an accident. MONDO decided it deliberately, and the decision is on the
    record in <a href="https://github.com/monarch-initiative/mondo/issues/745" target="_blank"
    rel="noopener">mondo#745</a>.</strong> ClinGen's Epilepsy Gene Curation Expert Panel asked for
    the split, and gave the reason: <em>"Dravet syndrome is a clinical diagnosis &mdash; most
    individuals do indeed have variants in SCN1A, but not all, which is one of the reasons why we do
    not want it exclusively tied to SCN1A, as is implied by using the term EIEE6."</em> MONDO agreed
    and split the class. <code>MONDO:0100079</code> now carries the note that <em>"in Mondo, DEE6A is
    treated as a distinct class from Dravet syndrome (contrary to OMIM), as not every case of Dravet
    syndrome is caused by a variation in SCN1A"</em>. The axis is <strong>clinical syndrome (Dravet)
    versus SCN1A-caused entity (DEE 6A)</strong> &mdash; not two names for one thing. The same ticket
    also removed DEE 4 (<em>STXBP1</em>) and DEE 19 (<em>GABRA1</em>) from underneath Dravet.
    <br><br>
    <strong>And the OMIM mapping is deliberate too.</strong> It is tempting to read
    <code>OMIM:607208</code> &mdash; whose title includes "DRAVET SYNDROME" &mdash; landing on DEE 6A
    as a mapping bug. It is the opposite: MONDO's lead curator ruled in that thread that
    <em>"OMIM:607208 should be equivalent to MONDO:0100079 (EIEE6), not to MONDO:0011794"</em>, and
    the change was made. EIEE6 was kept as a <em>related</em>, explicitly not exact, synonym of Dravet.
    So the collision this ingest detects is not an error to repair &mdash; it is two vocabularies
    drawing the boundary in different places, on purpose. <em>Orphanet</em> is the one still out of
    step: it treats the two as equivalent, which is why <code>Orphanet:33069</code> lands on the
    obsolete <code>MONDO:0011794</code>, and the thread records an unresolved request to sync.
    <br><br>
    <strong><code>MONDO:0800491</code> is a separate distinction, on different grounds.</strong> That
    term is early-infantile DEE &mdash; onset &le;3 months, abnormal interictal EEG &mdash; and it
    merges <em>two</em> Orphanet concepts (<code>Orphanet:1934</code> and <code>1935</code>). MONDO's
    own definition puts Dravet's onset in the first year, typically 4&ndash;5 months; the Dravet
    Syndrome European Federation put a mean of 5.2 months on the record in the same thread. The
    {shared_headline} variants the two terms share are SCN1A's phenotypic spectrum, not duplicate
    terminology.
    <br><br>
    <strong>So neither pair should be merged, and rolling evidence up a shared ancestor is the wrong
    fix in both cases.</strong> Two terms being hierarchically close is not evidence that they are
    redundant &mdash; here it is the residue of curators drawing a line carefully, at a domain expert
    panel's request. What MONDO does offer is guidance on which side to land on: the DEE 6A note
    recommends <em>"describing the phenotype, ie Dravet syndrome, or the genetic etiology, ie
    SCN1A"</em>, which points this ingest at <code>MONDO:0100135</code> plus the gene edge it already
    emits.
  </p>
  <p>
    <strong>The hierarchy explains it.</strong> <code>MONDO:0800491</code> (early-infantile DEE) and
    <code>MONDO:0100135</code> (Dravet syndrome) are <em>siblings</em> &mdash; both are direct children of
    <code>MONDO:0800490</code> (neonatal/infantile-onset epilepsy syndrome with DEE), and both descend from
    <code>MONDO:0100620</code> (developmental and epileptic encephalopathy). Neither is an ancestor of the
    other, so no amount of exact-id matching will ever unify them &mdash; and per the note above, they
    should not be unified. The hierarchy explains why the evidence splits; it does not license repairing
    the split by merging.
  </p>
  <p>
    <strong>Why ClinVar lands where it does.</strong> The term a variant gets is decided entirely by which
    MedGen concept the submitting lab chose. Labcorp/Invitae uses
    <code>C0393706:Early-infantile DEE</code> &rarr; <code>MONDO:0800491</code>, and is the <em>only</em>
    lab that does &mdash; hence a pooled submitter count of 1 across all its variants. The other ~87 labs use
    <code>C0751122:Severe myoclonic epilepsy in infancy</code> (the former name for Dravet syndrome)
    &rarr; <code>MONDO:0100135</code>. One gene's clinical spectrum, two vocabularies, two MONDO terms,
    and neither side can corroborate the other. The tempting move is to pool them at
    <code>MONDO:0100620</code>, where SCN1A's DEE evidence becomes overwhelming and multi-lab &mdash; but
    that pools terms MONDO separates on clinical grounds. The measurable problem is real; the rollup is
    not the remedy (see recommendation 2 in section {S.ingest_recommendation}).
  </p>
  <p>
    <strong>Monarch's side has the mirror-image problem.</strong> Its 13 SCN1A associations come from three
    curators that also disagree on terms: ClinGen asserts <code>MONDO:0100135</code> (Dravet), OMIM asserts
    the numbered series (<code>MONDO:0030268</code>, <code>MONDO:0100079</code>, <code>MONDO:0011461</code>,
    <code>MONDO:0012320</code>), and Orphanet asserts seven more including
    <code>MONDO:0011794</code> &mdash; which is <strong>flagged deprecated in MONDO</strong>
    ("obsolete Dravet syndrome") yet still carries a live gene-disease edge in the KG. Neither
    <code>MONDO:0800491</code> nor <code>MONDO:0100620</code> &mdash; ClinVar's two largest SCN1A buckets
    &mdash; appears anywhere in Monarch's SCN1A set.
  </p>
</div>

<h3>ClinVar pairs pointing at deprecated MONDO terms</h3>
<p class="subtitle">
  A related, independently fixable problem: <strong>{n_deprecated:,}</strong> of this ingest's
  {mc['n_clinvar']:,} gene-disease pairs resolve to a MONDO term the ontology marks
  <code>deprecated</code>. These come through <code>mondo.sssom.tsv</code> / MedGen mappings that still
  point at obsoleted ids, and nothing in <code>clinvar_helpers.py</code> filters them. Top 25 by variant
  count:
</p>
<div class="table-wrap">
<table>
<thead><tr><th>Gene</th><th>MONDO</th><th>Disease name</th><th># variants</th><th># submitters</th></tr></thead>
<tbody>{deprecated_rows_html}</tbody>
</table>
</div>
"""

    # --- consolidated ingest recommendation ----------------------------------
    all_support_table = support_table(mc["all_xtab"], TIER_LABELS)
    proj2, proj3 = mc["projected"][2], mc["projected"][3]
    gain3 = proj3 - mc["ingested_today"]
    gain2 = proj2 - mc["ingested_today"]

    ga = gene_attribution
    ga_art_pct_before = 100 * ga["geneinfo_artifact"] / max(ga["geneinfo_attributions"], 1)
    ga_art_pct_after = 100 * ga["vs_artifact"] / max(ga["vs_attributions"], 1)
    ga_artifact_reduction = 100 * (1 - ga["vs_artifact"] / max(ga["geneinfo_artifact"], 1))
    ga_removed = ga["geneinfo_attributions"] - ga["vs_attributions"]
    ga_removed_pct = 100 * ga_removed / max(ga["geneinfo_attributions"], 1)
    ga_multi_after = ga["multi_gene_variants"] - ga["multi_resolved_to_one"]
    ga_minus1_pct = 100 * ga["vs_minus1"] / max(ga["vs_attributions"], 1)

    recommendation_html = f"""{section_heading("ingest-recommendation")}
<p class="subtitle">
  Everything above is measurement. This section is the actionable summary: concrete changes to
  <code>src/clinvar_helpers.py</code>, ordered by value-per-unit-risk, each tied to the evidence that
  motivates it.
  <br><br>
  <strong>Two have been applied</strong> and are described at the end of this section: the
  <code>review_star_map</code> star-value correction, and the gene-attribution fix &mdash; the ingest no
  longer derives genes from the VCF's positional <code>GENEINFO</code> list. Every gene-disease pair in
  this report reflects that corrected attribution. <strong>Recommendations 1, 2 and 4 below are not
  implemented</strong>; they remain proposals.
</p>

<div class="summary-box">
  <strong>The core diagnosis.</strong> The stated goal is to capture gene&ndash;disease pairs with
  significant support from ClinVar. The current filter cannot express that, because both of its tests
  &mdash; the star threshold and the concordance rescue &mdash; run <em>inside</em>
  <code>variant_records_to_disease(record_list, ...)</code>, where <code>record_list</code> is a
  <strong>single variant's</strong> submission records. A pair is admitted when some <em>one</em> of its
  variants clears the bar alone, so evidence spread thinly across many variants &mdash; which is what a
  well-studied gene&ndash;disease relationship actually looks like &mdash; never accumulates. The metric
  that does express the goal, <code>pair_submitters</code>, is already computed in
  <code>compute_star_data()</code> and currently feeds only the ClinGen-coverage table.
</div>

<h3>Evidence: pair-level support vs. what the current filter admits</h3>
<p class="subtitle">
  Every ClinVar gene-disease pair by current evidence tier (rows) and pooled distinct submitters
  (columns). The filter is <strong>precise</strong> &mdash; only {mc['thin_ingested']:,} of the
  {mc['ingested_today']:,} pairs it admits rest on a single submitter. It is <strong>not
  sensitive</strong> &mdash; it drops {mc['dropped_ge3']:,} pairs backed by &ge;3 independent labs and
  {mc['dropped_ge2']:,} backed by &ge;2.
</p>
{all_support_table}

<h3>Recommendation 1 &mdash; decide inclusion per pair, not per variant</h3>
<p class="subtitle">
  Replace the per-variant test with: <strong>admit a (gene, disease) pair if any variant carries raw
  &ge;{PRODUCTION_STAR_MIN}&#9733; evidence, OR if &ge;N distinct submitters have made a
  Pathogenic/Likely-pathogenic assertion about that pair across all its variants.</strong> Keeping the
  &ge;{PRODUCTION_STAR_MIN}&#9733; clause unconditional means an expert-panel review is never discarded
  for lack of headcount.
</p>
<div class="table-wrap">
<table>
<thead><tr><th>Rule</th><th>Pairs captured</th><th>vs today</th></tr></thead>
<tbody>
<tr><td>current per-variant rule</td><td class="num">{mc['ingested_today']:,}</td><td class="num">&mdash;</td></tr>
<tr class="prod"><td>&ge;{PRODUCTION_STAR_MIN}&#9733; <strong>or</strong> N = 3 pooled submitters <strong>(recommended)</strong></td><td class="num">{proj3:,}</td><td class="num">+{gain3:,}</td></tr>
<tr><td>&ge;{PRODUCTION_STAR_MIN}&#9733; <strong>or</strong> N = 2 pooled submitters</td><td class="num">{proj2:,}</td><td class="num">+{gain2:,}</td></tr>
</tbody>
</table>
</div>
<p class="subtitle" style="font-size:12.5px;">
  N = 3 is the defensible "significant support" bar and costs nothing in precision. Implementation note:
  this decision cannot live inside <code>variant_records_to_disease()</code>, which only ever sees one
  variant. Pair-level support has to be accumulated across variants first (a pre-pass, exactly as
  <code>compute_star_data()</code> already does) and the admit/reject decision applied in the transform.
</p>

<h3>Recommendation 2 &mdash; report term collisions for curation; do not roll evidence up the hierarchy</h3>
<p class="subtitle">
  <code>concordance_groups()</code> keys agreement on <code>(mondo_id, ClinicalSignificance)</code>, so
  two labs describing one gene's disease under two MONDO terms never corroborate each other. Section
  {S.monarch_kg} measured the scale: <strong>{reconcilable:,} of {kinship_total:,}
  ({reconcilable_pct:.1f}%)</strong> ClinVar-only pairs involve a gene Monarch already links to a
  hierarchically related term. The obvious fix is to count support over the MONDO ancestor closure.
  <strong>That fix is wrong, and SCN1A is the reason.</strong>
</p>
<div class="summary-box" style="border-left:4px solid #b91c1c;">
  <strong>&#9888; The terms that fragment SCN1A's evidence are separated deliberately, at a domain
  expert panel's request.</strong>
  <a href="https://github.com/monarch-initiative/mondo/issues/745" target="_blank"
  rel="noopener">mondo#745</a> records ClinGen's Epilepsy Gene Curation Expert Panel asking MONDO to
  split Dravet syndrome from DEE 6A, because <em>"Dravet syndrome is a clinical diagnosis &mdash; most
  individuals do indeed have variants in SCN1A, but not all"</em>. MONDO agreed and split the class;
  <code>MONDO:0100079</code> still carries the note. <code>MONDO:0800491</code> (early-infantile DEE,
  onset &le;3 months) is separately distinct from Dravet on onset. Pooling at
  <code>MONDO:0100620</code> would silently undo both decisions and make the resulting evidence counts
  look <em>better</em> while meaning less.
</div>
<p class="subtitle">
  <strong>Emit the collisions instead.</strong> Scan <code>mondo.sssom.tsv</code> for source labels
  claimed by more than one MONDO class and publish the list as a QC artifact for curation. This ingest
  already reads that file, so the cost is one pass, and it finds
  <strong>{mc['sssom']['n_collisions']:,} colliding labels</strong> out of {mc['sssom']['n_labels']:,}.
  "dravet syndrome" is one of them, surfaced without being looked for: {sssom_dravet_html}.
  <br><br>
  <strong>Flag, never auto-merge &mdash; and do not assume a collision is a bug.</strong> The SCN1A
  collision is the cautionary case. <code>OMIM:607208</code> is titled "Dravet syndrome" and
  <code>skos:exactMatch</code>es <code>MONDO:0100079</code>, which reads like an obvious mapping error
  &mdash; but MONDO's lead curator ruled in that same thread that <em>"OMIM:607208 should be equivalent
  to MONDO:0100079 (EIEE6), not to MONDO:0011794"</em>. The mapping is intentional; the two
  vocabularies simply draw the boundary in different places. Merging on <code>exactMatch</code> closure
  would therefore have collapsed a distinction MONDO and ClinGen built on purpose.
  <br><br>
  The mappings cannot support more than flagging in any case: {mc['sssom']['exact_pct']:.1f}% of the
  {mc['sssom']['total']:,} are <code>skos:exactMatch</code> and
  <strong>{mc['sssom']['unspecified_pct']:.0f}% carry <code>semapv:UnspecifiedMatching</code></strong>,
  so no mapping records how it was made and there is no provenance to adjudicate a collision with. A
  human reading the ticket is what resolved this one.
  <br><br>
  A smaller, independent fix in the same function: concordance also requires the <em>exact same</em>
  <code>ClinicalSignificance</code> string, so <em>Pathogenic</em> + <em>Likely pathogenic</em> on one
  variant and disease do not currently corroborate each other. Grouping by the P/LP family, which
  <code>predicate_map</code> already defines, would fix that on its own.
</p>

<h3>Gene attribution</h3>
<div class="summary-box">
  <strong>The gene comes from ClinVar's own attribution, not the VCF's <code>GENEINFO</code>.</strong>
  <code>GENEINFO</code> is populated <em>positionally</em> &mdash; it lists every gene whose span
  covers the variant, which is what the field is for. Treating it as a causal claim gives antisense
  transcripts, readthrough fusions, locus control regions and NCBI <code>LOC</code> placeholders
  their own gene edges, from which they inherit the real gene's entire disease and submitter
  roster. <code>variant_summary.txt.gz</code> carries the gene ClinVar assigns, one per variant,
  and that is what <code>make_variant_gene_map()</code> reads.
  <br><br>
  Where ClinVar declines to attribute a gene (<code>GeneID == -1</code>) no gene edge is emitted
  and there is deliberately no <code>GENEINFO</code> fallback &mdash; an unasserted gene is left
  unasserted, and the variant keeps its disease edges. <code>HGNC_ID</code> in the same file is the
  id space the Monarch KG keys genes on, so emitting it would remove the HGNC&harr;Entrez crosswalk
  section {S.monarch_kg} builds.
</div>

<h3>Recommendation 4 &mdash; drop deprecated MONDO terms</h3>
<p class="subtitle">
  <strong>{n_deprecated:,}</strong> of this ingest's {mc['n_clinvar']:,} pairs resolve to a MONDO term the
  ontology marks <code>deprecated</code> (e.g. <code>MONDO:0010086</code> "obsolete sudden infant death
  syndrome"). They arrive through <code>mondo.sssom.tsv</code> / MedGen mappings that still point at
  obsoleted ids, and nothing in <code>clinvar_helpers.py</code> filters them. Smallest and least
  contentious change on this list. Worth noting the same class of problem exists upstream in the KG
  itself &mdash; Orphanet still carries a live SCN1A edge to <code>MONDO:0011794</code>
  ("obsolete Dravet syndrome").
</p>

<h3>Already applied: <code>review_star_map</code> correction</h3>
<p class="subtitle">
  One change from this analysis <em>has</em> been made, because it was an outright error rather than a
  design choice: <code>no_assertion_criteria_provided</code>, <code>no_classification_provided</code> and
  <code>no_classifications_from_unflagged_records</code> were scored 1&#9733; where
  <a href="https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/" target="_blank" rel="noopener">ClinVar
  documents them as 0&#9733;</a> &mdash; contradicting this report's own curation table. Corrected in
  <code>clinvar_helpers.py</code>. <strong>Production output is unchanged</strong>
  (<code>star_min={PRODUCTION_STAR_MIN}</code>, so 0 and 1 both fail identically, and the concordance
  rescue ignores stars); it moves 112,154 variants from 1&#9733; to 0&#9733; in section {S.crossfilter}'s crossfilter
  and fixes the star column in section {S.clinvar_curation}'s table.
</p>
"""

    crossfilter_html = f"""{section_heading("crossfilter")}
<p class="subtitle">
  Built from all {total_variants:,} variants in <code>clinvar.vcf.gz</code>, using ClinVar's own
  aggregate, variant-level fields (CLNSIG / CLNREVSTAT / CLNVC) &mdash; not the per-submission
  evidence used in section {S.star_cutoff}. <code>CLNREVSTAT</code> here is ClinVar's real aggregate
  review status, so it genuinely includes the "2 star" tier that's absent per-record in section {S.star_cutoff}
  &mdash; <strong>two different star systems share the word "stars" in this report</strong>, and
  the star panel below is the variant-level one. Because of that, the star panel alone can't
  express production's inclusion rule; the <strong>"In production ingest"</strong> panel does,
  re-running section {S.star_cutoff}'s exact criterion
  (&ge;{PRODUCTION_STAR_MIN}&#9733; on some individual submission record, or
  &ge;{MIN_CONCORDANT_SUBMITTERS} concordant submitters) per variant, so every panel here can be
  crossfiltered against what the ingest actually keeps.
  Toggle checkboxes in any panel; the other panels re-filter to match, and the totals update live.
  Each panel's own bars stay visible when toggled off (dimmed) so you can see what you're excluding.
  Gene-disease pairs are only ever produced from Pathogenic/Likely-pathogenic evidence (same as
  section {S.star_cutoff}) &mdash; filtering to VUS/Benign/Conflicting/Other/Not-classified alone will
  always show 0 pairs, which is expected: production never creates disease edges for those.
  "Has literature submission" and "Multiple concordant submitters" are per-variant flags computed
  from <code>submission_summary.txt</code> (same source as section {S.star_cutoff}, not the VCF fields
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
    <h3>Aggregate review status &mdash; <code>CLNREVSTAT</code> (stars)</h3>
    <div id="xf-star-rows"></div>
    <div class="xf-actions">
      <button data-dim="star" data-action="all">All</button>
      <button data-dim="star" data-action="none">None</button>
    </div>
    <p style="font-size:11.5px; color:#64748b; margin:0.6rem 0 0;">
      ClinVar's <em>variant-level</em> star rating, scored by the same
      <code>review_star_map</code> section {S.star_cutoff} use. The 2&#9733; tier is real here (and
      large) because it's an aggregate across a variant's submitters &mdash; it never appears on
      an individual submission record, so it can't be selected per-record in section {S.star_cutoff}.
      Ticking &ge;3&#9733; here is <strong>not</strong> the same as production's filter; use the
      "In production ingest" panel for that.
    </p>
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
  <div class="xf-panel">
    <h3>In production ingest (&ge;{PRODUCTION_STAR_MIN}&#9733; per-record or &ge;{MIN_CONCORDANT_SUBMITTERS} concordant)</h3>
    <div id="xf-production-rows"></div>
    <div class="xf-actions">
      <button data-dim="production" data-action="all">All</button>
      <button data-dim="production" data-action="none">None</button>
    </div>
    <p style="font-size:11.5px; color:#64748b; margin:0.6rem 0 0;">
      The actual section {S.star_cutoff} criterion, re-run per variant:
      <code>variant_records_to_disease(star_min={PRODUCTION_STAR_MIN},
      rescue_min_submitters={MIN_CONCORDANT_SUBMITTERS})</code> yields &ge;1 MONDO disease.
      Computed from per-submission <code>ReviewStatus</code>, not from the <code>CLNREVSTAT</code>
      stars in the panel above &mdash; so this is the dimension to tick when you want "what
      production actually ingests" rather than "what ClinVar rates highly". Like sections
      1&ndash;{S.illustrative_examples} it ignores <code>process_row()</code>'s separate HPO-overlap requirement.
      Selecting "Yes" alone reproduces section {S.star_cutoff}'s "2&#9733; (computed)" <em>variant</em> count exactly.
      The gene-disease pair total will read slightly higher: pairs in this section are enumerated from
      every Pathogenic/Likely-pathogenic record on a variant (any star), so a variant that passes on one
      disease also contributes its other diseases here &mdash; section {S.star_cutoff}'s pair count applies the
      threshold per (variant, disease) instead.
    </p>
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
  Genes are the single ClinVar-asserted gene per variant, not the VCF's positional
  <code>GENEINFO</code> list &mdash; see section {S.ingest_recommendation}.
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
"""

    production_pairs_today = results[PRODUCTION_STAR_MIN]["gene_disease_pairs"]
    biolink_proposal_html = f"""{section_heading("biolink-proposal")}
<p class="subtitle">
  Everything above measures how well production's <em>existing</em> model performs. This closing section is
  different in kind: not an analysis of the downloaded data, but a proposal for discussion, informed by what
  that analysis found. It describes the Biolink model this ingest emits today, the gaps sections {S.phenotype_terms}&ndash;{S.ingest_recommendation}
  expose (structural variants excluded entirely; no way to represent per-gene or functional-mechanism
  differences), and a concrete reframing that addresses both &mdash; grounded in Biolink Model classes and
  qualifier slots that already exist, not invented from scratch.
</p>

<h3>Current state</h3>
<p class="subtitle">
  Today, <code>process_row()</code> in <code>src/clinvar_helpers.py</code> emits exactly one node class and
  three association classes, built entirely from <code>clinvar.vcf.gz</code>:
</p>
<p class="subtitle">
  <strong>Every ClinVar node lands on <code>biolink:SequenceVariant</code>, never the more specific
  <code>biolink:Snv</code> subclass</strong> &mdash; even for records that are, molecularly, plain
  single-nucleotide variants. This isn't just true of the ingest code (<code>SequenceVariant(...)</code> is
  the only node class it ever constructs); it's confirmed live in the production Monarch graph: querying
  <a href="https://api-v3.monarchinitiative.org/v3/api/entity/CLINVAR:55630" target="_blank" rel="noopener">
  api-v3.monarchinitiative.org's entity endpoint for CLINVAR:55630</a> (the BRCA1 SNV example above) returns
  <code>"category":"biolink:SequenceVariant"</code>. The distinction between "this is a SNV" and "this is a
  CNV" today lives only in the free-text <code>type</code> slot (Sequence Ontology terms), not in the node's
  Biolink category.
</p>
<div class="table-wrap">
{diagram_sequencevariant_hierarchy()}
</div>
<div class="table-wrap" style="margin-top:1rem;">
{diagram_current_model()}
</div>
<div class="table-wrap" style="margin-top:1rem;">
<table>
<tr><th>Slot</th><th>Populated from</th><th>What it actually captures</th></tr>
<tr><td><code>SequenceVariant.id</code></td><td><code>CLINVAR:{{ClinVar VariationID}}</code></td><td>One node per ClinVar variant passing the production filter (&ge;{PRODUCTION_STAR_MIN}&#9733; or &ge;{MIN_CONCORDANT_SUBMITTERS} concordant submitters, see section {S.star_cutoff})</td></tr>
<tr><td><code>SequenceVariant.type</code></td><td>VCF <code>MC</code> field</td><td>Sequence Ontology <strong>molecular consequence</strong> terms (e.g. missense_variant) &mdash; one list, for the whole variant, regardless of how many genes it touches</td></tr>
<tr><td><code>VariantToGeneAssociation</code></td><td>VCF <code>GENEINFO</code> field</td><td>One identical <code>is_sequence_variant_of</code> edge per gene in <code>GENEINFO</code> &mdash; no per-gene distinction at all</td></tr>
<tr><td><code>VariantToDiseaseAssociation</code></td><td><code>submission_summary.txt</code></td><td><code>causes</code> / <code>associated_with_increased_likelihood_of</code>, qualified only by <code>CLNREVSTAT</code></td></tr>
<tr><td><code>VariantToPhenotypicFeatureAssociation</code></td><td>VCF INFO groups</td><td><code>contributes_to</code>, gated on a disease edge already existing</td></tr>
</table>
</div>
<p class="subtitle">
  With today's production settings this produces <strong>{production_pairs_today:,}</strong> gene-disease
  pairs (star_min&ge;{PRODUCTION_STAR_MIN} column, section {S.star_cutoff}) &mdash; and, as section {S.structural_variants} showed, misses at
  least <strong>{sv_only_pairs_count:,}</strong> more pairs whose only evidence is a CNV, inversion, or other
  structural variant that <code>clinvar.vcf.gz</code> cannot represent at all.
</p>

<h3>The problem, precisely</h3>
<ul style="color:#334155; font-size:14px; line-height:1.7;">
  <li><strong>Variant class and molecular consequence are conflated.</strong> <code>type</code> today holds
  SO terms describing what the variant does at the sequence level (missense, frameshift, ...) &mdash; but
  extending this ingest to CNVs/SVs means <code>type</code> would ALSO need to hold SO terms describing what
  <em>kind</em> of variant it structurally is (SO:0001742 copy_number_gain, SO:0001743 copy_number_loss,
  SO:1000036 inversion, SO:0000289 microsatellite, ...). Both are legitimate uses of the same multivalued
  slot today, which makes "what does <code>type</code> mean" ambiguous the moment a CNV shows up.</li>
  <li><strong>No per-gene differentiation.</strong> A CNV spanning 5 genes gets 5 identical
  <code>is_sequence_variant_of</code> edges. There's no way today to say gene A is completely deleted, gene B
  is truncated at exon 5, and gene C merely sits inside the deleted region's flanking regulatory sequence
  without losing any coding sequence &mdash; three very different biological claims, one indistinguishable
  edge shape.</li>
  <li><strong>No functional/allelic mechanism at all.</strong> The model has no slot for whether a variant
  destroys gene function, reduces it, increases it, actively opposes the wild-type product, or gives the
  gene a new function &mdash; despite this report's own LMNA and SMCHD1 examples showing that mechanism, not
  just "which gene, which disease," is often what actually explains why the same gene produces such
  different phenotypes.</li>
</ul>

<h3>Proposed reframing: three orthogonal axes</h3>
<p class="subtitle">
  Split what <code>type</code> is doing today into three separate questions, asked at the right place in the
  graph:
</p>
<ol style="color:#334155; font-size:14px; line-height:1.8;">
  <li><strong>Variant class</strong> (what kind of DNA change is this?) &mdash; a property of the
  <code>SequenceVariant</code> node itself, independent of any gene. SO term(s): SNV, insertion, deletion,
  microsatellite/STR expansion, copy_number_gain, copy_number_loss, inversion, complex rearrangement.
  This is exactly the <code>Type</code>/<code>CLNVC</code> value ClinVar already assigns every variant,
  SNV or CNV alike &mdash; production just never reads it for anything but VCF rows today.</li>
  <li><strong>Molecular consequence</strong> (what does it do to <em>this specific</em> overlapping gene?)
  &mdash; moves from the node onto each <code>VariantToGeneAssociation</code> <strong>edge</strong>, because
  it can legitimately differ per gene for the same variant (a CNV's breakpoint can be exonic in one gene and
  purely intronic in its neighbor). Still SO terms: missense_variant, stop_gained, frameshift_variant,
  exon_loss_variant, ...</li>
  <li><strong>Functional morph</strong> (what does it do to the gene <em>product's function</em>, relative to
  wild-type?) &mdash; also an edge-level qualifier, using
  <a href="https://en.wikipedia.org/wiki/Muller%27s_morphs" target="_blank" rel="noopener">Muller's 1932
  classification</a>: <strong>amorph</strong> (null/complete loss), <strong>hypomorph</strong> (partial loss),
  <strong>hypermorph</strong> (excess of normal function), <strong>antimorph</strong> (dominant-negative,
  actively opposes the wild-type product), <strong>neomorph</strong> (a genuinely new function). This is an
  assertion about mechanism, not sequence &mdash; it usually has to be curated or inferred, not computed
  directly from HGVS notation.</li>
</ol>
<div class="table-wrap">
{diagram_proposed_model()}
</div>

<h3>Grounding in Biolink Model as it exists today</h3>
<p class="subtitle">
  None of this requires inventing new infrastructure from scratch &mdash; Biolink Model already ships most of
  the pieces, just not wired together for variants this way:
</p>
<div class="table-wrap">
<table>
<tr><th>Proposed slot</th><th>Existing Biolink support</th><th>Gap</th></tr>
<tr>
  <td>Variant class</td>
  <td><code>SequenceVariant.type</code> (already multivalued, already SO-typed)</td>
  <td>None &mdash; just needs a documented convention that this slot means <em>class</em>, not consequence</td>
</tr>
<tr>
  <td>Molecular consequence, per gene</td>
  <td><code>VariantToGeneAssociation.qualifiers</code> (free-form list, already used elsewhere in this ingest
  for <code>VariantToDiseaseAssociation.qualifiers=[CLNREVSTAT]</code>)</td>
  <td>None structurally, but an unconstrained string list doesn't stop two curators from writing the same
  concept two different ways</td>
</tr>
<tr>
  <td>Functional morph: hypomorph / hypermorph</td>
  <td><code>EntityToFeatureOrVariantQualifiersMixin.object_aspect_qualifier</code> +
  <code>object_direction_qualifier</code> &mdash; <code>aspect_qualifier=activity_or_abundance</code> with
  <code>direction_qualifier=decreased</code> (hypomorph) or <code>increased</code> (hypermorph) already
  exists and is used elsewhere in Biolink (e.g. chemical/gene regulatory associations)</td>
  <td>None &mdash; this half of Muller's classification maps onto existing Biolink qualifiers directly</td>
</tr>
<tr>
  <td>Functional morph: antimorph / neomorph</td>
  <td><code>DirectionQualifierEnum</code> today only defines <code>increased</code> / <code>upregulated</code>
  / <code>decreased</code> / <code>downregulated</code></td>
  <td>Real gap. Dominant-negative and novel-function mechanisms aren't a "direction" of the same activity at
  all &mdash; representing them needs either a local enum extension or an upstream Biolink Model proposal</td>
</tr>
<tr>
  <td>Per-gene impact category (whole/partial/regulatory)</td>
  <td><code>GenomicSequenceLocalization</code> already exists as an association class carrying
  <code>start_interbase_coordinate</code> / <code>end_interbase_coordinate</code> / <code>genome_build</code>
  / <code>strand</code> &mdash; i.e. Biolink already has a coordinate-bearing association shape</td>
  <td>Not wired to genes specifically; would need a controlled vocabulary for the impact category itself
  (whole_gene / partial_gene / regulatory_only / unknown) since Biolink doesn't define one</td>
</tr>
</table>
</div>

<h3>Multi-gene variants: worked example from this report</h3>
<p class="subtitle">
  The DMD examples in the mutation-type gallery above are the cleanest illustration: a single deletion or
  duplication of exons 4&ndash;{S.star_cutoff} is <em>whole_gene</em>-scale relative to nothing else nearby, but the exact
  same physical event, if it happened to span a neighboring gene's promoter, would need to be represented as
  <em>regulatory_only</em> for that second gene &mdash; a completely different claim from "this gene lost
  coding sequence," using the identical variant. Section 12's CNV/SV-only pairs table already had to restrict
  itself to single-gene rows precisely because <code>variant_summary.txt</code>'s own <code>GeneID</code>
  column degrades to an unusable <code>-1</code> sentinel the moment more than one gene is involved &mdash;
  ClinVar's own source data doesn't give us clean per-gene impact today, which is exactly the gap
  <code>GenomicSequenceLocalization</code>-style coordinate overlap (computed against a real gene model, the
  same way this report fetched LMNA's and SMCHD1's exon structure from Ensembl) would need to fill.
</p>

<h3>Discussion questions</h3>
<ol style="color:#334155; font-size:14px; line-height:1.8;">
  <li>Should variant <em>class</em> and molecular <em>consequence</em> really share one <code>type</code>
  slot going forward (as they do today), or should class move to a distinct slot even though both are SO
  terms &mdash; and if we split them, which existing consumers of <code>SequenceVariant.type</code> break?</li>
  <li>Molecular consequence and functional morph are proposed as <strong>edge</strong> qualifiers (per gene),
  not node properties &mdash; does that hold up for SNVs too, or only matters once CNVs are in scope? An SNV
  only ever touches one gene today, so would this be over-engineering for the common case?</li>
  <li>Biolink's <code>object_aspect_qualifier</code> + <code>object_direction_qualifier</code> already covers
  hypomorph/hypermorph. For antimorph and neomorph: extend locally with our own enum, or take this upstream
  as a Biolink Model change? What breaks for other consumers of <code>DirectionQualifierEnum</code> if we do?</li>
  <li>Who asserts <code>functional_morph</code>, and how confidently? Is this exactly the kind of call a
  ClinGen VCEP is positioned to make (see the ClinGen-coverage section above) &mdash; or does it need a
  dedicated functional-genomics source (e.g. MGI, FlyBase, or primary literature) that ClinVar itself
  doesn't carry at all?</li>
  <li><code>impact_category</code> (whole_gene / partial_gene / regulatory_only) is proposed as a new,
  locally-invented controlled vocabulary. Does an existing standard already cover this &mdash; ClinGen's own
  CNV interpretation framework, or the ACMG/ClinGen technical standards for constitutional CNV
  classification &mdash; that we should align to instead of inventing our own terms?</li>
  <li>Practically: does ingesting CNVs/SVs need a new node category (e.g. a dedicated structural-variant
  class) alongside <code>SequenceVariant</code>, or can one node class serve both, differentiated only by
  the variant-class value in <code>type</code>?</li>
  <li>Since ClinVar's own gene attribution for multi-gene CNVs isn't usable (the <code>GeneID</code>
  <code>-1</code> problem from section {S.structural_variants}), whose job is computing real per-gene coordinate overlap &mdash;
  this ingest, at transform time, using a fetched gene model per variant (expensive, as the Ensembl calls in
  this report already show); or a separate upstream annotation step this ingest would just consume?</li>
</ol>

"""

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
  html {{ scroll-behavior: smooth; }}
  body {{ font-family: -apple-system, sans-serif; max-width: 1000px; margin: 40px auto; color: #0f172a; }}
  body {{ padding: 0 16px; }}
  h1 {{ font-size: 1.6rem; font-weight: 700; margin-bottom: 0.25rem; }}
  h2 {{ font-size: 1.4rem; font-weight: 700; margin: 4rem 0 1.1rem; padding-bottom: 0.6rem; }}
  h2 {{ border-bottom: 2px solid #e2e8f0; color: #0f172a; scroll-margin-top: 64px; }}
  h3 {{ font-size: 1.05rem; font-weight: 700; margin: 2.75rem 0 0.9rem; color: #1e293b; }}
  h3 {{ padding-left: 12px; border-left: 3px solid #2563eb; scroll-margin-top: 64px; }}
  .xf-panel h3 {{ font-size: 1rem; font-weight: 600; margin: 0 0 8px; color: #334155; }}
  .xf-panel h3 {{ padding-left: 0; border-left: none; }}
  .section-nav {{ position: sticky; top: 0; z-index: 20; background: rgba(255,255,255,0.96); }}
  .section-nav {{ backdrop-filter: blur(4px); display: flex; flex-wrap: wrap; gap: 6px; }}
  .section-nav {{ padding: 10px 0; margin: 0 0 1rem; border-bottom: 1px solid #e2e8f0; }}
  .section-nav a {{ padding: 5px 12px; border-radius: 999px; background: #f1f5f9; }}
  .section-nav a {{ color: #334155; font-size: 12.5px; font-weight: 600; text-decoration: none; white-space: nowrap; }}
  .section-nav a:hover {{ background: #dbeafe; color: #1d4ed8; }}
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
{section_nav()}

{section_heading("purpose")}
<p class="subtitle">
  The purpose of this interactive document is to provide an overview of the ClinVar data ingest into the
  Monarch Knowledge Graph (KG). ClinVar, a database of pathogenic human variation, is an open source
  NCBI-hosted repository that the KG leverages to support disease:gene pairings, and, specifically,
  variant:disease pairings. Within this exploration we detail what's in ClinVar, what we currently ingest
  into the KG, and how we can expand the ingest and the Biolink model for variant representation to better
  capture human pathogenic variation.
</p>

<div class="summary-box">
  <strong><a href="https://www.ncbi.nlm.nih.gov/clinvar/" target="_blank" rel="noopener">ClinVar</a></strong>
  is a free, public archive run by NCBI (part of NIH) that aggregates reported relationships between human
  genetic variants and health conditions. Clinical testing laboratories, research groups, and expert panels
  submit their own interpretation of a variant &mdash; e.g. <em>Pathogenic</em>, <em>Likely pathogenic</em>,
  <em>Uncertain significance</em>, <em>Likely benign</em>, or <em>Benign</em> &mdash; along with the disease
  it's associated with and the evidence behind that call. Because submitters vary widely in rigor, every
  submission carries a <strong>review status</strong> ("star rating", 0&ndash;{S.input_files}) reflecting how it was
  assessed: from a single submitter with no assessed criteria (0&#9733;) up to a formal expert panel or
  practice guideline (4&#9733;). Multiple independent submitters can (and often do) submit interpretations
  for the same variant, sometimes agreeing and sometimes conflicting. This ingest, and the analyses below,
  are about deciding which of that submitted evidence is trustworthy enough to promote into the knowledge
  graph as a variant-disease association.
</div>

{section_heading("illustrative-examples")}
<p class="subtitle">
  ClinVar covers several fundamentally different kinds of variation, from a single altered base up to
  megabase-scale rearrangements. Each card below is a simplified schematic (not to scale, and not a real
  genomic-coordinate plot -- see the gene models further down for that) of the underlying mechanism, paired
  with a real example pulled from the actual downloaded release and linked to its ClinVar record. The two
  "large recurrent" CNVs both come from the same 22q11.21 region (the classic DiGeorge/velocardiofacial
  locus) to show a deletion and its reciprocal duplication side by side; the two "single-gene" CNVs both
  come from DMD, whose exon-level deletions and duplications are the textbook example of this category.
  Structural types (CNVs, inversion) aren't in <code>clinvar.vcf.gz</code> at all -- see section {S.structural_variants}.
</p>
<div class="xf-panels" style="grid-template-columns: repeat(3, 1fr);">
{example_cards_html}
</div>

{gene_model_sections_html}

{section_heading("clinvar-curation")}
<p class="subtitle">
  ClinVar itself does not adjudicate whether a variant is truly pathogenic &mdash; it's a clearinghouse that
  aggregates independently-submitted classifications from clinical testing laboratories, research groups,
  and expert panels, each working from their own evidence. The "review status" used throughout this report
  (and gating production's own inclusion filter) is a proxy for how rigorously a classification was
  <em>reviewed</em>, not a guarantee that it's correct.
</p>
<div class="summary-box">
  <strong>Submitters</strong> report a variant, a classification (most commonly using the 2015 ACMG/AMP
  five-tier framework: Pathogenic / Likely pathogenic / Uncertain significance / Likely benign / Benign),
  the condition it's associated with, and &mdash; if they used a standardized framework &mdash; the specific
  evidence criteria they applied (e.g. PVS1, PM2, BP4). ClinVar accepts and displays these largely as
  submitted; it doesn't independently re-derive them itself.
  <br><br>
  <strong>Review status</strong> ("star rating") describes the <em>process</em> a classification went
  through, not ClinVar's own opinion of it:
</div>
<div class="table-wrap">
<table>
<tr><th>Stars</th><th>Review status</th><th>What it means</th></tr>
<tr>
  <td>0</td>
  <td>no assertion criteria provided / no classification provided / no classifications from unflagged
  records / no classification for the single variant / flagged submission</td>
  <td>A classification with no stated criteria, no classification at all, or one ClinVar has flagged as
  problematic.</td>
</tr>
<tr>
  <td>1</td><td>criteria provided, single submitter (or conflicting classifications)</td>
  <td>One submitter applied a standardized framework; or multiple submitters applied criteria but disagree
  with each other.</td>
</tr>
<tr>
  <td>2</td><td>criteria provided, multiple submitters, no conflicts</td>
  <td>&ge;2 independent submitters used criteria and agree &mdash; an aggregate ClinVar computes
  <em>across</em> a variant's submitters, which is why it can never appear on any single submission record
  (see section {S.star_cutoff}'s own docstring for why this matters to production's inclusion filter).</td>
</tr>
<tr>
  <td>3</td><td>reviewed by expert panel</td>
  <td>A <a href="https://clinicalgenome.org/affiliation/" target="_blank" rel="noopener">ClinGen Variant
  Curation Expert Panel</a> (VCEP) &mdash; a standing group of gene/disease specialists that develops
  refined, gene-specific ACMG/AMP rules &mdash; issued the classification directly.</td>
</tr>
<tr>
  <td>4</td><td>practice guideline</td>
  <td>Endorsed by a professional body's published guideline (e.g. pharmacogenomic dosing guidance) rather
  than a single variant-curation exercise.</td>
</tr>
</table>
</div>

<h3>What's actually in the downloaded data</h3>
<p class="subtitle">
  That's the conceptual picture. Empirically, <code>submission_summary.txt</code>'s per-record
  <code>ReviewStatus</code> field only ever takes {len([r for r in review_status_rows if r['count'] > 0]):,}
  distinct values out of the {len(review_status_rows):,} <code>review_star_map</code> knows how to score
  &mdash; the rest, greyed out below, are defined but never actually observed on any individual submission
  record in this release. Note in particular that the 2&#9733; row is empty: as discussed above, it's an
  aggregate ClinVar computes across a variant's submitters and only ever shows up in the VCF's
  <code>CLNREVSTAT</code> (section {S.crossfilter}), never here.
  The Stars column is scored by <code>clinvar_helpers.review_star_map</code> &mdash; the same single
  mapping that gates the production ingest, drives section {S.star_cutoff}'s <code>star_min</code> sweep, and
  labels section {S.crossfilter}'s <code>CLNREVSTAT</code> panel &mdash; and it now matches
  <a href="https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/" target="_blank" rel="noopener">ClinVar's
  published star table</a> value-for-value, including the four "no assertion / no classification" statuses
  that belong at 0&#9733; rather than 1&#9733;.
</p>
<div class="table-wrap">
<table>
<tr><th>Stars</th><th>ReviewStatus (as it appears in the file)</th><th># records</th><th>% of records</th></tr>
{review_status_table_rows}
</table>
</div>
<p class="subtitle">
  A few consequences worth knowing before the star-cutoff analyses below: (1) review status says nothing
  about how much evidence exists, only how it was reviewed &mdash; a variant seen once and thoroughly
  assessed by one lab can outrank, in stars, a well-established variant whose many submitters simply never
  used standardized criteria; (2) classifications are periodically revisited and can change as evidence
  accumulates, so a snapshot download like the one behind this report is already slightly out of date the
  moment it's downloaded; (3) ClinGen's expert panels cover only a small, deliberately prioritized slice of
  genes &mdash; most genes have no VCEP at all, so 3&ndash;{S.input_files}&#9733; evidence is systematically concentrated
  in a handful of well-studied disease genes (like the LMNA and SMCHD1 examples above), and its absence
  elsewhere reflects a lack of formal panel review, not necessarily a lack of confidence. Section 6
  below quantify exactly how much of production's evidence lives below the expert-panel tier, and test
  whether independent multi-submitter agreement is a reasonable proxy for curation quality when a formal
  panel review hasn't happened.
</p>

<p class="subtitle">
  Sections 5&ndash;{S.evidence_tiers} are the analyses over the downloaded ClinVar release: (4) review-star cutoff
  impact on the production disease-mapping filter, (5) reconciliation against the Monarch KG's
  curated gene-disease associations, and (6) a crossfilter over ClinVar's own variant-level
  classification and review-status fields. Section 10 turns those findings into concrete ingest
  changes; section {S.structural_variants} covers the structural variants the VCF omits entirely.
</p>

{input_files_html}

{scope_decisions_html}

{phenotype_terms_html}

{two_star_html}

{pairing_html}

{section_heading("star-cutoff")}
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
<p class="subtitle" style="font-size:12.5px; margin-top:0.5rem;">
  <span style="color:#6d28d9; font-weight:700;">*</span> <strong>2&#9733; (computed)</strong> is not a raw
  <code>ReviewStatus</code> filter &mdash; as the table in "How ClinVar curation works" above shows, no
  individual submission record ever actually carries the aggregate 2-star status, so a literal
  <code>star_min=2</code> filter over per-record data is identical to <code>star_min=3</code> (compare the
  two rows above). This row instead <em>reconstructs</em> what a true &ge;2&#9733; population would be:
  the raw &ge;2&#9733; pairs, plus any pair rescued by &ge;{MIN_CONCORDANT_SUBMITTERS} independent concordant
  submitters &mdash; the same proxy mechanism section {S.star_cutoff} below uses, applied here directly to the headline
  counts rather than shown as a separate rescue analysis.
</p>

<h3>Gene &rarr; MONDO disease pairs, by evidence tier</h3>
<p class="subtitle">
  Every gene-disease pair with &ge;1 Pathogenic/Likely-pathogenic submission record splits into exactly
  <strong>one</strong> of the three tiers below &mdash; each pair appears in exactly one table, never more
  than one. Tier 1 (raw &ge;3&#9733;) is checked first; anything left over is checked against 2&#9733;
  concordance; anything <em>still</em> left over falls to tier 3. Gene attribution throughout is the single gene ClinVar
  asserts for each variant (<code>variant_summary.txt.gz</code>), not the VCF's positional
  <code>GENEINFO</code> list &mdash; see section {S.ingest_recommendation} for what that changed.
</p>

<h3 style="margin-top:1.75rem;">Tier 1: star_min &ge;3 (raw)</h3>
<p class="subtitle">
  Real, unambiguous evidence: at least one submission record for this pair actually reaches 3&#9733; or
  4&#9733; on its own (raw star_min=2 and =3 are identical here, see script docstring, so this is also
  every pair a literal &ge;2&#9733; filter would find). The &ge;4&#9733; column highlights which of these
  also survive the strictest cutoff.
</p>
<div class="controls">
  <input id="pairs-ge3-search" class="search-box" type="text" placeholder="Filter by gene or MONDO id...">
  <label style="font-size:13px; color:#475569; display:flex; align-items:center; gap:4px; white-space:nowrap;">
    <input id="pairs-ge3-only4" type="checkbox"> only &ge;4&#9733; survivors
  </label>
  <span id="pairs-ge3-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th data-sort-for="pairs-ge3-rows" data-sort-key="gene">Gene &#8597;</th>
  <th data-sort-for="pairs-ge3-rows" data-sort-key="gene_id">Gene ID &#8597;</th>
  <th data-sort-for="pairs-ge3-rows" data-sort-key="mondo">MONDO disease &#8597;</th>
  <th data-sort-for="pairs-ge3-rows" data-sort-key="disease_name">Disease name &#8597;</th>
  <th data-sort-for="pairs-ge3-rows" data-sort-key="v2">Variants (&ge;3&#9733;) &#8597;</th>
  <th data-sort-for="pairs-ge3-rows" data-sort-key="v4">Variants (&ge;4&#9733;) &#8597;</th>
  <th data-sort-for="pairs-ge3-rows" data-sort-key="n_submitters"># submitters (pooled) &#8597;</th>
</tr>
</thead>
<tbody id="pairs-ge3-rows"></tbody>
</table>
</div>

<h3 style="margin-top:1.75rem;">Tier 2: 2&#9733; concordance (computed)</h3>
<p class="subtitle">
  Pairs with <strong>no</strong> raw &ge;3&#9733; evidence anywhere, reconstructed instead via
  &ge;{MIN_CONCORDANT_SUBMITTERS} independent concordant submitters &mdash; the same proxy mechanism
  section {S.star_cutoff} below uses for its own, more detailed rescue analysis (submitter counts, literature flags,
  etc. per pair). Mutually exclusive of Tier 1 by construction: any pair already in Tier 1 is excluded
  here even if it would also separately qualify for concordance.
</p>
<div class="controls">
  <input id="pairs-2star-search" class="search-box" type="text" placeholder="Filter by gene or MONDO id...">
  <span id="pairs-2star-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th data-sort-for="pairs-2star-rows" data-sort-key="gene">Gene &#8597;</th>
  <th data-sort-for="pairs-2star-rows" data-sort-key="gene_id">Gene ID &#8597;</th>
  <th data-sort-for="pairs-2star-rows" data-sort-key="mondo">MONDO disease &#8597;</th>
  <th data-sort-for="pairs-2star-rows" data-sort-key="disease_name">Disease name &#8597;</th>
  <th data-sort-for="pairs-2star-rows" data-sort-key="n_variants"># variants &#8597;</th>
  <th data-sort-for="pairs-2star-rows" data-sort-key="n_submitters"># submitters (pooled) &#8597;</th>
</tr>
</thead>
<tbody id="pairs-2star-rows"></tbody>
</table>
</div>

<h3 style="margin-top:1.75rem;">Tier 3: remaining (below concordance threshold)</h3>
<p class="subtitle">
  Everything left over: &ge;1 Pathogenic/Likely-pathogenic submission record, but neither raw &ge;3&#9733;
  evidence nor &ge;{MIN_CONCORDANT_SUBMITTERS}-concordant-submitter rescue. This is the largest tier by far
  &mdash; mostly single-submitter, 0/1&#9733; calls with no independent corroboration yet.
  <strong>Read the two count columns separately.</strong> Both tier assignment tests are
  <em>per-variant</em>: a pair reaches tier 1 or 2 when some <em>single</em> variant of that pair
  clears &ge;3&#9733; or draws &ge;{MIN_CONCORDANT_SUBMITTERS} submitters <em>by itself</em>.
  "# submitters (pooled)" is the different, pair-level question &mdash; how many distinct submitters
  have said anything P/LP about this gene-disease relationship across <em>all</em> its variants.
  Sorting this table by that column (the default) surfaces the pairs the per-variant test discards
  despite broad submitter support.
</p>
<div class="summary-box" style="margin-bottom:1rem;">
  <strong>Worked example &mdash; why "# variants" is not an evidence score.</strong>
  <code>SCN1A</code> &rarr; <code>MONDO:0800491</code> (early-infantile DEE) sits in this tier with
  <strong>the most variants of any of its disease terms</strong>, yet a pooled submitter count of
  <strong>1</strong> &mdash; see the table above. Every one of
  those variants came from Labcorp/Invitae, which is the only submitter using MedGen
  <code>C0393706:Early-infantile DEE</code>. The other ~87 labs describing the same SCN1A biology
  submit MedGen <code>C0751122:Severe myoclonic epilepsy in infancy</code>, which maps to a
  <em>different</em> MONDO id &mdash; <code>MONDO:0100135</code> (Dravet syndrome), sitting in tier 1
  on just 4 variants. Because concordance is keyed on <code>(mondo_id, ClinicalSignificance)</code>,
  submitter agreement expressed through a different condition name is invisible to it: one gene's
  clinical spectrum splits across two terms and neither side can corroborate the other. That variant
  count measures one lab's submission volume, not corroboration. (Those two terms are genuinely
  distinct diseases in MONDO, so the fix is not to merge them &mdash; see section {S.monarch_kg}.) Two further leaks compound this: agreement also
  requires the <em>exact same</em> <code>ClinicalSignificance</code> string, so
  <em>Pathogenic</em> + <em>Likely pathogenic</em> on the same variant and disease does not count;
  and high-volume labs frequently submit the placeholder <code>C3661900:not provided</code>
  (GeneDx contributes 535 such SCN1A P/LP records, CeGaT 148), which maps to no disease at all and
  so can never support any pair.
</div>
<div class="controls">
  <input id="pairs-remaining-search" class="search-box" type="text" placeholder="Filter by gene or MONDO id...">
  <span id="pairs-remaining-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th data-sort-for="pairs-remaining-rows" data-sort-key="gene">Gene &#8597;</th>
  <th data-sort-for="pairs-remaining-rows" data-sort-key="gene_id">Gene ID &#8597;</th>
  <th data-sort-for="pairs-remaining-rows" data-sort-key="mondo">MONDO disease &#8597;</th>
  <th data-sort-for="pairs-remaining-rows" data-sort-key="disease_name">Disease name &#8597;</th>
  <th data-sort-for="pairs-remaining-rows" data-sort-key="n_variants"># variants &#8597;</th>
  <th data-sort-for="pairs-remaining-rows" data-sort-key="n_submitters"># submitters (pooled) &#8597;</th>
</tr>
</thead>
<tbody id="pairs-remaining-rows"></tbody>
</table>
</div>

<h3 style="margin-top:1.75rem;">ClinGen-only single-submitter pairs</h3>
<p class="subtitle">
  A related question: among gene-disease pairs whose entire evidence comes from exactly one submitter (and
  which can therefore never be concordance-rescued above &mdash; that needs &ge;2 independent submitters),
  how many of those lone submitters are <a href="https://clinicalgenome.org/" target="_blank"
  rel="noopener">ClinGen</a>-affiliated (a Variant Curation Expert Panel, the ACMG Interpretation Working
  Group, etc. &mdash; detected as a "ClinGen" substring in the Submitter field, <code>submission_summary.txt</code>'s
  only proxy for submitter identity)? And more broadly: for how many pairs is ClinGen the <em>only</em>
  source of Pathogenic/Likely-pathogenic evidence at all &mdash; i.e. the pair would have zero P/LP
  evidence, from anyone, if ClinGen's own submissions were removed?
</p>
<div class="summary-box">
  <strong>{clingen_summary["single_submitter_clingen_pairs"]:,}</strong> of
  <strong>{clingen_summary["single_submitter_pairs"]:,}</strong> single-submitter gene-disease pairs
  ({clingen_single_pct:.1f}%) have ClinGen as that lone submitter &mdash; these can only ever reach
  production today if that one ClinGen submission itself is a 3&#9733; expert-panel review.<br>
  <strong>{clingen_summary["clingen_only_pairs"]:,}</strong> of
  <strong>{clingen_summary["total_pairs"]:,}</strong> gene-disease pairs overall
  ({clingen_only_pct:.1f}%) have ClinGen as their <em>only</em> source of evidence (one or more
  ClinGen-affiliated submitters, and no independent non-ClinGen submitter at any star level).
</div>
<div class="controls">
  <input id="clingen-table-search" class="search-box" type="text"
    placeholder="Filter by gene, MONDO id, or disease name...">
  <span id="clingen-table-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th data-sort-for="clingen-table-rows" data-sort-key="gene">Gene &#8597;</th>
  <th data-sort-for="clingen-table-rows" data-sort-key="gene_id">Gene ID &#8597;</th>
  <th data-sort-for="clingen-table-rows" data-sort-key="mondo">MONDO disease &#8597;</th>
  <th data-sort-for="clingen-table-rows" data-sort-key="disease_name">Disease name &#8597;</th>
  <th data-sort-for="clingen-table-rows" data-sort-key="n_variants"># variants &#8597;</th>
  <th data-sort-for="clingen-table-rows" data-sort-key="n_clingen_submitters">ClinGen submitters &#8597;</th>
  <th>Submitter(s)</th>
</tr>
</thead>
<tbody id="clingen-table-rows"></tbody>
</table>
</div>
<div class="pagination">
  <button id="clingen-table-prev">&larr; Prev</button>
  <span id="clingen-table-page-info" class="page-info"></span>
  <button id="clingen-table-next">Next &rarr;</button>
</div>

{evidence_tiers_html}

{monarch_section_html}

{recommendation_html}

{ingest_compare_html}

{crossfilter_html}

{section_heading("structural-variants")}
<p class="subtitle">
  Sections 5&ndash;{S.evidence_tiers} above are built entirely from <code>clinvar.vcf.gz</code> and
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
  One row per variant across copy number gain/loss, Translocation, and Complex, any
  ClinicalSignificance &mdash; GRCh38 where ClinVar provides it, GRCh37 otherwise, never both.
  <strong>That build choice dominates this section.</strong> Of the
  {sv_b38['count'] + sv_b37['count']:,} SV/CNV variants placed on a usable build,
  <strong>{sv_b38['count']:,} are represented by a GRCh38 row and only
  {100 * sv_b38['resolved'] / max(sv_b38['count'], 1):.0f}% of those resolve to a disease id</strong>
  &mdash; they are mostly ClinGen-style dosage-sensitivity regions with no single named condition,
  whose <code>PhenotypeList</code> is the literal placeholder <code>"See cases"</code>. The
  <strong>{sv_b37['count']:,} that ClinVar has only ever placed on GRCh37 are
  {100 * sv_b37['resolved'] / max(sv_b37['count'], 1):.0f}% resolved.</strong>
  ({sv_summary['unplaced']:,} more are on NCBI36 or no build at all and cannot be placed here.)
  An earlier version of this report
  filtered to GRCh38 and consequently reported a ~6&ndash;{S.ingest_recommendation}% resolution rate as though it were a
  property of ClinVar's CNV curation; it was an artifact of dropping the older, better-curated
  build. Uncheck "SV resolved to disease" to browse the unresolved subset. Note the
  per-variant gene attribution used elsewhere in this report has no equivalent here &mdash;
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

<h3>Gene-disease pairs invisible to production: CNV/SV-only evidence</h3>
<p class="subtitle">
  Production's ingest reads <code>clinvar.vcf.gz</code> exclusively &mdash; a fixed-REF/ALT format that
  structurally cannot represent a CNV, translocation, or inversion at all (see above). That means any
  gene-disease pair whose <em>entire</em> body of evidence is structural has zero chance of ever producing a
  <code>VariantToGeneAssociation</code>/<code>VariantToDiseaseAssociation</code> edge today, no matter how
  well-reviewed that structural evidence is. This table finds exactly those pairs: Pathogenic/Likely-pathogenic
  CNV/SV rows (from <code>variant_summary.txt.gz</code>) unambiguously about one gene (multi-gene and
  "covers N genes" regional entries are excluded &mdash; there's no reliable per-gene id to pair with a
  disease for those, see <code>is_single_clean_gene_symbol()</code>), resolved to a MONDO disease, whose
  (gene, disease) pair <strong>never once appears</strong> anywhere in the star_min=0 SNV/indel pair set
  from section {S.star_cutoff} above (any review status, any evidence at all).
</p>
<div class="summary-box">
  <strong>{sv_only_pairs_count:,}</strong> of <strong>{sv_gene_disease_pairs_count:,}</strong> single-gene
  CNV/SV-derived gene-disease pairs ({sv_only_pct:.1f}%) have no corresponding SNV/indel evidence anywhere
  in <code>clinvar.vcf.gz</code>.
</div>
<div class="controls">
  <input id="sv-only-table-search" class="search-box" type="text"
    placeholder="Filter by gene, MONDO id, or disease name...">
  <span id="sv-only-table-count" class="count-label"></span>
</div>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th data-sort-for="sv-only-table-rows" data-sort-key="gene">Gene &#8597;</th>
  <th data-sort-for="sv-only-table-rows" data-sort-key="gene_id">Gene ID &#8597;</th>
  <th data-sort-for="sv-only-table-rows" data-sort-key="mondo">MONDO disease &#8597;</th>
  <th data-sort-for="sv-only-table-rows" data-sort-key="disease_name">Disease name &#8597;</th>
  <th>SV type(s)</th>
  <th data-sort-for="sv-only-table-rows" data-sort-key="n_variants"># variants &#8597;</th>
  <th>ClinVar IDs</th>
</tr>
</thead>
<tbody id="sv-only-table-rows"></tbody>
</table>
</div>
<div class="pagination">
  <button id="sv-only-table-prev">&larr; Prev</button>
  <span id="sv-only-table-page-info" class="page-info"></span>
  <button id="sv-only-table-next">Next &rarr;</button>
</div>

{biolink_proposal_html}

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

// Shared by all three evidence-tier pairs tables in section {S.star_cutoff} (and reusable for similar
// gene/mondo pairs tables elsewhere): plain sort+search over an in-memory array, no
// pagination -- matches this report's existing "show everything, let search narrow it
// down" pattern for pairs-scale tables (hundreds to tens of thousands of rows).
function setupPairsTable(config) {{
  // config: {{ data, tbodyId, countElId, searchElId, defaultSortKey, extraCheckboxElId,
  //   extraCheckboxPredicate, rowHtml(r) }}
  var state = {{ sortKey: config.defaultSortKey, sortDir: -1, query: "" }};
  var tbody = document.getElementById(config.tbodyId);
  var countEl = document.getElementById(config.countElId);
  var searchEl = document.getElementById(config.searchElId);
  var checkboxEl = config.extraCheckboxElId ? document.getElementById(config.extraCheckboxElId) : null;

  function render() {{
    var rows = config.data.filter(function(r) {{
      if (checkboxEl && checkboxEl.checked && config.extraCheckboxPredicate && !config.extraCheckboxPredicate(r)) return false;
      if (!state.query) return true;
      var q = state.query.toLowerCase();
      return r.gene.toLowerCase().indexOf(q) !== -1 ||
             r.gene_id.toLowerCase().indexOf(q) !== -1 ||
             r.mondo.toLowerCase().indexOf(q) !== -1 ||
             (r.disease_name || "").toLowerCase().indexOf(q) !== -1;
    }});

    rows.sort(function(a, b) {{
      var av = a[state.sortKey], bv = b[state.sortKey];
      if (av === null || av === undefined) av = -1;
      if (bv === null || bv === undefined) bv = -1;
      if (typeof av === "string") return av.localeCompare(bv) * state.sortDir;
      return (av - bv) * state.sortDir;
    }});

    countEl.textContent = rows.length + " of " + config.data.length + " pairs";
    tbody.innerHTML = rows.map(config.rowHtml).join("");
  }}

  searchEl.addEventListener("input", function(e) {{
    state.query = e.target.value;
    render();
  }});
  if (checkboxEl) checkboxEl.addEventListener("change", render);
  document.querySelectorAll("[data-sort-for='" + config.tbodyId + "']").forEach(function(th) {{
    th.addEventListener("click", function() {{
      var k = th.getAttribute("data-sort-key");
      if (state.sortKey === k) {{
        state.sortDir *= -1;
      }} else {{
        state.sortKey = k;
        state.sortDir = (k === "gene" || k === "gene_id" || k === "mondo" || k === "disease_name") ? 1 : -1;
      }}
      render();
    }});
  }});

  render();
}}

(function() {{
  setupPairsTable({{
    data: {pairs_ge3_json},
    tbodyId: "pairs-ge3-rows",
    countElId: "pairs-ge3-count",
    searchElId: "pairs-ge3-search",
    defaultSortKey: "v2",
    extraCheckboxElId: "pairs-ge3-only4",
    extraCheckboxPredicate: function(r) {{ return r.v4 !== null; }},
    rowHtml: function(r) {{
      var v4cell = r.v4 === null ? "<span class='no'>&mdash;</span>" : "<span class='yes'>" + r.v4 + "</span>";
      return "<tr>" +
        "<td>" + r.gene + "</td>" +
        "<td class='mono'>" + r.gene_id + "</td>" +
        "<td class='mono'>" + r.mondo + "</td>" +
        "<td>" + (r.disease_name || "<span class='no'>&mdash;</span>") + "</td>" +
        "<td class='num'>" + r.v2 + "</td>" +
        "<td class='num'>" + v4cell + "</td>" +
        "<td class='num'>" + r.n_submitters + "</td>" +
      "</tr>";
    }},
  }});
}})();

(function() {{
  setupPairsTable({{
    data: {pairs_2star_json},
    tbodyId: "pairs-2star-rows",
    countElId: "pairs-2star-count",
    searchElId: "pairs-2star-search",
    defaultSortKey: "n_variants",
    rowHtml: function(r) {{
      return "<tr>" +
        "<td>" + r.gene + "</td>" +
        "<td class='mono'>" + r.gene_id + "</td>" +
        "<td class='mono'>" + r.mondo + "</td>" +
        "<td>" + (r.disease_name || "<span class='no'>&mdash;</span>") + "</td>" +
        "<td class='num'>" + r.n_variants + "</td>" +
        "<td class='num'>" + r.n_submitters + "</td>" +
      "</tr>";
    }},
  }});
}})();

(function() {{
  setupPairsTable({{
    data: {pairs_remaining_json},
    tbodyId: "pairs-remaining-rows",
    countElId: "pairs-remaining-count",
    searchElId: "pairs-remaining-search",
    defaultSortKey: "n_submitters",
    rowHtml: function(r) {{
      return "<tr>" +
        "<td>" + r.gene + "</td>" +
        "<td class='mono'>" + r.gene_id + "</td>" +
        "<td class='mono'>" + r.mondo + "</td>" +
        "<td>" + (r.disease_name || "<span class='no'>&mdash;</span>") + "</td>" +
        "<td class='num'>" + r.n_variants + "</td>" +
        "<td class='num'>" + r.n_submitters + "</td>" +
      "</tr>";
    }},
  }});
}})();


(function() {{
  function bind(id, target) {{
    var cb = document.getElementById(id);
    if (!cb) return;
    cb.addEventListener("change", function() {{
      (target || []).forEach(function(g) {{
        var el = document.getElementById(g);
        if (el) el.style.display = cb.checked ? "" : "none";
      }});
    }});
  }}
  bind("scn1a-toggle-shares", ["scn1a-shares", "scn1a-shares-labels"]);
  bind("scn1a-toggle-xrefs", ["scn1a-xrefs"]);
}})();

(function() {{
  var cube = {filter_cube_json};
  var STAR_OPTS = [
    {{v: 99, label: "Off"}}, {{v: 4, label: "&ge;4&#9733;"}}, {{v: 3, label: "&ge;3&#9733; (shipping)"}},
    {{v: 1, label: "&ge;1&#9733;"}}, {{v: 0, label: "&ge;0&#9733; (any record)"}}
  ];
  var CONC_OPTS = [
    {{v: 99, label: "Off"}}, {{v: 5, label: "&ge;5 submitters"}}, {{v: 3, label: "&ge;3 submitters"}},
    {{v: 2, label: "&ge;2 submitters (shipping)"}}
  ];
  var AGG_OPTS = [
    {{v: 99, label: "Off"}}, {{v: 4, label: "&ge;4&#9733;"}}, {{v: 3, label: "&ge;3&#9733;"}},
    {{v: 2, label: "&ge;2&#9733; (shipping)"}}
  ];
  var state = {{star: 3, conc: 2, agg: 2, geneOnly: false}};

  function count(starMin, concMin, aggMin, geneOnly) {{
    var v = 0, e = 0;
    cube.variants.forEach(function(r) {{
      if (geneOnly && !r.g) return;
      if (r.s >= starMin || r.c >= concMin || r.a >= aggMin) v += r.n;
    }});
    cube.edges.forEach(function(r) {{
      if (r.s >= starMin || r.c >= concMin || r.a >= aggMin) e += r.n;
    }});
    return {{variants: v, edges: e}};
  }}

  function renderOpts(containerId, opts, key) {{
    document.getElementById(containerId).innerHTML = opts.map(function(o) {{
      var on = state[key] === o.v;
      return "<div class='xf-row' data-fp='" + key + "' data-v='" + o.v + "'>" +
        "<input type='radio'" + (on ? " checked" : "") + " data-fp='" + key + "' data-v='" + o.v + "'>" +
        "<span style='font-size:12.5px; color:#334155;'>" + o.label + "</span></div>";
    }}).join("");
  }}

  var PRESETS = [
    {{name: "Per-record + concordance only", star: 3, conc: 2, agg: 99}},
    {{name: "All three paths (shipping)", star: 3, conc: 2, agg: 2}},
    {{name: "Expert panel only", star: 3, conc: 99, agg: 99}},
    {{name: "Aggregate 2&#9733; alone", star: 99, conc: 99, agg: 2}},
    {{name: "Concordance alone", star: 99, conc: 2, agg: 99}},
    {{name: "Everything P/LP", star: 0, conc: 99, agg: 99}}
  ];

  function fmtStar(v) {{ return v === 99 ? "off" : "&ge;" + v + "&#9733;"; }}
  function fmtConc(v) {{ return v === 99 ? "off" : "&ge;" + v; }}

  function render() {{
    renderOpts("fp-star-rows", STAR_OPTS, "star");
    renderOpts("fp-conc-rows", CONC_OPTS, "conc");
    renderOpts("fp-agg-rows", AGG_OPTS, "agg");
    var r = count(state.star, state.conc, state.agg, state.geneOnly);
    document.getElementById("fp-variants").textContent = r.variants.toLocaleString();
    document.getElementById("fp-edges").textContent = r.edges.toLocaleString();
    document.getElementById("fp-preset-rows").innerHTML = PRESETS.map(function(p) {{
      var c = count(p.star, p.conc, p.agg, false);
      var isCur = p.star === state.star && p.conc === state.conc && p.agg === state.agg;
      return "<tr" + (isCur ? " class='prod'" : "") + "><td>" + p.name + "</td>" +
        "<td>" + fmtStar(p.star) + "</td><td>" + fmtConc(p.conc) + "</td><td>" + fmtStar(p.agg) + "</td>" +
        "<td class='num'>" + c.variants.toLocaleString() + "</td>" +
        "<td class='num'>" + c.edges.toLocaleString() + "</td></tr>";
    }}).join("");
    var sc = count(3, 2, 2, false).variants;
    document.getElementById("fp-selfcheck").textContent = sc.toLocaleString() + " variants";
  }}

  document.body.addEventListener("change", function(e) {{
    if (e.target.id === "fp-gene-only") {{ state.geneOnly = e.target.checked; render(); return; }}
    if (e.target.tagName !== "INPUT" || !e.target.hasAttribute("data-fp")) return;
    state[e.target.getAttribute("data-fp")] = parseInt(e.target.getAttribute("data-v"), 10);
    render();
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
  // CLNREVSTAT text behind each star level, so the panel says which aggregate review status is
  // being selected rather than a bare number (see review_star_map in clinvar_helpers.py). Short
  // form goes in the label; full form is the hover tooltip.
  var starShortText = {{
    0: "no criteria / no classification",
    1: "single submitter, or conflicting",
    2: "multiple submitters, no conflicts",
    3: "expert panel",
    4: "practice guideline"
  }};
  var starFullText = {{
    0: "no assertion criteria provided / no classification provided / no classifications from unflagged records / no classification for the single variant / flagged submission / .",
    1: "criteria provided, single submitter | criteria provided, conflicting classifications",
    2: "criteria provided, multiple submitters, no conflicts",
    3: "reviewed by expert panel",
    4: "practice guideline"
  }};

  var clnsigColors = {{
    "P": "#b91c1c", "LP": "#ea580c", "P/LP": "#c2410c",
    "VUS": "#a16207", "LB": "#16a34a", "B": "#15803d", "B/LB": "#166534",
    "Conflicting": "#7c3aed", "Other": "#64748b", "Not classified": "#94a3b8"
  }};

  var boolValues = [true, false];
  var ALL_DIMS = ["clnsig", "star", "clnvc", "size", "literature", "concordant", "strchive", "production"];

  var state = {{
    clnsig: new Set(clnsigBuckets),
    star: new Set(starLevels),
    clnvc: new Set(clnvcTypes),
    size: new Set(sizeBuckets),
    literature: new Set(boolValues),
    concordant: new Set(boolValues),
    strchive: new Set(boolValues),
    production: new Set(boolValues)
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

  function renderPanel(dim, values, otherDims, containerId, labelFn, colorFn, titleFn) {{
    var totals = breakdown(dim, otherDims);
    var maxVal = Math.max.apply(null, values.map(function(v) {{ return totals[v] || 0; }})) || 1;
    var container = document.getElementById(containerId);
    container.innerHTML = values.map(function(v) {{
      var n = totals[v] || 0;
      var on = state[dim].has(v);
      var pct = (n / maxVal) * 100;
      var color = colorFn ? colorFn(v) : "#2563eb";
      var title = titleFn ? " title='" + titleFn(v).replace(/'/g, "&#39;") + "'" : "";
      return "<div class='xf-row' data-dim='" + dim + "' data-value='" + v + "'" + title + ">" +
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
      function(v) {{
        return "<strong>" + v + "&#9733;</strong>" +
          "<span style='color:#94a3b8;'> &mdash; " + starShortText[v] + "</span>";
      }},
      function() {{ return "#2563eb"; }},
      function(v) {{ return "CLNREVSTAT: " + starFullText[v]; }});
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
    renderPanel("production", boolValues, otherDims("production"), "xf-production-rows",
      function(v) {{ return v ? "Yes &mdash; ingested today" : "No &mdash; filtered out"; }},
      function(v) {{ return v ? "#f59e0b" : "#2563eb"; }});
  }}

  document.body.addEventListener("change", function(e) {{
    if (e.target.tagName !== "INPUT" || !e.target.hasAttribute("data-dim")) return;
    var dim = e.target.getAttribute("data-dim");
    if (!(dim in state)) return;
    var raw = e.target.getAttribute("data-value");
    var value;
    if (dim === "star") {{
      value = parseInt(raw, 10);
    }} else if (dim === "literature" || dim === "concordant" || dim === "strchive" || dim === "production") {{
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

(function() {{
  var clingenTable = createPaginatedTable({{
    tbodyId: "clingen-table-rows",
    countElId: "clingen-table-count",
    searchElId: "clingen-table-search",
    prevBtnId: "clingen-table-prev",
    nextBtnId: "clingen-table-next",
    pageInfoElId: "clingen-table-page-info",
    pageSize: 25,
    defaultSortKey: "n_variants",
    defaultSortDir: -1,
    searchFields: ["gene", "gene_id", "mondo", "disease_name"],
    columns: [
      {{ key: "gene" }},
      {{ key: "gene_id", className: "mono" }},
      {{ key: "mondo", className: "mono" }},
      {{ key: "disease_name", format: function(v) {{ return v || "<span class='no'>&mdash;</span>"; }} }},
      {{ key: "n_variants", className: "num" }},
      {{ key: "n_clingen_submitters", className: "num" }},
      {{ key: "clingen_submitters", format: function(v) {{ return v.join(", "); }} }},
    ],
  }});
  clingenTable.setData({clingen_only_json});
}})();

(function() {{
  function clinvarLink(id) {{
    var url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/" + id + "/";
    return "<a href='" + url + "' target='_blank' rel='noopener'>" + id + "</a>";
  }}
  var svOnlyTable = createPaginatedTable({{
    tbodyId: "sv-only-table-rows",
    countElId: "sv-only-table-count",
    searchElId: "sv-only-table-search",
    prevBtnId: "sv-only-table-prev",
    nextBtnId: "sv-only-table-next",
    pageInfoElId: "sv-only-table-page-info",
    pageSize: 25,
    defaultSortKey: "n_variants",
    defaultSortDir: -1,
    searchFields: ["gene", "gene_id", "mondo", "disease_name"],
    columns: [
      {{ key: "gene" }},
      {{ key: "gene_id", className: "mono" }},
      {{ key: "mondo", className: "mono" }},
      {{ key: "disease_name", format: function(v) {{ return v || "<span class='no'>&mdash;</span>"; }} }},
      {{ key: "types", format: function(v) {{ return v.join(", "); }} }},
      {{ key: "n_variants", className: "num" }},
      {{ key: "variant_sample", format: function(v) {{ return v.map(clinvarLink).join(", "); }} }},
    ],
  }});
  svOnlyTable.setData({sv_only_json});
}})();

(function() {{
  var TIER_TEXT = {{1: "1 (>=3 stars)", 2: "2 (concordance)", 3: "3 (not ingested)"}};
  var KINSHIP_TEXT = {{
    "ancestor in Monarch": "<span style='color:#b45309;'>ancestor in Monarch</span>",
    "descendant in Monarch": "<span style='color:#b45309;'>descendant in Monarch</span>",
    "gene in Monarch, unrelated term": "gene in Monarch, unrelated term",
    "gene absent from Monarch": "<span style='color:#15803d;'>gene absent from Monarch</span>"
  }};

  var monarchOnlyTable = createPaginatedTable({{
    tbodyId: "monarch-only-rows",
    countElId: "monarch-only-count",
    searchElId: "monarch-only-search",
    prevBtnId: "monarch-only-prev",
    nextBtnId: "monarch-only-next",
    pageInfoElId: "monarch-only-page-info",
    pageSize: 25,
    defaultSortKey: "n_submitters",
    defaultSortDir: -1,
    searchFields: ["gene", "gene_id", "mondo", "disease_name", "kinship"],
    columns: [
      {{ key: "gene" }},
      {{ key: "mondo", className: "mono" }},
      {{ key: "disease_name", format: function(v) {{ return v || "<span class='no'>&mdash;</span>"; }} }},
      {{ key: "n_variants", className: "num" }},
      {{ key: "n_submitters", className: "num" }},
      {{ key: "tier", format: function(v) {{ return TIER_TEXT[v] || v; }} }},
      {{ key: "kinship", format: function(v) {{ return KINSHIP_TEXT[v] || v; }} }},
    ],
  }});
  monarchOnlyTable.setData({monarch_only_json});

  var monarchSupportedTable = createPaginatedTable({{
    tbodyId: "monarch-supported-rows",
    countElId: "monarch-supported-count",
    searchElId: "monarch-supported-search",
    prevBtnId: "monarch-supported-prev",
    nextBtnId: "monarch-supported-next",
    pageInfoElId: "monarch-supported-page-info",
    pageSize: 25,
    defaultSortKey: "n_submitters",
    defaultSortDir: -1,
    searchFields: ["gene", "gene_id", "mondo", "disease_name"],
    columns: [
      {{ key: "gene" }},
      {{ key: "mondo", className: "mono" }},
      {{ key: "disease_name", format: function(v) {{ return v || "<span class='no'>&mdash;</span>"; }} }},
      {{ key: "n_variants", className: "num" }},
      {{ key: "n_submitters", className: "num" }},
      {{ key: "tier", format: function(v) {{ return TIER_TEXT[v] || v; }} }},
    ],
  }});
  monarchSupportedTable.setData({monarch_supported_json});

  var MONDO_STATUS_TEXT = {{
    "current in MONDO": "<span style='color:#15803d;'>current</span>",
    "deprecated in MONDO": "<span style='color:#b91c1c;'>deprecated</span>",
    "not in MONDO": "<span style='color:#b91c1c;'>not in MONDO</span>"
  }};

  var monarchKgOnlyTable = createPaginatedTable({{
    tbodyId: "monarch-kg-only-rows",
    countElId: "monarch-kg-only-count",
    searchElId: "monarch-kg-only-search",
    prevBtnId: "monarch-kg-only-prev",
    nextBtnId: "monarch-kg-only-next",
    pageInfoElId: "monarch-kg-only-page-info",
    checkboxFilterElId: "monarch-kg-only-notcurrent",
    checkboxFilterPredicate: function(r) {{ return r.mondo_status !== "current in MONDO"; }},
    pageSize: 25,
    defaultSortKey: "gene",
    defaultSortDir: 1,
    searchFields: ["gene", "gene_id", "mondo", "disease_name", "sources", "mondo_status"],
    columns: [
      {{ key: "gene", format: function(v) {{ return v || "<span class='no'>&mdash;</span>"; }} }},
      {{ key: "gene_id", className: "mono" }},
      {{ key: "mondo", className: "mono" }},
      {{ key: "disease_name", format: function(v) {{ return v || "<span class='no'>&mdash;</span>"; }} }},
      {{ key: "sources" }},
      {{ key: "mondo_status", format: function(v) {{ return MONDO_STATUS_TEXT[v] || v; }} }},
    ],
  }});
  monarchKgOnlyTable.setData({monarch_kg_only_json});
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

    # Same curated attribution the production transform uses -- see variant_genes_for()
    ensure_variant_summary_downloaded(args.data_dir)
    variant_genes = make_variant_gene_map(args.data_dir / "variant_summary.txt.gz")
    print(f"Gene attribution: {len(variant_genes):,} variants with a ClinVar-asserted gene")

    submission_profile = summarize_submission_file(var_records)
    print(
        f"submission_summary: {submission_profile['n_records']:,} records, "
        f"{submission_profile['n_submitters']:,} submitters, "
        f"{submission_profile['placeholder_phenotype']:,} with a 'not provided' phenotype"
    )

    review_status_rows = summarize_review_status(var_records)
    print("\nReviewStatus distribution in submission_summary.txt (per-record):")
    for r in review_status_rows:
        stars = r["stars"] if r["stars"] is not None else "?"
        print(f"{stars:>2}★ | {r['count']:>9,} | {r['status']}")

    (
        results,
        pairs_ge3,
        pairs_2star_concordance,
        pairs_remaining,
        pair_submitters,
        snv_pair_set,
        scn1a_sets,
    ) = compute_star_data(clinvar_tsv, var_records, map_to_mondo, mondo_labels, variant_genes)
    print(
        f"\nGene-disease pair tiers: {len(pairs_ge3):,} >=3-star, "
        f"{len(pairs_2star_concordance):,} 2-star-concordance-only, "
        f"{len(pairs_remaining):,} remaining (unrescued 0/1-star-only)"
    )

    print(f"{'star_min':>8} | {'variants':>10} | {'gene_disease_pairs':>19}")
    for star in STAR_LEVELS:
        r = results[star]
        print(f"{star:>8} | {r['variants']:>10,} | {r['gene_disease_pairs']:>19,}")

    clingen_summary, clingen_only_rows = summarize_clingen_coverage(pair_submitters, mondo_labels)
    print(
        f"\nSingle-submitter pairs: {clingen_summary['single_submitter_pairs']:,} total, "
        f"{clingen_summary['single_submitter_clingen_pairs']:,} submitted solely by ClinGen"
    )
    print(f"Pairs with ClinGen as the only submitter(s) at all: {clingen_summary['clingen_only_pairs']:,}")

    monarch = load_monarch_gene_disease(args.data_dir)
    monarch_comparison = build_monarch_comparison(
        pairs_ge3, pairs_2star_concordance, pairs_remaining, monarch, pair_submitters
    )
    monarch_comparison["sssom"] = sssom_label_collisions(args.data_dir)
    print(
        f"  SSSOM label collisions: {monarch_comparison['sssom']['n_collisions']:,} of "
        f"{monarch_comparison['sssom']['n_labels']:,} source labels claimed by >1 MONDO class"
    )
    print(
        f"\nMonarch KG reconciliation: {monarch_comparison['n_monarch']:,} curated gene-disease pairs, "
        f"{monarch_comparison['n_clinvar']:,} ClinVar-derived pairs"
    )
    print(
        f"  in both: {monarch_comparison['n_both']:,} | "
        f"ClinVar-only: {monarch_comparison['n_clinvar_only']:,} | "
        f"Monarch-only: {monarch_comparison['n_monarch_only']:,}"
    )
    for kind, n in monarch_comparison["kinship_counts"].most_common():
        print(f"    ClinVar-only, {kind}: {n:,}")
    print(f"  ClinVar pairs on deprecated MONDO terms: {len(monarch_comparison['deprecated_rows']):,}")

    gene_attribution = compare_gene_attribution(clinvar_tsv, args.data_dir)
    print(
        f"\nGene attribution: GENEINFO {gene_attribution['geneinfo_attributions']:,} attributions "
        f"({gene_attribution['geneinfo_artifact']:,} antisense/ORF) -> variant_summary "
        f"{gene_attribution['vs_attributions']:,} ({gene_attribution['vs_artifact']:,} antisense/ORF); "
        f"{gene_attribution['multi_resolved_to_one']:,}/{gene_attribution['multi_gene_variants']:,} "
        f"multi-gene variants collapse to one; GeneID=-1: {gene_attribution['vs_minus1']:,}"
    )

    # the report runs from analysis/, so the KGX artifacts sit alongside data/
    scn1a_subgraph = load_scn1a_subgraph(args.data_dir)
    emitted = load_emitted_summary(args.data_dir.parent / "output")
    if emitted["available"]:
        print(f"Emitted artifacts: {emitted['nodes']:,} nodes, {emitted['edges']:,} edges")
    else:
        print("Emitted artifacts not found -- run `just transform` first for section counts")

    hp_labels = load_hp_labels(args.data_dir)
    phenotype_profile = build_phenotype_profile(clinvar_tsv, var_records, map_to_mondo, mondo_labels, hp_labels)
    pp = phenotype_profile
    print(
        f"Phenotype: {pp['n_distinct_hp']:,} distinct HPO terms; "
        f"{pp['hist_all'][0]:,} VCF variants with none; "
        f"{len(pp['rows']):,} diseases with >=2 HPO-bearing variants; "
        f"{pp['n_obsolete']:,} obsolete HPO terms in use"
    )
    consist = [r["consistency"] for r in pp["rows"]]
    if consist:
        print(f"  disease HPO consistency: mean {sum(consist)/len(consist):.3f}, "
              f"perfect (1.0) on {sum(1 for c in consist if c == 1.0):,}/{len(consist):,} diseases")

    pubmed_variants = load_pubmed_variants(args.data_dir)
    evidence_tiers = build_evidence_tiers(clinvar_tsv, var_records, map_to_mondo, variant_genes, pubmed_variants)
    print(
        f"Evidence tiers: stage3 {evidence_tiers['stage_variants'][3]:,} -> "
        f"stage4 {evidence_tiers['stage_variants'][4]:,} -> stage5 {evidence_tiers['stage_variants'][5]:,} variants; "
        f"{evidence_tiers['pairs_causes']:,} pairs @causes, {evidence_tiers['pairs_assoc']:,} @associated_with"
    )

    filter_cube = build_filter_cube(clinvar_tsv, var_records, map_to_mondo, variant_genes)
    print(
        f"Filter cube: {len(filter_cube['variants']):,} variant cells, "
        f"{len(filter_cube['edges']):,} edge cells"
    )

    strchive_intervals = load_strchive_intervals(args.data_dir)
    cube, clnvc_types, pair_list = build_clnsig_cube(
        clinvar_tsv, var_records, map_to_mondo, mondo_labels, strchive_intervals, variant_genes
    )
    total = sum(r["count"] for r in cube)
    by_clnsig: Counter = Counter()
    for r in cube:
        by_clnsig[r["clnsig"]] += r["count"]
    print(f"\nCLNSIG distribution (all {total:,} variants):")
    for bucket in CLNSIG_BUCKETS:
        n = by_clnsig[bucket]
        print(f"{bucket:>15} | {n:>10,} | {100 * n / total:>5.1f}%")
    by_star: Counter = Counter()
    for r in cube:
        by_star[r["star"]] += r["count"]
    print("\nCLNREVSTAT star distribution (variant-level, all variants):")
    for level in STAR_LEVELS:
        n = by_star[level]
        print(f"{level:>15}★ | {n:>10,} | {100 * n / total:>5.1f}%")
    in_prod = sum(r["count"] for r in cube if r["production"])
    print(
        f"\nVariants passing the production filter (>={PRODUCTION_STAR_MIN}★ per-record or "
        f">={MIN_CONCORDANT_SUBMITTERS} concordant): {in_prod:,} "
        f"(section {S.star_cutoff} star_min={PRODUCTION_STAR_MIN} variants: {results[PRODUCTION_STAR_MIN]['variants']:,}; "
        f"section {S.star_cutoff} '2★ computed' variants: {results['2c']['variants']:,})"
    )
    all_pairs = {p["i"] for r in cube for p in r["pairs"]}
    print(f"\nGene-disease pairs across all filters: {len(all_pairs):,}")
    n_strchive_loci = sum(len(v) for v in strchive_intervals.values())
    n_strchive_variants = sum(r["count"] for r in cube if r["strchive"])
    print(f"\nSTRchive loci loaded: {n_strchive_loci:,}; variants overlapping a locus: {n_strchive_variants:,}")

    sv_summary = build_sv_summary(args.data_dir, map_to_mondo)
    print(f"\nStructural variant types (variant_summary.txt, GRCh38, {sv_summary['total_variants']:,} rows):")
    for r in sv_summary["sv_type_rows"]:
        if r["count"]:
            print(f"{r['type']:>20} | {r['count']:>7,} | resolved={r['resolved']:>6,} ({r['resolved_pct']:.1f}%)")
        else:
            print(f"{r['type']:>20} | {r['count']:>7,}")
    resolved_count = sum(1 for r in sv_summary["sv_rows"] if r["resolved"])
    sv_only_pairs = sv_only_gene_disease_pairs(sv_summary["gene_disease_pairs"], snv_pair_set, mondo_labels)
    print(
        f"\nSingle-gene CNV/SV gene-disease pairs: {len(sv_summary['gene_disease_pairs']):,} total, "
        f"{len(sv_only_pairs):,} with no SNV/indel evidence anywhere in clinvar.vcf.gz"
    )
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
            clinvar_tsv, var_records, map_to_mondo, mondo_labels, variant_genes, gene, example.get("diseases")
        )
        gene_model_colors, gene_model_legend = assign_gene_model_colors(gene_model_variants)
        gene_model_svg = render_gene_model_svg(gene, gene_model, gene_model_variants, gene_model_colors)
        print(f"\n{gene} illustrative gene model: {len(gene_model_variants):,} SNVs, {len(gene_model_legend):,} disease colors")
        domain_info = GENE_DOMAIN_DIAGRAMS.get(gene)
        gene_model_results.append(
            {
                "heading": example["heading"],
                "blurb": example["blurb"],
                "svg": gene_model_svg,
                "legend": gene_model_legend,
                "count": len(gene_model_variants),
                "domain_svg": domain_info["svg_fn"]() if domain_info else None,
                "domain_uniprot_id": domain_info["uniprot_id"] if domain_info else None,
                "domain_caption": domain_info["caption"] if domain_info else None,
                "phenotype_svg": diagram_bams_fshd2_venn() if gene == "SMCHD1" else None,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        render_html(
            results,
            pairs_ge3,
            pairs_2star_concordance,
            pairs_remaining,
            cube,
            clnvc_types,
            pair_list,
            sv_summary,
            example_variants,
            gene_model_results,
            clingen_summary,
            clingen_only_rows,
            sv_only_pairs,
            review_status_rows,
            monarch_comparison,
            gene_attribution,
            variant_genes,
            submission_profile,
            filter_cube,
            phenotype_profile,
            evidence_tiers,
            emitted,
            scn1a_subgraph,
        )
    )
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
