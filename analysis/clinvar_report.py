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
   significance / star / variant-type filters and see live variant counts.

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
import json
from collections import Counter
from pathlib import Path

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

CLNSIG_BUCKETS = ["P", "LP", "P/LP", "VUS", "LB", "B", "B/LB", "Conflicting", "Other", "Not classified"]
NOT_CLASSIFIED_VALUES = {
    ".",
    "not_provided",
    "no_classification_provided",
    "no_classification_for_the_single_variant",
    "no_classifications_from_unflagged_records",
}


# ---------------------------------------------------------------------------
# Section 1 & 2: review-star cutoff + multi-submitter concordance rescue
# ---------------------------------------------------------------------------


def load_maps(data_dir: Path):
    var_records = make_variant_record_map(str(data_dir / "submission_summary.txt.gz"))
    map_to_mondo = make_mondo_map(str(data_dir / "mondo.sssom.tsv"))
    medgen_to_mondo = make_medgen_to_mondo_map(str(data_dir / "MedGenIDMappings.txt.gz"))
    map_to_mondo.update(medgen_to_mondo)
    return var_records, map_to_mondo


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


def compute_star_data(clinvar_tsv: Path, var_records: dict, map_to_mondo: dict):
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

            gene_ids, gene_symbols = make_genes_from_row(row["GENEINFO"])
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


def build_clnsig_cube(clinvar_tsv: Path, var_records: dict, map_to_mondo: dict) -> tuple[list[dict], list[str]]:
    """Full pass over clinvar.tsv (every row, not just those with a
    submission_summary match) tallying (star, clnsig_bucket, clnvc), plus --
    for rows with a submission_summary match -- the set of (gene, MONDO
    disease) pairs implied by that variant (any Pathogenic/Likely-pathogenic
    submission record, any star; same disease-mapping logic as section 1/2's
    star_min=0 pass, since P/LP-family classification is a prerequisite for
    any disease mapping regardless of the variant's own aggregate CLNSIG/star
    shown here). Non-Pathogenic/Likely-pathogenic cells (VUS, B, LB,
    Conflicting, Other, Not classified) will always have an empty pair set --
    production never creates disease edges for those, so that's the expected
    (not a bug) result of toggling those buckets on in the crossfilter.
    """
    counts: Counter = Counter()
    clnvc_types: set = set()
    cell_pairs: dict = {}
    pair_index: dict = {}

    with open(clinvar_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            star = review_star_map.get(row["CLNREVSTAT"])
            if star is None:
                star = 0
            clnsig = classify_clnsig(row["CLNSIG"])
            clnvc = row["CLNVC"] if row["CLNVC"] != "." else "unspecified"
            clnvc_types.add(clnvc)
            cell_key = (star, clnsig, clnvc)
            counts[cell_key] += 1

            varid = row["ID"]
            if varid not in var_records:
                continue
            dis, _, _ = variant_records_to_disease(var_records[varid], map_to_mondo, star_min=0)
            if not dis:
                continue
            gene_ids, _gene_symbols = make_genes_from_row(row["GENEINFO"])
            cell_set = cell_pairs.setdefault(cell_key, set())
            for gene_id in gene_ids:
                for mondo_id in dis:
                    # keyed by (gene_id, disease) only, matching section 1/2's canonical pair
                    # identity -- NCBIGene id is what production actually emits; a handful of
                    # ids (135 observed) carry >1 symbol spelling across records (gene renames),
                    # which would otherwise inflate this count vs section 1's headline number.
                    pkey = (gene_id, mondo_id)
                    idx = pair_index.setdefault(pkey, len(pair_index))
                    cell_set.add(idx)

    cube = [
        {
            "star": star,
            "clnsig": clnsig,
            "clnvc": clnvc,
            "count": n,
            "pairs": sorted(cell_pairs.get((star, clnsig, clnvc), set())),
        }
        for (star, clnsig, clnvc), n in sorted(counts.items())
    ]
    return cube, sorted(clnvc_types)


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

    pairs_json = json.dumps(pairs)
    rescue_json = json.dumps(rescue_rows)
    cube_json = json.dumps(cube)
    clnsig_buckets_json = json.dumps(CLNSIG_BUCKETS)
    clnvc_types_json = json.dumps(clnvc_types)
    total_variants = sum(r["count"] for r in cube)

    rescue_pct = (
        100 * rescue_summary["sub_expert_pairs_rescued"] / rescue_summary["sub_expert_pairs_total"]
        if rescue_summary["sub_expert_pairs_total"]
        else 0
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
</style>
</head>
<body>
<h1>ClinVar exploration report</h1>
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
  <input id="rescue-search" class="search-box" type="text" placeholder="Filter by gene or MONDO id...">
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
</div>

<script>
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
             r.mondo.toLowerCase().indexOf(q) !== -1;
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
  var clnsigBuckets = {clnsig_buckets_json};
  var clnvcTypes = {clnvc_types_json};
  var starLevels = [0, 1, 2, 3, 4];

  var clnsigColors = {{
    "P": "#b91c1c", "LP": "#ea580c", "P/LP": "#c2410c",
    "VUS": "#a16207", "LB": "#16a34a", "B": "#15803d", "B/LB": "#166534",
    "Conflicting": "#7c3aed", "Other": "#64748b", "Not classified": "#94a3b8"
  }};

  var state = {{
    clnsig: new Set(clnsigBuckets),
    star: new Set(starLevels),
    clnvc: new Set(clnvcTypes)
  }};

  function matches(row, dims) {{
    for (var i = 0; i < dims.length; i++) {{
      var d = dims[i];
      if (!state[d].has(row[d])) return false;
    }}
    return true;
  }}

  function total() {{
    var sum = 0;
    cube.forEach(function(row) {{ if (matches(row, ["clnsig", "star", "clnvc"])) sum += row.count; }});
    return sum;
  }}

  function pairsTotal() {{
    var seen = new Set();
    cube.forEach(function(row) {{
      if (!matches(row, ["clnsig", "star", "clnvc"])) return;
      row.pairs.forEach(function(i) {{ seen.add(i); }});
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
    renderPanel("clnsig", clnsigBuckets, ["star", "clnvc"], "xf-clnsig-rows",
      function(v) {{ return v; }}, function(v) {{ return clnsigColors[v] || "#2563eb"; }});
    renderPanel("star", starLevels, ["clnsig", "clnvc"], "xf-star-rows",
      function(v) {{ return v + (v === 1 ? " star" : " stars"); }}, function() {{ return "#2563eb"; }});
    renderPanel("clnvc", clnvcTypes, ["clnsig", "star"], "xf-clnvc-rows",
      function(v) {{ return v.replace(/_/g, " "); }}, function() {{ return "#2563eb"; }});
  }}

  document.body.addEventListener("change", function(e) {{
    if (e.target.tagName !== "INPUT" || !e.target.hasAttribute("data-dim")) return;
    var dim = e.target.getAttribute("data-dim");
    if (!(dim in state)) return;
    var raw = e.target.getAttribute("data-value");
    var value = (dim === "star") ? parseInt(raw, 10) : raw;
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
    var values = dim === "clnsig" ? clnsigBuckets : dim === "star" ? starLevels : clnvcTypes;
    state[dim] = action === "all" ? new Set(values) : new Set();
    renderAll();
  }});

  renderAll();
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
    clinvar_tsv = args.data_dir / "clinvar.tsv"
    results, pairs, rescue_summary, rescue_rows = compute_star_data(clinvar_tsv, var_records, map_to_mondo)

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

    cube, clnvc_types = build_clnsig_cube(clinvar_tsv, var_records, map_to_mondo)
    total = sum(r["count"] for r in cube)
    by_clnsig: Counter = Counter()
    for r in cube:
        by_clnsig[r["clnsig"]] += r["count"]
    print(f"\nCLNSIG distribution (all {total:,} variants):")
    for bucket in CLNSIG_BUCKETS:
        n = by_clnsig[bucket]
        print(f"{bucket:>15} | {n:>10,} | {100 * n / total:>5.1f}%")
    all_pairs = {i for r in cube for i in r["pairs"]}
    print(f"\nGene-disease pairs across all filters: {len(all_pairs):,}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_html(results, pairs, rescue_summary, rescue_rows, cube, clnvc_types))
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
