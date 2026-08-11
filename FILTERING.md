# Filtering approach

What this ingest keeps, what it discards, which input file each decision reads, and which
component in this repo implements it. Every claim below is traceable to a named function.

## Input files

| File | Grain | Downloaded by | Read by |
|---|---|---|---|
| `clinvar.vcf.gz` → `data/clinvar.tsv` | one variant | `download.yaml` → `scripts/vcf_to_tsv.py` (`just preprocess`) | `src/transform.yaml` reader |
| `submission_summary.txt.gz` | one submission (variant × lab) | `download.yaml` | `make_variant_record_map()` |
| `variant_summary.txt.gz` | one variant × genome build | `download.yaml` | `make_variant_gene_map()` |
| `mondo.sssom.tsv` | one xref mapping | `download.yaml` | `make_mondo_map()` |
| `MedGenIDMappings.txt.gz` | one MedGen mapping | `download.yaml` | `make_medgen_to_mondo_map()` |

## Components

| Component | Role |
|---|---|
| `src/transform.yaml` | Koza config: declares the reader over `clinvar.tsv`, and which node/edge properties reach the output TSVs |
| `src/transform.py` | Koza entry point. Loads the four auxiliary maps once at module import, then delegates every row to `process_row()` |
| `src/clinvar_helpers.py` | **All filtering logic lives here.** Constants, mapping builders, and `process_row()` |
| `scripts/vcf_to_tsv.py` | Flattens the VCF INFO column into 43 TSV columns |
| `src/versions.py` / `scripts/write_metadata.py` | Emit `output/release-metadata.yaml`; no filtering |

## Controlling constants

All in `src/clinvar_helpers.py`:

```python
var2disease_star_min      = 3               # per-submission review-status floor
min_concordant_submitters = 2               # independent-agreement rescue threshold
aggregate_star_min        = 2               # ClinVar's own variant-level CLNREVSTAT floor
ASSEMBLY_PREFERENCE       = ("GRCh38", "GRCh37")
predicate_map             = {...}           # which ClinicalSignificance values survive at all
review_star_map           = {...}           # ReviewStatus text -> 0-4 stars
```

## The pipeline, in the order filters actually apply

Koza streams `clinvar.tsv` one variant at a time into `process_row()`. A variant either yields a
`SequenceVariant` node plus its edges, or yields nothing.

### Stage 1 — variant must have submission records
**File:** `submission_summary.txt.gz` · **Code:** `process_row()`

```python
if varid not in var_records:
    return []
```

A variant present in the VCF but absent from `submission_summary.txt` is dropped whole. No node.

### Stage 2 — classification must be pathogenic-family
**File:** `submission_summary.txt.gz` (`ClinicalSignificance`) · **Code:** `variant_records_to_disease()`

Each submission record is discarded unless its `ClinicalSignificance` is a key of `predicate_map`:

| ClinicalSignificance | Predicate emitted |
|---|---|
| `Pathogenic`, `Pathogenic, low penetrance`, `Pathogenic/Likely pathogenic`, `Likely pathogenic`, `Likely pathogenic, low penetrance` | `biolink:causes` |

All five assert causation. "Likely pathogenic" expresses the curator's confidence in that
causal claim rather than a weaker kind of relationship, so it is not downgraded to
`associated_with_increased_likelihood_of`. Exactly one predicate is emitted per
(variant, disease); the full set of submitted classifications is preserved in
`original_predicate`.

Everything else — VUS, Benign, Likely benign, Conflicting, drug response — contributes nothing.
Matching is on the **exact string**.

### Stage 3 — disease must resolve to a MONDO id
**Files:** `submission_summary.txt.gz` + both mapping files · **Code:** `variant_records_to_disease()`

Per record, in order:

1. `ReportedPhenotypeInfo` → `MedGen:<CUI>` → MONDO, via `make_medgen_to_mondo_map()`.
2. Only if that yields nothing for the record: `SubmittedPhenotypeInfo` → MONDO, via
   `make_mondo_map()` (OMIM / Orphanet / MeSH / direct MONDO).

A record whose phenotype resolves to no MONDO id contributes nothing. This silently removes a
large fraction of the file: `ReportedPhenotypeInfo` is the placeholder `C3661900:not provided` on
about 28% of all submission records, and that concept maps to no MONDO term.

### Stage 4 — review status OR concordance OR ClinVar's aggregate
**File:** `submission_summary.txt.gz` (`ReviewStatus`, `Submitter`) · **Code:**
`variant_records_to_disease()` + `concordant_disease_pairs()`

This is the core inclusion rule. For each (record, disease):

```python
if not accept_on_aggregate and stars < star_min and (disease, clinsig) not in concordant_pairs:
    continue
```

A (variant, disease) survives if **any** of three paths holds:

- **Path A — per-record review status.** Some record for it reaches
  `review_star_map[ReviewStatus] >= 3`, i.e. `reviewed by expert panel` (3) or
  `practice guideline` (4).
- **Path B — concordance rescue.** `concordant_disease_pairs()` finds **≥2 distinct
  `Submitter` values** asserting the *same* MONDO disease with the *exact same*
  `ClinicalSignificance` string, at any star level.
- **Path C — ClinVar's aggregate review status.** The variant's `CLNREVSTAT` (from the VCF,
  not the submission file) reaches `>= aggregate_star_min`, i.e. ClinVar itself scored it
  *"criteria provided, multiple submitters, no conflicts"* or better. This is the only place
  the 2-star tier exists — it is a statement about agreement *between* records, so no single
  record can carry it, and paths A and B structurally cannot see it.

Path C dominates: it accounts for roughly 95k of the ~102k variants the ingest emits. Paths A
and B together contribute only ~7k that C would not have found on its own.

Two consequences that are easy to misread:

- **Concordance is computed per variant, not per gene-disease pair.** `concordant_disease_pairs()`
  receives one variant's records. Evidence spread thinly across many variants of the same gene
  never accumulates.
- **`Pathogenic` + `Likely pathogenic` do not corroborate each other**, because the group key is
  the raw string. Two labs agreeing the variant is disease-causing fail the test if they picked
  different tiers.

Note also that per-record `ReviewStatus` never carries ClinVar's 2-star tier — that value is an
aggregate computed *across* submitters and appears only in variant-level files. A per-record
`star_min=2` filter is therefore identical to `star_min=3`.

### Stage 5 — disease must also appear in the VCF's `CLNDISDB`
**File:** `clinvar.tsv` (`CLNDISDB`) · **Code:** `parse_CLNDISDB()` → `map_CLNDISDB_to_mondo()` →
`map_mondo_to_hp()`

```python
corroborated = map_mondo_to_hp(diss_info, disease_ids)
if len(corroborated) == 0:
    return []
```

**Despite the function name this is not an "HPO overlap" requirement.** `map_mondo_to_hp()`
creates an entry for any disease that survived stage 4 and also appears in a `CLNDISDB` group —
whether or not that group carries a single HPO term. So the real requirement is: *at least one
surviving disease must be echoed in the variant's `CLNDISDB`*. It is a cross-check between the
two input files, and it outlives the removal of phenotype edges.

It is also **all-or-nothing per variant**. If no surviving disease is echoed, the variant is
dropped entirely — node, gene edge and every disease edge. If at least one is, *all* surviving
disease edges are emitted, including ones not echoed in `CLNDISDB`.

In practice this gate is nearly inert: of 102,074 variants that clear stage 4, 101,909 emit. It
gates hard in principle and almost never in fact.

### Stage 5b — `SequenceVariant.type`
**File:** `clinvar.tsv` (`CLNVCSO`)

`type` carries the variant **class** — `SO:0001483` SNV, `SO:0000159` deletion, `SO:1000035`
duplication — not the molecular consequence from `MC` (missense, frameshift). Consequence is a
property of variant × transcript, not of the variant. `MC` is no longer read by the ingest.

### Stage 6 — gene attribution
**File:** `variant_summary.txt.gz` · **Code:** `make_variant_gene_map()`

Not a filter on the variant — a filter on which gene edge is emitted.

- One row per variant per genome build. `ASSEMBLY_PREFERENCE = ("GRCh38", "GRCh37")` selects the
  most-preferred build present per `VariationID`, never more than one. Filtering to GRCh38 alone
  would discard ~56k variants ClinVar has only ever placed on GRCh37, over 99% of them structural.
- `GeneID == -1` means ClinVar declines to attribute a gene. **No `GENEINFO` fallback** — no gene
  edge is emitted, and the variant keeps its disease and phenotype edges.

> **The VCF's `GENEINFO` field is deliberately not used for this.** It lists every locus whose span
> covers the variant — a positional annotation, not a causal claim. Reading it as causal gave
> antisense transcripts, readthrough fusions and `LOC` placeholders their own gene edges, from
> which they inherited the real gene's diseases and submitter roster (`CFTR-AS1` → cystic fibrosis
> with 76 submitters behind it). Switching to the curated attribution removed 6,570 such edges and
> left disease/phenotype edge counts untouched.

### Stage 7 — phenotype edges: **removed**

The ingest no longer emits `VariantToPhenotypicFeatureAssociation`.

ClinVar's HPO ids appear only inside `CLNDISDB` groups, alongside the MONDO / MedGen / OMIM /
Orphanet ids for the *same* condition — cross-references naming that condition in HPO's
vocabulary. All 214,880 HPO-bearing groups also carry a disease id; none is HPO-only. They are
not supplied by submitters either (12 of 6,374,759 records put an HPO id in
`ReportedPhenotypeInfo`): a submitter names a condition, ClinVar normalises it to a MedGen
concept, and MedGen's cross-reference table supplies the HPO id.

Measured consequence: across 585 diseases with ≥2 HPO-bearing variants, mean consistency is
**0.9994**, and **584 of 585** show zero variation — every variant of a disease carries an
identical HPO set, because the set is a property of the disease concept. A `contributes_to`
edge built from them restated the variant's disease edge in a second vocabulary rather than
asserting a distinct phenotype.

## What comes out

| Output | Predicate | Gated by |
|---|---|---|
| `SequenceVariant` node | — | stages 1–5 |
| `VariantToGeneAssociation` | `biolink:is_sequence_variant_of` | stages 1–5, plus stage 6 |
| `VariantToDiseaseAssociation` | `biolink:causes` | stages 1–5 |

Current release: 101,909 nodes; 101,661 gene edges and 220,824 disease edges
(140,324 `causes` + 80,500 `associated_with_increased_likelihood_of`). No phenotype edges.

## Known limitations of this approach

1. **Inclusion is decided per variant, not per gene-disease pair** (stage 4). A pair enters the KG
   when some single variant clears the bar alone; pooled evidence across variants never counts.
   Path C partly sidesteps this by using ClinVar's own cross-submitter aggregate, but only within
   a variant — never across a gene-disease pair's variants.
2. **Disease terms are matched by exact MONDO id.** MONDO carries co-existing terms for
   overlapping entities, so labs describing the same condition under synonymous terms cannot
   corroborate each other.
3. **~28% of submission records are invisible** because their phenotype is `not provided`
   (stage 3), which biases every submitter count downward — worst for the highest-volume labs.
4. **Deprecated MONDO terms are not filtered**; the mapping files still point at some obsoleted ids.

See `analysis/clinvar_report.py` (sections 5–8 of the generated report) for the measurements behind
each of these.
