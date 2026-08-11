# ClinVar

ClinVar aggregates information about genomic variation and its relationship to human health. 

## Ingest files and how nodes and edges are generated

Three files downloaded from ClinVar are leveraged in this ingest:

- **clinvar.vcf** — A single line per ClinVar variant, with the variant's associated terms in the INFO column, grouped by which submission record(s) they originated from. Defines which variants the ingest considers at all. https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38 (hg38)
- **submission_summary.txt** — A single line per submission record: one lab's assertion about one variant. Multiple records usually exist per variant. This is the evidence the inclusion filter runs over. https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
- **variant_summary.txt** — A single line per variant per genome build, carrying ClinVar's own curated `GeneID`/`GeneSymbol`/`HGNC_ID` for the variant. Used for gene attribution; see the Variant to Gene section below for why the VCF's `GENEINFO` is not. https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

See [FILTERING.md](FILTERING.md) for the full stage-by-stage filtering pipeline.

### Variant nodes (SequenceVariant)

SequenceVariant nodes are created from ClinVar variants deemed Pathogenic or Likely Pathogenic. A (variant, disease) is accepted if **any** of three paths holds:

1. **Per-record review status** — some submission record reaches 3 or more stars (4 maximum), i.e. reviewed by expert panel or practice guideline.
2. **Multi-submitter concordance** — at least 2 independent submitters agree on the same disease with the same classification, regardless of each submitter's own review status.
3. **ClinVar's aggregate review status** — the variant's `CLNREVSTAT` in the VCF reaches 2 or more stars ("criteria provided, multiple submitters, no conflicts").

Path 3 exists because the 2-star tier is a statement about agreement *between* records, so no individual submission record can carry it; paths 1 and 2 structurally cannot see it. In the current release path 3 accounts for the large majority of accepted variants.

The `type` slot on the node carries the variant **class** (SNV, deletion, duplication — from `CLNVCSO`), not the molecular consequence (missense, frameshift), which is a property of variant × transcript rather than of the variant.

### Variant to Disease edges (VariantToDiseaseAssociation)

Disease IDs are derived from the ReportedPhenotypeInfo column within the submission_summary.txt file. This column consists of MedGen IDs that we then map to a Mondo ID. Alternatively, if a Mondo ID cannot be found, then the SubmittedPhenotypeInfo column will be used instead. If neither column maps to a Mondo ID then no edge will be made.

Predicates are derived from the ClinicalSignificance column within the submission_summary.txt file. Pathogenic, Pathogenic/Likely pathogenic and Likely pathogenic (including the low-penetrance variants of each) all map to **"causes"** — "Likely pathogenic" expresses the curator's confidence in a causal claim rather than a weaker kind of relationship. Exactly one predicate is emitted per (variant, disease); every submitted classification is preserved in `original_predicate`. A gene-disease pair is included if it meets the review-status or multi-submitter concordance rule described above.

### Variant to Phenotype edges — not emitted

This ingest does **not** produce `VariantToPhenotypicFeatureAssociation` edges.

ClinVar's HPO ids appear only inside `CLNDISDB` groups, alongside the MONDO/MedGen/OMIM/Orphanet ids for the *same* condition — they are cross-references naming that condition in HPO's vocabulary, not observed phenotypes. They are supplied by MedGen's cross-reference table rather than by submitters (12 of 6.37M submission records carry an HPO id directly). Measured across 585 diseases with 2 or more HPO-bearing variants, mean consistency is 0.9994 and 584 of 585 show *zero* variation between variants — the term set is a property of the disease, not the variant. A `contributes_to` edge built from them restated the disease edge in a second vocabulary rather than asserting a distinct phenotype.

### Variant to Gene edges (VariantToGeneAssociation)

These edges are created only if a Variant to Disease edge can be made. The gene comes from **ClinVar's own curated attribution** in variant_summary.txt (`GeneID`), not from the VCF's `GENEINFO` field.

`GENEINFO` is populated *positionally* — it lists every locus whose span covers the variant, which is what the field is documented to do. Building gene edges from it gave antisense transcripts, readthrough fusions and NCBI `LOC` placeholders their own edges, from which they inherited the causal gene's diseases and submitters (`CFTR-AS1` appeared as a cystic fibrosis gene with 76 submitters behind it). Switching to the curated attribution removed 6,570 such edges while leaving disease edge counts untouched.

Where ClinVar declines to attribute a gene (`GeneID == -1`) no gene edge is emitted and there is deliberately no `GENEINFO` fallback; the variant keeps its disease edges.

The predicate "is_sequence_variant_of" is used.

