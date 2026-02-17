# ClinVar

ClinVar aggregates information about genomic variation and its relationship to human health. 

## Ingest files and how nodes and edges are generated

Two files downloaded from ClinVar are leveraged in this ingest:

- **clinvar.vcf** — Contains a single line per ClinVar variant with each variant's associated terms reported in the INFO column and grouped by which submission record(s) they originated from. https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38 (hg38 genome version)
- **submission_summary.txt** — Contains a single line per variant record. These records contain in-depth information about the variant in question that we can leverage in the ingest process. Multiple records often exist per one variant. https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

### Variant nodes (SequenceVariant)

SequenceVariant nodes are created from ClinVar variants that are deemed Pathogenic or Likely Pathogenic. Additionally, only variants that have a ClinVar review status of 3 or more stars (4 maximum) will be included. This subset corresponds to the most credible set of variants that are pathogenic within the ClinVar dataset.

### Variant to Disease edges (VariantToDiseaseAssociation)

Disease IDs are derived from the ReportedPhenotypeInfo column within the submission_summary.txt file. This column consists of MedGen IDs that we then map to a Mondo ID. Alternatively, if a Mondo ID cannot be found, then the SubmittedPhenotypeInfo column will be used instead. If neither column maps to a Mondo ID then no edge will be made.

Predicates are derived from the ClinicalSignificance column within the submission_summary.txt file. Currently only Pathogenic and Likely pathogenic are included as "causes" and "associated_with_increased_likelihood_of" respectively.

### Variant to Phenotype edges (VariantToPhenotypicFeatureAssociation)

These edges are only created if a Variant to Disease edge can be made. The phenotype terms themselves are derived from the INFO column of the clinvar.vcf file. The information within this column is reported as groups of terms that map back to an individual record(s) from the submission_summary.txt file. Any group of terms that contains the disease ID from the Variant to Disease edge will have Variant to Phenotype edges created for all reported Human Phenotype Ontology terms reported within the group.

The predicate for these edges is "contributes_to".

### Variant to Gene edges (VariantToGeneAssociation)

These edges are created only if a Variant to Disease edge can be made. Gene symbols are derived from the INFO column within the clinvar.vcf file and gene symbols are mapped to NCBI genes.

The predicate "is_sequence_variant_of" is used. Sequence Ontology terms are also reported within the INFO column pertaining to the variant's "molecular consequence" (MC subfield within INFO). These terms are recorded in the "type" slot for the SequenceVariant node that is created.

