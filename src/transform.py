import koza

from clinvar_helpers import (
    make_medgen_to_mondo_map,
    make_mondo_map,
    load_pubmed_variants,
    make_variant_gene_map,
    make_variant_record_map,
    process_row,
)

# File paths to accessory data
sub_path = "./data/submission_summary.txt.gz"
sssom_path = "./data/mondo.sssom.tsv"
medgen_path = "./data/MedGenIDMappings.txt.gz"
variant_summary_path = "./data/variant_summary.txt.gz"
citations_path = "./data/var_citations.txt"

@koza.on_data_begin()
def load_auxiliary_data(koza_transform):
    """Load auxiliary data files into koza_transform.state for use during transform."""
    # File paths to acessory data
    sub_path = "./data/submission_summary.txt.gz"
    sssom_path = "./data/mondo.sssom.tsv"
    medgen_path = "./data/MedGenIDMappings.txt.gz"

# Map records to each clinvar variant id
var_records = make_variant_record_map(sub_path)

# Make general map back to mondo terms
map_to_mondo = make_mondo_map(sssom_path)

# Make medgen to mondo map and merge
medgen_to_mondo = make_medgen_to_mondo_map(medgen_path)
map_to_mondo.update(medgen_to_mondo)

# ClinVar's own per-variant gene attribution -- see make_variant_gene_map for why this
# replaces the VCF's positional GENEINFO field as the source of variant-gene edges
variant_genes = make_variant_gene_map(variant_summary_path)

# VariationIDs with a PubMed citation -- the evidence behind the <=1-star
# associated_with tier (see publication_star_max in clinvar_helpers)
pubmed_variants = load_pubmed_variants(citations_path)


@koza.transform_record()
def transform(koza_transform, row):
    entities = process_row(row, var_records, map_to_mondo, variant_genes, pubmed_variants)
    if entities:
        koza_transform.write(*entities)
