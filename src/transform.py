import koza

from clinvar_helpers import (
    make_medgen_to_mondo_map,
    make_mondo_map,
    make_variant_record_map,
    process_row,
)

# File paths to accessory data
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


@koza.transform_record()
def transform(koza_transform, row):
    entities = process_row(row, var_records, map_to_mondo)
    if entities:
        koza_transform.write(*entities)
