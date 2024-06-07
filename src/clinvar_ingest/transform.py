import uuid  # For generating UUIDs for associations

from biolink_model.datamodel.pydanticmodel_v2 import *  # Replace * with any necessary data classes from the Biolink Model
from koza.cli_utils import get_koza_app

koza_app = get_koza_app("clinvar_variant")

CONTRIBUTES_TO = "biolink:contributes_to"
HAS_PHENOTYPE = "biolink:has_phenotype"
IS_SEQUENCE_VARIANT_OF = 'biolink:is_sequence_variant_of'

# def parse_vcf_info_column(info_column):
#     rel_columns = {"GENEINFO":'',      # Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)
#                    "ORIGIN":'',       # Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">

#                    "CLNDISDB":'',     # Tag-value pairs of disease database name and identifier submitted for germline classifications
#                    "CLNREVSTAT":'',   # ClinVar review status of germline classification for the Variation ID

#                    "ONCDISDB":'',     # Tag-value pairs of disease database name and identifier submitted for oncogenicity classifications, e.g. MedGen:NNNNNN">
#                    "ONCREVSTAT":'',   # ClinVar review status of oncogenicity classification for the Variation ID">"

#                    "SCIDISDB":'',     # Tag-value pairs of disease database name and identifier submitted for somatic clinial impact classifications, e.g. MedGen:NNNNNN
#                    "SCIREVSTAT":'',   # ClinVar review status of somatic clinical impact for the Variation ID
#                    "CLNHGVS":'',      # Variant info / name

# Only ontology terms we are interested in are HPO and MONDO
#'AF_ESP=0.00008;AF_EXAC=0.00002;ALLELEID=980658;CLNDISDB=MedGen:C3661900|Human_Phenotype_Ontology:HP:0000730,Human_Phenotype_Ontology:HP:0001249,Human_Phenotype_Ontology:HP:0001267,Human_Phenotype_Ontology:HP:0001286,Human_Phenotype_Ontology:HP:0002122,Human_Phenotype_Ontology:HP:0002192,Human_Phenotype_Ontology:HP:0002316,Human_Phenotype_Ontology:HP:0002382,Human_Phenotype_Ontology:HP:0002386,Human_Phenotype_Ontology:HP:0002402,Human_Phenotype_Ontology:HP:0002458,Human_Phenotype_Ontology:HP:0002482,Human_Phenotype_Ontology:HP:0002499,Human_Phenotype_Ontology:HP:0002543,Human_Phenotype_Ontology:HP:0003767,Human_Phenotype_Ontology:HP:0006833,Human_Phenotype_Ontology:HP:0007154,Human_Phenotype_Ontology:HP:0007176,Human_Phenotype_Ontology:HP:0007180,MONDO:MONDO:0001071,MeSH:D008607,MedGen:C3714756|Human_Phenotype_Ontology:HP:0001380,Human_Phenotype_Ontology:HP:0001383,Human_Phenotype_Ontology:HP:0001388,Human_Phenotype_Ontology:HP:0002771,MedGen:C0086437;CLNDN=not_provided|Intellectual_disability|Joint_laxity;CLNHGVS=NC_000001.11:g.2304067C>T;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=SKI:6497;MC=SO:0001583|missense_variant;ORIGIN=1;RS=367916348'}


# TO DO (how to link to existing gene_ids already in database)
def make_genes_from_row(gene_list):
    if gene_list == '.':
        return []
    gene_ids = []
    genes = gene_list.split('|')
    for gene in genes:
        values = gene.split(':')[1:]  # sometimes there's more than one Entrez ID after the symbol
        for value in values:
            gene_ids.append("NCBIGene:{}".format(value))
    return gene_ids

    # ZNRF2:223082|LOC105375218:105375218|LOC129998188:12999818


def extract_ids(prefix, value):
    ids = []
    groups = value.split('|')
    for group in groups:
        items = group.split(',')
        for item in items:
            # replace the MONDO banana if necessary
            item.replace('MONDO:MONDO:', 'MONDO:')
            # deal with HP DB name + prefix
            item.replace('Human_Phenotype_Ontology:HP:', 'HP:')
            if item.startswith(prefix):
                ids.append(item)
    return ids


while (row := koza_app.get_row()) is not None:
    # Code to transform each row of data
    # For more information, see https://koza.monarchinitiative.org/Ingests/transform

    entities = []

    gene_ids = make_genes_from_row(row["GENEINFO"])

    seq_var = SequenceVariant(
        id="CLINVAR:{}".format(row["ID"]),
        name=row["CLNHGVS"],
        xref=["DBSNP:{}".format(row["RS"])],
        has_gene=gene_ids,
        in_taxon=["NCBITaxon:9606"],
        in_taxon_label="Homo sapiens",
        # type? could this be a SO term?
        # has_bioligical_sequence  do we want it? not so sure
    )

    entities.append(seq_var)

    for gene_id in gene_ids:
        entities.append(
            VariantToGeneAssociation(
                id=str(uuid.uuid4()),
                subject=seq_var.id,
                predicate=IS_SEQUENCE_VARIANT_OF,  # TODO: more specific predicates might be possible, is_missense_variant_of etc
                object=gene_id,
                primary_knowledge_source="infores:clinvar",
                aggregator_knowledge_source=["infores:monarchinitiative"],
                knowledge_level=KnowledgeLevelEnum.knowledge_assertion,  # TODO: we should confirm this
                agent_type=AgentTypeEnum.manual_agent,  # TODO: we should confirm this
            )
        )

    mondo_ids = extract_ids('MONDO', row['CLNDISDB'])
    if len(mondo_ids) == 0:
        continue  # exit without writing the sequence variant entity if there's no disease association

    for mondo_id in mondo_ids:
        entities.append(
            VariantToDiseaseAssociation(
                id=str(uuid.uuid4()),
                subject=seq_var.id,
                predicate=CONTRIBUTES_TO,
                object=mondo_id,
                primary_knowledge_source="infores:clinvar",
                aggregator_knowledge_source=["infores:monarchinitiative"],
                knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
                agent_type=AgentTypeEnum.manual_agent,
            )
        )

    for hp_id in extract_ids('Human_Phenotype_Ontology:HP:', row['CLNDISDB']):
        entities.append(
            VariantToPhenotypicFeatureAssociation(
                id=str(uuid.uuid4()),
                subject=seq_var.id,
                predicate=CONTRIBUTES_TO,
                object=hp_id,
                primary_knowledge_source="infores:clinvar",
                aggregator_knowledge_source=["infores:monarchinitiative"],
                knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
                agent_type=AgentTypeEnum.manual_agent,
            )
        )

    koza_app.write(*entities)
