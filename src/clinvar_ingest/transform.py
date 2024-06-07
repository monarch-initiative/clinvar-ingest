import uuid  # For generating UUIDs for associations

from biolink_model.datamodel.pydanticmodel_v2 import *  # Replace * with any necessary data classes from the Biolink Model
from koza.cli_utils import get_koza_app

koza_app = get_koza_app("clinvar_variant")


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
def make_genes_from_row(gene_row):
        gene_objs = []
        if "|" in gene_row:
             gene_rows = gene_row.split('|')
        else:
             gene_rows = [gene_row]
        
        for gene_info in gene_rows:
            symbol, idnum = gene_info.split(':')[0], gene_info.split(':')[1]
            gene_obj = Gene(id=idnum, 
                            symbol=symbol, 
                            in_taxon="Homo sapians"
                            )
        
        return gene_objs

    #ZNRF2:223082|LOC105375218:105375218|LOC129998188:12999818



while (row := koza_app.get_row()) is not None:
    # Code to transform each row of data
    # For more information, see https://koza.monarchinitiative.org/Ingests/transform
    
    #rel_fields = row["INFO"].split('|')
    
    #print(row)

    # Buggy still
    gen_objs = make_genes_from_row(row["GENEINFO"])

    # Do we need to create a gene object 
    seq_var = SequenceVariant(id="CLINVAR:{}".format(row["ID"]),
                              name=row["CLNHGVS"],
                              xref=row["RS"],
                              has_gene=gen_objs[0])
    
    #print(seq_var)
    
                              #has_gene)

                              #has_attribute
                                # has_biological_sequence
                                # has_gene

                                # in_taxon
                                # in_taxon_label   ##Homo sapiens


                                # iri

                                # provided_by
                                # synonym
                                # type
    
    
    
    
    
    
    
    # entity_a = Entity(
    #     id=f"XMPL:00000{row['example_column_1'].split('_')[-1]}",
    #     name=row["example_column_1"],
    #     category=["biolink:Entity"],
    # )
    # entity_b = Entity(
    #     id=f"XMPL:00000{row['example_column_2'].split('_')[-1]}",
    #     name=row["example_column_2"],
    #     category=["biolink:Entity"],
    # )
    # association = Association(
    #     id=str(uuid.uuid1()),
    #     subject=row["example_column_1"],
    #     predicate=row["example_column_3"],
    #     object=row["example_column_2"],
    #     subject_category="SUBJ",
    #     object_category="OBJ",
    #     category=["biolink:Association"],
    #     knowledge_level="not_provided",
    #     agent_type="not_provided",
    # )
    # koza_app.write(entity_a, entity_b, association)
