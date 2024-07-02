"""
An example test file for the transform script.

It uses pytest fixtures to define the input data and the mock koza transform.
The test_example function then tests the output of the transform script.

See the Koza documentation for more information on testing transforms:
https://koza.monarchinitiative.org/Usage/testing/
"""

import pytest
from koza.utils.testing_utils import mock_koza
from biolink_model.datamodel.pydanticmodel_v2 import (
    SequenceVariant,
    VariantToGeneAssociation,
    VariantToPhenotypicFeatureAssociation,
)

# Define the ingest name and transform script path
INGEST_NAME = "clinvar_variant"
TRANSFORM_SCRIPT = "./src/clinvar_ingest/transform.py"


#################################################################
### Rows taken from vcf file that have our desired test cases ###

# Single record, Single mondoID, No HPO, Single Gene
@pytest.fixture
def test_case1_row():
    row = {
        'CHROM':'1',
        'POS':'11128107',
        'ID':'1296989',
        'REF':'G',
        'ALT':'C',
        'QUAL':'.',
        'FILTER':'.',
        'AF_ESP':'.',
        'AF_EXAC':'.',
        'AF_TGP':'.',
        'ALLELEID':'1286779',
        'CLNDN':'Overgrowth_syndrome_and/or_cerebral_malformations_due_to_abnormalities_in_MTOR_pathway_genes',
        'CLNDNINCL':'.',
        'CLNDISDB':'MONDO:MONDO:0100283,MedGen:CN300503',
        'CLNDISDBINCL':'.',
        'CLNHGVS':'NC_000001.11:g.11128107G>C',
        'CLNREVSTAT':'reviewed_by_expert_panel',
        'CLNSIG':'Pathogenic',
        'CLNSIGCONF':'.',
        'CLNSIGINCL':'.',
        'CLNVC':'single_nucleotide_variant',
        'CLNVCSO':'SO:0001483',
        'CLNVI':'.',
        'DBVARID':'.',
        'GENEINFO':'MTOR:2475',
        'MC':'SO:0001583|missense_variant',
        'ONCDN':'.',
        'ONCDNINCL':'.',
        'ONCDISDB':'.',
        'ONCDISDBINCL':'.',
        'ONC':'.',
        'ONCINCL':'.',
        'ONCREVSTAT':'.',
        'ONCCONF':'.',
        'ORIGIN':'1',
        'RS':'587777893',
        'SCIDN':'.',
        'SCIDNINCL':'.',
        'SCIDISDB':'.',
        'SCIDISDBINCL':'.',
        'SCIREVSTAT':'.',
        'SCI':'.',
        'SCIINCL':'.',
        }
    return row

# Single record, Single mondoID, 1 HPO term, Single gene
@pytest.fixture
def test_case2_row():
    row = {
        'CHROM':'1',
        'POS':'171636092',
        'ID':'2505295',
        'REF':'T',
        'ALT':'A',
        'QUAL':'.',
        'FILTER':'.',
        'AF_ESP':'.',
        'AF_EXAC':'.',
        'AF_TGP':'.',
        'ALLELEID':'2670262',
        'CLNDN':'Glaucoma_of_childhood',
        'CLNDNINCL':'.',
        'CLNDISDB':'Human_Phenotype_Ontology:HP:0001087,MONDO:MONDO:0020367,MedGen:C2981140,Orphanet:98977',
        'CLNDISDBINCL':'.',
        'CLNHGVS':'NC_000001.11:g.171636092T>A',
        'CLNREVSTAT':'reviewed_by_expert_panel',
        'CLNSIG':'Likely_pathogenic',
        'CLNSIGCONF':'.',
        'CLNSIGINCL':'.',
        'CLNVC':'single_nucleotide_variant',
        'CLNVCSO':'SO:0001483',
        'CLNVI':'.',
        'DBVARID':'.',
        'GENEINFO':'MYOC:4653',
        'MC':'SO:0001583|missense_variant',
        'ONCDN':'.',
        'ONCDNINCL':'.',
        'ONCDISDB':'.',
        'ONCDISDBINCL':'.',
        'ONC':'.',
        'ONCINCL':'.',
        'ONCREVSTAT':'.',
        'ONCCONF':'.',
        'ORIGIN':'1',
        'RS':'.',
        'SCIDN':'.',
        'SCIDNINCL':'.',
        'SCIDISDB':'.',
        'SCIDISDBINCL':'.',
        'SCIREVSTAT':'.',
        'SCI':'.',
        'SCIINCL':'.',
            }
    return row

# Two records, Two mondoIDs, 2 HPO terms, Single gene
@pytest.fixture
def test_case3_row():
    row = {
        'CHROM':'2',
        'POS':'202464950',
        'ID':'8797',
        'REF':'C',
        'ALT':'G',
        'QUAL':'.',
        'FILTER':'.',
        'AF_ESP':'.',
        'AF_EXAC':'.',
        'AF_TGP':'.',
        'ALLELEID':'23836',
        'CLNDN':'Pulmonary_arterial_hypertension|Pulmonary_hypertension,_primary,_1',
        'CLNDNINCL':'.',
        'CLNDISDB':'Human_Phenotype_Ontology:HP:0002092,Human_Phenotype_Ontology:HP:0006546,MONDO:MONDO:0015924,MeSH:D000081029,MedGen:C2973725,Orphanet:182090|MONDO:MONDO:0024533,MedGen:C4552070,OMIM:178600,Orphanet:422',
        'CLNDISDBINCL':'.',
        'CLNHGVS':'NC_000002.12:g.202464950C>G',
        'CLNREVSTAT':'reviewed_by_expert_panel',
        'CLNSIG':'Pathogenic',
        'CLNSIGCONF':'.',
        'CLNSIGINCL':'.',
        'CLNVC':'single_nucleotide_variant',
        'CLNVCSO':'SO:0001483',
        'CLNVI':'ClinGen:CA278072|OMIM:600799.0003',
        'DBVARID':'.',
        'GENEINFO':'BMPR2:659',
        'MC':'SO:0001587|nonsense',
        'ONCDN':'.',
        'ONCDNINCL':'.',
        'ONCDISDB':'.',
        'ONCDISDBINCL':'.',
        'ONC':'.',
        'ONCINCL':'.',
        'ONCREVSTAT':'.',
        'ONCCONF':'.',
        'ORIGIN':'1',
        'RS':'137852742',
        'SCIDN':'.',
        'SCIDNINCL':'.',
        'SCIDISDB':'.',
        'SCIDISDBINCL':'.',
        'SCIREVSTAT':'.',
        'SCI':'.',
        'SCIINCL':'.',
            }
    return row

# Multiple records, Multiple mondoID, No HPO, Single gene
@pytest.fixture
def test_case4_row():
    row = {
    'CHROM':'1',
    'POS':'11128107',
    'ID':'156702',
    'REF':'G',
    'ALT':'T',
    'QUAL':'.',
    'FILTER':'.',
    'AF_ESP':'.',
    'AF_EXAC':'.',
    'AF_TGP':'.',
    'ALLELEID':'166562',
    'CLNDN':'not_provided|Overgrowth_syndrome_and/or_cerebral_malformations_due_to_abnormalities_in_MTOR_pathway_genes|Macrocephaly-intellectual_disability-neurodevelopmental_disorder-small_thorax_syndrome',
    'CLNDNINCL':'.',
    'CLNDISDB':'MedGen:C3661900|MONDO:MONDO:0100283,MedGen:CN300503|MONDO:MONDO:0014716,MedGen:C4225259,OMIM:616638,Orphanet:457485',
    'CLNDISDBINCL':'.',
    'CLNHGVS':'NC_000001.11:g.11128107G>T',
    'CLNREVSTAT':'reviewed_by_expert_panel',
    'CLNSIG':'Pathogenic',
    'CLNSIGCONF':'.',
    'CLNSIGINCL':'.',
    'CLNVC':'single_nucleotide_variant',
    'CLNVCSO':'SO:0001483',
    'CLNVI':'ClinGen:CA248390',
    'DBVARID':'.',
    'GENEINFO':'MTOR:2475',
    'MC':'SO:0001583|missense_variant',
    'ONCDN':'.',
    'ONCDNINCL':'.',
    'ONCDISDB':'.',
    'ONCDISDBINCL':'.',
    'ONC':'.',
    'ONCINCL':'.',
    'ONCREVSTAT':'.',
    'ONCCONF':'.',
    'ORIGIN':'1',
    'RS':'587777893',
    'SCIDN':'.',
    'SCIDNINCL':'.',
    'SCIDISDB':'.',
    'SCIDISDBINCL':'.',
    'SCIREVSTAT':'.',
    'SCI':'.',
    'SCIINCL':'.',
    }
    return row

# Multiple records, Multiple mondoID, Multiple HPO, Sigle gene
@pytest.fixture
def test_case5_row():
    row = {
    'CHROM':'2',
    'POS':'202464950',
    'ID':'8797',
    'REF':'C',
    'ALT':'G',
    'QUAL':'.',
    'FILTER':'.',
    'AF_ESP':'.',
    'AF_EXAC':'.',
    'AF_TGP':'.',
    'ALLELEID':'23836',
    'CLNDN':'Pulmonary_arterial_hypertension|Pulmonary_hypertension,_primary,_1',
    'CLNDNINCL':'.',
    'CLNDISDB':'Human_Phenotype_Ontology:HP:0002092,Human_Phenotype_Ontology:HP:0006546,MONDO:MONDO:0015924,MeSH:D000081029,MedGen:C2973725,Orphanet:182090|MONDO:MONDO:0024533,MedGen:C4552070,OMIM:178600,Orphanet:422',
    'CLNDISDBINCL':'.',
    'CLNHGVS':'NC_000002.12:g.202464950C>G',
    'CLNREVSTAT':'reviewed_by_expert_panel',
    'CLNSIG':'Pathogenic',
    'CLNSIGCONF':'.',
    'CLNSIGINCL':'.',
    'CLNVC':'single_nucleotide_variant',
    'CLNVCSO':'SO:0001483',
    'CLNVI':'ClinGen:CA278072|OMIM:600799.0003',
    'DBVARID':'.',
    'GENEINFO':'BMPR2:659',
    'MC':'SO:0001587|nonsense',
    'ONCDN':'.',
    'ONCDNINCL':'.',
    'ONCDISDB':'.',
    'ONCDISDBINCL':'.',
    'ONC':'.',
    'ONCINCL':'.',
    'ONCREVSTAT':'.',
    'ONCCONF':'.',
    'ORIGIN':'1',
    'RS':'137852742',
    'SCIDN':'.',
    'SCIDNINCL':'.',
    'SCIDISDB':'.',
    'SCIDISDBINCL':'.',
    'SCIREVSTAT':'.',
    'SCI':'.',
    'SCIINCL':'.',
    }
    return row
    #{'MONDO:0015924': ['HP:0002092', 'HP:0006546']}

# Multiple records, Multiple mondoID, Multiple HPO, Multiple gene
@pytest.fixture
def test_case6_row():
    row = {
    'CHROM':'1',
    'POS':'216084853',
    'ID':'179773',
    'REF':'C',
    'ALT':'T',
    'QUAL':'.',
    'FILTER':'.',
    'AF_ESP':'.',
    'AF_EXAC':'6e-05',
    'AF_TGP':'.',
    'ALLELEID':'172388',
    'CLNDN':'not_specified|not_provided|Retinitis_pigmentosa|Retinitis_pigmentosa_39|Usher_syndrome|USH2A-related_disorder|Retinal_dystrophy',
    'CLNDNINCL':'Retinitis_pigmentosa',
    'CLNDISDB':'MedGen:CN169374|MedGen:C3661900|Human_Phenotype_Ontology:HP:0000547,MONDO:MONDO:0019200,MeSH:D012174,MedGen:C0035334,OMIM:268000,OMIM:PS268000,Orphanet:791|MONDO:MONDO:0013436,MedGen:C3151138,OMIM:613809,Orphanet:791|MONDO:MONDO:0019501,MeSH:D052245,MedGen:C0271097,OMIM:PS276900,Orphanet:886|MedGen:CN239332|Human_Phenotype_Ontology:HP:0000556,Human_Phenotype_Ontology:HP:0007736,Human_Phenotype_Ontology:HP:0007910,Human_Phenotype_Ontology:HP:0007974,Human_Phenotype_Ontology:HP:0007982,MONDO:MONDO:0019118,MeSH:D058499,MedGen:C0854723,Orphanet:71862',
    'CLNDISDBINCL':'Human_Phenotype_Ontology:HP:0000547,MONDO:MONDO:0019200,MeSH:D012174,MedGen:C0035334,OMIM:268000,OMIM:PS268000,Orphanet:791',
    'CLNHGVS':'NC_000001.11:g.216084853C>T',
    'CLNREVSTAT':'reviewed_by_expert_panel',
    'CLNSIG':'Pathogenic',
    'CLNSIGCONF':'.',
    'CLNSIGINCL':'812445:Likely_pathogenic',
    'CLNVC':'single_nucleotide_variant',
    'CLNVCSO':'SO:0001483',
    'CLNVI':'ClinGen:CA185105',
    'DBVARID':'.',
    'GENEINFO':'USH2A:7399|USH2A-AS2:102723833',
    'MC':'SO:0001583|missense_variant',
    'ONCDN':'.',
    'ONCDNINCL':'.',
    'ONCDISDB':'.',
    'ONCDISDBINCL':'.',
    'ONC':'.',
    'ONCINCL':'.',
    'ONCREVSTAT':'.',
    'ONCCONF':'.',
    'ORIGIN':'1',
    'RS':'727505116',
    'SCIDN':'.',
    'SCIDNINCL':'.',
    'SCIDISDB':'.',
    'SCIDISDBINCL':'.',
    'SCIREVSTAT':'.',
    'SCI':'.',
    'SCIINCL':'.',
    }
    return row
    #{'MONDO:0019118': ['HP:0000556', 'HP:0007736', 'HP:0007910', 'HP:0007974', 'HP:0007982']}

# Multiple records, Single HPO, Single gene, Multiple MC info
@pytest.fixture
def test_case7_row():
    row = {
    'CHROM':'1',
    'POS':'173911974',
    'ID':'654211',
    'REF':'T',
    'ALT':'G',
    'QUAL':'.',
    'FILTER':'.',
    'AF_ESP':'.',
    'AF_EXAC':'1e-05',
    'AF_TGP':'.',
    'ALLELEID':'627126',
    'CLNDN':'not_provided|Hereditary_antithrombin_deficiency',
    'CLNDNINCL':'.',
    'CLNDISDB':'MedGen:C3661900|Human_Phenotype_Ontology:HP:0001976,MONDO:MONDO:0013144,MedGen:C0272375,OMIM:613118,Orphanet:82',
    'CLNDISDBINCL':'.',
    'CLNHGVS':'NC_000001.11:g.173911974T>G',
    'CLNREVSTAT':'reviewed_by_expert_panel',
    'CLNSIG':'Likely_pathogenic',
    'CLNSIGCONF':'.',
    'CLNSIGINCL':'.',
    'CLNVC':'single_nucleotide_variant',
    'CLNVCSO':'SO:0001483',
    'CLNVI':'.',
    'DBVARID':'.',
    'GENEINFO':'SERPINC1:462',
    'MC':'SO:0001583|missense_variant,SO:0001627|intron_variant',
    'ONCDN':'.',
    'ONCDNINCL':'.',
    'ONCDISDB':'.',
    'ONCDISDBINCL':'.',
    'ONC':'.',
    'ONCINCL':'.',
    'ONCREVSTAT':'.',
    'ONCCONF':'.',
    'ORIGIN':'1',
    'RS':'765445413',
    'SCIDN':'.',
    'SCIDNINCL':'.',
    'SCIDISDB':'.',
    'SCIDISDBINCL':'.',
    'SCIREVSTAT':'.',
    'SCI':'.',
    'SCIINCL':'.',
    }
    return row


####################################################################
### Generates koza like output from the input test_case_row data ###
 
@pytest.fixture
def test_case1_entities(test_case1_row, mock_koza):
    return mock_koza(INGEST_NAME,
                     test_case1_row,
                     TRANSFORM_SCRIPT)

@pytest.fixture
def test_case2_entities(test_case2_row, mock_koza):
    return mock_koza(INGEST_NAME,
                     test_case2_row,
                     TRANSFORM_SCRIPT)

@pytest.fixture
def test_case3_entities(test_case3_row, mock_koza):
    return mock_koza(INGEST_NAME,
                     test_case3_row,
                     TRANSFORM_SCRIPT)

@pytest.fixture
def test_case4_entities(test_case4_row, mock_koza):
    return mock_koza(INGEST_NAME,
                     test_case4_row,
                     TRANSFORM_SCRIPT)

@pytest.fixture
def test_case5_entities(test_case5_row, mock_koza):
    return mock_koza(INGEST_NAME,
                     test_case5_row,
                     TRANSFORM_SCRIPT)

@pytest.fixture
def test_case6_entities(test_case6_row, mock_koza):
    return mock_koza(INGEST_NAME,
                     test_case6_row,
                     TRANSFORM_SCRIPT)

@pytest.fixture
def test_case7_entities(test_case7_row, mock_koza):
    return mock_koza(INGEST_NAME,
                     test_case7_row,
                     TRANSFORM_SCRIPT)


########################
### Our actual tests ###

def test_case1(test_case1_entities):
    assert len(test_case1_entities) == 3 # SequenceVariant, VariantToGene. VariantToDisease

def test_case2(test_case2_entities):
    assert len(test_case2_entities) == 4 # SequenceVariant, VariantToGene, VariantToDisease, VariantToPhenotype

def test_case3(test_case3_entities):
    assert len(test_case3_entities) == 5 # SequenceVariant,  VariantToGene, VariantToDisease, VariantToPhenotype, VariantToPhenotype
    assert len([association for association in test_case3_entities if isinstance(association, VariantToPhenotypicFeatureAssociation)]) == 2
    assert test_case3_entities[2].object == "MONDO:0015924" # Two potential mondo ids are listed in the records

def test_case4(test_case4_entities):
    assert len(test_case4_entities) == 3 # SequenceVariant, VariantToGene, VariantToDisease
    assert test_case4_entities[2].object == "MONDO:0100283" # Multiple mondoids are available and this is the one that should be chosen

def test_case5(test_case5_entities):
    assert len(test_case5_entities) == 5 # SequenceVariant, VariantToGene, VariantToDisease, VariantToPhenotype1, VariantToPhenotype2
    assert test_case5_entities[2].object == "MONDO:0015924" # Multiple mondoids are available and this is the one that should be chosen
    assert len([association for association in test_case5_entities if isinstance(association, VariantToPhenotypicFeatureAssociation)]) == 2

def test_case6(test_case6_entities):
    assert len(test_case6_entities) == 9 # SequenceVariant, VariantToGene1, VariantToGene2, VariantToDisease, VariantToPhenotype_x_5 
    assert test_case6_entities[3].object == "MONDO:0019118" # Multiple mondoids are available and this is the one that should be chosen
    assert len([association for association in test_case6_entities if isinstance(association, VariantToGeneAssociation)]) == 2
    assert len([association for association in test_case6_entities if isinstance(association, VariantToPhenotypicFeatureAssociation)]) == 5

def test_case7(test_case7_entities):
    assert len(test_case7_entities) == 4 # SequenceVariant, VariantToGene, VariantToDisease, VariantToPhenotype 
    assert test_case7_entities[2].object == "MONDO:0013144" # Multiple records are available and we want to make sure we choose the proper on
    assert test_case7_entities[0].type == ["SO:0001583", "SO:0001627"]


# TO DO: Tests for proper predicates? Test rows for variants that are not of review status 3 stars or above? 
# (Tricky because paramters / decisions about what information actually gets pulled and when can alter the results obtained for variants with a review status of less than 3 stars