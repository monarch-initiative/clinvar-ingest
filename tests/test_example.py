"""
An example test file for the transform script.

It uses pytest fixtures to define the input data and the mock koza transform.
The test_example function then tests the output of the transform script.

See the Koza documentation for more information on testing transforms:
https://koza.monarchinitiative.org/Usage/testing/
"""

import pytest
from koza.utils.testing_utils import mock_koza
from biolink_model.datamodel.pydanticmodel_v2 import VariantToPhenotypicFeatureAssociation, VariantToGeneAssociation, VariantToDiseaseAssociation


# Define the ingest name and transform script path
INGEST_NAME = "clinvar_variant"
TRANSFORM_SCRIPT = "./src/clinvar_ingest/transform.py"


# Define an example row to test (as a dictionary)
@pytest.fixture
def no_mondo_row():

    # First line in the current clinvar vcf
    row = {'CHROM': '1', 'POS': '69134', 'ID': '2205837', 'REF': 'A', 'ALT': 'G', 'QUAL': '.', 'FILTER': '.', 'INFO': 'ALLELEID=2193183;CLNDISDB=MedGen:CN169374;CLNDN=not_specified;CLNHGVS=NC_000001.11:g.69134A>G;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=OR4F5:79501;MC=SO:0001583|missense_variant;ORIGIN=1'}
    return row

@pytest.fixture
def mondo_and_hpo_row():
    row = {'CHROM': '1', 
           'POS': '2304067', 
           'ID': '992795', 
           'REF': 'C', 
           'ALT': 'T', 
           'QUAL': '.', 
           'FILTER': '.', 
           'INFO': 'AF_ESP=0.00008;AF_EXAC=0.00002;ALLELEID=980658;CLNDISDB=MedGen:C3661900|Human_Phenotype_Ontology:HP:0000730,Human_Phenotype_Ontology:HP:0001249,Human_Phenotype_Ontology:HP:0001267,Human_Phenotype_Ontology:HP:0001286,Human_Phenotype_Ontology:HP:0002122,Human_Phenotype_Ontology:HP:0002192,Human_Phenotype_Ontology:HP:0002316,Human_Phenotype_Ontology:HP:0002382,Human_Phenotype_Ontology:HP:0002386,Human_Phenotype_Ontology:HP:0002402,Human_Phenotype_Ontology:HP:0002458,Human_Phenotype_Ontology:HP:0002482,Human_Phenotype_Ontology:HP:0002499,Human_Phenotype_Ontology:HP:0002543,Human_Phenotype_Ontology:HP:0003767,Human_Phenotype_Ontology:HP:0006833,Human_Phenotype_Ontology:HP:0007154,Human_Phenotype_Ontology:HP:0007176,Human_Phenotype_Ontology:HP:0007180,MONDO:MONDO:0001071,MeSH:D008607,MedGen:C3714756|Human_Phenotype_Ontology:HP:0001380,Human_Phenotype_Ontology:HP:0001383,Human_Phenotype_Ontology:HP:0001388,Human_Phenotype_Ontology:HP:0002771,MedGen:C0086437;CLNDN=not_provided|Intellectual_disability|Joint_laxity;CLNHGVS=NC_000001.11:g.2304067C>T;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=SKI:6497;MC=SO:0001583|missense_variant;ORIGIN=1;RS=367916348'}
    return row


@pytest.fixture
def no_mondo_entities(no_mondo_row):
    return mock_koza(
        INGEST_NAME,
        no_mondo_row,
        TRANSFORM_SCRIPT,
    )

@pytest.fixture
def mondo_hpo_entities(mondo_and_hpo_row):
    return mock_koza(
        INGEST_NAME,
        mondo_and_hpo_row,
        TRANSFORM_SCRIPT,
    )


# Test the output of the transform
def test_no_mondo_row(no_mondo_entities):
    assert len(no_mondo_entities) == 0

def test_mondo_hpo_entities_hpo(mondo_hpo_entities):
    assert len([association for association in mondo_hpo_entities if isinstance(association, VariantToPhenotypicFeatureAssociation)]) == 1 # Or however many there actually are

def test_single_gene(mondo_hpo_entites):
    assert len([association for association in mondo_hpo_entities if isinstance(association, VariantToGeneAssociation)]) == 1
