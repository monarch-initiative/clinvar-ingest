"""
Tests for the ClinVar transform.

Calls process_row() directly with mock var_records and map_to_mondo data
rather than using the removed koza 1.x mock_koza pattern.
"""

import pytest
from biolink_model.datamodel.pydanticmodel_v2 import (
    VariantToGeneAssociation,
    VariantToPhenotypicFeatureAssociation,
)

from clinvar_helpers import process_row


def _make_record(clinsig="Pathogenic", medgen_cui="CN000000", disease_name="test_disease"):
    """Create a minimal submission record that passes the 3-star review filter."""
    return {
        "VariationID": "0",
        "ClinicalSignificance": clinsig,
        "DateLastEvaluated": "2024-01-01",
        "Description": "",
        "SubmittedPhenotypeInfo": ".",
        "ReportedPhenotypeInfo": "{}:{}".format(medgen_cui, disease_name),
        "ReviewStatus": "reviewed by expert panel",
        "CollectionMethod": "clinical testing",
        "OriginCounts": "",
        "Submitter": "TestLab",
        "SCV": "SCV000000001",
    }


# ---- Mock data: var_records keyed by variant ID ----
# Each test case's variant ID must have matching records here.

VAR_RECORDS = {
    # Case 1: ID=1296989 → single record → MONDO:0100283
    "1296989": [_make_record("Pathogenic", "CN300503", "MTOR_overgrowth")],
    # Case 2: ID=2505295 → single record → MONDO:0020367
    "2505295": [_make_record("Likely pathogenic", "C2981140", "Glaucoma_of_childhood")],
    # Case 3 & 5: ID=8797 → only first record passes 3-star filter → MONDO:0015924
    "8797": [
        _make_record("Pathogenic", "C2973725", "Pulmonary_arterial_hypertension"),
    ],
    # Case 4: ID=156702 → only first record passes 3-star filter → MONDO:0100283
    "156702": [
        _make_record("Pathogenic", "CN300503", "MTOR_overgrowth"),
    ],
    # Case 6: ID=179773 → records for MONDO:0019118 (retinal dystrophy)
    "179773": [
        _make_record("Pathogenic", "C0854723", "Retinal_dystrophy"),
    ],
    # Case 7: ID=654211 → records for MONDO:0013144
    "654211": [
        _make_record("Likely pathogenic", "C0272375", "Hereditary_antithrombin_deficiency"),
    ],
}

# ---- Mock data: map_to_mondo ----
# Maps MedGen CUIs, MONDO self-references, and other IDs to MONDO IDs.

MAP_TO_MONDO = {
    # MedGen → MONDO mappings (used by variant_records_to_disease via ReportedPhenotypeInfo)
    "MedGen:CN300503": {"MONDO:0100283": ""},
    "MedGen:C2981140": {"MONDO:0020367": ""},
    "MedGen:C2973725": {"MONDO:0015924": ""},
    "MedGen:C4552070": {"MONDO:0024533": ""},
    "MedGen:C4225259": {"MONDO:0014716": ""},
    "MedGen:C0854723": {"MONDO:0019118": ""},
    "MedGen:C0272375": {"MONDO:0013144": ""},
    "MedGen:C3661900": {"MONDO:0000001": ""},  # "not_provided" — doesn't match CLNDISDB
    "MedGen:C0035334": {"MONDO:0019200": ""},
    "MedGen:C3151138": {"MONDO:0013436": ""},
    "MedGen:C0271097": {"MONDO:0019501": ""},
    "MedGen:CN239332": {"MONDO:0099999": ""},
    "MedGen:CN169374": {"MONDO:0099998": ""},
    # MONDO self-mappings (used by parse_CLNDISDB → map_CLNDISDB_to_mondo)
    "MONDO:MONDO:0100283": {"MONDO:0100283": ""},
    "MONDO:0100283": {"MONDO:0100283": ""},
    "MONDO:MONDO:0020367": {"MONDO:0020367": ""},
    "MONDO:0020367": {"MONDO:0020367": ""},
    "MONDO:MONDO:0015924": {"MONDO:0015924": ""},
    "MONDO:0015924": {"MONDO:0015924": ""},
    "MONDO:MONDO:0024533": {"MONDO:0024533": ""},
    "MONDO:0024533": {"MONDO:0024533": ""},
    "MONDO:MONDO:0014716": {"MONDO:0014716": ""},
    "MONDO:0014716": {"MONDO:0014716": ""},
    "MONDO:MONDO:0019200": {"MONDO:0019200": ""},
    "MONDO:0019200": {"MONDO:0019200": ""},
    "MONDO:MONDO:0013436": {"MONDO:0013436": ""},
    "MONDO:0013436": {"MONDO:0013436": ""},
    "MONDO:MONDO:0019501": {"MONDO:0019501": ""},
    "MONDO:0019501": {"MONDO:0019501": ""},
    "MONDO:MONDO:0019118": {"MONDO:0019118": ""},
    "MONDO:0019118": {"MONDO:0019118": ""},
    "MONDO:MONDO:0013144": {"MONDO:0013144": ""},
    "MONDO:0013144": {"MONDO:0013144": ""},
    # OMIM → MONDO mappings (used by CLNDISDB mapping)
    "OMIM:178600": {"MONDO:0024533": ""},
    "OMIM:268000": {"MONDO:0019200": ""},
    "OMIM:613809": {"MONDO:0013436": ""},
    "OMIM:616638": {"MONDO:0014716": ""},
    "OMIM:613118": {"MONDO:0013144": ""},
    # Orphanet → MONDO mappings
    "Orphanet:182090": {"MONDO:0015924": ""},
    "Orphanet:422": {"MONDO:0024533": ""},
    "Orphanet:98977": {"MONDO:0020367": ""},
    "Orphanet:791": {"MONDO:0019200": ""},
    "Orphanet:886": {"MONDO:0019501": ""},
    "Orphanet:457485": {"MONDO:0014716": ""},
    "Orphanet:82": {"MONDO:0013144": ""},
    "Orphanet:71862": {"MONDO:0019118": ""},
    # MeSH → MONDO mappings
    "mesh:D000081029": {"MONDO:0015924": ""},
    "mesh:D012174": {"MONDO:0019200": ""},
    "mesh:D052245": {"MONDO:0019501": ""},
    "mesh:D058499": {"MONDO:0019118": ""},
}


#################################################################
### Rows taken from vcf file that have our desired test cases ###


# Single record, Single mondoID, No HPO, Single Gene
@pytest.fixture
def test_case1_row():
    return {
        "CHROM": "1",
        "POS": "11128107",
        "ID": "1296989",
        "REF": "G",
        "ALT": "C",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": ".",
        "AF_TGP": ".",
        "ALLELEID": "1286779",
        "CLNDN": "Overgrowth_syndrome_and/or_cerebral_malformations_due_to_abnormalities_in_MTOR_pathway_genes",
        "CLNDNINCL": ".",
        "CLNDISDB": "MONDO:MONDO:0100283,MedGen:CN300503",
        "CLNDISDBINCL": ".",
        "CLNHGVS": "NC_000001.11:g.11128107G>C",
        "CLNREVSTAT": "reviewed_by_expert_panel",
        "CLNSIG": "Pathogenic",
        "CLNSIGCONF": ".",
        "CLNSIGINCL": ".",
        "CLNVC": "single_nucleotide_variant",
        "CLNVCSO": "SO:0001483",
        "CLNVI": ".",
        "DBVARID": ".",
        "GENEINFO": "MTOR:2475",
        "MC": "SO:0001583|missense_variant",
        "ONCDN": ".",
        "ONCDNINCL": ".",
        "ONCDISDB": ".",
        "ONCDISDBINCL": ".",
        "ONC": ".",
        "ONCINCL": ".",
        "ONCREVSTAT": ".",
        "ONCCONF": ".",
        "ORIGIN": "1",
        "RS": "587777893",
        "SCIDN": ".",
        "SCIDNINCL": ".",
        "SCIDISDB": ".",
        "SCIDISDBINCL": ".",
        "SCIREVSTAT": ".",
        "SCI": ".",
        "SCIINCL": ".",
    }


# Single record, Single mondoID, 1 HPO term, Single gene
@pytest.fixture
def test_case2_row():
    return {
        "CHROM": "1",
        "POS": "171636092",
        "ID": "2505295",
        "REF": "T",
        "ALT": "A",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": ".",
        "AF_TGP": ".",
        "ALLELEID": "2670262",
        "CLNDN": "Glaucoma_of_childhood",
        "CLNDNINCL": ".",
        "CLNDISDB": "Human_Phenotype_Ontology:HP:0001087,MONDO:MONDO:0020367,MedGen:C2981140,Orphanet:98977",
        "CLNDISDBINCL": ".",
        "CLNHGVS": "NC_000001.11:g.171636092T>A",
        "CLNREVSTAT": "reviewed_by_expert_panel",
        "CLNSIG": "Likely_pathogenic",
        "CLNSIGCONF": ".",
        "CLNSIGINCL": ".",
        "CLNVC": "single_nucleotide_variant",
        "CLNVCSO": "SO:0001483",
        "CLNVI": ".",
        "DBVARID": ".",
        "GENEINFO": "MYOC:4653",
        "MC": "SO:0001583|missense_variant",
        "ONCDN": ".",
        "ONCDNINCL": ".",
        "ONCDISDB": ".",
        "ONCDISDBINCL": ".",
        "ONC": ".",
        "ONCINCL": ".",
        "ONCREVSTAT": ".",
        "ONCCONF": ".",
        "ORIGIN": "1",
        "RS": ".",
        "SCIDN": ".",
        "SCIDNINCL": ".",
        "SCIDISDB": ".",
        "SCIDISDBINCL": ".",
        "SCIREVSTAT": ".",
        "SCI": ".",
        "SCIINCL": ".",
    }


# Two records, Two mondoIDs, 2 HPO terms, Single gene
@pytest.fixture
def test_case3_row():
    return {
        "CHROM": "2",
        "POS": "202464950",
        "ID": "8797",
        "REF": "C",
        "ALT": "G",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": ".",
        "AF_TGP": ".",
        "ALLELEID": "23836",
        "CLNDN": "Pulmonary_arterial_hypertension|Pulmonary_hypertension,_primary,_1",
        "CLNDNINCL": ".",
        "CLNDISDB": "Human_Phenotype_Ontology:HP:0002092,Human_Phenotype_Ontology:HP:0006546,MONDO:MONDO:0015924,MeSH:D000081029,MedGen:C2973725,Orphanet:182090|MONDO:MONDO:0024533,MedGen:C4552070,OMIM:178600,Orphanet:422",
        "CLNDISDBINCL": ".",
        "CLNHGVS": "NC_000002.12:g.202464950C>G",
        "CLNREVSTAT": "reviewed_by_expert_panel",
        "CLNSIG": "Pathogenic",
        "CLNSIGCONF": ".",
        "CLNSIGINCL": ".",
        "CLNVC": "single_nucleotide_variant",
        "CLNVCSO": "SO:0001483",
        "CLNVI": "ClinGen:CA278072|OMIM:600799.0003",
        "DBVARID": ".",
        "GENEINFO": "BMPR2:659",
        "MC": "SO:0001587|nonsense",
        "ONCDN": ".",
        "ONCDNINCL": ".",
        "ONCDISDB": ".",
        "ONCDISDBINCL": ".",
        "ONC": ".",
        "ONCINCL": ".",
        "ONCREVSTAT": ".",
        "ONCCONF": ".",
        "ORIGIN": "1",
        "RS": "137852742",
        "SCIDN": ".",
        "SCIDNINCL": ".",
        "SCIDISDB": ".",
        "SCIDISDBINCL": ".",
        "SCIREVSTAT": ".",
        "SCI": ".",
        "SCIINCL": ".",
    }


# Multiple records, Multiple mondoID, No HPO, Single gene
@pytest.fixture
def test_case4_row():
    return {
        "CHROM": "1",
        "POS": "11128107",
        "ID": "156702",
        "REF": "G",
        "ALT": "T",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": ".",
        "AF_TGP": ".",
        "ALLELEID": "166562",
        "CLNDN": "not_provided|Overgrowth_syndrome_and/or_cerebral_malformations_due_to_abnormalities_in_MTOR_pathway_genes|Macrocephaly-intellectual_disability-neurodevelopmental_disorder-small_thorax_syndrome",
        "CLNDNINCL": ".",
        "CLNDISDB": "MedGen:C3661900|MONDO:MONDO:0100283,MedGen:CN300503|MONDO:MONDO:0014716,MedGen:C4225259,OMIM:616638,Orphanet:457485",
        "CLNDISDBINCL": ".",
        "CLNHGVS": "NC_000001.11:g.11128107G>T",
        "CLNREVSTAT": "reviewed_by_expert_panel",
        "CLNSIG": "Pathogenic",
        "CLNSIGCONF": ".",
        "CLNSIGINCL": ".",
        "CLNVC": "single_nucleotide_variant",
        "CLNVCSO": "SO:0001483",
        "CLNVI": "ClinGen:CA248390",
        "DBVARID": ".",
        "GENEINFO": "MTOR:2475",
        "MC": "SO:0001583|missense_variant",
        "ONCDN": ".",
        "ONCDNINCL": ".",
        "ONCDISDB": ".",
        "ONCDISDBINCL": ".",
        "ONC": ".",
        "ONCINCL": ".",
        "ONCREVSTAT": ".",
        "ONCCONF": ".",
        "ORIGIN": "1",
        "RS": "587777893",
        "SCIDN": ".",
        "SCIDNINCL": ".",
        "SCIDISDB": ".",
        "SCIDISDBINCL": ".",
        "SCIREVSTAT": ".",
        "SCI": ".",
        "SCIINCL": ".",
    }


# Multiple records, Multiple mondoID, Multiple HPO, Single gene
@pytest.fixture
def test_case5_row():
    return {
        "CHROM": "2",
        "POS": "202464950",
        "ID": "8797",
        "REF": "C",
        "ALT": "G",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": ".",
        "AF_TGP": ".",
        "ALLELEID": "23836",
        "CLNDN": "Pulmonary_arterial_hypertension|Pulmonary_hypertension,_primary,_1",
        "CLNDNINCL": ".",
        "CLNDISDB": "Human_Phenotype_Ontology:HP:0002092,Human_Phenotype_Ontology:HP:0006546,MONDO:MONDO:0015924,MeSH:D000081029,MedGen:C2973725,Orphanet:182090|MONDO:MONDO:0024533,MedGen:C4552070,OMIM:178600,Orphanet:422",
        "CLNDISDBINCL": ".",
        "CLNHGVS": "NC_000002.12:g.202464950C>G",
        "CLNREVSTAT": "reviewed_by_expert_panel",
        "CLNSIG": "Pathogenic",
        "CLNSIGCONF": ".",
        "CLNSIGINCL": ".",
        "CLNVC": "single_nucleotide_variant",
        "CLNVCSO": "SO:0001483",
        "CLNVI": "ClinGen:CA278072|OMIM:600799.0003",
        "DBVARID": ".",
        "GENEINFO": "BMPR2:659",
        "MC": "SO:0001587|nonsense",
        "ONCDN": ".",
        "ONCDNINCL": ".",
        "ONCDISDB": ".",
        "ONCDISDBINCL": ".",
        "ONC": ".",
        "ONCINCL": ".",
        "ONCREVSTAT": ".",
        "ONCCONF": ".",
        "ORIGIN": "1",
        "RS": "137852742",
        "SCIDN": ".",
        "SCIDNINCL": ".",
        "SCIDISDB": ".",
        "SCIDISDBINCL": ".",
        "SCIREVSTAT": ".",
        "SCI": ".",
        "SCIINCL": ".",
    }


# Multiple records, Multiple mondoID, Multiple HPO, Multiple gene
@pytest.fixture
def test_case6_row():
    return {
        "CHROM": "1",
        "POS": "216084853",
        "ID": "179773",
        "REF": "C",
        "ALT": "T",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": "6e-05",
        "AF_TGP": ".",
        "ALLELEID": "172388",
        "CLNDN": "not_specified|not_provided|Retinitis_pigmentosa|Retinitis_pigmentosa_39|Usher_syndrome|USH2A-related_disorder|Retinal_dystrophy",
        "CLNDNINCL": "Retinitis_pigmentosa",
        "CLNDISDB": "MedGen:CN169374|MedGen:C3661900|Human_Phenotype_Ontology:HP:0000547,MONDO:MONDO:0019200,MeSH:D012174,MedGen:C0035334,OMIM:268000,OMIM:PS268000,Orphanet:791|MONDO:MONDO:0013436,MedGen:C3151138,OMIM:613809,Orphanet:791|MONDO:MONDO:0019501,MeSH:D052245,MedGen:C0271097,OMIM:PS276900,Orphanet:886|MedGen:CN239332|Human_Phenotype_Ontology:HP:0000556,Human_Phenotype_Ontology:HP:0007736,Human_Phenotype_Ontology:HP:0007910,Human_Phenotype_Ontology:HP:0007974,Human_Phenotype_Ontology:HP:0007982,MONDO:MONDO:0019118,MeSH:D058499,MedGen:C0854723,Orphanet:71862",
        "CLNDISDBINCL": "Human_Phenotype_Ontology:HP:0000547,MONDO:MONDO:0019200,MeSH:D012174,MedGen:C0035334,OMIM:268000,OMIM:PS268000,Orphanet:791",
        "CLNHGVS": "NC_000001.11:g.216084853C>T",
        "CLNREVSTAT": "reviewed_by_expert_panel",
        "CLNSIG": "Pathogenic",
        "CLNSIGCONF": ".",
        "CLNSIGINCL": "812445:Likely_pathogenic",
        "CLNVC": "single_nucleotide_variant",
        "CLNVCSO": "SO:0001483",
        "CLNVI": "ClinGen:CA185105",
        "DBVARID": ".",
        "GENEINFO": "USH2A:7399|USH2A-AS2:102723833",
        "MC": "SO:0001583|missense_variant",
        "ONCDN": ".",
        "ONCDNINCL": ".",
        "ONCDISDB": ".",
        "ONCDISDBINCL": ".",
        "ONC": ".",
        "ONCINCL": ".",
        "ONCREVSTAT": ".",
        "ONCCONF": ".",
        "ORIGIN": "1",
        "RS": "727505116",
        "SCIDN": ".",
        "SCIDNINCL": ".",
        "SCIDISDB": ".",
        "SCIDISDBINCL": ".",
        "SCIREVSTAT": ".",
        "SCI": ".",
        "SCIINCL": ".",
    }


# Multiple records, Single HPO, Single gene, Multiple MC info
@pytest.fixture
def test_case7_row():
    return {
        "CHROM": "1",
        "POS": "173911974",
        "ID": "654211",
        "REF": "T",
        "ALT": "G",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": "1e-05",
        "AF_TGP": ".",
        "ALLELEID": "627126",
        "CLNDN": "not_provided|Hereditary_antithrombin_deficiency",
        "CLNDNINCL": ".",
        "CLNDISDB": "MedGen:C3661900|Human_Phenotype_Ontology:HP:0001976,MONDO:MONDO:0013144,MedGen:C0272375,OMIM:613118,Orphanet:82",
        "CLNDISDBINCL": ".",
        "CLNHGVS": "NC_000001.11:g.173911974T>G",
        "CLNREVSTAT": "reviewed_by_expert_panel",
        "CLNSIG": "Likely_pathogenic",
        "CLNSIGCONF": ".",
        "CLNSIGINCL": ".",
        "CLNVC": "single_nucleotide_variant",
        "CLNVCSO": "SO:0001483",
        "CLNVI": ".",
        "DBVARID": ".",
        "GENEINFO": "SERPINC1:462",
        "MC": "SO:0001583|missense_variant,SO:0001627|intron_variant",
        "ONCDN": ".",
        "ONCDNINCL": ".",
        "ONCDISDB": ".",
        "ONCDISDBINCL": ".",
        "ONC": ".",
        "ONCINCL": ".",
        "ONCREVSTAT": ".",
        "ONCCONF": ".",
        "ORIGIN": "1",
        "RS": "765445413",
        "SCIDN": ".",
        "SCIDNINCL": ".",
        "SCIDISDB": ".",
        "SCIDISDBINCL": ".",
        "SCIREVSTAT": ".",
        "SCI": ".",
        "SCIINCL": ".",
    }


####################################################################
### Entity fixtures using process_row directly ###


@pytest.fixture
def test_case1_entities(test_case1_row):
    return process_row(test_case1_row, VAR_RECORDS, MAP_TO_MONDO)


@pytest.fixture
def test_case2_entities(test_case2_row):
    return process_row(test_case2_row, VAR_RECORDS, MAP_TO_MONDO)


@pytest.fixture
def test_case3_entities(test_case3_row):
    return process_row(test_case3_row, VAR_RECORDS, MAP_TO_MONDO)


@pytest.fixture
def test_case4_entities(test_case4_row):
    return process_row(test_case4_row, VAR_RECORDS, MAP_TO_MONDO)


@pytest.fixture
def test_case5_entities(test_case5_row):
    return process_row(test_case5_row, VAR_RECORDS, MAP_TO_MONDO)


@pytest.fixture
def test_case6_entities(test_case6_row):
    return process_row(test_case6_row, VAR_RECORDS, MAP_TO_MONDO)


@pytest.fixture
def test_case7_entities(test_case7_row):
    return process_row(test_case7_row, VAR_RECORDS, MAP_TO_MONDO)


########################
### Our actual tests ###


def test_case1(test_case1_entities):
    assert len(test_case1_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease


def test_case2(test_case2_entities):
    assert len(test_case2_entities) == 4  # SequenceVariant, VariantToGene, VariantToDisease, VariantToPhenotype


def test_case3(test_case3_entities):
    assert len(test_case3_entities) == 5  # SequenceVariant, VariantToGene, VariantToDisease, VariantToPhenotype x2
    assert len([a for a in test_case3_entities if isinstance(a, VariantToPhenotypicFeatureAssociation)]) == 2
    assert test_case3_entities[2].object == "MONDO:0015924"


def test_case4(test_case4_entities):
    assert len(test_case4_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease
    assert test_case4_entities[2].object == "MONDO:0100283"


def test_case5(test_case5_entities):
    assert len(test_case5_entities) == 5  # SequenceVariant, VariantToGene, VariantToDisease, VariantToPhenotype x2
    assert test_case5_entities[2].object == "MONDO:0015924"
    assert len([a for a in test_case5_entities if isinstance(a, VariantToPhenotypicFeatureAssociation)]) == 2


def test_case6(test_case6_entities):
    assert len(test_case6_entities) == 9  # SequenceVariant, VariantToGene x2, VariantToDisease, VariantToPhenotype x5
    assert test_case6_entities[3].object == "MONDO:0019118"
    assert len([a for a in test_case6_entities if isinstance(a, VariantToGeneAssociation)]) == 2
    assert len([a for a in test_case6_entities if isinstance(a, VariantToPhenotypicFeatureAssociation)]) == 5


def test_case7(test_case7_entities):
    assert len(test_case7_entities) == 4  # SequenceVariant, VariantToGene, VariantToDisease, VariantToPhenotype
    assert test_case7_entities[2].object == "MONDO:0013144"
    assert test_case7_entities[0].type == ["SO:0001583", "SO:0001627"]
