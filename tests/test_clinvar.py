"""
Tests for the ClinVar transform.

Calls process_row() directly with mock var_records and map_to_mondo data
rather than using the removed koza 1.x mock_koza pattern.
"""

import pytest
from biolink_model.datamodel.pydanticmodel_v2 import (
    VariantToDiseaseAssociation,
    VariantToGeneAssociation,
)
from koza.runner import KozaTransform, PassthroughWriter

from clinvar_helpers import process_row


def _make_record(
    clinsig="Pathogenic",
    medgen_cui="CN000000",
    disease_name="test_disease",
    review_status="reviewed by expert panel",
    submitter="TestLab",
):
    """Create a minimal submission record. Defaults pass the 3-star review filter."""
    return {
        "VariationID": "0",
        "ClinicalSignificance": clinsig,
        "DateLastEvaluated": "2024-01-01",
        "Description": "",
        "SubmittedPhenotypeInfo": ".",
        "ReportedPhenotypeInfo": "{}:{}".format(medgen_cui, disease_name),
        "ReviewStatus": review_status,
        "CollectionMethod": "clinical testing",
        "OriginCounts": "",
        "Submitter": submitter,
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
    # Case 8: ID=9000001 → two 1-star submitters agreeing on Pathogenic + same disease
    # → rescued via concordance despite neither reaching star_min=3
    "9000001": [
        _make_record(
            "Pathogenic",
            "C2981140",
            "Glaucoma_of_childhood",
            review_status="criteria provided, single submitter",
            submitter="LabA",
        ),
        _make_record(
            "Pathogenic",
            "C2981140",
            "Glaucoma_of_childhood",
            review_status="criteria provided, single submitter",
            submitter="LabB",
        ),
    ],
    # Case 9: ID=9000002 → a single 1-star submitter alone → NOT rescued
    "9000002": [
        _make_record(
            "Pathogenic",
            "C2981140",
            "Glaucoma_of_childhood",
            review_status="criteria provided, single submitter",
            submitter="LabA",
        ),
    ],
    # Case 10: ID=9000003 → two 1-star submitters, same disease, but DIFFERENT
    # classifications → NOT rescued (concordance requires the exact same clinsig)
    "9000003": [
        _make_record(
            "Pathogenic",
            "C2981140",
            "Glaucoma_of_childhood",
            review_status="criteria provided, single submitter",
            submitter="LabA",
        ),
        _make_record(
            "Likely pathogenic",
            "C2981140",
            "Glaucoma_of_childhood",
            review_status="criteria provided, single submitter",
            submitter="LabB",
        ),
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


# ---- Mock data: variant_genes (ClinVar's curated per-variant gene attribution) ----
# Mirrors make_variant_gene_map()'s output: VariationID -> (gene curie, symbol), one
# gene per variant. Contrast with each row's GENEINFO field, which lists every
# overlapping locus -- case 6's GENEINFO carries USH2A-AS2 alongside USH2A, and the
# antisense transcript must NOT produce a gene association.
VARIANT_GENES = {
    "1296989": ("HGNC:3942", "MTOR"),
    "2505295": ("HGNC:7610", "MYOC"),
    "8797": ("HGNC:1078", "BMPR2"),
    "156702": ("HGNC:3942", "MTOR"),
    "179773": ("HGNC:12601", "USH2A"),
    "654211": ("HGNC:775", "SERPINC1"),
    "9000001": ("HGNC:7610", "MYOC"),
    "9000002": ("HGNC:7610", "MYOC"),
    "9000003": ("HGNC:7610", "MYOC"),
    # 9000004 deliberately absent: ClinVar GeneID == -1, i.e. no asserted gene
}


# Minimal VCF row shared by tests that construct their own cases. Mirrors the MYOC /
# MONDO:0020367 wiring used by the rescue fixtures.
_TEST_ROW_TEMPLATE = {
    "CHROM": "1", "POS": "171636092", "ID": "0", "REF": "T", "ALT": "A",
    "CLNDISDB": "Human_Phenotype_Ontology:HP:0001087,MONDO:MONDO:0020367,MedGen:C2981140,Orphanet:98977",
    "CLNHGVS": "NC_000001.11:g.171636092T>A",
    "CLNREVSTAT": "criteria_provided,_single_submitter",
    "CLNVC": "single_nucleotide_variant", "CLNVCSO": "SO:0001483",
    "GENEINFO": "MYOC:4653", "MC": "SO:0001583|missense_variant", "RS": ".",
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


# Rescue cases: reuse the MYOC/MONDO:0020367/HP:0001087 wiring from test_case2 to
# isolate the multi-submitter-concordance behavior without new mondo/HPO fixture data.

# Two 1-star submitters, same disease, same classification → rescued
@pytest.fixture
def test_case8_row():
    return {
        "CHROM": "1",
        "POS": "171636092",
        "ID": "9000001",
        "REF": "T",
        "ALT": "A",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": ".",
        "AF_TGP": ".",
        "ALLELEID": "2670263",
        "CLNDN": "Glaucoma_of_childhood",
        "CLNDNINCL": ".",
        "CLNDISDB": "Human_Phenotype_Ontology:HP:0001087,MONDO:MONDO:0020367,MedGen:C2981140,Orphanet:98977",
        "CLNDISDBINCL": ".",
        "CLNHGVS": "NC_000001.11:g.171636092T>A",
        "CLNREVSTAT": "criteria_provided,_multiple_submitters,_no_conflicts",
        "CLNSIG": "Pathogenic",
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


# Single 1-star submitter alone → NOT rescued
@pytest.fixture
def test_case9_row():
    return {
        "CHROM": "1",
        "POS": "171636092",
        "ID": "9000002",
        "REF": "T",
        "ALT": "A",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": ".",
        "AF_TGP": ".",
        "ALLELEID": "2670264",
        "CLNDN": "Glaucoma_of_childhood",
        "CLNDNINCL": ".",
        "CLNDISDB": "Human_Phenotype_Ontology:HP:0001087,MONDO:MONDO:0020367,MedGen:C2981140,Orphanet:98977",
        "CLNDISDBINCL": ".",
        "CLNHGVS": "NC_000001.11:g.171636092T>A",
        "CLNREVSTAT": "criteria_provided,_single_submitter",
        "CLNSIG": "Pathogenic",
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


# Two 1-star submitters, same disease, DIFFERENT classifications → NOT rescued
@pytest.fixture
def test_case10_row():
    return {
        "CHROM": "1",
        "POS": "171636092",
        "ID": "9000003",
        "REF": "T",
        "ALT": "A",
        "QUAL": ".",
        "FILTER": ".",
        "AF_ESP": ".",
        "AF_EXAC": ".",
        "AF_TGP": ".",
        "ALLELEID": "2670265",
        "CLNDN": "Glaucoma_of_childhood",
        "CLNDNINCL": ".",
        "CLNDISDB": "Human_Phenotype_Ontology:HP:0001087,MONDO:MONDO:0020367,MedGen:C2981140,Orphanet:98977",
        "CLNDISDBINCL": ".",
        "CLNHGVS": "NC_000001.11:g.171636092T>A",
        "CLNREVSTAT": "criteria_provided,_conflicting_classifications",
        "CLNSIG": "Conflicting_classifications_of_pathogenicity",
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


####################################################################
### Entity fixtures using process_row directly ###


@pytest.fixture
def test_case1_entities(test_case1_row):
    return process_row(test_case1_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


@pytest.fixture
def test_case2_entities(test_case2_row):
    return process_row(test_case2_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


@pytest.fixture
def test_case3_entities(test_case3_row):
    return process_row(test_case3_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


@pytest.fixture
def test_case4_entities(test_case4_row):
    return process_row(test_case4_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


@pytest.fixture
def test_case5_entities(test_case5_row):
    return process_row(test_case5_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


@pytest.fixture
def test_case6_entities(test_case6_row):
    return process_row(test_case6_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


@pytest.fixture
def test_case7_entities(test_case7_row):
    return process_row(test_case7_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


@pytest.fixture
def test_case8_entities(test_case8_row):
    return process_row(test_case8_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


@pytest.fixture
def test_case9_entities(test_case9_row):
    return process_row(test_case9_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


@pytest.fixture
def test_case10_entities(test_case10_row):
    return process_row(test_case10_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)


########################
### Our actual tests ###


def test_case1(test_case1_entities):
    assert len(test_case1_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease


def test_case2(test_case2_entities):
    assert len(test_case2_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease


def test_case3(test_case3_entities):
    assert len(test_case3_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease
    assert test_case3_entities[2].object == "MONDO:0015924"


def test_case4(test_case4_entities):
    assert len(test_case4_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease
    assert test_case4_entities[2].object == "MONDO:0100283"


def test_case5(test_case5_entities):
    assert len(test_case5_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease
    assert test_case5_entities[2].object == "MONDO:0015924"


def test_case6_only_the_asserted_gene(test_case6_entities):
    """GENEINFO for this row is "USH2A:7399|USH2A-AS2:102723833", but ClinVar attributes
    the variant to USH2A alone. The antisense transcript must not gain a gene association
    (and so must not inherit USH2A's diseases)."""
    assert len(test_case6_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease
    gene_assocs = [a for a in test_case6_entities if isinstance(a, VariantToGeneAssociation)]
    assert len(gene_assocs) == 1
    assert gene_assocs[0].object == "HGNC:12601"
    assert test_case6_entities[2].object == "MONDO:0019118"


def test_case7_type_is_variant_class_not_consequence(test_case7_entities):
    """SequenceVariant.type carries the variant CLASS from CLNVCSO, not the molecular
    consequence from MC. This row's MC lists SO:0001583 (missense) and SO:0001627
    (intron) -- neither belongs in `type`; the variant is a single nucleotide variant,
    SO:0001483. Consequence is a property of variant x transcript, not of the variant."""
    assert len(test_case7_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease
    assert test_case7_entities[2].object == "MONDO:0013144"
    assert test_case7_entities[0].type == ["SO:0001483"]


def test_case8_rescued_by_concordant_submitters(test_case8_entities):
    """Two 1-star submitters independently agreeing on the same disease and the
    same classification should be rescued despite neither reaching star_min=3."""
    assert len(test_case8_entities) == 3  # SequenceVariant, VariantToGene, VariantToDisease
    diseases = [e for e in test_case8_entities if isinstance(e, VariantToDiseaseAssociation)]
    assert len(diseases) == 1
    assert diseases[0].object == "MONDO:0020367"


def test_case9_not_rescued_single_low_star_submitter(test_case9_entities):
    """A single 1-star submitter alone (no concordant corroboration) must NOT be rescued."""
    assert test_case9_entities == []


def test_no_dbsnp_xref_when_rs_absent(test_case2_row, test_case1_row):
    """RS is "." on variants ClinVar has not mapped to dbSNP; emitting "DBSNP:." would be
    a junk CURIE, so xref is left empty instead."""
    no_rs = process_row(test_case2_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)
    assert test_case2_row["RS"] == "."
    assert no_rs[0].xref == []
    with_rs = process_row(test_case1_row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES)
    assert with_rs[0].xref == ["DBSNP:587777893"]


def test_unasserted_gene_yields_no_gene_association(test_case1_row):
    """When ClinVar declines to attribute a variant to a gene (GeneID == -1, so the
    variant is absent from variant_genes), no VariantToGeneAssociation is emitted and
    there is deliberately NO fallback to GENEINFO -- but the variant's disease edges
    still stand."""
    entities = process_row(test_case1_row, VAR_RECORDS, MAP_TO_MONDO, {})
    assert [a for a in entities if isinstance(a, VariantToGeneAssociation)] == []
    assert entities[0].has_gene == []
    diseases = [a for a in entities if isinstance(a, VariantToDiseaseAssociation)]
    assert len(diseases) == 1
    assert diseases[0].object == "MONDO:0100283"


def test_case10_not_rescued_differing_classifications(test_case10_entities):
    """Two 1-star submitters on the same disease but different classifications must
    NOT be rescued -- concordance requires the exact same ClinicalSignificance."""
    assert test_case10_entities == []


def test_no_phenotype_edges_are_emitted(test_case3_entities, test_case6_entities):
    """ClinVar's HPO ids are cross-references naming the same condition as the disease id
    in the same CLNDISDB group -- not observed phenotypes. Emitting
    VariantToPhenotypicFeatureAssociation from them restated the disease edge in a second
    vocabulary, so no phenotype edges are produced. Case 6's row carries five such HPO
    terms and must still yield only variant, gene and disease entities."""
    for entities in (test_case3_entities, test_case6_entities):
        assert [type(e).__name__ for e in entities] == [
            "SequenceVariant",
            "VariantToGeneAssociation",
            "VariantToDiseaseAssociation",
        ]


def test_one_predicate_per_variant_disease():
    """A variant often carries both Pathogenic and Likely-pathogenic records for the same
    disease. Emitting an edge per predicate asserted both `causes` and
    `associated_with_increased_likelihood_of` for the same (variant, disease), which is
    contradictory. Only the strongest is emitted; the submitted classifications remain
    visible in original_predicate."""
    row = {**_TEST_ROW_TEMPLATE, "ID": "9000005"}
    records = [
        _make_record("Pathogenic", "C2981140", "Glaucoma_of_childhood"),
        _make_record("Likely pathogenic", "C2981140", "Glaucoma_of_childhood"),
    ]
    entities = process_row(
        row, {"9000005": records}, MAP_TO_MONDO, {"9000005": ("HGNC:7610", "MYOC")}
    )
    disease_edges = [e for e in entities if isinstance(e, VariantToDiseaseAssociation)]
    assert len(disease_edges) == 1
    assert disease_edges[0].predicate == "biolink:causes"
    assert disease_edges[0].original_predicate == "Likely pathogenic:Pathogenic"


def test_likely_pathogenic_also_causes():
    """Likely pathogenic expresses the curator's confidence in a causal claim, not a
    weaker kind of relationship, so it maps to `causes` alongside Pathogenic and
    Pathogenic/Likely pathogenic. associated_with_increased_likelihood_of is not emitted."""
    row = {**_TEST_ROW_TEMPLATE, "ID": "9000006"}
    records = [_make_record("Likely pathogenic", "C2981140", "Glaucoma_of_childhood")]
    entities = process_row(
        row, {"9000006": records}, MAP_TO_MONDO, {"9000006": ("HGNC:7610", "MYOC")}
    )
    disease_edges = [e for e in entities if isinstance(e, VariantToDiseaseAssociation)]
    assert len(disease_edges) == 1
    assert disease_edges[0].predicate == "biolink:causes"
    assert disease_edges[0].original_predicate == "Likely pathogenic"


def test_non_snv_indel_classes_are_pruned():
    """Only SNVs and indels are ingested. Microsatellites, inversions and the catch-all
    "Variation" class are dropped up front -- a repeat expansion or inversion is not well
    represented by a fixed REF/ALT."""
    for clnvc in ("Microsatellite", "Inversion", "Variation"):
        row = {**_TEST_ROW_TEMPLATE, "ID": "1296989", "CLNVC": clnvc}
        assert process_row(row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES) == [], clnvc
    kept = {**_TEST_ROW_TEMPLATE, "ID": "1296989", "CLNVC": "single_nucleotide_variant",
            "CLNDISDB": "MONDO:MONDO:0100283,MedGen:CN300503"}
    assert process_row(kept, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES) != []


def test_low_star_from_literature_is_associated_with():
    """A <=1-star call has no multi-submitter corroboration, but a classification the lab
    recorded as coming from the literature enters under the weaker predicate. Without that,
    it is dropped -- "the variant appears in a paper" is not the same claim."""
    row = {**_TEST_ROW_TEMPLATE, "ID": "9000007",
           "CLNREVSTAT": "criteria_provided,_single_submitter"}
    records = [_make_record("Pathogenic", "C2981140", "Glaucoma_of_childhood",
                            review_status="criteria provided, single submitter")]
    genes = {"9000007": ("HGNC:7610", "MYOC")}

    without = process_row(row, {"9000007": records}, MAP_TO_MONDO, genes, frozenset())
    assert [e for e in without if isinstance(e, VariantToDiseaseAssociation)] == []

    with_lit = process_row(row, {"9000007": records}, MAP_TO_MONDO, genes, {"9000007"})
    edges = [e for e in with_lit if isinstance(e, VariantToDiseaseAssociation)]
    assert len(edges) == 1
    assert edges[0].predicate == "biolink:associated_with_increased_likelihood_of"


def test_literature_only_variants_reads_collection_method():
    """The publication signal is CollectionMethod on a P/LP record, not the mere existence
    of a citation for the variant."""
    from clinvar_helpers import literature_only_variants

    lit = _make_record("Pathogenic", "C2981140", "Glaucoma_of_childhood")
    lit["CollectionMethod"] = "literature only"
    clinical = _make_record("Pathogenic", "C2981140", "Glaucoma_of_childhood")
    benign_lit = _make_record("Benign", "C2981140", "Glaucoma_of_childhood")
    benign_lit["CollectionMethod"] = "literature only"

    found = literature_only_variants({"a": [lit], "b": [clinical], "c": [benign_lit]})
    assert found == {"a"}, "only a literature-derived P/LP record counts"


def test_pair_below_variant_threshold_is_dropped():
    """A gene-disease pair supported by fewer than min_variants_per_pair distinct variants
    is dropped, even when the variant itself qualifies on review status."""
    row = {**_TEST_ROW_TEMPLATE, "ID": "1296989",
           "CLNDISDB": "MONDO:MONDO:0100283,MedGen:CN300503"}
    pair = ("HGNC:3942", "MONDO:0100283")

    below = process_row(row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES,
                        frozenset(), {pair: 1})
    assert [e for e in below if isinstance(e, VariantToDiseaseAssociation)] == []

    at = process_row(row, VAR_RECORDS, MAP_TO_MONDO, VARIANT_GENES,
                     frozenset(), {pair: 2})
    assert len([e for e in at if isinstance(e, VariantToDiseaseAssociation)]) == 1


def test_strong_evidence_never_downgraded_by_publication_tier():
    """A variant that already qualifies at >=2 stars keeps `causes` for that disease even
    if it also has a citation -- the tiers cannot both claim the same (variant, disease)."""
    row = {**_TEST_ROW_TEMPLATE, "ID": "9000008",
           "CLNREVSTAT": "criteria_provided,_multiple_submitters,_no_conflicts"}
    records = [_make_record("Pathogenic", "C2981140", "Glaucoma_of_childhood",
                            review_status="criteria provided, single submitter")]
    edges = [e for e in process_row(row, {"9000008": records}, MAP_TO_MONDO,
                                    {"9000008": ("HGNC:7610", "MYOC")}, {"9000008"})
             if isinstance(e, VariantToDiseaseAssociation)]
    assert len(edges) == 1
    assert edges[0].predicate == "biolink:causes"


def test_zero_star_variant_emits_nothing():
    """A 0-star aggregate is the absence of stated criteria, not weak evidence, so it is
    excluded by every path -- including the publication tier, which used to accept
    anything at or below `publication_star_max`."""
    row = {**_TEST_ROW_TEMPLATE, "ID": "9000010",
           "CLNREVSTAT": "no_assertion_criteria_provided"}
    records = [_make_record("Pathogenic", "C2981140", "Glaucoma_of_childhood",
                            review_status="reviewed by expert panel")]
    edges = process_row(row, {"9000010": records}, MAP_TO_MONDO,
                        {"9000010": ("HGNC:7610", "MYOC")}, {"9000010"})
    assert [e for e in edges if isinstance(e, VariantToDiseaseAssociation)] == []


def test_zero_star_records_do_not_count_toward_concordance():
    """The concordance rescue deliberately ignores how *high* each submitter's review
    status is, but two submitters who both state no assertion criteria are not
    independent corroboration, so neither may count."""
    row = {**_TEST_ROW_TEMPLATE, "ID": "9000011",
           "CLNREVSTAT": "criteria_provided,_single_submitter"}
    zero = [_make_record("Pathogenic", "C2981140", "Glaucoma_of_childhood",
                         review_status="no assertion criteria provided", submitter=s)
            for s in ("LabA", "LabB")]
    assert [e for e in process_row(row, {"9000011": zero}, MAP_TO_MONDO,
                                   {"9000011": ("HGNC:7610", "MYOC")})
            if isinstance(e, VariantToDiseaseAssociation)] == []

    # the same two submitters, each stating criteria, do corroborate
    one = [_make_record("Pathogenic", "C2981140", "Glaucoma_of_childhood",
                        review_status="criteria provided, single submitter", submitter=s)
           for s in ("LabA", "LabB")]
    edges = [e for e in process_row(row, {"9000011": one}, MAP_TO_MONDO,
                                    {"9000011": ("HGNC:7610", "MYOC")})
             if isinstance(e, VariantToDiseaseAssociation)]
    assert len(edges) == 1
    assert edges[0].predicate == "biolink:causes"


def test_gene_edges_use_hgnc_ids():
    """The Monarch KG keys human genes on HGNC -- its NCBIGene nodes are other species --
    so a VariantToGeneAssociation subjected on NCBIGene resolves to nothing on merge."""
    row = {**_TEST_ROW_TEMPLATE, "ID": "9000012",
           "CLNREVSTAT": "criteria_provided,_multiple_submitters,_no_conflicts"}
    records = [_make_record("Pathogenic", "C2981140", "Glaucoma_of_childhood")]
    gene_assocs = [e for e in process_row(row, {"9000012": records}, MAP_TO_MONDO,
                                          {"9000012": ("HGNC:7610", "MYOC")})
                   if isinstance(e, VariantToGeneAssociation)]
    assert len(gene_assocs) == 1
    assert gene_assocs[0].object.startswith("HGNC:")
