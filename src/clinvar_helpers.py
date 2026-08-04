import gzip
import uuid

from biolink_model.datamodel.pydanticmodel_v2 import (
    AgentTypeEnum,
    KnowledgeLevelEnum,
    SequenceVariant,
    VariantToDiseaseAssociation,
    VariantToGeneAssociation,
    VariantToPhenotypicFeatureAssociation,
)

# Variant to gene predicate
IS_SEQUENCE_VARIANT_OF = "biolink:is_sequence_variant_of"

# Variant to disease
CAUSES = "biolink:causes"
ASSOCIATED_WITH = "biolink:associated_with_increased_likelihood_of"

# Variant to phenotype
CONTRIBUTES_TO = "biolink:contributes_to"

pred_to_negated = {CAUSES: False, ASSOCIATED_WITH: False, CONTRIBUTES_TO: False}

# Star ratings exactly as ClinVar documents them:
# https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
# Note "criteria_provided,_multiple_submitters,_no_conflicts" (2 stars) is an
# aggregate ClinVar computes *across* a variant's submitters -- it appears in the
# VCF's variant-level CLNREVSTAT but never on an individual submission record in
# submission_summary.txt, so a per-record star_min=2 filter is identical to
# star_min=3. See analysis/clinvar_report.py sections 1-3.
review_star_map = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_classifications_from_unflagged_records": 0,
    "no_classification_for_the_single_variant": 0,
    "flagged_submission": 0,
    ".": 0,
}

var2disease_star_min = 3
min_concordant_submitters = 2

predicate_map = {
    "Pathogenic": CAUSES,
    "Pathogenic, low penetrance": CAUSES,
    "Pathogenic/Likely pathogenic": CAUSES,
    "Likely pathogenic": ASSOCIATED_WITH,
    "Likely pathogenic, low penetrance": ASSOCIATED_WITH,
}


def make_variant_record_map(submission_path):
    var_records = {}
    with gzip.open(submission_path, "rt") as infile:
        for line in infile:
            line = line.strip("\r").strip("\n")
            if line[0] == "#":
                header = line.split("\t")
                header[0] = header[0][1:]
                hcols = {k: i for i, k in enumerate(header)}
            else:
                cols = line.split("\t")
                varid = cols[hcols["VariationID"]]
                if varid not in var_records:
                    var_records[varid] = []
                rec = {k: cols[hcols[k]] for k in hcols}
                var_records[varid].append(rec)
    return var_records


def make_mondo_map(sssom_path):
    hcount = 0
    dups = 0
    map_to_mondo = {}
    with open(sssom_path, "r") as infile:
        for line in infile:
            line = line.strip("\r").strip("\n")
            if not line:
                continue
            if line[0] == "#":
                continue
            else:
                hcount += 1
                if hcount == 1:
                    header = line.split("\t")
                    hcols = {k: i for i, k in enumerate(header)}
                    continue
                cols = line.split("\t")
                obj_id, subj_id = cols[hcols["object_id"]], cols[hcols["subject_id"]]
                if obj_id not in map_to_mondo:
                    map_to_mondo[obj_id] = {}
                else:
                    dups += 1
                map_to_mondo[obj_id][subj_id] = ""

    # Add mondo_id mapping to self
    mondo_set = {k.split(":")[-1]: "" for kv in map_to_mondo for k in map_to_mondo[kv]}
    for k in mondo_set:
        v = "MONDO:{}".format(k)
        kv = "MONDO:MONDO:{}".format(k)
        map_to_mondo[v] = {v: ""}
        map_to_mondo[kv] = {v: ""}

    return map_to_mondo


def make_medgen_to_mondo_map(medgen_path):
    map_to_mondo = {}
    with gzip.open(medgen_path, "rt") as infile:
        for line in infile:
            line = line.strip("\r").strip("\n")
            if line[0] == "#":
                continue
            else:
                info = line.split("|")
                dis_id = info[2]
                mdg_id = "MedGen:{}".format(info[0])
                if "MONDO:" in dis_id:
                    if mdg_id not in map_to_mondo:
                        map_to_mondo[mdg_id] = {}
                    map_to_mondo[mdg_id][dis_id] = ""
    return map_to_mondo


def make_variant_gene_map(variant_summary_path, assembly="GRCh38"):
    """ClinVar's own per-variant gene attribution, from variant_summary.txt.gz.

    The VCF's GENEINFO field is populated *positionally* -- it lists every gene
    whose span covers the variant, so antisense transcripts, divergent
    transcripts and neighbouring loci ride along with the causal gene. Building
    VariantToGeneAssociation edges from it asserts non-viable gene-disease
    relationships: a CFTR variant also tags CFTR-AS1, which then inherits
    cystic fibrosis and CFTR's entire submitter roster.

    variant_summary.txt.gz carries the gene ClinVar itself assigns to each
    variant record. Measured over the current release, it makes 4,438,486
    attributions where GENEINFO makes 4,982,691, and reduces the
    antisense/uncharacterized-ORF share from 2.97% to 0.09%. Every one of the
    479,140 variants GENEINFO assigns to >1 gene collapses to exactly one
    curated GeneID.

    GeneID == -1 means ClinVar declines to attribute the variant to a gene
    (4,109 variants, 0.09%). Those are omitted here and deliberately get NO
    fallback to GENEINFO -- an unasserted gene is left unasserted, and the
    variant still contributes its disease and phenotype edges.
    """
    gene_map = {}
    symbol_pool = {}
    hcols = None
    with gzip.open(variant_summary_path, "rt") as infile:
        for line in infile:
            line = line.strip("\r").strip("\n")
            if not line:
                continue
            if line[0] == "#":
                header = line.split("\t")
                header[0] = header[0][1:]
                hcols = {k: i for i, k in enumerate(header)}
                continue
            cols = line.split("\t")
            if cols[hcols["Assembly"]] != assembly:
                continue
            gene_id = cols[hcols["GeneID"]]
            if gene_id == "-1" or not gene_id:
                continue
            symbol = cols[hcols["GeneSymbol"]]
            # gene symbols repeat across millions of rows -- intern them so the map
            # holds one string per gene rather than one per variant
            symbol = symbol_pool.setdefault(symbol, symbol)
            gene_map[cols[hcols["VariationID"]]] = ("NCBIGene:{}".format(gene_id), symbol)
    return gene_map


def make_genes_from_row(gene_list):
    if gene_list == ".":
        return [], []
    gene_ids, gene_symbols = [], []
    genes = gene_list.split("|")
    for gene in genes:
        values = gene.split(":")[1:]
        gene_sym = gene.split(":")[0]
        for value in values:
            gene_ids.append("NCBIGene:{}".format(value))
            gene_symbols.append(gene_sym)
    return gene_ids, gene_symbols


def format_id_to_map(info):
    idnum = info.split(":")[-1]
    if "MONDO:" in info:
        idname = "MONDO:{}".format(idnum)
    elif "HP:" in info:
        idname = "HP:{}".format(idnum)
    elif "MeSH:" in info:
        idname = "mesh:{}".format(idnum)
    elif "." == info:
        idname = None
    else:
        idname = info
    return idname


def concordant_disease_pairs(record_list, map_to_mondo, min_submitters):
    """(mondo_id, ClinicalSignificance) pairs supported by >=min_submitters
    distinct Submitters across record_list, regardless of each submitter's
    own review status. Used by variant_records_to_disease() to rescue
    gene-disease pairs whose individual submission records are all below
    star_min but where multiple independent submitters agree on the same
    disease and the same classification.
    """
    groups = {}
    for rec in record_list:
        clinsig = rec["ClinicalSignificance"]
        if clinsig not in predicate_map:
            continue

        mondo_ids = set()
        for mg_mapping in rec["ReportedPhenotypeInfo"].split(";"):
            mg_map = "MedGen:{}".format(mg_mapping.split(":")[0])
            if mg_map in map_to_mondo:
                mondo_ids.update(map_to_mondo[mg_map].keys())

        if not mondo_ids:
            for dis_id in rec["SubmittedPhenotypeInfo"].split(";"):
                dis_id = format_id_to_map(dis_id)
                if dis_id in map_to_mondo:
                    mondo_ids.update(map_to_mondo[dis_id].keys())
                elif dis_id is not None and "MONDO:" in dis_id:
                    mondo_ids.add(dis_id)

        for mondo_id in mondo_ids:
            groups.setdefault((mondo_id, clinsig), set()).add(rec["Submitter"])

    return {key for key, submitters in groups.items() if len(submitters) >= min_submitters}


def variant_records_to_disease(record_list, map_to_mondo, star_min=3, rescue_min_submitters=None):
    concordant_pairs = (
        concordant_disease_pairs(record_list, map_to_mondo, rescue_min_submitters) if rescue_min_submitters else set()
    )

    dis = {}
    preds = {}
    org_preds = {}
    for rec in record_list:
        stars = int(review_star_map[rec["ReviewStatus"].replace(" ", "_")])
        clinsig = rec["ClinicalSignificance"]

        if clinsig not in predicate_map:
            continue

        mapped_predicate = predicate_map[clinsig]
        org_predicate = clinsig

        mapped_terms = 0
        for mg_mapping in rec["ReportedPhenotypeInfo"].split(";"):
            mg_map = "MedGen:{}".format(mg_mapping.split(":")[0])
            mondo_ids = []

            if mg_map in map_to_mondo:
                mondo_ids = list(map_to_mondo[mg_map].keys())

            if len(mondo_ids) == 0:
                continue

            for d in mondo_ids:
                mapped_terms += 1
                if stars < star_min and (d, clinsig) not in concordant_pairs:
                    continue
                dis[d] = ""
                if d not in preds:
                    preds[d] = {}
                    org_preds[d] = {}
                preds[d][mapped_predicate] = ""
                org_preds[d][org_predicate] = ""

        if mapped_terms == 0:
            for dis_id in rec["SubmittedPhenotypeInfo"].split(";"):
                dis_id = format_id_to_map(dis_id)
                mondo_ids = []

                if dis_id in map_to_mondo:
                    mondo_ids = list(map_to_mondo[dis_id].keys())
                elif dis_id is not None and "MONDO:" in dis_id:
                    mondo_ids = [dis_id]

                if len(mondo_ids) == 0:
                    continue

                for d in mondo_ids:
                    mapped_terms += 1
                    if stars < star_min and (d, clinsig) not in concordant_pairs:
                        continue
                    dis[d] = ""
                    if d not in preds:
                        preds[d] = {}
                        org_preds[d] = {}
                    preds[d][mapped_predicate] = ""
                    org_preds[d][org_predicate] = ""

    return dis, preds, org_preds


def parse_CLNDISDB(column):
    diss = []
    for group_info in column.split("|"):
        default = {"MAP_TERMS": [], "HP": []}
        for info in group_info.split(","):
            idname = format_id_to_map(info)
            if idname is None:
                continue
            if "HP:" in idname:
                default["HP"].append(idname)
            else:
                default["MAP_TERMS"].append(idname)
        diss.append(default)
    diss = [d for d in diss if len(d["HP"]) > 0 or len(d["MAP_TERMS"]) > 0]
    return diss


def map_CLNDISDB_to_mondo(parse_results, map_to_mondo, map_stats=None):
    if map_stats is None:
        map_stats = {"MONDO": 0, "mesh": 0, "OMIM": 0, "Orphanet": 0}
    for i, d in enumerate(parse_results):
        map_terms = d["MAP_TERMS"]
        mondo_ids = []
        for gterm in map_terms:
            if "MONDO:" in gterm:
                mondo_ids.append(gterm)
                map_stats["MONDO"] += 1
            elif gterm in map_to_mondo:
                mondo_ids += list(map_to_mondo[gterm].keys())
                for k in map_stats:
                    if k in gterm:
                        map_stats[k] += 1
            else:
                unknown = "unkown_source_{}".format("".join(gterm.split(":")[:-1]))
                if unknown not in map_stats:
                    map_stats[unknown] = 0
                map_stats[unknown] += 1
        parse_results[i]["MAP_TERMS"] = mondo_ids
    return parse_results, map_stats


def map_mondo_to_hp(group_info, disease_ids):
    mondo_to_hp = {}
    for g in group_info:
        for d in set(g["MAP_TERMS"]):
            if d in disease_ids:
                if d not in mondo_to_hp:
                    mondo_to_hp[d] = []
                mondo_to_hp[d] += g["HP"]
    return mondo_to_hp


def process_row(row, var_records, map_to_mondo, variant_genes):
    """Process a single row from the ClinVar VCF and return a list of biolink entities.

    Returns an empty list if the row should be skipped (no records, no associations).

    Gene attribution comes from variant_genes (see make_variant_gene_map), NOT from
    the row's GENEINFO field -- GENEINFO lists every locus overlapping the variant's
    position, which would mint gene-disease associations for antisense transcripts and
    neighbours. A variant ClinVar declines to attribute to a gene simply produces no
    VariantToGeneAssociation; its disease and phenotype edges are unaffected.
    """
    entities = []

    varid = str(row["ID"])
    raw_diss_info = row["CLNDISDB"]
    so_info = [v.split("|")[0] for v in row["MC"].split(",") if "SO:" in v]

    if varid not in var_records:
        return []

    gene_entry = variant_genes.get(varid)
    gene_ids = [gene_entry[0]] if gene_entry else []

    disease_ids, disease_predicates, org_predicates = variant_records_to_disease(
        var_records[varid],
        map_to_mondo,
        star_min=var2disease_star_min,
        rescue_min_submitters=min_concordant_submitters,
    )

    diss_info = parse_CLNDISDB(raw_diss_info)
    diss_info, _ = map_CLNDISDB_to_mondo(diss_info, map_to_mondo)
    mondo_to_hp = map_mondo_to_hp(diss_info, disease_ids)

    if len(mondo_to_hp) == 0:
        return []

    seq_var = SequenceVariant(
        id="CLINVAR:{}".format(row["ID"]),
        name=row["CLNHGVS"],
        xref=["DBSNP:{}".format(row["RS"])],
        has_gene=gene_ids,
        in_taxon=["NCBITaxon:9606"],
        in_taxon_label="Homo sapiens",
        type=so_info,
    )
    entities.append(seq_var)

    for gene_id in gene_ids:
        entities.append(
            VariantToGeneAssociation(
                id=str(uuid.uuid4()),
                subject=seq_var.id,
                predicate=IS_SEQUENCE_VARIANT_OF,
                object=gene_id,
                primary_knowledge_source="infores:clinvar",
                aggregator_knowledge_source=["infores:monarchinitiative"],
                knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
                agent_type=AgentTypeEnum.manual_agent,
            )
        )

    for dis_id, predicate in disease_predicates.items():
        og_preds = sorted(list(org_predicates[dis_id].keys()))
        for pred in list(predicate.keys()):
            negated = pred_to_negated[pred]
            entities.append(
                VariantToDiseaseAssociation(
                    id=str(uuid.uuid4()),
                    subject=seq_var.id,
                    predicate=pred,
                    qualifiers=[row["CLNREVSTAT"]],
                    object=dis_id,
                    negated=negated,
                    original_predicate=":".join(og_preds),
                    primary_knowledge_source="infores:clinvar",
                    aggregator_knowledge_source=["infores:monarchinitiative"],
                    knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
                    agent_type=AgentTypeEnum.manual_agent,
                )
            )

    for mondo_id, hp_terms in mondo_to_hp.items():
        for hp_id in hp_terms:
            entities.append(
                VariantToPhenotypicFeatureAssociation(
                    id=str(uuid.uuid4()),
                    subject=seq_var.id,
                    predicate=CONTRIBUTES_TO,
                    object=hp_id,
                    primary_knowledge_source="infores:clinvar",
                    aggregator_knowledge_source=["infores:monarchinitiative"],
                    knowledge_level=KnowledgeLevelEnum.observation,
                    agent_type=AgentTypeEnum.manual_agent,
                )
            )

    return entities
