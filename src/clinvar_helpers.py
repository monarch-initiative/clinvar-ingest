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

# https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
review_star_map = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "no_classifications_from_unflagged_records": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 1,
    "no_classification_provided": 1,
    "no_classification_for_the_single_variant": 0,
    "flagged_submission": 0,
    ".": 0,
}

var2disease_star_min = 3

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


def variant_records_to_disease(record_list, map_to_mondo, star_min=3):
    dis = {}
    preds = {}
    org_preds = {}
    for rec in record_list:
        stars = int(review_star_map[rec["ReviewStatus"].replace(" ", "_")])
        clinsig = rec["ClinicalSignificance"]

        if stars < star_min:
            continue
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
                dis[d] = ""
                if d not in preds:
                    preds[d] = {}
                    org_preds[d] = {}
                preds[d][mapped_predicate] = ""
                org_preds[d][org_predicate] = ""
                mapped_terms += 1

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
                    dis[d] = ""
                    if d not in preds:
                        preds[d] = {}
                        org_preds[d] = {}
                    preds[d][mapped_predicate] = ""
                    org_preds[d][org_predicate] = ""
                    mapped_terms += 1

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


def process_row(row, var_records, map_to_mondo):
    """Process a single row from the ClinVar VCF and return a list of biolink entities.

    Returns an empty list if the row should be skipped (no records, no associations).
    """
    entities = []

    varid = str(row["ID"])
    raw_diss_info = row["CLNDISDB"]
    ginfo = row["GENEINFO"]
    so_info = [v.split("|")[0] for v in row["MC"].split(",") if "SO:" in v]

    if varid not in var_records:
        return []

    gene_ids, gene_symbols = make_genes_from_row(ginfo)

    disease_ids, disease_predicates, org_predicates = variant_records_to_disease(
        var_records[varid], map_to_mondo, star_min=var2disease_star_min
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
