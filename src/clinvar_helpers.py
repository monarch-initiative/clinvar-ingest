import gzip
import uuid

from biolink_model.datamodel.pydanticmodel_v2 import (
    AgentTypeEnum,
    KnowledgeLevelEnum,
    SequenceVariant,
    VariantToDiseaseAssociation,
    VariantToGeneAssociation,
)

# Variant to gene predicate
IS_SEQUENCE_VARIANT_OF = "biolink:is_sequence_variant_of"

# Variant to disease
CAUSES = "biolink:causes"
ASSOCIATED_WITH = "biolink:associated_with_increased_likelihood_of"

pred_to_negated = {CAUSES: False, ASSOCIATED_WITH: False}

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
# ClinVar's own variant-level aggregate review status (the VCF's CLNREVSTAT). Its 2-star
# tier -- "criteria provided, multiple submitters, no conflicts" -- is computed ACROSS a
# variant's submitters and never appears on an individual submission record, so the
# per-record star filter above can never see it. Accepting it as a third inclusion path
# uses ClinVar's own multi-submitter inference rather than re-deriving one.
aggregate_star_min = 2

# Variant classes this ingest keeps (the VCF's CLNVC field). SNVs and indels only --
# Microsatellite (38,744 rows), Inversion (1,519) and the catch-all "Variation" (459)
# are excluded: repeat expansions and inversions are not well represented by a fixed
# REF/ALT, and "Variation" carries no class at all. Adjust this set to change scope.
KEPT_VARIANT_CLASSES = frozenset(
    {
        "single_nucleotide_variant",
        "Deletion",
        "Duplication",
        "Indel",
        "Insertion",
    }
)

# Predicate by evidence strength, not by ClinicalSignificance. Every Pathogenic-family
# classification asserts causation; what differs is how well corroborated it is.
#   >=aggregate_star_min (ClinVar's own cross-submitter aggregate)  -> causes
#   exactly 1 star, but derived from published evidence             -> associated_with
# A 1-star call with no published basis is not emitted at all, and 0 stars is excluded
# outright by min_review_stars below -- so "<=publication_star_max" means exactly 1.
#
# "Published evidence" means the submitting lab recorded CollectionMethod == "literature
# only" on its Pathogenic/Likely-pathogenic record -- i.e. the classification itself came
# from the literature. An earlier version used "the variant has a PubMed citation" from
# var_citations.txt; that admitted 2,438,423 variants (55% of the VCF), because large
# cohort papers cite thousands at once and a citation asserts nothing about pathogenicity.
# This signal is attached to the assertion rather than the variant, and needs no extra file.
publication_star_max = 1

# Hard floor: nothing with a 0-star review status enters the graph, by any path.
#
# 0 stars is not "weak evidence", it is the absence of an assertion: the aggregate values
# scored 0 in review_star_map are no_assertion_criteria_provided, no_classification_provided,
# no_classifications_from_unflagged_records, no_classification_for_the_single_variant and
# flagged_submission. A classification with no stated criteria cannot be weighed, and a
# published basis does not repair that -- so the publication tier does not get to admit it
# either, and neither does the multi-submitter concordance rescue (several submitters
# agreeing without criteria is still nobody stating criteria).
#
# Applied on BOTH star axes, which are different measurements:
#   row["CLNREVSTAT"]     ClinVar's aggregate ACROSS a variant's submitters (0-4)
#   rec["ReviewStatus"]   one submission record's own status (0, 1, 3 or 4 -- never 2)
min_review_stars = 1
LITERATURE_ONLY = "literature only"

# A gene-disease pair must be supported by at least this many distinct variants. Inclusion
# is otherwise decided per variant, so a pair could enter the graph on a single variant's
# evidence; requiring corroboration across variants is a different axis from ClinVar's
# per-variant review status and disproportionately removes pairs resting only on the
# weaker publication tier.
min_variants_per_pair = 2

# Which ClinicalSignificance values are eligible at all. "Likely pathogenic" expresses
# the curator's confidence in a causal claim, not a weaker kind of relationship, so it is
# not downgraded here. The predicate is chosen by
# evidence tier (see above), not by which of these matched.
predicate_map = {
    "Pathogenic": CAUSES,
    "Pathogenic, low penetrance": CAUSES,
    "Pathogenic/Likely pathogenic": CAUSES,
    "Likely pathogenic": CAUSES,
    "Likely pathogenic, low penetrance": CAUSES,
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


# variant_summary.txt.gz carries one row per variant per genome build, so a variant
# typically appears on both GRCh37 and GRCh38. Take the most-preferred build present for
# each VariationID and never more than one -- filtering to GRCh38 alone would silently
# drop variants ClinVar has only ever placed on GRCh37, while taking every row would
# count the same variant twice.
ASSEMBLY_PREFERENCE = ("GRCh38", "GRCh37")


def make_variant_gene_map(variant_summary_path, assembly_preference=ASSEMBLY_PREFERENCE):
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
    variant still contributes its disease edges. The same applies
    across builds: if the preferred build's row says -1, that is ClinVar's
    answer and a lower-preference build is not consulted to override it.

    The emitted id is the HGNC id, not NCBIGene, because that is the id space the
    Monarch KG keys HUMAN genes on. Its NCBIGene nodes are other species (dog, rat,
    ...), so NCBIGene-subjected VariantToGeneAssociation edges resolve to nothing on
    merge -- they dangle rather than error, which is why this went unnoticed. A row
    with no HGNC_ID is treated exactly like GeneID == -1: no gene edge, disease edges
    unaffected. variant_summary.txt.gz carries HGNC_ID alongside GeneID, so this costs
    no extra input.
    """
    gene_map = {}
    best_rank = {}
    symbol_pool = {}
    rank_of = {name: i for i, name in enumerate(assembly_preference)}
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
            rank = rank_of.get(cols[hcols["Assembly"]])
            if rank is None:
                continue
            varid = cols[hcols["VariationID"]]
            seen = best_rank.get(varid)
            if seen is not None and seen <= rank:
                continue
            best_rank[varid] = rank
            gene_id = cols[hcols["GeneID"]]
            hgnc_id = cols[hcols["HGNC_ID"]].strip()
            # Both must be present. GeneID == -1 is ClinVar declining to attribute the
            # variant at all; a blank or "-" HGNC_ID means the gene it attributes has no
            # HGNC record (mostly LOC placeholders and non-coding loci), and an id the KG
            # cannot resolve is worse than no edge -- it dangles silently.
            if gene_id == "-1" or not gene_id or not hgnc_id or hgnc_id == "-":
                # the preferred build declines to attribute a gene -- drop any attribution
                # picked up from a lower-preference build rather than letting it stand
                gene_map.pop(varid, None)
                continue
            symbol = cols[hcols["GeneSymbol"]]
            # gene symbols repeat across millions of rows -- intern them so the map
            # holds one string per gene rather than one per variant
            symbol = symbol_pool.setdefault(symbol, symbol)
            # HGNC ids repeat as often as symbols do (one per gene across millions of
            # rows), so intern them through the same pool rather than holding ~4.5M
            # separate equal strings.
            hgnc_id = symbol_pool.setdefault(hgnc_id, hgnc_id)
            gene_map[varid] = (hgnc_id, symbol)
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
    distinct Submitters across record_list. Used by variant_records_to_disease()
    to rescue gene-disease pairs whose individual submission records are all
    below star_min but where multiple independent submitters agree on the same
    disease and the same classification.

    Submitters below min_review_stars do not count toward the total. The rescue
    deliberately ignores how *high* each submitter's own review status is -- that is
    the point of it -- but a 0-star record states no assertion criteria at all, and
    several submitters stating no criteria is not independent corroboration.
    """
    groups = {}
    for rec in record_list:
        clinsig = rec["ClinicalSignificance"]
        if clinsig not in predicate_map:
            continue
        if review_star_map[rec["ReviewStatus"].replace(" ", "_")] < min_review_stars:
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


def variant_records_to_disease(record_list, map_to_mondo, star_min=3, rescue_min_submitters=None, aggregate_stars=0):
    concordant_pairs = (
        concordant_disease_pairs(record_list, map_to_mondo, rescue_min_submitters) if rescue_min_submitters else set()
    )
    # ClinVar has already aggregated this variant to >=2 stars across its submitters, which
    # is the same signal the concordance rescue reconstructs -- accept every P/LP mapping
    # without re-testing per-record stars or exact-string agreement.
    accept_on_aggregate = aggregate_stars >= aggregate_star_min

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
                if not accept_on_aggregate and stars < star_min and (d, clinsig) not in concordant_pairs:
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


def _disease_edge(subject, dis_id, predicate, og_preds, row):
    """One VariantToDiseaseAssociation. Exactly one predicate is emitted per
    (variant, disease): the strong tier claims a disease first, and the publication
    tier only sees what is left, so the two can never contradict each other."""
    return VariantToDiseaseAssociation(
        id=str(uuid.uuid4()),
        subject=subject,
        predicate=predicate,
        qualifiers=[row["CLNREVSTAT"]],
        object=dis_id,
        negated=pred_to_negated[predicate],
        original_predicate=":".join(og_preds),
        primary_knowledge_source="infores:clinvar",
        aggregator_knowledge_source=["infores:monarchinitiative"],
        knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
        agent_type=AgentTypeEnum.manual_agent,
    )


def literature_only_variants(var_records):
    """VariationIDs with a Pathogenic/Likely-pathogenic submission whose CollectionMethod
    is "literature only" -- the lab's pathogenic call was derived from published evidence.

    Computed from var_records, so it costs no extra I/O and no extra download.
    """
    return {
        varid
        for varid, records in var_records.items()
        if any(
            rec["CollectionMethod"] == LITERATURE_ONLY
            and rec["ClinicalSignificance"] in predicate_map
            for rec in records
        )
    }


def qualifying_diseases(row, records, map_to_mondo, lit_only_variants):
    """Which diseases this variant qualifies for, and under which predicate.

    Returns (causes, causes_org, assoc, assoc_org). The two tiers are disjoint: the
    strong tier claims a disease first and the publication tier only sees what is left,
    so no (variant, disease) can carry both predicates.

    Shared by build_pair_variant_counts() and process_row() so the pre-pass that counts
    supporting variants per pair and the pass that emits edges cannot disagree.
    """
    aggregate_stars = review_star_map.get(row["CLNREVSTAT"], 0)

    # Nothing with a 0-star aggregate qualifies by any path -- see min_review_stars.
    if aggregate_stars < min_review_stars:
        return {}, {}, {}, {}

    causes, _preds, causes_org = variant_records_to_disease(
        records,
        map_to_mondo,
        star_min=var2disease_star_min,
        rescue_min_submitters=min_concordant_submitters,
        aggregate_stars=aggregate_stars,
    )

    assoc: dict = {}
    assoc_org: dict = {}
    if aggregate_stars <= publication_star_max and row["ID"] in lit_only_variants:
        # star_min=min_review_stars, not 0: a 0-star submission record cannot supply the
        # disease term either, or the floor would leak back in on the record axis.
        all_ids, _all_preds, all_org = variant_records_to_disease(
            records, map_to_mondo, star_min=min_review_stars
        )
        for d in all_ids:
            if d not in causes:
                assoc[d] = ""
                assoc_org[d] = all_org[d]

    return causes, causes_org, assoc, assoc_org


def build_pair_variant_counts(clinvar_tsv, var_records, map_to_mondo, variant_genes, lit_only_variants):
    """(gene, disease) -> how many distinct variants qualify for it.

    A pre-pass over the same rows the transform will stream, applying the same tier logic,
    so process_row() can drop pairs supported by fewer than min_variants_per_pair variants.
    Inclusion is otherwise a per-variant decision and a pair can enter the graph on one
    variant alone.
    """
    import csv as _csv

    counts = {}
    with open(clinvar_tsv, newline="") as fh:
        for row in _csv.DictReader(fh, delimiter="\t"):
            varid = row["ID"]
            records = var_records.get(varid)
            if records is None or row["CLNVC"] not in KEPT_VARIANT_CLASSES:
                continue
            gene_entry = variant_genes.get(varid)
            if gene_entry is None:
                continue
            causes, _co, assoc, _ao = qualifying_diseases(row, records, map_to_mondo, lit_only_variants)
            for d in list(causes) + list(assoc):
                key = (gene_entry[0], d)
                counts[key] = counts.get(key, 0) + 1
    return counts


def process_row(
    row,
    var_records,
    map_to_mondo,
    variant_genes,
    lit_only_variants=frozenset(),
    pair_variant_counts=None,
):
    """Process a single row from the ClinVar VCF and return a list of biolink entities.

    Returns an empty list if the row should be skipped (no records, no associations).

    Gene attribution comes from variant_genes (see make_variant_gene_map), NOT from
    the row's GENEINFO field -- GENEINFO lists every locus overlapping the variant's
    position, which would mint gene-disease associations for antisense transcripts and
    neighbours. A variant ClinVar declines to attribute to a gene simply produces no
    VariantToGeneAssociation; its disease edges are unaffected.

    No phenotype edges are emitted. ClinVar's HPO ids live inside CLNDISDB groups
    alongside the MONDO/MedGen/OMIM ids for the SAME condition -- they are
    cross-references naming that condition in HPO's vocabulary, supplied by MedGen
    rather than by the submitter, and measured at 0.9994 mean consistency across a
    disease's variants. A VariantToPhenotypicFeatureAssociation built from them
    restated the disease edge in a second vocabulary instead of adding phenotype
    information, so they were removed.
    """
    entities = []

    varid = str(row["ID"])
    raw_diss_info = row["CLNDISDB"]
    # SequenceVariant.type carries the variant CLASS (SNV / deletion / duplication ...) from
    # CLNVCSO, not the molecular consequence from MC (missense / frameshift ...) -- the class
    # is a property of the variant, the consequence is a property of variant x transcript.
    variant_class = [row["CLNVCSO"]] if row["CLNVCSO"].startswith("SO:") else []

    if varid not in var_records:
        return []

    # Prune to SNVs and indels up front -- see KEPT_VARIANT_CLASSES
    if row["CLNVC"] not in KEPT_VARIANT_CLASSES:
        return []

    gene_entry = variant_genes.get(varid)
    gene_ids = [gene_entry[0]] if gene_entry else []

    records = var_records[varid]
    disease_ids, org_predicates, assoc_ids, assoc_org = qualifying_diseases(
        row, records, map_to_mondo, lit_only_variants
    )

    # A pair supported by fewer than min_variants_per_pair distinct variants is dropped.
    # Without this, a gene-disease pair enters the graph on one variant's evidence.
    if pair_variant_counts is not None and gene_ids:
        gene_id = gene_ids[0]
        disease_ids = {
            d: v for d, v in disease_ids.items()
            if pair_variant_counts.get((gene_id, d), 0) >= min_variants_per_pair
        }
        assoc_ids = {
            d: v for d, v in assoc_ids.items()
            if pair_variant_counts.get((gene_id, d), 0) >= min_variants_per_pair
        }

    # Corroboration check: at least one disease derived from the submission records must
    # also appear in the VCF's own aggregate CLNDISDB list for this variant. This gate
    # predates the removal of phenotype edges (it used to double as "are there HPO terms
    # to emit?"), but it stands on its own as a cross-check between the two input files.
    # All-or-nothing per variant, and very nearly inert in practice: it drops 2 of 56,270.
    diss_info = parse_CLNDISDB(raw_diss_info)
    diss_info, _ = map_CLNDISDB_to_mondo(diss_info, map_to_mondo)
    # Both tiers are subject to the gate -- a publication-tier disease must be echoed in
    # CLNDISDB just as a causes-tier one must, so the union is what gets checked.
    corroborated = map_mondo_to_hp(diss_info, {**disease_ids, **assoc_ids})

    if len(corroborated) == 0:
        return []

    seq_var = SequenceVariant(
        id="CLINVAR:{}".format(row["ID"]),
        name=row["CLNHGVS"],
        # RS is "." on variants with no dbSNP mapping -- emitting "DBSNP:." would be a junk CURIE
        xref=["DBSNP:{}".format(row["RS"])] if row["RS"] != "." else [],
        has_gene=gene_ids,
        in_taxon=["NCBITaxon:9606"],
        in_taxon_label="Homo sapiens",
        type=variant_class,
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

    for dis_id in disease_ids:
        og_preds = sorted(list(org_predicates[dis_id].keys()))
        entities.append(_disease_edge(seq_var.id, dis_id, CAUSES, og_preds, row))

    # <=1 star with published support: a weaker claim, so a weaker predicate
    for dis_id in assoc_ids:
        og_preds = sorted(list(assoc_org[dis_id].keys()))
        entities.append(_disease_edge(seq_var.id, dis_id, ASSOCIATED_WITH, og_preds, row))

    return entities
