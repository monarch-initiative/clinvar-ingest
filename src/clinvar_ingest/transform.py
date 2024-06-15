import uuid  # For generating UUIDs for associations
import gzip
import pandas as pd
#import matplotlib.pyplot as plt
from collections import Counter

from biolink_model.datamodel.pydanticmodel_v2 import *  # Replace * with any necessary data classes from the Biolink Model
from koza.cli_utils import get_koza_app


### Create variant--> record mapping {variant_id:[record,...]}
def make_variant_record_map(submission_filepath):
    
    # Last line with "#" character in it will be read as the header
    var_records = {}
    rec_count = 0
    with gzip.open(submission_filepath, 'rt') as infile:
        
        for line in infile:
            line = line.strip('\r').strip('\n')
            if line[0] == "#":
                header = line.split('\t')
                
                # Removes leading "#" character
                header[0] = header[0][1:]
                
                # Links column name to column index
                hcols = {k:i for i, k in enumerate(header)}
            
            else:
                cols = line.split('\t')
                varid = cols[hcols["VariationID"]]
                if varid not in var_records:
                    var_records.update({varid:[]})
                
                # Make our record and append
                rec = {k:cols[hcols[k]] for k in hcols}
                var_records[varid].append(rec)
                rec_count += 1
    
    print("- {} variants read into memory across {} records".format(format(len(var_records), ','), 
                                                                    format(rec_count, ',')))
    return var_records


### Create generalized disease_id to mondo_id map
def make_mondo_map(sssom_path):
    
    # Last line with "#" character in it will be read as the header
    var_records = {}
    rec_count = 0
    dups = 0
    hcount = 0
    map_to_mondo = {}
    with open(sssom_path, 'r') as infile:

        for line in infile:
            line = line.strip('\r').strip('\n')

            if not line:
                continue

            if line[0] == "#":
                continue
            else:
                hcount += 1

                # Make header to data
                if hcount == 1:
                    header = line.split('\t')

                    # Links column name to column index
                    hcols = {k:i for i, k in enumerate(header)}
                    continue

                cols = line.split('\t')
                obj_id, subj_id = cols[hcols["object_id"]], cols[hcols["subject_id"]]

                if obj_id not in map_to_mondo:
                    map_to_mondo.update({obj_id:{}})
                else:
                    dups += 1
                    #print(obj_id, map_to_mondo[obj_id])
                    #print(obj_id, line)
                    #print("###")

                # A small amount of objects will map to more than one mondo id (hence the dict structure)
                map_to_mondo[obj_id].update({subj_id:''})

    print("- Total mappings produced {}".format(len(map_to_mondo)))
    print("- Multi mapping objects found {}".format(format(dups, ',')))
    
    # Lastly, add in bannana and mondo_id mapping to self
    mondo_set = {k.split(":")[-1]:'' for kv in map_to_mondo for k in map_to_mondo[kv]}
    for k in mondo_set:
        v = "MONDO:{}".format(k)
        kv = "MONDO:MONDO:{}".format(k)
        
        map_to_mondo.update({v:{v:''}})
        map_to_mondo.update({kv:{v:''}})
    
    return map_to_mondo
    

def make_genes_from_row(gene_list):
    if gene_list == '.':
        return [], []
    gene_ids, gene_symbols = [], []
    genes = gene_list.split('|')
    for gene in genes:
        values = gene.split(':')[1:]  # sometimes there's more than one Entrez ID after the symbol
        gene_sym = gene.split(':')[0]
        for value in values:
            gene_ids.append("NCBIGene:{}".format(value))
            gene_symbols.append(gene_sym)
    return gene_ids, gene_symbols


def format_id_to_map(info):
    
    idnum = info.split(':')[-1]     
    if "MONDO:" in info:
        idname = "MONDO:{}".format(idnum)

    elif "HP" in info:
        idname = "HP:{}".format(idnum)

    elif "MeSH" in info:
        idname = "mesh:{}".format(idnum)
        
    elif "." == info:
        idname = None
        
    else:
        idname = info
    
    return idname
  

def get_variant_to_disease_predicate(term, pathogenicity_lookup):
    
    if term in pathogenic_lookup:
        predicate = "biolink:causes"
    else:
        predicate = "biolink:contributes_to"
    
    return predicate    


def variant_records_to_disease(record_list, review_star_map, map_to_mondo, pathogenicity_lookup, star_min=3):
    
    # Each record is a dictionary. We have a list of records that we can loop through
    # keys = 
    # VariationID, ClinicalSignificance, DateLastEvaluated, 
    # Description, SubmittedPhenotypeInfo, ReportedPhenotypeInfo, 
    # ReviewStatus, CollectionMethod, OriginCounts, Submitter, SCV
    dis = {}
    preds = {}
    org_preds = {}
    for rec in record_list:
        stars = int(review_star_map[rec["ReviewStatus"].replace(" ", "_")])
        if stars < star_min:
            continue
        
        clinsig = rec["ClinicalSignificance"]
        mapped_predicate = get_variant_to_disease_predicate(clinsig, pathogenicity_lookup)
        org_predicate = clinsig
        
        for dis_id in rec["SubmittedPhenotypeInfo"].split(';'):
            
            # Convert to MONDO_ID if we can
            dis_id = format_id_to_map(dis_id)
            mondo_ids = []
            if dis_id in map_to_mondo:
                mondo_ids = list(map_to_mondo[dis_id].keys())
            
            elif "MONDO:" in dis_id:
                mondo_ids = [dis_id]
            
            # Can't map this one back
            if len(mondo_ids) == 0:
                continue
            
            for d in mondo_ids:
                dis.update({d:''})
                if d not in preds:
                    preds.update({d:{}})
                    org_preds.update({d:{}})
                preds[d].update({mapped_predicate:''})
                org_preds[d].update({org_predicate:''})
                
    return dis, preds, org_preds
   

def parse_CLNDISDB(column):
    
    # Output datastructure [{"MAP_TERMS":[], "HP":[]}, {...]
    diss = []

    # Loop through each disease the variant is associated with
    for group_info in column.split('|'):

        # Each one can have multiple ids associated with it
        default = {"MAP_TERMS":[],
                   "HP":[]}
        
        # Can have multiple hp terms associated with a single disease (can tune formatting for each term)
        for info in group_info.split(','):
            
            idname = format_id_to_map(info)
            if idname == None:
                continue
            
            # Separate by hp terms and disease ids (that should be mapped back to mondo)
            if "HP:" in idname:
                default["HP"].append(idname)
            else:
                default["MAP_TERMS"].append(idname)
                
        # Add our parsed information to return datastructure
        diss.append(default)
    
    # Filter for things that have at least 1 or more hp or disease terms
    diss = [d for d in diss if len(d["HP"]) > 0 or len(d["MAP_TERMS"]) > 0]
    
    return diss


def map_CLNDISDB_to_mondo(parse_results, map_to_mondo, map_stats={"MONDO":0, "mesh":0, "OMIM":0, "Orphanet":0}):
    
    for i, d in enumerate(parse_results):
        
        map_terms = d["MAP_TERMS"]
        mondo_ids = []
        for gterm in map_terms:

            if "MONDO:" in gterm:
                mondo_ids.append(gterm)
                map_stats["MONDO"] += 1
                
            elif gterm in map_to_mondo:
                mondo_ids += list(map_to_mondo[gterm].keys())
                
                ### STATS KEEPING
                for k in map_stats:
                    if k in gterm:
                        map_stats[k] += 1
            
            else:
                unknown = "unkown_source_{}".format(''.join(gterm.split(":")[:-1]))
                if unknown not in map_stats:
                    map_stats.update({unknown:0})
                map_stats[unknown] += 1
        
        parse_results[i]["MAP_TERMS"] = mondo_ids
    
    return parse_results, map_stats


def map_mondo_to_hp(group_info, disease_ids):
    mondo_to_hp = {}
    for g in group_info:
        
        # Loop through each disease_id term and see if it matches any disease this variant is associated with
        for d in g["MAP_TERMS"]:
            if d in disease_ids:
                
                # Check if exists or not currently. We can pull in multiple hp terms from records that
                # have the same disease id as the "cononical" disease for this set of records
                if d not in mondo_to_hp:
                    mondo_to_hp.update({d:[]})
                    
                mondo_to_hp[d] += g["HP"]

    return mondo_to_hp


def map_records_to_disease_hp(record_list, 
                              clndisdb_column, 
                              review_star_map, 
                              map_to_mondo, 
                              pathogenicity_lookup, 
                              star_min=3):
    
    # Disease --> HP modeling pipeline
    
    # Variant records --> MONDO:ID based on review status of ACMG. Ratings < star_min will not have an association
    disease_ids, disease_predicates, org_predicates  = variant_records_to_disease(record_list, 
                                                                                   review_star_map, 
                                                                                   map_to_mondo, 
                                                                                   pathogenicity_lookup,
                                                                                   star_min)
    
    # Pull and format disease_ids and HP terms
    diss_info = parse_CLNDISDB(clndisdb_column)
    
    # Map inormation to mondo id if possible and discard otherwise
    diss_info = map_CLNDISDB_to_mondo(diss_info, map_to_mondo)
    
    # Map each disease to HP terms if possible
    mondo_to_hp = map_mondo_to_hp(diss_info, disease_ids)

    return mondo_to_hp, disease_predicates, org_predicates
    

########################################################################################
### Pipeline is to first create a variant to record(s) map {variantid:[{},{},{}...]} ###

### Then second is to create a disease id to mondo id map in the form of...
### map_to_mondo --> {"OMIM:123":"MONDO:123", ...}.
### where the key (disease id) is from the sources found within the sssom file

### The map_to_mondo is able to handle a few different types of data
### Orphanet, OMIM, Mondo, and some MeSH for this data set 

### Then we attempt to map all disease ids found within any gene back to this map (map_to_mondo).
### If MONDO is found within the id name and it isn't found in the map, then it is included anyways.
### This results in only making variantToDisease associations if the disease id is found in our map / is a mondo id.
### We are currently lacking MedGen which seems to be a large portion of the variant annotaitions within the vcf, but
### doesn't necessarily mean they all contribute to connections we want to make

### Map disease id --> mondo id
### Review status term of record --> stars determines if variantToDiesease association is made, AND HP associations (if possible)
### Can use 0,1,2,3,4 star mapping min to get different results (see review_star_map below)


koza_app = get_koza_app("clinvar_variant")

CONTRIBUTES_TO = "biolink:contributes_to"
CAUSES = "biolink:causes"
HAS_PHENOTYPE = "biolink:has_phenotype"
IS_SEQUENCE_VARIANT_OF = 'biolink:is_sequence_variant_of'

# Manually curated terms derived from files for data modeling purposes
review_star_map = {"practice_guideline":4,
                   "reviewed_by_expert_panel":3,
                   "criteria_provided,_multiple_submitters,_no_conflicts":2,
                   "criteria_provided,_conflicting_classifications":2,
                   "no_classifications_from_unflagged_records":1,
                   "criteria_provided,_single_submitter":1,
                   "no_assertion_criteria_provided":1,
                   "no_classification_provided":1,
                   "no_classification_for_the_single_variant":0,
                   "flagged_submission":0, # ??? Conflicting information is what the really means...
                   ".":0} #Means that there is no data submitted for germline classification"

var2disease_star_min = 3

# Manually curated to help determine predicate
pathogenic_enums = ["Likely pathogenic", 
                    "Likely pathogenic, low penetrance", 
                    "Pathogenic", "Pathogenic, low penetrance",
                    "Pathogenic/Likely pathogenic"]

pathogenic_lookup = {k:'' for k in pathogenic_enums}

# File paths to acessory data
sub_path = "./data/submission_summary.txt.gz"
vcf_tsv_path = "./data/clinvar.tsv"
vcf_path = "./data/clinvar.vcf"
sssom_path = "./data/mondo.sssom.tsv"

# Map records to each clinvar variant id
var_records = make_variant_record_map(sub_path)
print("- Var records read in {}".format(format(len(var_records), ',')))

# Make general map back to mondo terms for things that are in our vcf / submission_summary files
map_to_mondo = make_mondo_map(sssom_path)
print("- mondo sssom read in {}".format(format(len(map_to_mondo), ',')))

# Koza app loop through each line of the file as a dictionary 
no_record = 0
with_record = 0
var_to_diss = 0
dis_hp_counts = {i:0 for i in range(0, 50)}
map_stats = {"MONDO":0, "mesh":0, "OMIM":0, "Orphanet":0}

vars_added = 0
var2gene_added = 0
var2dis_added = 0
var2hp_added = 0

tot_count = 0
while (row := koza_app.get_row()) is not None:
    # Code to transform each row of data
    # For more information, see https://koza.monarchinitiative.org/Ingests/transform

    tot_count += 1
    if tot_count % 100000 == 0:
        print("- Processed {}".format(format(tot_count, ',')))

    # Graph level objects we need koza to write
    entities = []

    # Values we need pull
    varid = str(row["ID"])
    clinical_significance = row["CLNSIG"]
    crev = row["CLNREVSTAT"]
    ginfo = row["GENEINFO"]
    raw_diss_info = row["CLNDISDB"]
    
    # Initial acmg filtering here
    # This is also performed from variant_records_to_disease function so it is left out
    ###stars = review_star_map[crev]
    ###if stars < var2disease_star_min:
    ###    continue
    
    if varid not in var_records:
        no_record += 1
        continue
    else:
        with_record += 1
    
    # Make SequenceVariant (must first find genes that are associated with it to pass in)
    gene_ids, gene_symbols = make_genes_from_row(ginfo)

#     # This should do four steps in one fell swoop
#     # Disease --> Disease + HP pipeline
#     # Finds uniq mondo_ids, and linked hp terms (if possible)
#     # Filters for records < star_min (ACMG review status "exper panel and practice guidlines")
#     mondo_to_hp, d_preds, org_preds = map_records_to_disease_hp(record_list=var_records[varid],
#                                                                 clndisdb_column=raw_diss_info, 
#                                                                 review_star_map=review_star_map, 
#                                                                 map_to_mondo=map_to_mondo,
#                                                                 pathogenicity_lookup=pathogenic_lookup,
#                                                                 star_min=3)
    
    # Variant records --> MONDO:ID based on review status of ACMG. Ratings < star_min will not have an association
    disease_ids, disease_predicates, org_predicates  = variant_records_to_disease(var_records[varid], 
                                                                                  review_star_map, 
                                                                                  map_to_mondo, 
                                                                                  pathogenic_lookup,
                                                                                  star_min=var2disease_star_min)
    
    # Pull and format disease_ids and HP terms
    diss_info = parse_CLNDISDB(raw_diss_info)
    
    # Map inormation to mondo id if possible and discard otherwise
    diss_info, map_stats = map_CLNDISDB_to_mondo(diss_info, map_to_mondo, map_stats)
    
    # Map each disease to HP terms if possible
    # {MONDO:123:[HP:123, HP:234, ...], ...}
    mondo_to_hp = map_mondo_to_hp(diss_info, disease_ids) #--> {disease_id:[HP:123, HP:234, ...], }

    # TO DO: Fill out more?
    # Start creating graph data starting with the variant itself
    seq_var = SequenceVariant(
        id="CLINVAR:{}".format(row["ID"]),
        name=row["CLNHGVS"],
        xref=["DBSNP:{}".format(row["RS"])],
        has_gene=gene_ids,
        in_taxon=["NCBITaxon:9606"],
        in_taxon_label="Homo sapiens")
        # type? could this be a SO term?
        # has_bioligical_sequence  do we want it? not so sure
# #     )

    entities.append(seq_var)
    vars_added += 1
    
    # TO DO: Fill out more? MC (molecular consequence value(s)) How to ascribe multiple terms with the gene info...
    # Make Gene Associations
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
        var2gene_added += 1
    

    # TO DO: Predicate when and what still needs work 
    # Create VariantToDiseaseAssociations 
    for dis_id, predicate in disease_predicates.items():
        
        # We can curate predicate from anywhere here (original record file, or what the vcf provides)
        # Can be able to process multiple predicates (or we can remove this)
        #for og_pred in o/g_preds:
        #    if og_pred in pathogenic_lookup:
        #        predicate = CAUSES
        #        negated = False
        #    else:
        #        predicate = CONTRIBUTES_TO
        #        negated = True
        
        ### predicate is a dictionary of possible predicate values (in case a variant has multiple status's? (not sure possible...))
        for pred in list(predicate.keys()):
            entities.append(
                VariantToDiseaseAssociation(
                    id=str(uuid.uuid4()),
                    subject=seq_var.id,
                    predicate=pred,
                    qualifiers=[row["CLNREVSTAT"]],
                    object=dis_id,
                    original_predicate=row["CLNSIG"], # This can also be pulled from submission_summary records...,
                    primary_knowledge_source="infores:clinvar",
                    aggregator_knowledge_source=["infores:monarchinitiative"],
                    knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
                    agent_type=AgentTypeEnum.manual_agent,
                )
            )
            var2dis_added += 1
    

    # TO DO: Predicate when and what still needs work
    # Create Variant to HP assocations (Currently dependent on an existing VariantToDisease association)
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
                    knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
                    agent_type=AgentTypeEnum.manual_agent,
                )
            )
            var2hp_added += 1
    
    # Stats here
    #dis_counts[len(hp_terms)] += 1
    
    # Test case with 11 hp terms
    # for dis, hp_terms in mondo_to_hp.items():
    #     if len(hp_terms) == 11:
    #         print(dis, hp_terms)
    #         for rr in var_records[varid]:
    #             print(rr)
    
    # Record keeping
    if len(mondo_to_hp) > 0:
        var_to_diss += len(mondo_to_hp)
    
    for m_id, hps in mondo_to_hp.items():
        dis_hp_counts[len(hps)] += 1
    

    koza_app.write(*entities)


# This won't actually will be printed, but would be nice to keep track of in a stats file or something   
print("- No records {}".format(no_record))
print("- With records {}".format(with_record))
print("- VariantsToDisease {}".format(format(var_to_diss, ',')))
print("")
print("-Variants added {}".format(format(vars_added, ',')))
print("-Variant to gene edges {}".format(format(var2gene_added, ',')))
print("-Variant to disease edges {}".format(format(var2dis_added, ',')))
print("-Variant to hp edges {}".format(format(var2hp_added, ',')))
for k, v in map_stats.items():
    print("- {} --> mondo_id , {}".format(k, format(v, ',')))
    
for k, v in dis_hp_counts.items():
    print("- HP Associations {}, Number of cases {}".format(k, v))
