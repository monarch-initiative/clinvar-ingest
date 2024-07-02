import uuid  # For generating UUIDs for associations
import gzip
import pandas as pd
#import matplotlib.pyplot as plt
from collections import Counter

from biolink_model.datamodel.pydanticmodel_v2 import *  # Replace * with any necessary data classes from the Biolink Model
from koza.cli_utils import get_koza_app


### Create variant--> record mapping {variant_id:[record,...]}
def make_variant_record_map(submission_path):
    
    # Last line with "#" character in it will be read as the header
    var_records = {}
    rec_count = 0
    with gzip.open(submission_path, 'rt') as infile:
        
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
    
    ##print("- {} variants read into memory across {} records".format(format(len(var_records), ','), 
    ##                                                                format(rec_count, ',')))
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

    ##print("- Total mappings produced {}".format(len(map_to_mondo)))
    ##print("- Multi mapping objects found {}".format(format(dups, ',')))
    
    # Lastly, add in bannana and mondo_id mapping to self
    mondo_set = {k.split(":")[-1]:'' for kv in map_to_mondo for k in map_to_mondo[kv]}
    for k in mondo_set:
        v = "MONDO:{}".format(k)
        kv = "MONDO:MONDO:{}".format(k)
        
        map_to_mondo.update({v:{v:''}})
        map_to_mondo.update({kv:{v:''}})
    
    return map_to_mondo


def make_medgen_to_mondo_map(medgen_path):
    
    # Last line with "#" character in it will be read as the header
    map_to_mondo = {}
    rec_count = 0
    with gzip.open(medgen_path, 'rt') as infile:
        
        for line in infile:
            line = line.strip('\r').strip('\n')
            if line[0] == "#":
                continue
            else:
                info = line.split("|")
                dis_id = info[2]
                mdg_id = "MedGen:{}".format(info[0])
                
                if "MONDO:" in dis_id:
                    if mdg_id not in map_to_mondo:
                        map_to_mondo.update({mdg_id:{}})
                    
                    # Will account for multimapping terms if there are any
                    map_to_mondo[mdg_id].update({dis_id:''})
    
    multi_mappers = len([k for k, v in map_to_mondo.items() if len(v) > 1])
    ##print("- Total mappings found {}".format(format(len(map_to_mondo), ',')))
    ##print("- MedGen ids that map to multiple mondo {}".format(format(multi_mappers, ',')))
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

    elif "HP:" in info:
        idname = "HP:{}".format(idnum)

    elif "MeSH:" in info:
        idname = "mesh:{}".format(idnum)
        
    elif "." == info:
        idname = None
        
    else:
        idname = info
    
    return idname
  
 
def variant_records_to_disease(record_list, review_star_map, map_to_mondo, predicate_map, star_min=3):
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
        clinsig = rec["ClinicalSignificance"]
        
        if stars < star_min:
            continue
        
        if clinsig not in predicate_map:
            continue
            
        mapped_predicate = predicate_map[clinsig]
        org_predicate = clinsig
        
        #if len(rec["SubmittedPhenotypeInfo"].split(';')) != len(rec["ReportedPhenotypeInfo"].split(';')):
        #    print("- ERROR, Submitted vs Reported PhenotypeInfo is off kilter...")
        #    print(rec["SubmittedPhenotypeInfo"], rec["ReportedPhenotypeInfo"])
        
        # Submitted and Reported "phenotype info" columns can have a different number of submissions.
        # So we need to loop through these terms one by one first for the reported terms, and then if no
        # mondo id is found, then we try and pull from the supported terms
        
        mapped_terms = 0
        for mg_mapping in rec["ReportedPhenotypeInfo"].split(';'):
            
            # Convert to MONDO_ID if we can
            mg_map = "MedGen:{}".format(mg_mapping.split(":")[0]) # number: name and description of disease
            mondo_ids = []
            
            # Use reported medgen number from clinver First if we can...
            if mg_map in map_to_mondo:
                mondo_ids = list(map_to_mondo[mg_map].keys())
            
            # Can't map this one back
            if len(mondo_ids) == 0:
                continue
            
            # Create disease id to predicate mapping (for ingest and original)
            for d in mondo_ids:
                dis.update({d:''})
                if d not in preds:
                    preds.update({d:{}})
                    org_preds.update({d:{}})
                preds[d].update({mapped_predicate:''})
                org_preds[d].update({org_predicate:''})
                mapped_terms += 1
    
        # Try to query the Submitted info for mondo id if none is found for Reported info
        if mapped_terms == 0:
            for dis_id in rec["SubmittedPhenotypeInfo"].split(';'):
            
                # Convert to MONDO_ID if we can
                dis_id = format_id_to_map(dis_id)

                # Otherwise we default to querying the sssom for a mondo map
                if dis_id in map_to_mondo:
                    mondo_ids = list(map_to_mondo[dis_id].keys())

                elif "MONDO:" in dis_id:
                    mondo_ids = [dis_id]

                # Can't map this one back
                if len(mondo_ids) == 0:
                    continue
                
                # Create disease id to predicate mapping (for ingest and original)
                for d in mondo_ids:
                    dis.update({d:''})
                    if d not in preds:
                        preds.update({d:{}})
                        org_preds.update({d:{}})
                    preds[d].update({mapped_predicate:''})
                    org_preds[d].update({org_predicate:''})
                    mapped_terms += 1
                    ##print(rec["SubmittedPhenotypeInfo"], rec["ReportedPhenotypeInfo"])
    
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
                
                ### STATS KEEPING for different submission ids
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
        for d in set(g["MAP_TERMS"]):
            if d in disease_ids:
                
                # Check if exists or not currently. We can pull in multiple hp terms from records that
                # have the same disease id as the "cononical" disease for this set of records
                if d not in mondo_to_hp:
                    mondo_to_hp.update({d:[]})
                    
                mondo_to_hp[d] += g["HP"]

    return mondo_to_hp


###############################
### Ingest pipeline / steps ###

### First, create a variant to record(s) map {variantid:[{},{},{}...]}

### Second, create a disease id to mondo id map in the form of... map_to_mondo --> {"OMIM:123":"MONDO:123", ...}
### where the key (disease id) is from the sources found within the sssom and medgen files
### The map_to_mondo is able to handle a few different types of data
### Orphanet, OMIM, Mondo, MedGen, and some MeSH for this data set 

### Third, loop through each variant within the vcf...
###  - Pull all variant records
###  - Pull out records that pass review star minimum criteria 
###  - Map relevant record(s) medGen disease id(s) back to mondo id(s)
###  - Match reported vcf disease/hpo group(s) back to relevant records via disease id
###  - Make SequenceVariant node, variant-->gene edges, variant-->disease edges, variant-->hpo edges
###  - Write biolink certified nodes/edges via koza app

koza_app = get_koza_app("clinvar_variant")

# Variant to gene predicate
IS_SEQUENCE_VARIANT_OF = "biolink:is_sequence_variant_of"

# Variant to disease 
CAUSES = "biolink:causes"
ASSOCIATED_WITH = "biolink:associated_with_increased_likelihood_of"

# Variant to phenotype
CONTRIBUTES_TO = "biolink:contributes_to" 

pred_to_negated = {CAUSES:False,
                   ASSOCIATED_WITH:False,
                   CONTRIBUTES_TO:False}


# Manually curated terms derived from files for data modeling purposes
# https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
review_star_map = {"practice_guideline":4,
                   "reviewed_by_expert_panel":3,
                   "criteria_provided,_multiple_submitters,_no_conflicts":2,
                   "criteria_provided,_conflicting_classifications":1,
                   "no_classifications_from_unflagged_records":1,
                   "criteria_provided,_single_submitter":1,
                   "no_assertion_criteria_provided":1,
                   "no_classification_provided":1,
                   "no_classification_for_the_single_variant":0,
                   "flagged_submission":0, # ??? Conflicting information is what the really means...
                   ".":0} #Means that there is no data submitted for germline classification"

var2disease_star_min = 3 ### 3 is reviewed by expert panel and above

# Manually curated to help determine predicate
# For using submission_summary file (pulling terms from the clinicial significance column)
predicate_map = {"Pathogenic":CAUSES,
                 "Pathogenic, low penetrance":CAUSES,
                 "Pathogenic/Likely pathogenic":CAUSES,
                
                 "Likely pathogenic":ASSOCIATED_WITH,
                 "Likely pathogenic, low penetrance":ASSOCIATED_WITH}

                #  # If we want to include more terms we can do that here...
                #  "Benign":RELATED_TO,
                #  "Benign/Likely benign":RELATED_TO,
                #  "Likely benign":RELATED_TO,

                #  "Affects":RELATED_TO,
                #  "Established risk allele":RELATED_TO,
                #  "Likely risk allele":RELATED_TO,
                #  "Uncertain risk allele":RELATED_TO,
                #  "Uncertain significance":RELATED_TO,
                #  "association":RELATED_TO,
                #  "association not found":RELATED_TO, # use negate
                #  "confers sensitivity":RELATED_TO,
                #  "conflicting data from submitters":RELATED_TO,
                #  "drug response":RELATED_TO, # Can also probably change this for target
                #  "not provided":RELATED_TO,
                #  "protective":RELATED_TO,
                #  "risk factor":RELATED_TO} # PREDESPOSES_TO_CONDITION (Can alter predicate here)

## Allows us to convert to hgnc (before merge step. Only used for develpoment purposes)
##hgnc_df = pd.read_csv("./data/hgnc_complete_set.txt", sep='\t', header=0)
##symbol_to_hgnc = {sym:hgnc_id for sym, hgnc_id in zip(list(hgnc_df["symbol"]), list(hgnc_df["hgnc_id"]))}
##name_to_hgnc = {name:hgnc_id for name, hgnc_id in zip(list(hgnc_df["name"]), list(hgnc_df["hgnc_id"]))}

##print("- Symbol to hgnc map {}".format(format(len(symbol_to_hgnc), ',')))
##print("- Name to hgnc map {}".format(format(len(name_to_hgnc), ',')))


# File paths to acessory data
sub_path = "./data/submission_summary.txt.gz"
sssom_path = "./data/mondo.sssom.tsv"
medgen_path = "./data/MedGenIDMappings.txt.gz"

# Map records to each clinvar variant id
var_records = make_variant_record_map(sub_path)

# Make general map back to mondo terms for things that are in our vcf / submission_summary files
map_to_mondo = make_mondo_map(sssom_path)

# Make medgen to mondo map
medgen_to_mondo = make_medgen_to_mondo_map(medgen_path)

# Merge the two maps we made into one by updatating one dictionary with the other
map_to_mondo.update(medgen_to_mondo)

# Record keeping variables 
no_record = 0
with_record = 0
map_stats = {"MONDO":0, "mesh":0, "OMIM":0, "Orphanet":0, "MedGen":0}

vars_added = 0
var2gene_added = 0
var2dis_added = 0
var2hp_added = 0
tot_count = 0

while (row := koza_app.get_row()) is not None:
    # Code to transform each row of data
    # For more information, see https://koza.monarchinitiative.org/Ingests/transform

    # Progress tracking
    ##tot_count += 1
    ##if tot_count % 100000 == 0:
    ##    print("- Processed {}".format(format(tot_count, ',')))

    # Graph level objects we need koza to write
    entities = []

    # Values we need pull
    varid = str(row["ID"])
    clinical_significance = row["CLNSIG"]
    crev = row["CLNREVSTAT"]
    ginfo = row["GENEINFO"]
    raw_diss_info = row["CLNDISDB"]
    so_info = [v.split("|")[0] for v in row["MC"].split(",") if "SO:" in v] # Pull out sequence ontology term(s) ### Example MC column SO:0001575|splice_donor_variant,SO:0001587|nonsense
    ### Note, that the Sequence ontology term could be derived from the "CLNVCSO" vcf column as well, however the terms listed in that column are much less specific and are not actually particularly useful (Too broad)
    ### The terms listed within the MC column are far more specific to the effect(s) a variant will have on any given gene it overlaps, thus making it the preffered choice. 

    # No record info means we don't want to include
    if varid not in var_records:
        no_record += 1
        continue
    else:
        with_record += 1
    
    # Make SequenceVariant (must first find genes that are associated with it to pass in)
    gene_ids, gene_symbols = make_genes_from_row(ginfo)

    # Variant records --> MONDO:ID based on review status of ACMG. Ratings < star_min will not have an association
    disease_ids, disease_predicates, org_predicates  = variant_records_to_disease(var_records[varid], 
                                                                                  review_star_map, 
                                                                                  map_to_mondo, 
                                                                                  predicate_map,
                                                                                  star_min=var2disease_star_min)
    
    # Pull and format disease_ids and HP terms
    diss_info = parse_CLNDISDB(raw_diss_info)
    
    # Map inormation to mondo id if possible and discard otherwise
    diss_info, map_stats = map_CLNDISDB_to_mondo(diss_info, map_to_mondo, map_stats)
    
    # Map each disease to HP terms if possible
    # {MONDO:123:[HP:123, HP:234, ...], ...}
    mondo_to_hp = map_mondo_to_hp(diss_info, disease_ids) #--> {disease_id:[HP:123, HP:234, ...], }

    # This means we were not able to make a variant --> disease association. 
    # Therefore we do not want to add any information to the graph
    if len(mondo_to_hp) == 0:
        continue

    # Start creating graph data starting with the variant itself
    seq_var = SequenceVariant(
                    id="CLINVAR:{}".format(row["ID"]),
                    name=row["CLNHGVS"],
                    xref=["DBSNP:{}".format(row["RS"])],
                    has_gene=gene_ids,
                    in_taxon=["NCBITaxon:9606"],
                    in_taxon_label="Homo sapiens",
                    type=so_info)

    entities.append(seq_var)
    vars_added += 1
    
    # # Make Gene Associations (If we want to pre-convert ncbi geneIds to hgnc geneIds... )
    # # This is done at the merge step so not necessary here, but a good initial sanity check to ensure majority of genes are being converted
    # for gene_id, gene_symbol in zip(gene_ids, gene_symbols):
    #     if gene_symbol in symbol_to_hgnc:
    #         gene_id = symbol_to_hgnc[gene_symbol]
    #     elif gene_symbol in name_to_hgnc:
    #         gene_id = name_to_hgnc[gene_symbol]

    # Make variant to gene associations
    for gene_id in gene_ids:
        entities.append(
            VariantToGeneAssociation(
                id=str(uuid.uuid4()),
                subject=seq_var.id,
                predicate=IS_SEQUENCE_VARIANT_OF, # Can't be more specific than this for variants with multiple gene overlaps / annotations
                object=gene_id,
                primary_knowledge_source="infores:clinvar",
                aggregator_knowledge_source=["infores:monarchinitiative"],
                knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
                agent_type=AgentTypeEnum.manual_agent,
            )
        )
        var2gene_added += 1
    
    # Make variant to disease associations
    for dis_id, predicate in disease_predicates.items():
        
        ### predicate is a dictionary of possible predicate values (in case a variant has multiple status's? (not sure possible...))
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
                    original_predicate=":".join(og_preds), # Pulled from the submission_summary file
                    primary_knowledge_source="infores:clinvar",
                    aggregator_knowledge_source=["infores:monarchinitiative"],
                    knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
                    agent_type=AgentTypeEnum.manual_agent,
                )
            )
            var2dis_added += 1
    
    # Make variant to HPO assocations (Currently dependent on an existing VariantToDisease association)
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
            var2hp_added += 1
    
    koza_app.write(*entities)