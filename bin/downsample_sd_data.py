#!/usr/bin/env python3

import pandas as pd
import re
import yaml
from Bio import SeqIO
from unidecode import unidecode
import gzip
import tqdm
import time
import copy
import sys
import os

repo_path = sys.argv[1]
lineage_report_path = os.path.join(repo_path, "lineage_report.csv")
metadata_path = os.path.join(repo_path, "metadata.csv")
sequences_dir = os.path.join(repo_path, "consensus_sequences")
lineages_yml_file = sys.argv[2]
out_dir = sys.argv[3]
regions_mapping = sys.argv[4]
seed = int(sys.argv[5]) if len(sys.argv) > 4 else 112358

# lineage_report_path = "../../../../HCoV-19-Genomics/lineage_report.csv"
# metadata_path = "../../../../HCoV-19-Genomics/metadata.csv"
# sequences_path = "../../../../HCoV-19-Genomics/2021-11-27_seq.fasta.gz"
# lineages_yml_file = "../data/lineages.yml"
# seed = 112358

metadata = pd.read_csv(metadata_path)
lineages = pd.read_csv(lineage_report_path)

# Structure of IDS: SEARCH
def generate_mapping(id):
    mapped_id = ""
    get_id = lambda id, re_exp: re.findall(re_exp, id)[0] if len(re.findall(re_exp, id)) > 0 else id
    if "SEARCH" in id:
        mapped_id = get_id(id, "SEARCH-[0-9]+")
    elif "CA-SDCPHL" in id:
        mapped_id = get_id(id, "CA-SDCPHL-[0-9A-za-z]+")
    elif "STM" in id:
        mapped_id = get_id(id, "STM-[\-A-Za-z0-9]+")
    else:
        mapped_id = id
    return mapped_id

metadata["mapping_id"] = metadata["ID"].apply(generate_mapping)

lineages["mapping_id"] = lineages["taxon"].apply(generate_mapping)

merged_df = pd.merge(metadata, lineages, how="inner", on="mapping_id")

# Include sequences with known zipcode
merged_df["zipcode"] = merged_df["zipcode"].astype(str).apply(lambda x: x.split(".")[0])
sd_regions = pd.read_csv(regions_mapping, dtype=str)

merged_df = pd.merge(merged_df, sd_regions, left_on="zipcode", right_on="zip")

merged_df = merged_df[(~merged_df["region"].isna())]



# Read in lineages
lineages = [
    [["B.1.2", False]],              # name, recursive
    [["B.1.617.2", True]],
    [["B.1.1.7", True]],
    [["P.1", True]],
    [["B.1.427", True], ["B.1.429", True]],
    [["B.1", False]]
]

# Read in lineage
lineage_f = open(lineages_yml_file)
lineages_list = yaml.load(lineage_f, Loader=yaml.FullLoader)
lineage_f.close()

lineage_queries = []

for lins in lineages:
    flines = []
    for lin in lins:
        lin_info = next(i for i in lineages_list if i["name"] == lin[0])
        flines.append(lin[0])
        if lin[1]:          # If resursive
            flines.extend(lin_info["children"])
    lineage_queries.append(flines)


def rename_sequence(attrs):
    new_name = "|".join(attrs)
    new_name = new_name.replace(" ", "_")
    new_name = unidecode(new_name)
    return new_name

lineage_metadata_dfs = []
lineage_names = ["_".join([j[0] for j in i]) for i in lineages]
for lineage_query, lineage_name in zip(lineage_queries, lineage_names):
    print("\nLineage: {}\n".format(lineage_name))

    lineage_metadata = merged_df[merged_df["lineage"].str.upper().isin(lineage_query)]
    if lineage_metadata.shape[0] > 2500:
        lineage_metadata = lineage_metadata.sample(n=2500, random_state = seed)
    print("Selected {} sequence from location  ... ".format(lineage_metadata.shape[0]))

    lineage_metadata["country"] = "USA"
    lineage_metadata["division"] = "California"
    lineage_metadata["location"] = lineage_metadata[["zip", "region"]].apply(lambda x: "San_Diego_County_zip-{}_region-{}".format(x["zip"], x["region"]), axis = 1)
    
    lineage_metadata["new_name"] = lineage_metadata[["mapping_id", "collection_date", "country", "division", "location"]].fillna("None").apply(rename_sequence, axis = 1)
    lineage_metadata = lineage_metadata.set_index("mapping_id")

    #lineage_metadata.to_csv("../2021.11.27_results/sd/{}_{}_{}_downsampled_metadata.tsv.gz".format(lineage_name, seed, 2500), sep="\t", compression="gzip")
    lineage_metadata_dfs.append(lineage_metadata)

loc_accession_ids = [i.index.tolist() for i in lineage_metadata_dfs]

# Read in fasta
seqs = []
for i in os.listdir(sequences_dir):
    if i[-6:] != ".fasta":
        continue
    seq = SeqIO.read(os.path.join(sequences_dir, i), format="fasta")
    seqs.append(seq)

start_time = time.time()
new_seqs = dict([[i, []] for i in lineage_names])
for rec in tqdm.tqdm(seqs):
    for loc_accession_id, lineage_name, lineage_loc_metadata in zip(loc_accession_ids, lineage_names, lineage_metadata_dfs):
        seq_id = generate_mapping(rec.name)
        if seq_id in loc_accession_id:
            new_rec = copy.deepcopy(rec)
            n = lineage_loc_metadata.loc[seq_id, "new_name"]
            new_rec.name = n
            new_rec.id = n
            new_rec.description = n
            new_seqs[lineage_name].append(new_rec)

for i,j in zip(lineage_names, lineage_metadata_dfs):
    if len(new_seqs[i]) != j.shape[0]:
        print("{}: Expected {} seqs got {} seqs ".format(i, j.shape[0], len(new_seqs[i])))
        print(j.index[~j.index.isin([k.name.split("|")[0] for k in new_seqs[i]])])
        print("\n")

for name, new_seq in new_seqs.items():
        if len(new_seq) > 0:
            out_name = os.path.join(out_dir, "downsampled-{}-{}.fa.gz".format(name, seed))
            print("Writing {} sequences to {}".format(len(new_seq), out_name))
            fout = gzip.open(out_name, "wt")
            SeqIO.write(new_seq, fout, 'fasta')
            fout.close()
            print("Finished writing.")
        else:
            print("Found {} sequences for {}".format(len(new_seq), name))

sys.exit(0)