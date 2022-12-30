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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import tol_colors as tc
import numpy as np
from scipy import stats

plt.rc('axes', prop_cycle=plt.cycler('color', list(tc.tol_cset('bright'))))

cset = tc.tol_cset('bright')


repo_path = sys.argv[1]
lineage_report_path = os.path.join(repo_path, "lineage_report.csv")
metadata_path = os.path.join(repo_path, "metadata.csv")
sequences_dir = os.path.join(repo_path, "consensus_sequences")
lineages_yml_file = sys.argv[2]
out_dir = sys.argv[3]
regions_mapping = sys.argv[4]
cases_path = sys.argv[5]
first_last_dates_path = sys.argv[6]
seed = int(sys.argv[7]) if len(sys.argv) > 6 else 112358

# repo_path = "../../../HCoV-19-Genomics/"
# lineage_report_path = os.path.join(repo_path, "lineage_report.csv")
# metadata_path = os.path.join(repo_path, "metadata.csv")
# sequences_dir = os.path.join(repo_path, "consensus_sequences")

# lineages_yml_file = "../2021-12-22/data/lineages.yml"
# regions_mapping = "../2021-12-22/data/SanDiegoZIP_region.csv"
# cases_path = "../../cases/cases_over_time.csv"
# seed = 112358
# out_dir = "./tmp_out/"
# first_last_dates_path = "../../cases/first_last_dates.csv"

metadata = pd.read_csv(metadata_path)
lineages = pd.read_csv(lineage_report_path)
cases = pd.read_csv(cases_path, dtype={"ziptext": str, "zipcode_zip": str})
first_last_dates = pd.read_csv(first_last_dates_path)
first_last_dates["min_date"] =  pd.to_datetime(first_last_dates["min_date"], format="%Y-%m-%d")
first_last_dates["max_date"] =  pd.to_datetime(first_last_dates["max_date"], format="%Y-%m-%d")

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

# Get cases by region
cases = pd.merge(cases, sd_regions, left_on="ziptext", right_on="zip").assign(updatedate = lambda x: pd.to_datetime(x["updatedate"], format="%Y-%m-%d"))
cases = cases.groupby(["region", "updatedate"]).sum().reset_index()

# Read in lineages
lineages = [
    # [["B.1.2", False]],              # name, recursive
    # [["B.1.617.2", True]],
    # [["B.1.1.7", True]],
    # [["P.1", True]],
    # [["B.1.427", True], ["B.1.429", True]],
    # [["B.1", False]],
    [["BA.1", True]],
    [["BA.2", True]],
    [["BA.2.12.1", True]],
    [["BA.4", True]],
    [["BA.5", True]]
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
        if lin[1]:          # If recursive
            flines.extend(lin_info["children"])
	if lin[0] == "BA.2": # Exclude BA.2.12.1
            ba_2_12_1 = next(i for i in lineages_list if i["name"] == "BA.2.12.1")
            flines = [i for i in flines if i not in ba_2_12_1["children"]]
    lineage_queries.append(flines)

def rename_sequence(attrs):
    new_name = "|".join(attrs)
    new_name = new_name.replace(" ", "_")
    new_name = unidecode(new_name)
    return new_name

def get_cumulative_cases_by_region(x):
    x = x.sort_values("updatedate")
    row = x.iloc[-1]
    row["case_count"] = x["case_count"].iloc[-1] - x["case_count"].iloc[0]
    return row

lineage_metadata_dfs = []
lineage_names = ["_".join([j[0] for j in i]) for i in lineages]
for lineage_query, lineage_name in zip(lineage_queries, lineage_names):
    print("\nLineage: {}\n".format(lineage_name))

    lineage_metadata = merged_df[merged_df["lineage"].str.upper().isin(lineage_query)].assign(collection_date = lambda x: pd.to_datetime(x["collection_date"], format="%Y-%m-%d"))
    min_max_date = first_last_dates[first_last_dates["lineage"] == lineage_name]
    if min_max_date.shape[0] == 0:
        print("For {} does not exist in first_last_dates".format(lineage_name))
    lineage_cases = cases[(cases["updatedate"] <= min_max_date["max_date"].iloc[0]) & (cases["updatedate"] >= min_max_date["min_date"].iloc[0])]
    cumulative_lineage_cases = lineage_cases.groupby("region").apply(get_cumulative_cases_by_region)
    cumulative_lineage_cases = cumulative_lineage_cases.assign(cases_prop = lambda x: x["case_count"]/x["pop"]).reset_index(drop=True)

    # Plot cases
    out_name = os.path.join(out_dir, "{}_{}_cases.png".format(lineage_name, seed))
    pivoted_lineage_cases = lineage_cases.pivot(index="updatedate", columns = "region", values="case_count")
    f,ax = plt.subplots(figsize=(7.5, 5))
    pivoted_lineage_cases.plot(ax = ax)
    plt.tight_layout()
    plt.savefig(out_name, dpi = 300)
    plt.close()

    # Sample based on infections 
    sampling = pd.merge(cumulative_lineage_cases, lineage_metadata.groupby("region").size().reset_index().rename(columns={0:"nseqs"}), on="region")
    sampling["nseqs_prop"] = sampling["cases_prop"]/sampling["cases_prop"].max()
    most_nseqs = sampling[sampling["nseqs_prop"] == 1]["nseqs"].max()
    sampling["nseqs_to_sample"] = np.ceil(sampling["nseqs_prop"] * most_nseqs)
    sampled_lineage_metadata = []
    for n,grp in lineage_metadata.groupby("region"):
        sample_n = int(sampling[sampling["region"] == n]["nseqs_to_sample"].iloc[0])
        region_lineage_df = grp
        if sample_n < grp.shape[0]:
            region_lineage_df = grp.sample(n=sample_n, random_state = seed)
        sampled_lineage_metadata.append(region_lineage_df)
    
    sampled_lineage_metadata = pd.concat(sampled_lineage_metadata)

    if sampled_lineage_metadata.shape[0] > 2500:
        sampled_lineage_metadata = sampled_lineage_metadata.sample(n=2500, random_state = seed)
    
    final_sampling = pd.merge(cumulative_lineage_cases, sampled_lineage_metadata.groupby("region").size().reset_index().rename(columns = {0: "nseqs"}), on="region").sort_values("cases_prop")
    res = stats.linregress(final_sampling["cases_prop"], final_sampling["nseqs"])
    final_sampling["predicted"] = final_sampling["cases_prop"] * res.slope + res.intercept
    final_sampling["residuals"] = final_sampling["predicted"] - final_sampling["nseqs"]
    final_sampling[["region", "updatedate", "case_count", "pop", "cases_prop", "nseqs", "predicted", "residuals"]].to_csv(os.path.join(out_dir, "{}_{}_final_sampling.csv".format(lineage_name, seed)))

    f,ax = plt.subplots(figsize=(7.5, 5))
    ax.scatter(final_sampling["cases_prop"], final_sampling["nseqs"], color = cset[:6])
    ax.plot(final_sampling["cases_prop"], res.intercept + res.slope*final_sampling["cases_prop"], 'k', ls="--")
    #ax.text(final_sampling["cases_prop"], final_sampling["nseqs"], final_sampling["region"])
    ax.set_title("{}_{} ".format(lineage_name, seed)+ f"R-squared: {res.rvalue**2:.6f}" + " nseqs: {}".format(sampled_lineage_metadata.shape[0]))

    handles = [mpatches.Patch(color=c, label=r) for r,c in zip(final_sampling["region"].tolist(), cset[:6])]
    ax.legend(handles=handles)
    plt.tight_layout()
    #out_name = os.path.join(out_dir, "{}_{}_sampling_vs_cases.png".format(lineage_name, seed))
    plt.show()
    #plt.savefig(out_name, dpi = 300)
    plt.close()
    
    # Sampling randomly
    # if lineage_metadata.shape[0] > 2500:
    #     lineage_metadata = lineage_metadata.sample(n=2500, random_state = seed)
    print("Selected {} sequence from location  ... ".format(sampled_lineage_metadata.shape[0]))

    sampled_lineage_metadata["country"] = "USA"
    sampled_lineage_metadata["division"] = "California"
    sampled_lineage_metadata["location"] = sampled_lineage_metadata[["zip", "region"]].apply(lambda x: "San_Diego_County_zip-{}_region-{}".format(x["zip"], x["region"]), axis = 1)
    
    sampled_lineage_metadata["new_name"] = sampled_lineage_metadata[["mapping_id", "collection_date", "country", "division", "location"]].fillna("None").assign(collection_date = lambda x: x["collection_date"].dt.strftime("%Y-%m-%d")).apply(rename_sequence, axis = 1)
    sampled_lineage_metadata = sampled_lineage_metadata.set_index("mapping_id")

    #lineage_metadata.to_csv("../2021.11.27_results/sd/{}_{}_{}_downsampled_metadata.tsv.gz".format(lineage_name, seed, 2500), sep="\t", compression="gzip")
    lineage_metadata_dfs.append(sampled_lineage_metadata)

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
