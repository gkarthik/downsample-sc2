#!/usr/bin/env python3

import pandas as pd
import sys
import yaml
from unidecode import unidecode

# Create new fasta hdr names in dataframe
def rename_sequence(attrs):
    new_name = "|".join(attrs)
    new_name = new_name.replace(" ", "_")
    new_name = unidecode(new_name)
    return new_name


if __name__ == "__main__":
    metadata_file = sys.argv[1]
    lineages_yml_file = sys.argv[2]
    out_dir = sys.argv[3]
    seed = int(sys.argv[4]) if len(sys.argv) > 3 else 112358
    total_seqs = 2500

    # metadata_file = "../data/metadata_2021-10-15_23-20.tsv.gz"
    # lineages_yml_file = "../data/lineages.yml"

    lineages = [
        [["B.1.2", False]],              # name, recursive
        [["B.1.617.2", True]],
        [["B.1.1.7", True]],
        # [["P.1", True]],
        [["B.1.427", True], ["B.1.429", True]],
        [["B.1", False]],
        [["BA.1", True]],
        [["BA.2", True]],
        [["BA.2.12.1", True]],
        [["BA.4", True]],
        [["BA.5", True]]
    ]
    # lineages = [
    #     # [["B.1.427", True], ["B.1.429", True]] # 2021-11-02
    # ]

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
            if lin[0] == "BA.2": # Exclude BA.2.12.1
                ba_2_12_1 = next(i for i in lineages_list if i["name"] == "BA.2.12.1")
                flines = [i for i in flines if i not in ba_2_12_1["children"]]
        lineage_queries.append(flines)

    # Read in metadata
    print("Reading metadata from {} ... ".format(metadata_file))
    metadata = pd.read_csv(metadata_file, sep="\t", low_memory = False, compression="gzip")
    print("Filtering out incomplete dates ... ")
    metadata = metadata[metadata["date"].apply(lambda x: len(x.split("-")) == 3)]

    # Get reference sequence metadata
    ref_metadata = metadata[metadata["gisaid_epi_isl"] == "EPI_ISL_402124"]

    # Filter location metadata
    print("Getting location metadata ... ")
    loc_metadata = metadata[(~metadata["location"].isna()) & (metadata["location"].str.contains("Diego"))]
    loc_metadata = loc_metadata[loc_metadata["submitting_lab"] == "Andersen lab at Scripps Research"]

    # Downsample background
    print("Getting background metadata ... ")
    background_metadata = metadata[metadata["location"].apply(lambda x: False if (not pd.isna(x)) and "diego" in x.lower() else True)]

    lineage_loc_metadata_dfs = []
    lineage_names = ["_".join([j[0] for j in i]) for i in lineages]
    for lineage_query, lineage_name in zip(lineage_queries, lineage_names):
        print("\nLineage: {}\n".format(lineage_name))
        # Select lineage
        # loc_lineage_metadata = loc_metadata[loc_metadata["pango_lineage"].str.upper().isin(lineage_query)]
        # if loc_lineage_metadata.shape[0] > total_seqs/2:
        #     loc_lineage_metadata = loc_lineage_metadata.sample(n=int(total_seqs/2), random_state = seed)
        # print("Selected {} sequence from location  ... ".format(loc_lineage_metadata.shape[0]))

        # n_bg = total_seqs - loc_lineage_metadata.shape[0] # Max sequences set to 40000 by default
        n_bg = 2500
        background_lineage_metadata = background_metadata[background_metadata["pango_lineage"].str.upper().isin(lineage_query)]
        n_bg = n_bg if n_bg <= background_lineage_metadata.shape[0] else background_lineage_metadata.shape[0] # Check if n_bg > number of available sequences
        background_lineage_metadata = background_lineage_metadata.sample(n = n_bg, random_state = 112313) # Fix seed
        print("Selected {} sequences from background".format(n_bg))

        # Final sequence dataset
        # loc_lineage_metadata = pd.concat([loc_lineage_metadata, background_lineage_metadata, ref_metadata])
        loc_lineage_metadata = pd.concat([background_lineage_metadata, ref_metadata])
        print("Final metadata has {} sequences".format(loc_lineage_metadata.shape[0]))

        # Create new names
        print("Creating new names ... ")
        loc_lineage_metadata["new_name"] = loc_lineage_metadata[["gisaid_epi_isl", "date", "country", "division", "location"]].fillna("None").apply(rename_sequence, axis = 1)
        loc_lineage_metadata = loc_lineage_metadata.set_index("gisaid_epi_isl")
        loc_lineage_metadata.to_csv("{}/{}-{}-{}-downsampled_metadata.tsv.gz".format(out_dir, lineage_name, seed, total_seqs), sep="\t", compression="gzip")

    sys.exit(0)
