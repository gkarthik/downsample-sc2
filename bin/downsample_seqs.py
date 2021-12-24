#!/usr/bin/env python3

import pandas as pd
import sys
from Bio import SeqIO
import gzip
import yaml
from unidecode import unidecode
from multiprocessing import Pool
from itertools import repeat
import time
import itertools
import tqdm
import functools
import os
import copy

if __name__ == "__main__":
    out_file = sys.argv[1]
    metadata_file_paths = sys.argv[2: len(sys.argv)]

    # metadata_files_dir = "../2021.10.26_results"

    print("Reading metadata files: {}".format("\n".format(metadata_file_paths)))
    lineage_loc_metadata_dfs = [pd.read_csv(i, index_col = 0, sep="\t", compression = "gzip") for i in metadata_file_paths]
    
    # Metadata is in format: {name}_downsampled_metadata.tsv.gz
    lineage_names = [i.replace("_downsampled_metadata.tsv.gz", "") for i in metadata_file_paths]

    # Get list of accession ids
    loc_accession_ids = [i.index.tolist() for i in lineage_loc_metadata_dfs]

    # Read GISAID full fasta
    print("Reading FASTA ... ")

    seqs = SeqIO.parse(sys.stdin, "fasta")

    start_time = time.time()
    new_seqs = dict([[i, []] for i in lineage_names])
    for rec in tqdm.tqdm(seqs):
        for loc_accession_id, lineage_name, lineage_loc_metadata in zip(loc_accession_ids, lineage_names, lineage_loc_metadata_dfs):
            if rec.name in loc_accession_id:
                new_rec = copy.deepcopy(rec)
                n = lineage_loc_metadata.loc[new_rec.name, "new_name"]
                new_rec.name = n
                new_rec.id = n
                new_rec.description = n
                new_seqs[lineage_name].append(new_rec)

    print("Obtained {} sequences for final alignment".format(len(new_seqs)))
    print("Took {} seconds".format(time.time() - start_time))

    for name, new_seq in new_seqs.items():
        if len(new_seq) > 0:
            out_name = "{}-{}.fa.gz".format(out_file, name)
            print("Writing {} sequences to {}".format(len(new_seq), out_name))
            fout = gzip.open(out_name, "wt")
            SeqIO.write(new_seq, fout, 'fasta')
            fout.close()
            print("Finished writing.")
        else:
            print("Found {} sequences for {}".format(len(new_seq), name))

    print("#FINSIHED")

