#!/usr/bin/env python3

from Bio import Phylo
from Bio import SeqIO
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import baltic as bt
from scipy.optimize import fsolve
import numpy as np
import tqdm
import os
import gzip

tree_file = sys.argv[1]
alignment_file = sys.argv[2]
# metadata_file = sys.argv[3]

#tree_file = "../2021.11.02_results/aligned_masked_downsampled_B_1_427_B_1_429_112358.fa.gz.treefile"
#alignment_file = "../2021.11.02_results/aligned_masked_downsampled_B_1_427_B_1_429_112358.fa.gz"
#metadata_file = "../2021.11.02_results/B_1_427_B_1_429_112358_downsampled_metadata.tsv.gz"

alignment_out = sys.argv[3]
#metadata_out = os.path.join(out_dir, os.path.basename(metadata_file))
rtt_plot = sys.argv[4]

print("Reading tree from {} ...".format(tree_file))
tree = Phylo.read(tree_file, "newick")
print("Reading alignment from {} ...".format(alignment_file))
alignment = pd.read_csv(alignment_file, index_col = 0, sep="\t")
print("Reading metadata from {} ...".format(alignment_file))
#metadata = pd.read_csv(metadata_file, sep="\t")

root_name = "EPI_ISL_402124|2019-12-30|China|Hubei|Wuhan"
print("Rerooting tree at {}".format(root_name))
root = next(i for i in tree.get_terminals() if i.name == root_name)
tree.root_with_outgroup(root)

print("Ladderizing ...")
tree.ladderize()

print("Getting root to tip distance ... ")
data = []
for i in tqdm.tqdm(tree.get_terminals()):
    dist = tree.distance(i)
    date = bt.decimalDate(i.name.split("|")[1])
    data.append({
        "name": i.name,
        "dist": dist,
        "date": date
    })

rtt = pd.DataFrame(data)
rtt["gisaid_epi_isl"] = rtt["name"].apply(lambda x: x.split("|")[0])

clock_rate = 0.0008

def fixed_clock_rate(i):
    return ((rtt["dist"] - (clock_rate*rtt["date"] + i))**2).sum()

fixed_clock_intercept = fsolve(fixed_clock_rate, x0=0)[0]

_x = np.linspace(rtt["date"].min(), rtt["date"].max(), 100)

rtt["residual"] = rtt.apply(lambda x: x["dist"] - (clock_rate * x["date"] + fixed_clock_intercept), axis = 1)

residual_iqd = (rtt["residual"].quantile(0.75) - rtt["residual"].quantile(0.25))
rtt["ignore"] = rtt["residual"].apply(lambda x: abs(rtt["residual"].median() - x) >= (residual_iqd*3))

# Include root regardless
rtt.loc[rtt["name"] == root_name, "ignore"] = False

f,ax = plt.subplots(figsize=(10, 7))
rtt_colors = rtt["ignore"].apply(lambda x: "red" if x else "steelblue")
ax.scatter(rtt["date"].tolist(), rtt["dist"].tolist(), alpha = 0.5, c = rtt_colors)
ax.plot(_x, clock_rate*_x + fixed_clock_intercept, color="green", ls ="--", label="Fixed clock rate")
ax.set_ylim([rtt["dist"].min() - 0.001, rtt["dist"].max() + 0.001])
plt.tight_layout()
plt.savefig(rtt_plot, dpi = 300)
plt.close()

# Write filtered metadata
#filtered_metadata = pd.merge(metadata, rtt, on="gisaid_epi_isl")
#filtered_metadata = filtered_metadata[~filtered_metadata["ignore"]]
#filtered_metadata.to_csv(metadata_out, sep="\t", compression="gzip")
#print("Filtered metadata written to {}".format(metadata_out))

# Read, Filter and write alignment
print("Reading and filtering from {} ... ".format(alignment_file))
filtered_seqs = []
#filtered_ids = filtered_metadata["gisaid_epi_isl"].tolist()
filtered_ids = rtt[~rtt["ignore"]]["gisaid_epi_isl"].tolist()
with gzip.open(alignment_file, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id.split("|")[0] in filtered_ids:
            filtered_seqs.append(record)


fout = gzip.open(alignment_out, "wt")
SeqIO.write(filtered_seqs, fout, 'fasta')
fout.close()
print("Wrote {} filtered sequences to {}".format(len(filtered_seqs), alignment_out))
