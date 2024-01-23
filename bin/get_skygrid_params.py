#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import gzip
import baltic as bt
import math
import sys

name = sys.argv[1] # Format: <linegae_name>-<seed>
alignment_file = sys.argv[2]

# name = "B.1-112358"
# alignment_file = "../../2021-10-18/2021.11.27_results/clock_rate_filtered/B_1_112358_aligned_masked_cds.fa.gz"

f = gzip.open(alignment_file, "rt")
algn = SeqIO.parse(f, "fasta")

data = [i.name.split("|") for i in algn]
df = pd.DataFrame(data, columns = ["accession_id", "date", "country","division", "location"])
f.close()

approx_roots = {
    "B.1": "2019-10-01",
    "B.1.2": "2019-10-01",
    "B.1.1.7": "2020-09-01",
    "B.1.617.2": "2020-10-01",
    "P.1": "2020-10-01",
    "B.1.427_B.1.429": "2020-04-01",
    "BA.1": "2021-09-01",
    "BA.2": "2021-09-01",
    "BA.2.12.1": "2021-09-01",
    "BA.4": "2021-09-01",
    "BA.5": "2021-09-01"
}

approx_root_age = approx_roots[name.split("-")[0]]

# Exlcude earliest sequence (root)
df = df[df["accession_id"] != "EPI_ISL_402124"]

# Set grid point to every 2 weeks

grid_interval = 26
cutoff = round(df["date"].apply(lambda x: bt.decimalDate(x)).max() - bt.decimalDate(approx_root_age), 2)

ngrid_points = math.ceil(cutoff*grid_interval)

print("{} {}".format(cutoff, ngrid_points))
