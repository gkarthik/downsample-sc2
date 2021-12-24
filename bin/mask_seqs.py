#!/usr/bin/env python3

import pandas as pd
from Bio import AlignIO
import sys

in_path = sys.argv[1]
out_path = sys.argv[2]
mask_vcf = sys.argv[3]
algn_mask = pd.read_csv(mask_vcf,sep="\t", comment="#", names=["region", "pos", "ref", "alt", "x", "y", "mask", "comment"])
algn = AlignIO.read(in_path, "fasta")

for i in algn_mask[algn_mask["mask"] == "mask"]["pos"].tolist():
    pos = i-1
    for rec in algn:
        rec.seq = rec.seq[:pos] + "N" + rec.seq[pos+1:]

AlignIO.write(algn, out_path, "fasta")
