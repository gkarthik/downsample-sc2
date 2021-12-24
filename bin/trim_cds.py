#!/usr/bin/env python3

from Bio import AlignIO
import sys

fpath = sys.argv[1]
out_path = sys.argv[2]
# Aligned and renamed alignment file after clean_seqs.py
algn = AlignIO.read(fpath, "fasta")
# CDS: 265 to 29673
trimmed_algn = algn[:,265:29674]

AlignIO.write(trimmed_algn, out_path, "fasta")

