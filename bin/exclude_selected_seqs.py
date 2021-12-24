#!/usr/bin/env python3

from Bio import SeqIO
import sys
import gzip

exclude_ids_path = sys.argv[1]
alignment_file = sys.argv[2]
alignment_out = sys.argv[3]

f = open(exclude_ids_path)
exclude_ids = [i.strip() for i in f.readlines()]
f.close()

filtered_seqs = []
with gzip.open(alignment_file, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id.split("|")[0] not in exclude_ids:
            filtered_seqs.append(record)


fout = gzip.open(alignment_out, "wt")
SeqIO.write(filtered_seqs, fout, 'fasta')
fout.close()
print("Wrote {} filtered sequences to {}".format(len(filtered_seqs), alignment_out))
