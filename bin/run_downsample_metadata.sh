#!/bin/bash

mkdir -p /wd/results/bg /wd/results/sd

echo -e "112358\n112313\n112321" > seeds.txt

ls seeds.txt | xargs -n 1 -P 3 bash -c 'python downsample_metadata.py /wd/data/metadata_2021-12-10_16-25.tsv.gz /wd/data/lineages.yml /wd/results/bg $0'

# Downsample san diego sequences. Not parallelized since small number of sequences
ls seeds.txt | xargs -n 1 -P 3 bash -c 'python downsample_sd_data.py /hcov-genomics/ /wd/data/lineages.yml /wd/results/sd /wd/data/SanDiegoZIP_region.csv $0'