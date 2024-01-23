#!/bin/bash

# Mount HCoV-19-Genomics at /hcov-genomics/
# Mount working directory at /wd/
# Mount downsample_seqs code repo at /code/

docker run -v /kattegat/gk/hcov-19-genomics:/hcov-genomics -v /home/gk/code/downsample-sc2:/code -v /kattegat/gk/downsampling/2022-12-30_results:/wd -it ubuntu:seq_utils /bin/bash
