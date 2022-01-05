#!/bin/bash

# Mount HCoV-19-Genomics at /hcov-genomics/
# Mount working directory at /wd/
# Note: For run_persistence.nf mount /wd/results/beast-analysis at /wd
# Mount downsample_seqs code repo at /code/

docker run -v /Users/karthik/Documents/code/outbreak_db/hCoV/HCoV-19-Genomics:/hcov-genomics -v /Users/karthik/Documents/code/outbreak_db/hCoV/san_diego/timetree_attempts/downsample_seqs:/code -v /Users/karthik/Documents/code/outbreak_db/hCoV/san_diego/timetree_attempts/2021-12-22:/wd -it ubuntu:seq_utils /bin/bash