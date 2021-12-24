## Nextflow pipeline to downsample sequences from GISAID

* `./start_docker.sh` -> Mount docker container and associated volumes. See script for more details
* `chmod +x /code/bin/*` -> Make all external scripts executable
* Edit `bin/downsample_metadata.py` and `bin/downsample_sd_data.py` as per requirements
* Edit config_nf.txt as per requirements
* `nextflow -C config_nf.txt run generate_guide_tree.nf`

