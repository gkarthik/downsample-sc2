## Nextflow pipeline to downsample sequences from GISAID

* Download unmasked alignment and genomic epidemiology metadata from GISAID
* Edit `bin/downsample_metadata.py` and `bin/downsample_sd_data.py` as per requirements
* Edit config_nf.txt as per requirements
* `./start_docker.sh` -> Mount docker container and associated volumes. See script for more details
* `chmod +x /code/bin/*` -> Make all external scripts executable
* `nextflow -C config_nf.txt run generate_guide_tree.nf` to generate guide trees
* Nextflow workflows,
  * `generate_guide_tree.nf` to generate initial guide trees
  * `generate_guide_tree_after_check.nf` to generate rooted guide trees after manual checks
  * `generate_xml.nf` to generate XML from template from rooted guide trees
  * `generate_xml_discrete.nf` to get consensus trees, downsample trees from posterior and generate empirical discrete state XML from template
  * **TODO**: Write workflow to get consensus trees after discrete state analysis
  * `run_persistence.nf` to run persistence summarizer on trees from posterior annotated with complete jump history

Note: No workflow to run BEAST since versions might differ. 