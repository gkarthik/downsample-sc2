#!/bin/bash

# nextflow -C config_nf.txt run generate_guide_tree.nf -resume

nextflow -C config_nf.txt run generate_guide_tree.nf -with-report downsample_seqs_executaion_report.html -with-trace -with-timeline downsample_seqs_timeline.html -with-dag downsample_seqs_flowchart.dot 

# After manual checks run
nextflow -C config_nf.txt run generate_guide_tree_after_check.nf 

# Generate XML
nextflow -C config_nf.txt run generate_xml.nf

# Generate empirical discrete state xml
nextflow -C config_nf.txt run generate_xml_discrete.nf

