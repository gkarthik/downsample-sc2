Channel.fromPath(params.manual_exclude_ids+"/*.txt")
    .map {it -> [it.name.replace("_excluded.txt", ""), it]}
    .set {exclude_ids_ch}

Channel.fromPath(params.clock_rate_filter_dir+"/*.fa.gz")
    .map {it -> [it.name.replace("_aligned_masked_cds_filtered.fa.gz", ""), it]}
    .set {clock_rate_filter_ch}

ids_algn_ch = exclude_ids_ch.join(clock_rate_filter_ch)

process manual_exclude {
    input:
        tuple val(out_name), file(exclude_ids), file(clock_rate_filtered_algn) from ids_algn_ch
    output:
        tuple out_name, file("${out_name}_aligned_masked_cds_filtered_manual_exclude.fa.gz") into manual_exclude_ch
    publishDir params.manual_exclude_fasta_dir, mode: 'copy'
    script:
        """
        exclude_selected_seqs.py $exclude_ids $clock_rate_filtered_algn ${out_name}_aligned_masked_cds_filtered_manual_exclude.fa.gz
        """
    stub:
        """
        touch ${out_name}_aligned_masked_cds_filtered_manual_exclude.fa.gz
        """
}

// Optional. Required only sometimes. Ideall would just do iterative filter on first clock rate filter
process generate_ml_tree_after_manual_exclude {
    input:
        tuple out_name, file(cds_fa) from manual_exclude_ch
    output:
        tuple out_name, file(cds_fa), file("${out_name}_ml_manual_exclude.treefile") into ml_after_exclude_ch
    publishDir params.ml_manual_exclude, mode: 'copy'
    label 'ml'
    script:
        """
        iqtree2 -s $cds_fa -T 6 -m HKY -czb -fast --prefix ${out_name}_ml_manual_exclude
        """
    stub:
        """
        touch ${out_name}_ml_manual_exclude.treefile
        """
}

process rerun_clock_rate_filter {
    input:
        tuple out_name, file(cds_fa), file(ml_tree) from ml_after_exclude_ch
    output:
        tuple out_name, file("${out_name}_manual_check_filtered.fa.gz"), file("${out_name}_rtt.png") into clock_rate_filter_two_ch
    publishDir params.clock_rate_filter_two_dir, mode: 'copy'
    script:
        """
        clock_filter.py $ml_tree $cds_fa ${out_name}_manual_check_filtered.fa.gz ${out_name}_rtt.png 
        """
    stub:
        """
        touch ${out_name}_manual_check_filtered.fa.gz ${out_name}_rtt.png 
        """
}

process generate_guide_tree {
    input:
        tuple out_name, file(filtered_fasta), file(rtt_plot) from clock_rate_filter_two_ch
    output:
        tuple out_name, file("${out_name}_guide_tree.treefile"), file("${out_name}_guide_tree_rooted.treefile") into generate_guide_tree_ch
    publishDir params.final_ml_dir, mode: 'copy'
    label 'ml'
    script:
        """
        iqtree2 -s $filtered_fasta -T 6 -m HKY -czb -fast --prefix ${out_name}_guide_tree
        collapse_nz.R ${out_name}_guide_tree.treefile ${out_name}_guide_tree_rooted.treefile
        """
    stub:
        """
        touch ${out_name}_guide_tree.treefile ${out_name}_guide_tree_rooted.treefile
        """
}