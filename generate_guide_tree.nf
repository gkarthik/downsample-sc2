process downsample_metadata {
    input: 
        path gisaid_metadata from params.gisaid_metadata
        path lineages_yml from params.lineages_yml
        path out_dir from params.out_dir_downsample_metadata_bg
        val seed from params.seeds
    publishDir params.out_dir_downsample_metadata_bg
    output:
        file "*.tsv.gz" into downsampled_metadata_ch
    script:
        """
        downsample_metadata.py $gisaid_metadata $lineages_yml . $seed
        """
    stub:
        """
        touch B.1-${seed}-2500-downsampled_metadata.tsv.gz B.1.1.7-${seed}-2500-downsampled_metadata.tsv.gz B.1.2-${seed}-2500-downsampled_metadata.tsv.gz B.1.427_B.1.429-${seed}-2500-downsampled_metadata.tsv.gz B.1.617.2-${seed}-2500-downsampled_metadata.tsv.gz P.1-${seed}-2500-downsampled_metadata.tsv.gz
        """
}

process downsample_sd {
    input: 
        path seq_repo from params.seq_repo
        path lineages_yml from params.lineages_yml
        path out_dir from params.out_dir_downsample_metadata_sd
        path regions_mapping from params.regions_mapping
        val seed from params.seeds
    publishDir params.out_dir_downsample_metadata_sd
    output:
        file "*.fa.gz" into downsampled_sd_ch
    script:
        """
        downsample_sd_data.py $seq_repo $lineages_yml . $regions_mapping $seed
        """
    stub:
        """
        touch downsampled-B.1-${seed}.fa.gz downsampled-B.1.1.7-${seed}.fa.gz downsampled-B.1.2-${seed}.fa.gz downsampled-B.1.427_B.1.429-${seed}.fa.gz downsampled-B.1.617.2-${seed}.fa.gz downsampled-P.1-${seed}.fa.gz
        """
}

process downsample_gisaid_fasta {
    input:
        path downsampled_metadata_list from downsampled_metadata_ch.collect()
        path tmp from params.tmp_dir
        path gisaid_fasta from params.gisaid_fasta_file
        val threads from params.gisaid_downsample_threads
    publishDir params.downsampled_fa_bg_out_dir
    output:
        file "*.fa.gz" into downsample_gisaid_fasta_ch
    script:
        """
        mkdir -p tmp fa_out
        gunzip -c $gisaid_fasta | parallel --recend '\n' --recstart '>' --block 200M -j $threads --pipe --tmpdir tmp downsample_seqs.py $downsampled_metadata_list fa_out/downsampled_{#}
        cd fa_out/
        ls | cut -f 5- -d _ | sort | uniq | xargs -n 1 bash -c 'cat *\$0* > ../\$0'
        """
    stub:
        """
        echo $downsampled_metadata_list
        touch B.1-112358-2500.fa.gz B.1.1.7-112358-2500.fa.gz B.1.2-112358-2500.fa.gz B.1.427_B.1.429-112358-2500.fa.gz B.1.617.2-112358-2500.fa.gz P.1-112358-2500.fa.gz
        touch B.1-112313-2500.fa.gz B.1.1.7-112313-2500.fa.gz B.1.2-112313-2500.fa.gz B.1.427_B.1.429-112313-2500.fa.gz B.1.617.2-112313-2500.fa.gz P.1-112313-2500.fa.gz
        touch B.1-112321-2500.fa.gz B.1.1.7-112321-2500.fa.gz B.1.2-112321-2500.fa.gz B.1.427_B.1.429-112321-2500.fa.gz B.1.617.2-112321-2500.fa.gz P.1-112321-2500.fa.gz
        """
}

// File format for background metadata: <lineage_name>-<seed>-<total_seqs>.fa.gz
downsample_gisaid_fasta_ch
    .flatten()
    .map {it -> [it.name.split("-")[0] + "-" + it.name.split("-")[1], it]}
    .set{ key_fa_bg }

// File format for SD downsampeld fasta: downsampled-<lineage_name>-<seed>.fa.gz
downsampled_sd_ch
    .flatten()
    .map {it -> [it.name.split("-")[1] + "-" + it.name.split("-")[2].replace(".fa.gz", ""), it]}
    .set{ key_fa_sd }

key_fa_bg.join(key_fa_sd).set{key_fa_combined_ch}

process combine_fasta {
    input:
        tuple val(out_name), file(bg), file(sd) from key_fa_combined_ch
    output:
        tuple val(out_name), file("${out_name}.fa.gz") into combined_fa_ch
    publishDir params.combined_fasta_dir
    script:
        """
        cat $bg $sd > ${out_name}.fa.gz
        """
    stub:
        """
        touch ${out_name}.fa.gz
        """
}

process align_seqs {
    input:
        tuple val(out_name), file(fasta) from combined_fa_ch
    output:
        tuple out_name, file("${out_name}_aligned.fa") into aligned_fa_ch
    publishDir params.aligned_fasta_dir
    script:
        """
        mafft --auto --thread 6 --keeplength --addfragments $fasta $ref_genome > ${out_name}_aligned.fa
        """
    stub:
        """
        touch ${out_name}_aligned.fa
        """
}

process mask_seqs {
    input:
        tuple val(out_name), file(aligned_fasta) from aligned_fa_ch
        path mask_vcf from params.mask_vcf
    output:
        tuple out_name, file("${out_name}_aligned_masked.fa") into masked_fa_ch
    publishDir params.masked_fasta_dir
    script:
        """
        mask_seqs.py $aligned_fasta ${out_name}_aligned_masked.fa $mask_vcf
        """
    stub:
        """
        touch ${out_name}_aligned_masked.fa
        """
}

process trim_cds_seqs {
    input:
        tuple val(out_name), file(masked_fasta) from masked_fa_ch
    output:
        tuple out_name, file("${out_name}_aligned_masked_cds.fa.gz") into trim_cds_seqs_ch
    publishDir params.trim_cds_dir
    script:
        """
        trim_cds.py $masked_fasta ${out_name}_aligned_masked_cds.fa
        gzip ${out_name}_aligned_masked_cds.fa
        """
    stub:
        """
        touch ${out_name}_aligned_masked_cds.fa.gz
        """
}

process generate_ml_tree {
    input:
        tuple out_name, file(cds_fa) from trim_cds_seqs_ch
    output:
        tuple out_name, file(cds_fa), file("${out_name}_ml.treefile") into ml_tree_ch
    publishDir params.ml_dir
    script:
        """
        iqtree2 -s $cds_fa -T 6 -m HKY -czb -fast --prefix ${out_name}_ml
        """
    stub:
        """
        touch ${out_name}_ml.treefile
        """
}

process clock_rate_filter {
    input:
        tuple out_name, file(cds_fa), file(ml_tree) from ml_tree_ch
    output:
        tuple out_name, file("${out_name}_aligned_masked_cds_filtered.fa.gz"), file("${out_name}_rtt.png") into clock_rate_filter_ch
    publishDir params.clock_rate_filter_dir
    script:
        """
        clock_filter.py $ml_tree $cds_fa ${out_name}_aligned_masked_cds_filtered.fa.gz ${out_name}_rtt.png 
        """
    stub:
        """
        touch ${out_name}_aligned_masked_cds_filtered.fa.gz ${out_name}_rtt.png 
        """
}

process regenerate_ml_tree {
    input:
        tuple out_name, file(filtered_fasta), file(rtt_plot) from clock_rate_filter_ch
    output:
        tuple out_name, file("${out_name}_ml_filtered.treefile"), file("${out_name}_ml_filtered_rooted.treefile") into ml_filtered_tree_ch
    publishDir params.ml_filtered_dir
    script:
        """
        iqtree2 -s $filtered_fasta -T 6 -m HKY -czb -fast --prefix ${out_name}_ml_filtered
        collapse_nz.R ${out_name}_ml_filtered.treefile ${out_name}_ml_filtered_rooted.treefile
        """
    stub:
        """
        touch ${out_name}_ml_filtered.treefile ${out_name}_ml_filtered_rooted.treefile
        """
}
