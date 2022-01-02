Channel
    .fromPath(params.trees_dir+"/*.trees")
    .map {it -> [it.name.replace(".trees", ""), it]}
    .into {trees_ch1; trees_ch2}

process downsample_trees {
    input:
        tuple val(out_name), file(tree) from trees_ch1
        val burnin from params.burnin
        val resample from params.resample
    output:
        tuple out_name, file("${out_name}_downsampled.trees") into downsampled_trees_ch
    publishDir params.downsampled_trees, mode: 'copy'
    script:
        """
        logcombiner -burnin ${burnin} -trees -resample $resample ${tree} ${out_name}_downsampled.trees
        """
    stub:
        """
        touch ${out_name}_downsampled.trees
        """
}

process consensus_tree {
    input:
        tuple val(out_name), file(trees) from trees_ch2
        val burnin from params.burnin
    output:
        tuple out_name, file("${out_name}.nexus") into consensus_tree_ch
    publishDir params.consensus_trees, mode: 'copy'
    label 'treeannotator'
    script:
        """
        treeannotator -burnin ${burnin} -heights median ${trees} ${out_name}.nexus
        """
    stub:
        """
        touch ${out_name}.nexus
        """
}

Channel.fromPath(params.xml_dir+"/*_final.fa.gz")
    .map {it -> [it.name.replace("_final.fa.gz", ""), it]}
    .set {algn_ch}

process generate_traits {
    input:
        tuple out_name, file(algn) from algn_ch
    output:
        tuple out_name, file(algn), file("${out_name}.tsv") into traits_ch
    publishDir params.traits_dir, mode: 'copy'
    script:
        """
        gunzip -c $algn | grep "^>" | awk -F "|" 'BEGIN{print "taxon\tlocation"}{if(\$5 ~ /San_Diego_County/){split(\$5,n,"_region-");print substr(\$0,2)"\t""San_Diego_"n[2]} else {print substr(\$0,2)"\tOutside_San_Diego"}}' | grep -v "EPI_ISL_402124|2019-12-30|China|Hubei|Wuhan" > ${out_name}.tsv
        """
    stub:
        """
        touch ${out_name}.tsv
        """

}

discrete_xml_in_ch = traits_ch.join(downsampled_trees_ch)

process generate_discrete_xml {
    stageInMode 'copy'
    input:
        tuple out_name, file(algn), file(traits), file(downsampled_trees) from discrete_xml_in_ch
        path emp_discrete_xml_template from params.emp_discrete_xml_template
    output:
        tuple out_name, file("${out_name}_emp_discrete.xml") into generate_discrete_xml_ch
    publishDir params.discrete_xml, mode: 'copy'
    script:
        """
        exclude_selected_seqs.py <(echo "EPI_ISL_402124") $algn ${algn}_tmp
        gunzip -c ${algn}_tmp > ${algn.getName()}.tmp

        beastgen -D "name=./downsampled_trees/${downsampled_trees.baseName},outname=${out_name}_emp_discrete" -date_order 2 -date_prefix "|" -date_format "yyyy-MM-dd" -traits $traits $emp_discrete_xml_template ${algn.getName()}.tmp ${out_name}_emp_discrete.xml
        rm ${algn}_tmp ${algn.getName()}.tmp
        """
    stub:
        """
        touch ${out_name}_emp_discrete.xml
        """
}

// beastgen -D "name=./downsampled_1000/B.1-112358_downsampled.trees,outname=B.1-112358_emp_discrete" -date_order 2 -date_prefix "|" -date_format "yyyy-MM-dd" empirical_discrete.template B.1-112358_final.fa.gz.tmp B.1-112358_emp_discrete.xml