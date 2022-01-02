Channel
    .fromPath(params.trees_dir+"*.trees")
    .map {it -> [it.name.replace(".trees", ""), it]}
    .set {trees_ch}

process downsample_trees {
    input:
        tuple val(out_name), file(tree) from trees_ch
        val burnin from params.burnin
    output:
        tuple out_name, file("${out_name}_downsampled.trees") into downsampled_trees_ch
    publishDir params.downsampled_trees, mode: 'copy'
    script:
        """
        logcombiner -burnin ${burnin} -trees ${tree} ${out_name}_downsampled.trees
        """
    stub:
        """
        touch ${out_name}_downsampled.trees
        """
}

process consensus_tree {
    input:
        tuple val(out_name), file(trees) from trees_ch
        val burnin from params.burnin
    output:
        tuple out_name, file("${out_name}.nexus") into consensus_tree_ch
    publishDir params.consensus_trees, mode: 'copy'
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
        gunzip -c $algn | grep "^>" | awk -F "|" 'BEGIN{print "taxon\tlocation"}{if($5 ~ /San_Diego_County/){split($5,n,"_region-");print substr($0,2)"\t""San_Diego_"n[2]} else {print substr($0,2)"\tOutside_San_Diego"}}' > ${out_name}.tsv
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
        tuple out_name, file(algn), file(traits), file(downsampled_trees) from downsampled_trees_traits_ch
        path emp_discrete_xml_template from params.emp_discrete_xml_template
    output:
        tuple file("${out_name}_emp_discrete.xml") into generate_discrete_xml_ch
    publishDir params.discrete_xml, mode: 'copy'
    script:
        """
        gunzip -c $algn > ${out_name}_final.fa

        beastgen -D "name=./downsampled_1000/${out_name},outname=${out_name}_emp_discrete" -date_order 2 -date_prefix "|" -date_format "yyyy-MM-dd" $emp_discrete_xml_template ${out_name}_final.fa ${out_name}_emp_discrete.xml
        """
    stub:
        """
        touch ${out_name}_emp_discrete.xml
        """
}