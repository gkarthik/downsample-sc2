Channel
    .fromPath(params.downsampled_trees+"/*_downsampled.trees")
    .map {it -> [it.name.replace("_downsampled.trees", ""), it]}
    .set {downsampled_trees_ch}

Channel
    .fromPath(params.traits_dir+"/*.tsv")
    .map {it -> [it.name.replace(".tsv", ""), it]}
    .set {traits_ch}

Channel.fromPath(params.xml_dir+"/*_final.fa.gz")
    .map {it -> [it.name.replace("_final.fa.gz", ""), it]}
    .set {algn_ch}

Channel.fromPath(params.glm_predictors_dir+"/*_design_matrices.txt")
    .map {it -> [it.name.replace("_design_matrices.txt", ""), it]}
    .set {glm_predictors_ch}

generate_xml_glm_ch = glm_predictors_ch.join(downsampled_trees_ch).join(traits_ch).join(algn_ch)

process generate_xml_glm {
    stageInMode 'copy'
    input:
        tuple out_name, file(glm_predictors), file(downsampled_trees), file(traits), file(algn) from generate_xml_glm_ch
        path xml_template from params.emp_discrete_xml_template
    output:
        tuple out_name, file("${out_name}_emp_discrete_glm.xml") into xml_glm_ch
    publishDir params.glm_xml_dir, mode: 'copy'
    script:
        """
        exclude_selected_seqs.py <(echo "EPI_ISL_402124") $algn ${algn}_tmp
        gunzip -c ${algn}_tmp > ${algn.getName()}.tmp

        beastgen -D "name=./downsampled_trees/${downsampled_trees.baseName},outname=${out_name}_emp_discrete,designmatrix=${glm_predictors}" -date_order 2 -date_prefix "|" -date_format "yyyy-MM-dd" -traits $traits $xml_template ${algn.getName()}.tmp ${out_name}_emp_discrete_glm.xml
        rm ${algn}_tmp ${algn.getName()}.tmp
        """
    stub:
        """
        touch ${out_name}_emp_discrete_glm.xml
        """
}