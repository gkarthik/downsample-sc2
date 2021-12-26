Channel
    .fromPath(params.final_ml_dir+"*_guide_tree_rooted.treefile")
    .map {it -> [it.name.replace("_guide_tree_rooted.treefile", ""), it]}
    .set {guide_trees_ch}

Channel.fromPath(params.clock_rate_filter_two_dir+"/*_manual_check_filtered.fa.gz")
    .map {it -> [it.name.replace("_manual_check_filtered.fa.gz", ""), it]}
    .set {algn_ch}

ids_algn_ch = guide_trees_ch.join(algn_ch)

process generate_xml {
    stageInMode 'copy'
    input:
        tuple out_name, file(tree), file(algn) from ids_algn_ch
        path xml_template from params.xml_template
    output:
        tuple file("${out_name}.xml"), file("${out_name}_final.fa.gz"), file("${out_name}_final.treefile") into generate_xml_ch
    publishDir params.xml_dir, mode: 'copy'
    script:
        """
        exclude_selected_seqs.py <(echo "EPI_ISL_402124|2019-12-30|China|Hubei|Wuhan") $algn ${algn}_tmp
        gunzip -c ${algn}_tmp | sed "s/'/_/g" | sed 's/(/_/g' | sed 's/)/_/g' | sed 's/,/_/g' | gzip > ${out_name}_final.fa.gz
        rm ${algn}_tmp

        drop_root.R $tree ${out_name}_final.treefile

        a=\$(get_skygrid_params.py $out_name ${out_name}_final.fa.gz)
        b=(\$(echo \$a | tr ' ' '\n'))

        beastgen -D "name=${out_name},cutOff=\${b[0]},gridPoints=\${b[1]}" -date_order 2 -date_prefix "|" -date_format "yyyy-MM-dd" -tree ${out_name}_final.treefile $xml_template <(gunzip -c ${out_name}_final.fa.gz) ${out_name}.xml
        """
    stub:
        """
        touch ${out_name}.xml ${out_name}_final.fa.gz ${out_name}_final.treefile
        """
}