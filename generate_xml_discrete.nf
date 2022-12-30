Channel
    .fromPath(params.trees_dir+"/*.trees")
    .map {it -> [it.name.replace(".trees", ""), it]}
    .into {trees_ch1; trees_ch2}


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

