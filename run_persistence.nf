Channel
    .fromPath(params.discrete_consensus_trees+"/*_emp_discrete.nexus")
    .map {it -> [it.name.replace("_emp_discrete.nexus", ""), it]}
    .set {consesus_trees_ch}

Channel
    .fromPath(params.discrete_xml+"/*_emp_discrete.trees")
    .map {it -> [it.name.replace("_emp_discrete.trees", ""), it]}
    .set {trees_ch}

consesus_trees_ch.join(trees_ch).set {persistence_in_ch}

process run_persistence {
    input:
        tuple val(out_name), file(tree), file(trees_log) from persistence_in_ch
    output:
        tuple out_name, file("${out_name}_persistence_summarizer.csv") into run_persistence_ch
    publishDir params.persistence, mode: 'copy'
    script:
        """
        intervaltimes=\$(generate_evaluation_ancestral_times.py ${tree});
        java -cp /beast-mcmc/build/dist/beast.jar dr.app.tools.PersistenceSummarizer -nodeStateAnnotation location \$intervaltimes $trees_log ${out_name}_persistence_summarizer.csv
        """
    stub:
        """
        touch ${out_name}_persistence_summarizer.csv
        """
}
