Channel
    .fromPath(params.discrete_consensus_trees+"*_emp_discrete.nexus")
    .map {it -> [it.name.replace("_emp_discrete.nexus", ""), it]}
    .set {consesus_trees_ch}

process run_persistence {
    input:
        tuple val(out_name), file(tree) from consesus_trees_ch
    output:
        tuple out_name, file("${out_name}_persistence_summarizer.csv") into run_persistence_ch
    publishDir params.persistence_summarizer, mode: 'copy'
    script:
        """
        intervaltimes=\$(python3 generate_evaluation_ancestral_times.py ../consensus_trees/${0}_emp_discrete.nexus);
        java -cp /home/gk/code/beast-mcmc/build/dist/beast.jar dr.app.tools.PersistenceSummarizer -nodeStateAnnotation location \$intervaltimes $tree ${out_name}_persistence_summarizer.csv
        """
    stub:
        """
        touch ${out_name}_persistence_summarizer.csv
        """
}