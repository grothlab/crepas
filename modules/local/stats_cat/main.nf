process STATS_CAT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198ab15844726d095c90e454aa9b1cf91b7e3517dd62451791c596e2a2229082/data'
        : 'community.wave.seqera.io/library/coreutils_gawk_sed:e167f5dd848b5ae9'}"

    input:
    path t_stats

    output:
    path "*.tsv", emit: cat
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), emit: versions_awk, topic: versions
    tuple val("${task.process}"), val('coreutils'), eval("sort --version 2>&1 | head -1 | sed 's/^sort (GNU coreutils) //'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: 'samtools_stats_cat'
    """
    (head -n 1 ${t_stats[0]} && tail -q -n +2 ${t_stats.join(' ')}) \\
    | awk 'NR==1{print;next} {print \$0 | "sort -k1,1"}' \\
    > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "samtools_stats_cat"
    """
    touch ${prefix}.tsv
    """
}
