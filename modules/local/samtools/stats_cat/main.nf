process SAMTOOLS_STATS_CAT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:22.04'
        : 'nf-core/ubuntu:22.04'}"

    input:
    file t_stats

    output:
    path "*.tsv", emit: cat
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: 'samtools_stats_cat'
    """
    (head -n 1 ${t_stats[0]} && tail -q -n +2 ${t_stats.join(' ')}) \\
    | awk 'NR==1{print;next} {print \$0 | "sort -k1,1"}' \\
    > ${prefix}.tsv
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "samtools_stats_cat"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
