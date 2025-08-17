process SAMTOOLS_STATS_CAT {
    tag "$archive"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    file(t_stats)

    output:
    path "*.tsv"            , emit: cat
    path  "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'samtools_stats_cat'
    """
    (head -n 1 ${t_stats[0]} && tail -q -n +2 ${t_stats.join(' ')}) \\
    | awk 'NR==1{print;next} {print \$0 | "sort -k1,1"}' \\
    > ${prefix}.tsv
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        #samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "samtools_stats_cat"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        #samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
