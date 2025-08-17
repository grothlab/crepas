process PARTITION_OR_RFD_SMOOTH {
    tag "${meta.id}"
    label 'process_low_memory'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/perl:5.26.2'
        : 'quay.io/biocontainers/perl:5.26.2'}"

    input:
    tuple val(meta), val(partition_or_rfd), path(f_tab), path(r_tab)
    val radius
    val dradius
    val zradius

    output:
    tuple val(meta), path("*.tsv"), emit: rfd
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    partition_or_rfd_smooth.pl \\
        ${partition_or_rfd} \\
        ${args} \\
        ${f_tab} \\
        ${r_tab} \\
        ${radius} \\
        ${dradius} \\
        ${zradius} \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}
