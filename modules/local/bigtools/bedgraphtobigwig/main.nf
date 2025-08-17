process BIGTOOLS_BEDGRAPHTOBIGWIG {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bigtools:0.5.4--hc1c3326_1'
        : 'quay.io/biocontainers/bigtools:0.5.6--hc1c3326_0'}"

    input:
    tuple val(meta), path(bedgraph)
    tuple val(meta2), path(sizes)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bigtools bedgraphtobigwig \\
        ${args} \\
        --nthreads ${task.cpus} \\
        ${bedgraph} \\
        ${sizes} \\
        ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigtools-bedgraphtobigwig: \$( bigtools bedgraphtobigwig --version | sed 's/^bigtools-bedgraphtobigwig //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigtools-bedgraphtobigwig: \$( bigtools bedgraphtobigwig --version | sed 's/^bigtools-bedgraphtobigwig //' )
    END_VERSIONS
    """
}
