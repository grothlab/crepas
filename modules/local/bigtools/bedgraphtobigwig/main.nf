process BIGTOOLS_BEDGRAPHTOBIGWIG {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bigtools:0.5.4--hc1c3326_1'
        : 'quay.io/biocontainers/bigtools:0.5.4--hc1c3326_1'}"

    input:
    tuple val(meta), path(bedgraph)
    tuple val(meta2), path(sizes)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    tuple val("${task.process}"), val('bigtools'), eval("bigtools bedgraphtobigwig --version | sed 's/^bigtools-bedgraphtobigwig //'"), emit: versions_bigtools, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw
    """
}
