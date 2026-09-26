process FILE_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198ab15844726d095c90e454aa9b1cf91b7e3517dd62451791c596e2a2229082/data' :
        'community.wave.seqera.io/library/coreutils_gawk_sed:e167f5dd848b5ae9' }"

    input:
    tuple val(meta), path(file)
    val   extension

    output:
    tuple val(meta), path("*.${extension}"), emit: sorted
    tuple val("${task.process}"), val('coreutils'), eval("env sort --version | head -1 | sed 's/^sort (GNU coreutils) //'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def buffer  = task.memory ? "--buffer-size=${task.memory.toGiga().intdiv(2)}G" : ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    LC_COLLATE=C sort -T '.' \\
        $args \\
        $args2 \\
        --parallel=$task.cpus \\
        $buffer \\
        $file \\
        > ${prefix}.${extension}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.${extension}
    """
}
