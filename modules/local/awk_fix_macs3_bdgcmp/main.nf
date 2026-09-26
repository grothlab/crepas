process AWK_FIX_MACS3_BDGCMP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198ab15844726d095c90e454aa9b1cf91b7e3517dd62451791c596e2a2229082/data' :
        'community.wave.seqera.io/library/coreutils_gawk_sed:e167f5dd848b5ae9' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), emit: versions_awk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.fixed"

    """
    awk \\
        $args \\
        '{if (\$3 != -1) print \$0}' \\
        $bed \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.fixed"
    """
    touch  ${prefix}.bed
    """
}
