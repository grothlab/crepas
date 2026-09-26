process BEDGRAPH_NORMALIZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed4799c386bc676b388fb01df7e97792ed9dd5e84c8f20ef568feca5dc2cccdb/data' :
        'community.wave.seqera.io/library/gawk:5.4.0--0877e7a42fa88325' }"

    input:
    tuple val(meta), path(bdg)
    val extension

    output:
    tuple val(meta), path("*.${extension}") , emit: normalized
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*//'"), emit: versions_awk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def norm_factor = args2.contains('--norm_factor ') ? args2.split('--norm_factor ')[1].split(' ')[0] : '1'
    def norm_operation = args2.contains('--norm_operation ') ? args2.split('--norm_operation ')[1].split(' ')[0] : 'multiply'
    def pseudocount = args2.contains('--add_pseudocount') ? 1 : 0

    """
    awk \\
        ${args} \\
        -v norm_factor=${norm_factor} \\
        -v operation=${norm_operation} \\
        -v pseudocount=${pseudocount} \\
        'BEGIN { OFS="\\t" } { \\
        if (operation == "multiply") \$4 = (\$4 + pseudocount) * norm_factor; \\
        else if (operation == "divide") \$4 = (\$4 + pseudocount) / norm_factor; \\
        else if (operation == "add") \$4 = (\$4 + pseudocount) + norm_factor; \\
        else if (operation == "subtract") \$4 = (\$4 + pseudocount) - norm_factor; \\
        print }' \\
        ${bdg} \\
        > ${prefix}.${extension}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.${extension}
    """
}
