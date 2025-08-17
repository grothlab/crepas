process BEDGRAPH_NORMALIZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(bdg)

    output:
    tuple val(meta), path("*.bedgraph") , emit: bedgraph
    path  "versions.yml"                , emit: versions

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
        > ${prefix}.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
