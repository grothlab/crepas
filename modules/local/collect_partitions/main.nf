//
// Collect windows, bigwigaverageoverbed scores, normalized scores, and RFD scores into a partition file
// Generate also a bedGraph file from the partition file
//
process COLLECT_PARTITIONS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(windows), path(bwaob_fwd), path(bwaob_rev), path(norm_or_smi_fwd), path(norm_or_smi_rev), path(rfd)

    output:
    tuple val(meta), path("*.tsv") , emit: tsv
    tuple val(meta), path("*.bdg") , emit: bdg
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def buffer   = task.memory ? "--buffer-size=${task.memory.toGiga().intdiv(2)}G" : ''
    def sort_cmd = "| LC_ALL=C sort --parallel=$task.cpus $buffer -k1,1 -k2,2n"

    """
    # Check that all input files have the same number of lines
    n1=\$(wc -l < ${windows})
    n2=\$(wc -l < ${bwaob_fwd})
    n3=\$(wc -l < ${bwaob_rev})
    n4=\$(wc -l < ${norm_or_smi_fwd})
    n5=\$(wc -l < ${norm_or_smi_rev})
    n6=\$(wc -l < ${rfd})
    if [ \$n1 -ne \$n2 ] || [ \$n1 -ne \$n3 ] || [ \$n1 -ne \$n4 ] || [ \$n1 -ne \$n5 ] || [ \$n1 -ne \$n6 ]; then
        echo "Error: Input files have different number of lines" >&2
        exit 1
    fi

    # Paste columns and filter for bwaob_fwd > 0 and bwaob_rev > 0
    paste ${windows} ${bwaob_fwd} ${bwaob_rev} ${norm_or_smi_fwd} ${norm_or_smi_rev} ${rfd} \\
    | awk '\$8 > 0 && \$14 > 0' \\
    | cut -f -3,8,14,20,26,30- \\
    ${sort_cmd} \\
    > ${prefix}.tsv

    # Making bedGraph
    awk '{ printf "%s\\t%d\\t%d\\t%2.3f\\n", \$1, \$2, \$3, \$9 }' \\
    ${prefix}.tsv \\
    > ${prefix}.bdg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.tsv
    touch  ${prefix}.bdg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
