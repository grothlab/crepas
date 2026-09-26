//
// Collect windows, bigwigaverageoverbed scores, normalized scores, and RFD scores into a partition file
// Generate also a bedGraph file from the partition file
//
process COLLECT_PARTITIONS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198ab15844726d095c90e454aa9b1cf91b7e3517dd62451791c596e2a2229082/data'
        : 'community.wave.seqera.io/library/coreutils_gawk_sed:e167f5dd848b5ae9'}"

    input:
    tuple val(meta), path(windows), path(bwaob_fwd), path(bwaob_rev), path(norm_or_smi_fwd), path(norm_or_smi_rev), path(rfd)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.flT_by_counts.tsv"), emit: filtered_tsv
    tuple val(meta), path("${prefix}.flT_by_counts.bdg"), emit: filtered_bdg
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), topic: versions, emit: versions_awk
    tuple val("${task.process}"), val('coreutils'), eval("env paste --version | head -1 | sed 's/^paste (GNU coreutils) //'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.collect"
    def buffer = task.memory ? "--buffer-size=${task.memory.toGiga().intdiv(2)}G" : ''
    def sort_cmd = "| LC_ALL=C sort --parallel=${task.cpus} ${buffer} -k1,1 -k2,2n"

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

    # This file will contain the following 12 columns:
    #   1. chromosome
    #   2. start
    #   3. end
    #   4. bwaob_fwd_counts: Forward raw counts in bin
    #   5. bwaob_rev_counts: Reverse raw counts in bin
    #   6. bwaob_fwd_RPM: Forward RPMs in bin
    #   7. bwaob_rev_RPM: Reverse RPMs in bin
    #   8. RFD_raw: Raw partition or RFD score
    #   9. RFD_smooth: Smoothed partition or RFD score
    #  10. RFD_deriv: Value of the derivative of the partition/RFD at this bin
    #  11. score: not used
    #  12. zero_deriv: second derivative at this bin

    # Paste columns and filter
    paste ${args} ${windows} ${bwaob_fwd} ${bwaob_rev} ${norm_or_smi_fwd} ${norm_or_smi_rev} ${rfd} \\
    | cut -f -3,8,14,20,26,30- \\
    ${sort_cmd} \\
    > ${prefix}.tsv

    # Create filtered version with rows where either bwaob_fwd_counts or bwaob_rev_counts > 0
    awk '\$4 > 0 || \$5 > 0' ${prefix}.tsv \\
    > ${prefix}.flT_by_counts.tsv

    # The following creates a bedGraph file with the following 4 columns:
    #   1. chromosome
    #   2. start
    #   3. end
    #   4. RFD_smooth: Smoothed partition or RFD score

    awk ${args2} '{ printf "%s\\t%d\\t%d\\t%2.3f\\n", \$1, \$2, \$3, \$9 }' \\
    ${prefix}.flT_by_counts.tsv \\
    > ${prefix}.flT_by_counts.bdg
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.collect"
    """
    touch  ${prefix}.tsv
    touch  ${prefix}.flT_by_counts.tsv
    touch  ${prefix}.flT_by_counts.bdg
    """
}
