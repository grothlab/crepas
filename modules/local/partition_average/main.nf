process PARTITION_AVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198ab15844726d095c90e454aa9b1cf91b7e3517dd62451791c596e2a2229082/data' :
        'community.wave.seqera.io/library/coreutils_gawk_sed:e167f5dd848b5ae9' }"

    input:
    tuple val(meta) , path(partitions)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.flT_by_counts.tsv"), emit: filtered_tsv
    tuple val(meta), path("${prefix}.flT_by_counts.bdg"), emit: filtered_bdg
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), topic: versions, emit: versions_awk
    tuple val("${task.process}"), val('coreutils'), eval("env paste --version | head -1 | sed 's/^paste (GNU coreutils) //'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.average"
    """
    # This file will contain the following 12 columns, with 4-12 averaged across all input files:
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

    paste ${args} ${partitions.join(' ')} \\
        | awk -F'\\t' 'BEGIN {OFS="\\t"}
        {
            nfiles = NF / 12
            # Print first three columns once (chromosome, start, end)
            printf "%s\\t%s\\t%s", \$1, \$2, \$3

            # For each of the 9 columns to average (4..12)
            for (k = 4; k <= 12; k++) {
                sum = 0
                count = 0
                # Sum across all files for this column
                for (f = 0; f < nfiles; f++) {
                    col = 12*f + k
                    val = \$col
                    if (val != "NA") {
                        sum += val
                        count++
                    }
                }
                printf "\\t%s", (count > 0 ? sum/count : "NA")
            }
            printf "\\n"
        }' \\
        > ${prefix}.tsv

    # Create filtered version with rows where either bwaob_fwd_counts or bwaob_rev_counts > 0
    awk '\$4 > 0 || \$5 > 0' ${prefix}.tsv \\
    > ${prefix}.flT_by_counts.tsv

    # The following creates a bedGraph file with the 4th column averaged across all input files:
    #   1. chromosome
    #   2. start
    #   3. end
    #   4. RFD_smooth: Smoothed partition or RFD score

    awk ${args2} '{ printf "%s\\t%d\\t%d\\t%2.3f\\n", \$1, \$2, \$3, \$9 }' \\
    ${prefix}.flT_by_counts.tsv \\
    > ${prefix}.flT_by_counts.bdg
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.average"
    """
    touch ${prefix}.tsv
    touch  ${prefix}.flT_by_counts.tsv
    touch ${prefix}.flT_by_counts.bdg
    """
}
