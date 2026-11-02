//
// Convert epic2 output to diffusePeak (HOMER-compatible) and broadPeak (ENCODE) files
//
process EPIC2_TO_PEAK {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed4799c386bc676b388fb01df7e97792ed9dd5e84c8f20ef568feca5dc2cccdb/data'
        : 'community.wave.seqera.io/library/gawk:5.4.0--0877e7a42fa88325'}"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.diffusePeak"), emit: diffusepeak
    tuple val(meta), path("*.broadPeak"), emit: broadpeak
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*//'"), emit: versions_awk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Convert to diffusePeak (HOMER-compatible): chrom, start, end, name (peak_N), score, strand,
    # PValue, ChIPCount, InputCount, FDR, log2FoldChange
    awk \\
        -F'\\t' \\
        -v OFS='\\t' \\
        'NR==1 {next} \\
        { print \$1, \$2, \$3, "peak_" NR-1, \$5, \$6, \$4, \$7, \$8, \$9, \$10 }' \\
        ${bed} \\
        > ${prefix}.diffusePeak

    # Convert to broadPeak (ENCODE, as written by MACS3): chrom, start, end, name (peak_N), score,
    # strand, signalValue (fold enrichment, 2^log2FoldChange), pValue (-log10 PValue), qValue (-log10 FDR)
    awk \\
        -F'\\t' \\
        -v OFS='\\t' \\
        'function nlog10(x) { return x > 0 ? -log(x) / log(10) : 1000 } \\
        NR==1 {next} \\
        { print \$1, \$2, \$3, "peak_" NR-1, int(\$5), \$6, 2^\$10, nlog10(\$4), nlog10(\$9) }' \\
        ${bed} \\
        > ${prefix}.broadPeak
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.diffusePeak
    touch ${prefix}.broadPeak
    """
}
