//
// Filter LanceOtron candidate peaks by score and convert them to narrowPeak
//
process LANCEOTRON_TO_PEAK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed4799c386bc676b388fb01df7e97792ed9dd5e84c8f20ef568feca5dc2cccdb/data'
        : 'community.wave.seqera.io/library/gawk:5.4.0--0877e7a42fa88325'}"

    input:
    tuple val(meta), path(candidates)

    output:
    tuple val(meta), path("*.lanceotron.narrowPeak"), emit: narrowpeak
    tuple val(meta), path("*.lanceotron.tsv")       , emit: tsv
    tuple val(meta), path("*.lanceotron.scores")    , emit: scores
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*//'"), emit: versions_awk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Keep the candidates with overall_peak_score >= cutoff (-v cutoff=<score> in ext.args; default 0), as:
    # - the LanceOtron table, with its header
    # - narrowPeak: chrom, start, end, name (peak_N), score (overall_peak_score x 1000), strand,
    #   signalValue (enrichment_score), pValue (pvalue_input, or pvalue_chrom without a control),
    #   qValue (-1), summit (max_coverage_position, relative to start)
    # - headerless scores for the peak QC plots, with the same columns with or without a control: chrom,
    #   start, end, overall_peak_score, shape_score, enrichment_score, pValue
    # Columns are looked up by name: with a control, LanceOtron adds pvalue_input after enrichment_score.
    awk \\
        -F'\\t' \\
        -v OFS='\\t' \\
        -v cutoff=0 \\
        ${args} \\
        -v table=${prefix}.lanceotron.tsv \\
        -v scores=${prefix}.lanceotron.scores \\
        'NR == 1 {
            for (i = 1; i <= NF; i++) col[\$i] = i
            pval = ("pvalue_input" in col) ? col["pvalue_input"] : col["pvalue_chrom"]
            print > table
            next
        }
        \$col["overall_peak_score"] >= cutoff {
            n++
            print > table
            print \$1, \$2, \$3, "peak_" n, int(\$col["overall_peak_score"] * 1000 + 0.5), ".", \$col["enrichment_score"], \$pval, -1, \$col["max_coverage_position"] - \$2
            print \$1, \$2, \$3, \$col["overall_peak_score"], \$col["shape_score"], \$col["enrichment_score"], \$pval > scores
        }' \\
        ${candidates} \\
        > ${prefix}.lanceotron.narrowPeak

    # Empty outputs when there are no candidates or no peaks above the cutoff
    touch ${prefix}.lanceotron.tsv ${prefix}.lanceotron.scores
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lanceotron.narrowPeak
    touch ${prefix}.lanceotron.tsv
    touch ${prefix}.lanceotron.scores
    """
}
