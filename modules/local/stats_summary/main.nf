process STATS_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e467b0d8acd4e9f2efcc04855cd4f9796278de135ccd3a357246c8b800aa1790/data'
        : 'community.wave.seqera.io/library/r-argparse_r-dplyr_r-forcats_r-ggplot2_pruned:0969cc079210fc80'}"

    input:
    path table
    val endogenous_genome_name
    val exogenous_genome_name

    output:
    path "*.all.tsv", emit: tsv
    path "*.totals.tsv", emit: totals_tsv
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^R version //; s/ .*\$//'"), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'final_samtools_stats_summary'
    def exo = exogenous_genome_name ? "--exogenous_genome_name ${exogenous_genome_name}" : ''
    """
    process_stats_summary.R \\
        ${args} \\
        --summary_table ${table} \\
        --endogenous_genome_name ${endogenous_genome_name} \\
        ${exo} \\
        --prefix ${prefix} \\
        --outdir ./
    """

    stub:
    def prefix = task.ext.prefix ?: 'final_samtools_stats_summary'
    """
    touch ${prefix}.all.tsv
    touch ${prefix}.totals.tsv
    """
}
