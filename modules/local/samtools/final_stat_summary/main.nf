process FINAL_STAT_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/r-argparse_r-dplyr_r-forcats_r-ggplot2_pruned:e5fd467d162d4cbc'
        : 'community.wave.seqera.io/library/r-argparse_r-dplyr_r-forcats_r-ggplot2_pruned:0969cc079210fc80'}"

    input:
    path table
    val endogenous_genome_name
    val exogenous_genome_name

    output:
    path "*.tsv", emit: tsv
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
