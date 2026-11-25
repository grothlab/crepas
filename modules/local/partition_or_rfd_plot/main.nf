process PARTITION_OR_RFD_PLOT {
    tag "$meta.id"
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f4eebb916d1153809da29b40cbf56a5f78486e982a4a7373b3641cf6783236a/data' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_pip_python_pruned:6979014ed3536274' }"

    input:
    tuple val(meta), path(partition), path(strandedinput), path(minusinput), path(rfd_files), path(initiation_zones)
    tuple val(meta2), path(blacklist)
    tuple val(meta3), path(chrom_sizes)

    output:
    tuple val(meta), path("*_plot_smoothed.pdf"),       emit: plot_smoothed_pdf
    tuple val(meta), path("*_plot_smoothed.png"),       emit: plot_smoothed_png
    tuple val(meta), path("*_plot_smoothed.sd.pdf"),    emit: plot_smoothed_sd_pdf
    tuple val(meta), path("*_plot_smoothed.sd.png"),    emit: plot_smoothed_sd_png
    tuple val(meta), path("*_plot_smoothed.gam.pdf"),  emit: plot_smoothed_gam_pdf
    tuple val(meta), path("*_plot_smoothed.gam.png"),  emit: plot_smoothed_gam_png
    tuple val(meta), path("*_plot_raw.pdf"),        emit: plot_raw_pdf
    tuple val(meta), path("*_plot_raw.png"),        emit: plot_raw_png
    tuple val(meta), path("*_plot_raw.sd.pdf"),     emit: plot_raw_sd_pdf
    tuple val(meta), path("*_plot_raw.sd.png"),     emit: plot_raw_sd_png
    tuple val(meta), path("*_plot_smoothed.iqr.pdf"), emit: plot_smoothed_iqr_pdf
    tuple val(meta), path("*_plot_smoothed.iqr.png"), emit: plot_smoothed_iqr_png
    tuple val(meta), path("*_plot_raw.iqr.pdf"),      emit: plot_raw_iqr_pdf
    tuple val(meta), path("*_plot_raw.iqr.png"),      emit: plot_raw_iqr_png
    tuple val(meta), path("*_quantile_values.tsv"),   emit: quantile_values
    tuple val(meta), path("*_mean_values.tsv"),    emit: mean_values
    tuple val(meta), path("*_violin_plot_raw.pdf"),    emit: violin_pdf, optional:true
    tuple val(meta), path("*_violin_plot_raw.png"),    emit: violin_png, optional:true
    tuple val(meta), path("*_violin_values.tsv"),  emit: violin_values
    tuple val(meta), path("*_violin_stats.tsv"),   emit: violin_stats
    tuple val(meta), path("*_window_qc_plot.*.pdf"), emit: window_qc_pdf, optional:true
    tuple val(meta), path("*_window_qc_plot.*.png"), emit: window_qc_png, optional:true
    tuple val(meta), path("*_window_qc.tsv"),      emit: window_qc
    tuple val(meta), path("*_profile_difference_plot_raw.*.pdf"), emit: profile_difference_pdf, optional:true
    tuple val(meta), path("*_profile_difference_plot_raw.*.png"), emit: profile_difference_png, optional:true
    tuple val(meta), path("*_profile_differences.tsv"),     emit: profile_differences
    tuple val(meta), path("*_profile_tests.tsv"),           emit: profile_tests
    tuple val(meta), path("*.scatter_plot.pdf"),   emit: scatter_pdf, optional:true
    tuple val(meta), path("*.scatter_plot.png"),   emit: scatter_png, optional:true
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^R version //; s/ .*\$//'"), topic: versions, emit: versions_rbase

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def partition_arg      = partition ? "--partition_files ${partition.join(' ')}" : ''
    def minusinput_arg     = minusinput ? "--minusinput_partition_files ${minusinput.join(' ')}" : ''
    def strandedinput_arg  = strandedinput ? "--strandedinput_partition_files ${strandedinput.join(' ')}" : ''
    def rfd_arg            = rfd_files ? "--rfd_files ${rfd_files.join(' ')}" : ''
    def iz_arg             = initiation_zones ? "--initiation_zones ${initiation_zones}" : ''
    def blacklist_arg      = blacklist ? "--blacklist ${blacklist}" : ''

    """
    partition_or_rfd_plot.R \\
        ${partition_arg} \\
        ${minusinput_arg} \\
        ${strandedinput_arg} \\
        ${rfd_arg} \\
        ${iz_arg} \\
        ${blacklist_arg} \\
        --chrom_sizes ${chrom_sizes} \\
        --prefix ${prefix} \\
        --outdir ./ \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.scatter_plot.pdf
    touch  ${prefix}.scatter_plot.png
    touch  ${prefix}_plot_smoothed.pdf
    touch  ${prefix}_plot_smoothed.png
    touch  ${prefix}_plot_smoothed.sd.pdf
    touch  ${prefix}_plot_smoothed.sd.png
    touch  ${prefix}_plot_smoothed.gam.pdf
    touch  ${prefix}_plot_smoothed.gam.png
    touch  ${prefix}_plot_raw.pdf
    touch  ${prefix}_plot_raw.png
    touch  ${prefix}_plot_raw.sd.pdf
    touch  ${prefix}_plot_raw.sd.png
    touch  ${prefix}_mean_values.tsv
    touch  ${prefix}_plot_smoothed.iqr.pdf
    touch  ${prefix}_plot_smoothed.iqr.png
    touch  ${prefix}_plot_raw.iqr.pdf
    touch  ${prefix}_plot_raw.iqr.png
    touch  ${prefix}_quantile_values.tsv
    touch  ${prefix}_violin_plot_raw.pdf
    touch  ${prefix}_violin_plot_raw.png
    touch  ${prefix}_violin_values.tsv
    touch  ${prefix}_violin_stats.tsv
    touch  ${prefix}_window_qc_plot.fates.pdf
    touch  ${prefix}_window_qc_plot.coverage.pdf
    touch  ${prefix}_window_qc_plot.fates.png
    touch  ${prefix}_window_qc_plot.coverage.png
    touch  ${prefix}_window_qc.tsv
    touch  ${prefix}_profile_difference_plot_raw.partition.pdf
    touch  ${prefix}_profile_difference_plot_raw.partition.png
    touch  ${prefix}_profile_differences.tsv
    touch  ${prefix}_profile_tests.tsv
    """
}
