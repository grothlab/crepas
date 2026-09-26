process PARTITION_OR_RFD_PLOT {
    tag "$meta.id"
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/99b2a7149c1943265175ef013ca8c247c9847b14a0d2e802d0cbdda4a58458a5/data' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_r-argparse_r-ggpmisc_pruned:2c7c689513df97b4' }"

    input:
    tuple val(meta), path(partition), path(strandedinput), path(scarminusinput), path(okseq_rfd_file), path(initiation_zones)
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
    tuple val(meta), path("*_mean_values.tsv"),    emit: mean_values
    tuple val(meta), path("*.scatter_plot.pdf"),   emit: scatter_pdf, optional:true
    tuple val(meta), path("*.scatter_plot.png"),   emit: scatter_png, optional:true
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^R version //; s/ .*\$//'"), topic: versions, emit: versions_rbase

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def scar_partition_arg = partition ? "--scar_partition_file ${partition.join(' ')}" : ''
    def scarminusinput_arg = scarminusinput ? "--scarminusinput_partition_file ${scarminusinput.join(' ')}" : ''
    def strandedinput_arg  = strandedinput ? "--strandedinput_partition_file ${strandedinput.join(' ')}" : ''
    def okseq_rfd_arg      = okseq_rfd_file ? "--okseq_rfd_file ${okseq_rfd_file}" : ''
    def iz_arg             = initiation_zones ? "--initiation_zones ${initiation_zones}" : ''
    def blacklist_arg      = blacklist ? "--blacklist ${blacklist}" : ''

    """
    partition_or_rfd_plot.R \\
        ${scar_partition_arg} \\
        ${scarminusinput_arg} \\
        ${strandedinput_arg} \\
        ${okseq_rfd_arg} \\
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
    """
}
