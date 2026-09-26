process SPIKEIN_BARCODE_QC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d4b9eb3e405ac6ec68a83c8529e725b3fb318c29a34efd117e2f2c27602c8ab/data' :
        'community.wave.seqera.io/library/bioconductor-complexheatmap_r-argparse_r-circlize_r-tidyverse:a65ce78f73a78e8c' }"

    input:
    tuple val(meta), path(counts), val(uniq_totals)

    output:
    tuple val(meta), path("*.long.tsv")   , emit: long
    tuple val(meta), path("*.summary.tsv"), emit: summary
    tuple val(meta), path("*.pdf")        , emit: plot_pdf
    tuple val(meta), path("*.png")        , emit: plot_png
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^.*R version //; s/ .*\$//'"), emit: versions_rbase, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.spikein_barcode_summary"
    """
    spikein_barcode_qc.R \\
        ${args} \\
        --count_tables ${counts.join(' ')} \\
        --uniq_totals ${uniq_totals.join(' ')} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.spikein_barcode_summary"
    """
    touch ${prefix}.long.tsv
    touch ${prefix}.summary.tsv
    touch ${prefix}.pdf
    touch ${prefix}.png
    """
}
