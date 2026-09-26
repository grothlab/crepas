process PLOT_MACS3_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37ab8c363a0cdc55501790b81945b691ad907534510a78661f30ff1e6f65b8ec/data':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-complexheatmap_r-base_r-optparse_pruned:b699da7f1540e3bb' }"

    input:
    tuple val(meta), path(peaks)
    val is_narrow_peak

    output:
    tuple val(meta), path("*.txt")       , emit: txt
    tuple val(meta), path("*.pdf")       , emit: pdf
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^.*R version //; s/ .*\$//'"), emit: versions_rbase, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def peak_type = is_narrow_peak ? 'narrowPeak' : 'broadPeak'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_macs3_qc.r \\
        -i ${peaks.join(',')} \\
        -s ${peaks.join(',').replaceAll("_peaks.${peak_type}","")} \\
        -p $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    touch ${prefix}.pdf
    """
}
