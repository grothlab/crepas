process RFD_TO_IZ {
    tag "$meta.id"
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/99b2a7149c1943265175ef013ca8c247c9847b14a0d2e802d0cbdda4a58458a5/data' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_r-argparse_r-ggpmisc_pruned:2c7c689513df97b4' }"

    input:
    tuple val(meta), path(okseq_rfd_file)
    tuple val(meta2), path(blacklist)
    tuple val(meta3), path(chrom_sizes)

    output:
    tuple val(meta), path("*.prefiltered.bed"), emit: okseq_filtered_bed
    tuple val(meta), path("*.init_zones.bed"),  emit: iz_bed
    tuple val(meta), path("*.rm_overlaps.bed"), emit: iz_rm_overlaps_bed
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^R version //; s/ .*\$//'"), topic: versions, emit: versions_rbase

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def blacklist_arg      = blacklist ? "--blacklist $blacklist" : ''

    """
    rfd_to_iz.R \\
        --okseq_rfd_file ${okseq_rfd_file} \\
        ${blacklist_arg} \\
        --chrom_sizes ${chrom_sizes} \\
        --prefix ${prefix} \\
        --outdir ./ \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.prefiltered.bed
    touch  ${prefix}.init_zones.bed
    touch  ${prefix}.rm_overlaps.bed
    """
}
