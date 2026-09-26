process SPIKEIN_BARCODE_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d4b9eb3e405ac6ec68a83c8529e725b3fb318c29a34efd117e2f2c27602c8ab/data' :
        'community.wave.seqera.io/library/bioconductor-complexheatmap_r-argparse_r-circlize_r-tidyverse:a65ce78f73a78e8c' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.tsv"), emit: merged_counts
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^.*R version //; s/ .*\$//'"), emit: versions_rbase, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.spikein_barcode_merged"
    def count_files = counts.sort()
    if (count_files.size() > 1) {
        """
        spikein_barcode_merge.R \\
            ${args} \\
            --count_tables ${count_files.join(' ')} \\
            --prefix ${prefix}
        """
    } else {
        """
        ln -s ${count_files[0]} ${prefix}.tsv
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.spikein_barcode_merged"
    """
    touch ${prefix}.tsv
    """
}
