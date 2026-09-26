process MACE_PREPROCESSOR {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mace:1.2--py27he7e273a_2'
        : 'biocontainers/mace:1.2--py27he7e273a_2'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(chrom_sizes)

    output:
    tuple val(meta), path("*_Forward.wig"), emit: forward_wig
    tuple val(meta), path("*_Reverse.wig"), emit: reverse_wig
    tuple val("${task.process}"), val('mace'), eval("preprocessor.py --version 2>&1 | sed 's/^preprocessor\\.py //'"), topic: versions, emit: versions_mace

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def treatment = bam ? "--inputFile ${bam.join(',')}" : ""
    """
    preprocessor.py \\
        ${args} \\
        ${treatment} \\
        --chromSize ${chrom_sizes} \\
        --outPrefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_Forward.wig
    touch ${prefix}_Reverse.wig
    """
}
