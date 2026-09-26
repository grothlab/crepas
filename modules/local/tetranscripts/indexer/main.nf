process TETRANSCRIPTS_INDEXER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tetranscripts:2.2.3--pyh7cba7a3_0':
        'quay.io/biocontainers/tetranscripts:2.2.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(gtf)
    val index_type

    output:
    tuple val(meta), path("*.ind"),   emit: index
    tuple val("${task.process}"), val('TEtranscripts_indexer'), eval("TEtranscripts_indexer.py --version | sed 's/TEtranscripts_indexer //g'"), emit: versions_tetranscripts_indexer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args         ?: ""
    """
     TEtranscripts_indexer.py \\
        ${args} \\
        --afile ${gtf} \\
        --itype ${index_type}
    """

    stub:
    """
    touch TEtranscripts_index.ind
    """
}
