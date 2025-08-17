process TELOCAL_INDEXER {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://mhammelllab/telocal:latest'
        : 'docker://mhammelllab/telocal:latest'}"

    input:
    tuple val(meta), path(gtf)
    val index_type

    output:
    tuple val(meta), path("*.{ind,locInd}"), emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
     TElocal_indexer.py \\
        ${args} \\
        --afile ${gtf} \\
        --itype ${index_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TElocal_indexer: \$(TElocal_indexer.py --version | sed 's/TElocal_indexer //g')
    END_VERSIONS
    """

    stub:
    def suffix = index_type == 'gene' ? 'ind' : 'locInd'
    """
    touch TElocal_index.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TElocal_indexer: \$(TElocal_indexer.py --version | sed 's/TElocal_indexer //g')
    END_VERSIONS
    """
}
