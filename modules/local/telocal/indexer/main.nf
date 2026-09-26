process TELOCAL_INDEXER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a2/a2f1e953e388761ab988e1bbdeb3502aac651ca5ea5ee093e1e2fcc02a1f7773/data' :
        'community.wave.seqera.io/library/pip_telocal:5ebb38192e05596d' }"

    input:
    tuple val(meta), path(gtf)
    val index_type

    output:
    tuple val(meta), path("*.{ind,locInd}"), emit: index
    tuple val("${task.process}"), val('TElocal_indexer'), eval("TElocal_indexer.py --version | sed 's/TElocal_indexer //g'"), emit: versions_telocal_indexer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
     TElocal_indexer.py \\
        ${args} \\
        --afile ${gtf} \\
        --itype ${index_type}
    """

    stub:
    def suffix = index_type == 'gene' ? 'ind' : 'locInd'
    """
    touch TElocal_index.${suffix}
    """
}
