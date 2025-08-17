process TECOUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tetranscripts:2.2.3--pyh7cba7a3_0':
        'quay.io/biocontainers/tetranscripts:2.2.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(genic_gtf_or_index)
    tuple val(meta3), path(te_gtf_or_index)

    output:
    tuple val(meta), path("*.cntTable"),                        emit: counts
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args         ?: ""
    def prefix          = task.ext.prefix       ?: "${meta.id}"
    """
     TEcount \\
        ${args} \\
        --BAM ${bam} \\
        --GTF ${genic_gtf_or_index} \\
        --TE ${te_gtf_or_index} \\
        --project ${prefix} \\
        --outdir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TEcount: \$(TEcount --version | sed 's/TEcount //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cntTable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TEcount: \$(TEcount --version | sed 's/TEcount //g')
    END_VERSIONS
    """
}
