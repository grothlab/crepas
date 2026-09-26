
process MACS3_BDGCMP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macs3:3.0.1--py311h0152c62_3':
        'biocontainers/macs3:3.0.1--py311h0152c62_3' }"

    input:
    tuple val(meta), path(pileup_bdg), path(lambda_bdg)

    output:
    tuple val(meta), path("*.bedGraph"), emit: bdg
    tuple val("${task.process}"), val('macs3'), eval("macs3 --version | sed -e 's/macs3 //g'"), emit: versions_macs3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    macs3 \\
        bdgcmp \\
        ${args} \\
        --tfile ${pileup_bdg} \\
        --cfile ${lambda_bdg} \\
        --ofile ${prefix}.bedGraph
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bedGraph
    """
}
