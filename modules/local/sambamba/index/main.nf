process SAMBAMBA_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:1.0.1--h6f6fda4_0':
        'biocontainers/sambamba:1.0.1--h6f6fda4_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bai"), emit: bai
    tuple val("${task.process}"), val('sambamba'), eval("sambamba --version 2>&1 | grep -m1 sambamba | awk '{print \\\$2}'"), emit: versions_sambamba, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sambamba index \\
        ${args} \\
        --nthreads ${task.cpus} \\
        ${bam} \\
        ${prefix}.bai
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bai
    """
}
