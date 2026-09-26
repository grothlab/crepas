process SAMBAMBA_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:1.0.1--h6f6fda4_0':
        'biocontainers/sambamba:1.0.1--h6f6fda4_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('sambamba'), eval("sambamba --version 2>&1 | grep -m1 sambamba | awk '{print \\\$2}'"), emit: versions_sambamba, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // 0.1 GB is subtracted and rounded down to avoid issues with sambamba's memory allocation
    def memory = task.memory ? "--memory-limit ${(task.memory - 0.1.GB).toGiga()}GB" : ''
    """
    sambamba sort \\
        ${args} \\
        --nthreads ${task.cpus} \\
        ${memory} \\
        --tmpdir ./ \\
        --out ${prefix}.bam \\
        ${bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
