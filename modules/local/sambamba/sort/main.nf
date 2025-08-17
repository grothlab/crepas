process SAMBAMBA_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:1.0.1--h6f6fda4_0':
        'biocontainers/sambamba:1.0.1--h6f6fda4_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"             , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}

