process BIGTOOLS_BIGWIGAVERAGEOVERBED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bigtools:0.5.4--hc1c3326_1' :
        'quay.io/biocontainers/bigtools:0.5.6--hc1c3326_0' }"

    input:
    tuple val(meta), path(bigwig)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // From --help: "Note: for parts of the runtime, the actual usage may be nthreads+1"
    def nthreads = task.cpus > 1 ? task.cpus - 1 : task.cpus
    """
    bigtools bigwigaverageoverbed \\
        ${args} \\
        --nthreads ${nthreads} \\
        ${bigwig} \\
        ${bed} \\
        ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigtools-bigwigaverageoverbed: \$( bigtools bigwigaverageoverbed --version | sed 's/^bigtools-bigwigaverageoverbed //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigtools-bigwigaverageoverbed: \$( bigtools bigwigaverageoverbed --version | sed 's/^bigtools-bigwigaverageoverbed //' )
    END_VERSIONS
    """
}
