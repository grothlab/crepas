process MMR {
    tag "$meta.id"
    label 'process_high'

    // NOTE: there is no version information provided by the tool on the CLI
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/mmr:c5ce80a--0042edeccacdf2bc' :
        'community.wave.seqera.io/library/mmr:c5ce80a--1ec1f5037b507cc0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam")      , emit: bam
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def pair_usage_arg     = meta.single_end ? "" : "--pair-usage"

    """
    mmr \\
        $args \\
        $pair_usage_arg \\
        --threads $task.cpus \\
        --verbose \\
        -o ${prefix}.bam \\
        $bam 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmr: c5ce80a
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
                mmr: c5ce80a
    END_VERSIONS
    """
}
