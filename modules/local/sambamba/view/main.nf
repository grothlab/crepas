process SAMBAMBA_VIEW {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:1.0.1--h6f6fda4_0':
        'biocontainers/sambamba:1.0.1--h6f6fda4_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.bam"),         emit: bam,  optional: true
    tuple val(meta), path("*.sam"),         emit: sam, optional: true
    tuple val(meta), path("*.json"),        emit: json, optional: true
    tuple val(meta), path("*.msgpack"),     emit: msgpack, optional: true
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                = task.ext.args ?: ''
    def prefix              = task.ext.prefix ?: "${meta.id}"
    def whitelist_params    = bed ? "--regions $bed" : ''
    def extension           = args.contains("--format sam") ? "sam" :
                                args.contains("--format json") ? "json" :
                                    args.contains("--format msgpack") ? "msgpack" :
                                        "bam"
    """
    sambamba \\
        view \\
        $args \\
        $whitelist_params \\
        --nthreads $task.cpus \\
        $bam \\
        --output-filename ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}
