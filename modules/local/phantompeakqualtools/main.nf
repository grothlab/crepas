process PHANTOMPEAKQUALTOOLS {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a1cddfad5b503ee347cc5de17d172e1876c547fca00aa844559c9e764fb400f/data' :
        'community.wave.seqera.io/library/phantompeakqualtools:1.2.2--f8026fe2526a5e18' }"

    input:
    tuple val(meta), path(chip), path(input_control)

    output:
    tuple val(meta), path("*.ccscores")  , emit: ccscores
    tuple val(meta), path("*.pdf")  , emit: pdf
    tuple val(meta), path("*.Rdata"), emit: rdata
    tuple val(meta), path("*.narrowPeak.gz"), emit: narrowpeak, optional: true
    tuple val(meta), path("*.regionPeak.gz"), emit: regionpeak, optional: true
    tuple val("${task.process}"), val('phantompeakqualtools'), val('1.2.2'), topic: versions, emit: versions_phantompeakqualtools
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.spp"
    def input_control_arg = input_control ? "-i='${input_control}'" : ''
    def peaks_arg = input_control ? "-savn='${prefix}.narrowPeak' -savr='${prefix}.regionPeak'" : ""
    """
    RUN_SPP=`which run_spp.R`
    Rscript ${args} -e "library(caTools); source(\\"\$RUN_SPP\\")" \\
        -c="${chip}" \\
        ${input_control_arg} \\
        -savp="${prefix}.pdf" \\
        -savd="${prefix}.Rdata" \\
        -out="${prefix}.ccscores" \\
        ${args2} \\
        -p=${task.cpus} \\
        ${peaks_arg}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.spp"
    def peaks_arg = input_control ? "touch ${prefix}.narrowPeak; touch ${prefix}.regionPeak;" : ""
    """
    touch ${prefix}.pdf
    touch ${prefix}.Rdata
    touch ${prefix}.ccscores
    ${peaks_arg}
    """
}
