process GFF3SORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7e888288455c539735f35c165c8b4f914aa986d283c239056bbdcbe932349875/data':
        'community.wave.seqera.io/library/gff3sort:0.1.a1a2bc9--3bdc1d11d91dce42' }"

    input:
    tuple val(meta), path(gxf)

    output:
    tuple val(meta), path("*.${gxf.extension}"), emit: sorted
    // WARN: gff3sort.pl has no --version flag; report the perl interpreter version instead.
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | grep 'This is perl' | sed 's/.*(v\\(.*\\)) built.*/\\1/'"), topic: versions, emit: versions_perl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gxf.baseName}.sorted"
    """
    gff3sort.pl \\
        ${args} \\
        ${gxf} \\
        > ${prefix}.${gxf.extension}
    """

    stub:
    def prefix = task.ext.prefix ?: "${gxf.baseName}.sorted"
    """
    touch ${prefix}.${gxf.extension}
    """
}
