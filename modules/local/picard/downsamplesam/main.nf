process PICARD_DOWNSAMPLESAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    tuple val(meta), path(reads), path(index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam") , emit: bam,   optional: true
    tuple val(meta), path("*.bai") , emit: index, optional: true
    tuple val(meta), path("*.cram"), emit: cram,  optional: true
    tuple val(meta), path("*.metrics.txt"), emit: metrics, optional: true
    tuple val("${task.process}"), val('picard'), eval("picard DownsampleSam --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:"), emit: versions_picard, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix    ?: "${reads.getExtension()}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard DownsampleSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if ("$reads" == "${prefix}.${suffix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    if (meta.downsampling_prob < 1) {
        """
        picard \\
            -Xmx${avail_mem}g \\
            DownsampleSam \\
            $args \\
            --CREATE_INDEX \\
            --INPUT $reads \\
            --OUTPUT ${prefix}.${suffix} \\
            $reference \\
            --METRICS_FILE ${prefix}.DownsampleSam.metrics.txt
        """
    } else {
        """
        if [ -L ${prefix}.${suffix} ]; then
            rm ${prefix}.${suffix}
        fi
        ln -s $reads ${prefix}.${suffix}
        ln -s $index ${prefix}.bai
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix    ?: "${reads.getExtension()}"
    if ("$reads" == "${prefix}.${suffix}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${suffix}
    touch ${prefix}.bai
    touch ${prefix}.DownsampleSam.metrics.txt
    """
}
