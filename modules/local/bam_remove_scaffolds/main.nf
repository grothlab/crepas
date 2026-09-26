process BAM_REMOVE_SCAFFOLDS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/3388ed16c02ec833ab56da2cb4a6d1fbf2266460e1a04754692fe5c0716cf3e0/data' :
        'community.wave.seqera.io/library/htslib_samtools_gawk:235fd990e554cf33' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), emit: versions_awk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.flTsfds"
    // Split the number of extra threads between the two samtools commands
    def nthreads1 = Math.max(1, (task.cpus / 2) as int)
    def nthreads2 = Math.max(0, (task.cpus - nthreads1) as int)
    """
    samtools view \\
        --threads ${nthreads1} \\
        --with-header \\
        --no-PG \\
        ${bam} \\
        | awk 'BEGIN{OFS=FS="\t"} /^@CO/ {next} /^@PG/ {next} (\$2 ~ /SN:.*(\\.|_random)/) {next} (\$3 ~ /\\./ || \$3 ~ /_random/) {next} {print}' \\
        | samtools view \\
            --threads ${nthreads2} \\
            --bam \\
            - \\
            > ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.flTsfds"
    """
    touch ${prefix}.bam
    """
}
