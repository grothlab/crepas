process IDR_FILTER_THRESHOLD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198ab15844726d095c90e454aa9b1cf91b7e3517dd62451791c596e2a2229082/data' :
        'community.wave.seqera.io/library/coreutils_gawk_sed:e167f5dd848b5ae9' }"

    input:
    tuple val(meta), path(peak_file)
    val peak_type
    val idr_threshold

    output:
    tuple val(meta), path("${prefix}.${peak_type}"), emit: peaks
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), topic: versions, emit: versions_awk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.filtered_idr"
    def idr_threshold_arg = idr_threshold ?: 0.05
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    awk \\
        ${args} \\
        -v p=${idr_threshold_arg} 'BEGIN{OFS="\\t"; th=-log(p)/log(10)} \$12>=th {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' \\
        ${peak_file} \\
        | LC_COLLATE=C sort -T '.' \\
        | uniq \\
        | LC_COLLATE=C sort -T '.' -k7n,7n \\
        > ${prefix}.${peak_type}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.filtered_idr"
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    touch  ${prefix}.${peak_type}
    """
}
