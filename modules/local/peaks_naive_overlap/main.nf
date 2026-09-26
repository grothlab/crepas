process PEAKS_NAIVE_OVERLAP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4e1dbf5b874984146cb2b5ec7a90f3dd418a7b6015b5be1f172a1fb767b04002/data' :
        'community.wave.seqera.io/library/bedtools_gawk:996a038ce58c4d54' }"

    input:
    tuple val(meta), path(sample_peaks), path(peak_list), val(threshold)
    val peak_type

    output:
    tuple val(meta), path("*.${peak_type}"), emit: peak_overlap
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.naive_overlap"
    def threshold_arg = threshold ?: 0.05
    if (sample_peaks.toList().size < 2) {
        log.error "[ERROR] Naive peak overlapping needs at least two replicates only one provided."
    }
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    intersectBed ${args} -wo -a ${peak_list} -b ${sample_peaks[0]} | \\
    awk ${args2} 'BEGIN{FS="\\t";OFS="\\t"; th="${threshold_arg}"}{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= th) || (\$21/s2 >= th)) {print \$0}}' | cut -f 1-10 | LC_COLLATE=C sort -T '.' | uniq | \\
    intersectBed ${args3} -wo -a stdin -b ${sample_peaks[1]} | \\
    awk ${args4} 'BEGIN{FS="\\t";OFS="\\t"; th="${threshold_arg}"}{s1=\$3-\$2; s2=\$13-\$12; if ((\$21/s1 >= th) || (\$21/s2 >= th)) {print \$0}}' | cut -f 1-10 | LC_COLLATE=C sort -T '.' | uniq \\
    > ${prefix}.${peak_type}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.naive_overlap"
    if (sample_peaks.toList().size < 2) {
        log.error "[ERROR] Naive peak overlapping needs at least two replicates only one provided."
    }
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    touch "${prefix}.${peak_type}"
    """
}
