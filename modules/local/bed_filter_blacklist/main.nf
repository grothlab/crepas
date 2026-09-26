process BED_FILTER_BLACKLIST {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4e1dbf5b874984146cb2b5ec7a90f3dd418a7b6015b5be1f172a1fb767b04002/data' :
        'community.wave.seqera.io/library/bedtools_gawk:996a038ce58c4d54' }"

    input:
    tuple val(meta), path(peaks)
    tuple val(meta2), path(blacklist)
    val filter_chr
    val peak_type
    val max_score


    output:
    tuple val(meta), path("*.${peak_type}"), emit: peaks
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.flTbl"
    def zcat_arg = peaks.name.endsWith('.gz') ? 'zcat' : 'cat'
    def max_score_arg = max_score ?: 1000
    def filter_chr_arg = filter_chr ? "| { grep -P 'chr[\\dXY]+[ \\t]' || true; }" : ''
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    ${zcat_arg} ${peaks} \\
        | awk 'BEGIN{OFS="\\t"} {\$2=int(\$2); \$3=int(\$3); if (\$2<0) \$2=0; print \$0}' \\
         > ${prefix}.tmp.${peak_type}

    bedtools intersect \\
        ${args} \\
        -v \\
        -a ${prefix}.tmp.${peak_type} \\
        -b ${blacklist} \\
        ${filter_chr_arg} \\
        | awk ${args2} 'BEGIN{OFS="\\t"; ms="${max_score_arg}"} {if (\$5>ms) \$5=ms; print \$0}' \\
        > ${prefix}.${peak_type}

    rm ${prefix}.tmp.${peak_type}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.flTbl"
    def peak_types = ['narrowPeak', 'broadPeak', 'bed']
    if (!peak_types.contains(peak_type)) {
        log.error "[ERROR] Invalid option: '${peak_type}'. Valid options for 'peak_type': ${peak_types.join(', ')}."
    }
    """
    touch "${prefix}.${peak_type}"
    """
}
