process TAGALIGN_FRIP_SCORE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4e1dbf5b874984146cb2b5ec7a90f3dd418a7b6015b5be1f172a1fb767b04002/data' :
        'community.wave.seqera.io/library/bedtools_gawk:996a038ce58c4d54' }"

    input:
    tuple val(meta), path(tagalign), path(ccscores), path(peaks)
    tuple val(meta2), path(chromsizes)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.txt"), emit: frip
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.FRiP"
    def cat_arg = tagalign.name.endsWith('.gz') ? 'zcat' : 'cat'
    """
    # Get estimated fragment length from cross-correlation scores.
    FRAGLEN=\$(awk 'NR==1 {print int(\$3); exit}' ${ccscores})
    HALF_FRAGLEN=\$(( (FRAGLEN + 1) / 2 ))

    bedtools slop \\
        ${args} \\
        -i ${tagalign} \\
        -g ${chromsizes} \\
        -s \\
        -l -\${HALF_FRAGLEN} \\
        -r \${HALF_FRAGLEN} \\
        | awk 'BEGIN{OFS="\\t"} {if (\$2>=0 && \$3>=0 && \$2<=\$3) print \$0}' \\
        > ${prefix}.bed

    overlap_count=\$(bedtools intersect ${args2} -a ${prefix}.bed -b ${peaks} -wa -u | wc -l)
    total_count=\$(${cat_arg} ${tagalign} | wc -l)
    awk -v n="\${overlap_count}" -v d="\${total_count}" 'BEGIN{if (d==0) print 0; else print n/d}' > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.FRiP"
    """
    touch "${prefix}.bed"
    touch "${prefix}.txt"
    """
}
