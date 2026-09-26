process FRIP_SCORE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5cac195f18a1414900b9ee3c2830e8fc68b60684edbf3a5019c87a2594187a56/data':
        'community.wave.seqera.io/library/bedtools_samtools:458c8149e87aefd3' }"

    input:
    tuple val(meta), path(bam), path(peak)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | head -1 | sed 's/^samtools //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    READS_IN_PEAKS=\$(intersectBed -a $bam -b $peak $args | awk -F '\t' '{sum += \$NF} END {print sum}')
    samtools flagstat $bam > ${bam}.flagstat
    grep 'mapped (' ${bam}.flagstat | grep -v "primary" | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${prefix}", a/\$1}' > ${prefix}.FRiP.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.FRiP.txt
    """
}
