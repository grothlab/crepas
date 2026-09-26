/*
 * Prepare genome intervals for filtering by removing regions in blacklist file
 */
process GENOME_WHITELIST_REGIONS {
    tag "$sizes"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/2630bdd473cdd42149279090d3dd2a1c0e5d8a88af9346fff4c11ada3fc039ec/data':
        'community.wave.seqera.io/library/bedtools_gawk:3b83c7920e9b7f4a' }"

    input:
    tuple val(meta), path(sizes)
    tuple val(meta2), path(blacklist)

    output:
    tuple val(meta), path("*.bed")     , emit: bed
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${sizes.simpleName}.whitelist"
    if (blacklist) {
        """
        sortBed ${args} -i $blacklist -g $sizes | complementBed ${args2} -i stdin -g $sizes > ${prefix}.bed
        """
    } else {
        """
        awk ${args3} '{print \$1, '0' , \$2}' OFS='\t' $sizes > ${prefix}.bed
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${sizes.simpleName}.whitelist"
    """
    touch ${prefix}.bed
    """
}
