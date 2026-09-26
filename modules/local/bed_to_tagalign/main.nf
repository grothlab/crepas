//
// Converting BED to TAGALIGN format
// Based on: https://github.com/ENCODE-DCC/chip-seq-pipeline2; step 2a
//
process BED_TO_TAGALIGN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed4799c386bc676b388fb01df7e97792ed9dd5e84c8f20ef568feca5dc2cccdb/data' :
        'community.wave.seqera.io/library/gawk:5.4.0--0877e7a42fa88325' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.tagAlign"), emit: tagalign
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*//'"), emit: versions_awk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        awk \\
            ${args} \\
            'BEGIN{OFS="\\t"}{\$4="N";\$5="1000";print \$0}' \\
            ${bed} \\
            > ${prefix}.tagAlign
        """
    } else {
        """
        awk \\
            ${args} \\
            'BEGIN{OFS="\\t"}{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' \\
            ${bed} \\
            > ${prefix}.tagAlign
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.tagAlign
    """
}
