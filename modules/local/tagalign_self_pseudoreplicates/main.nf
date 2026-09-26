//
// Generating self-pseudoreplicates from TAGALIGN files
// Based on: https://github.com/ENCODE-DCC/chip-seq-pipeline2; step 2c
//
process TAGALIGN_SELF_PSEUDOREPLICATES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d6/d6696567851b54790cac3acd1e6744ca336f21d809e6fa3a1342cd3dae688198/data' :
        'community.wave.seqera.io/library/coreutils_gawk_openssl:8729727b42757e72' }"

    input:
    tuple val(meta), path(tagalign)

    output:
    tuple val(meta), path("*spr1.tagAlign"), emit: tagalign1
    tuple val(meta), path("*spr2.tagAlign"), emit: tagalign2
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), emit: versions_gawk, topic: versions
    tuple val("${task.process}"), val('coreutils'), eval("shuf --version | head -1 | sed 's/.*coreutils) //'"), emit: versions_coreutils, topic: versions
    tuple val("${task.process}"), val('openssl'), eval("openssl version | sed 's/^OpenSSL //; s/ .*//'"), emit: versions_openssl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        # Get total number of read pairs
        nlines=\$( wc -l < ${tagalign} )
        nlines=\$(( (nlines + 1) / 2 ))

        # Shuffle and split tagAlign file into 2 equal parts
        # Will produce prefix.shuf.split.tagAlign00 and prefix.shuf.split.tagAlign01
        shuf ${args} --random-source=<(openssl enc -aes-256-ctr -pass pass:\$(wc -c < ${tagalign}) -nosalt </dev/zero 2>/dev/null) ${tagalign} \\
            | split -d -l \${nlines} - ${prefix}.shuf.split.tagAlign

        # Convert reads into standard tagAlign file
        gzip -nc "${prefix}.shuf.split.tagAlign00" > ${prefix}.spr1.tagAlign
        rm "${prefix}.shuf.split.tagAlign00"
        gzip -nc "${prefix}.shuf.split.tagAlign01" > ${prefix}.spr2.tagAlign
        rm "${prefix}.shuf.split.tagAlign01"
        """
    } else {
        """
        # Make temporary fake BEDPE file from FINAL_TA_FILE
        sed 'N;s/\\n/\\t/' ${tagalign} > ${prefix}.tmp.bedpe

        # Get total number of read pairs
        nlines=\$( wc -l < ${prefix}.tmp.bedpe )
        nlines=\$(( (nlines + 1) / 2 ))

        # Shuffle and split BEDPE file into 2 equal parts
        shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:\$(wc -c < ${tagalign}) -nosalt </dev/zero 2>/dev/null) ${prefix}.tmp.bedpe \\
            | split -d -l \${nlines} - ${prefix}.shuf.split.tagAlign

        # Convert reads into standard tagAlign file
        awk 'BEGIN{OFS="\\t"}{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12}' "${prefix}.shuf.split.tagAlign00" \\
            > ${prefix}.spr1.tagAlign

        awk 'BEGIN{OFS="\\t"}{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12}' "${prefix}.shuf.split.tagAlign01" \\
            > ${prefix}.spr2.tagAlign

        rm "${prefix}.shuf.split.tagAlign00"
        rm "${prefix}.shuf.split.tagAlign01"
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.spr1.tagAlign
    touch  ${prefix}.spr2.tagAlign
    """
}
