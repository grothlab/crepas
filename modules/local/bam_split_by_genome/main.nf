/*
 * Split a BAM file by genome (in this case a string appended to chromosome names)
 */
process BAM_SPLIT_BY_GENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)
    val endo_genome   // e.g., 'hg38'
    val exo_genome    // e.g., 'dm6'
    val reads_to_keep // 'endo' or 'exo'

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${reads_to_keep}"
    if (reads_to_keep == 'endo') {
        """
        samtools view \\
            --threads ${task.cpus} \\
            --with-header \\
            --no-PG \\
            ${bam} \\
            | grep -v -e "_${exo_genome}" -e '^@CO' -e '^@PG' \\
            > ${prefix}.${endo_genome}.sam

        samtools view \\
            --threads ${task.cpus} \\
            --bam \\
            ${prefix}.${endo_genome}.sam \\
            > ${prefix}.${endo_genome}.bam

        rm ${prefix}.${endo_genome}.sam
        """
    } else if (reads_to_keep == 'exo') {
        """
        samtools view \\
            --threads ${task.cpus} \\
            --with-header \\
            --no-PG \\
            ${bam} \\
            | grep -E "^@RG|_${exo_genome}" | grep -v -e '^@CO' -e '^@PG' \\
            | sed "s/_${exo_genome}//g" \\
            > ${prefix}.${exo_genome}.sam

        samtools view \\
            --threads ${task.cpus} \\
            --no-PG \\
            --bam \\
            ${prefix}.${exo_genome}.sam \\
            > ${prefix}.${exo_genome}.bam

        rm ${prefix}.${exo_genome}.sam
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${reads_to_keep}"
    """
    touch ${prefix}.${reads_to_keep}.bam
    """
}
