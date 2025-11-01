/*
 * Split a BAM file by genome (in this case a string appended to chromosome names)
 */
process BAM_SPLIT_BY_GENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)
    val filter_genome_string
    val keep_genome_string
    val filter_out

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def grep_command     = filter_out ? "grep -v" : "grep"
    def reheader_command = filter_out ? "samtools reheader -c 'grep -v \"_${filter_genome_string}\" -e ^@CO -e ^@PG' ${prefix}.${keep_genome_string}.sam > ${prefix}.${keep_genome_string}.tmp.sam && mv ${prefix}.${keep_genome_string}.tmp.sam ${prefix}.${keep_genome_string}.sam" : ''
    """
    samtools view \\
        -h $bam | \\
            $grep_command "_${filter_genome_string}" > ${prefix}.${keep_genome_string}.sam

    # reaheader if filter_out is set
    $reheader_command

    samtools view -b ${prefix}.${keep_genome_string}.sam \\
        > ${prefix}.${keep_genome_string}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
