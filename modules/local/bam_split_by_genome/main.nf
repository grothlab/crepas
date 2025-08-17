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
    val endo_genome   // e.g., 'hg38'
    val exo_genome    // e.g., 'dm6'
    val reads_to_keep // 'endo' or 'exo'

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${reads_to_keep}"
    // Split the number of extra threads between the two samtools commands
    def nthreads1 = Math.max(1, (task.cpus / 2) as int)
    def nthreads2 = Math.max(0, (task.cpus - nthreads1) as int)
    if (reads_to_keep == 'endo') {
        """
        samtools view \\
            --threads ${nthreads1} \\
            --with-header \\
            --no-PG \\
            ${bam} \\
            | grep -v -e "_${exo_genome}" -e '^@CO' -e '^@PG' \\
            | samtools view \\
                --threads ${nthreads2} \\
                --bam \\
                - \\
                > ${prefix}.${endo_genome}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else if (reads_to_keep == 'exo') {
        """
        samtools view \\
            --threads ${nthreads1} \\
            --with-header \\
            --no-PG \\
            ${bam} \\
            | grep -E "^@RG|_${exo_genome}" | grep -v -e '^@CO' -e '^@PG' \\
            | samtools view \\
                --threads ${nthreads2} \\
                --no-PG \\
                --bam \\
                - \\
                > ${prefix}.${exo_genome}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${reads_to_keep}"
    """
    touch ${prefix}.${reads_to_keep}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
}