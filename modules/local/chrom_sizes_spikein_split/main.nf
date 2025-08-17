process CHROM_SIZES_SPIKEIN_SPLIT {
    tag "$sizes"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(sizes)
    val exo_genome_string
    val endo_genome_string

    output:
    tuple val(meta), path ("*.${endo_genome_string}.sizes") , emit: endo_sizes
    tuple val(meta), path ("*.${exo_genome_string}.sizes")  , emit: exo_sizes
    path  "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sizes.baseName}"
    """
    awk \\
        ${args} \\
        '!(\$1 ~ /_${exo_genome_string}\$/)' \\
        ${sizes} \\
        > ${prefix}.${endo_genome_string}.sizes

    awk \\
        ${args} \\
        '(\$1 ~ /_${exo_genome_string}\$/)' \\
        ${sizes} \\
        > ${prefix}.${exo_genome_string}.sizes


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${sizes.baseName}"
    """
    touch ${prefix}.${endo_genome_string}.sizes
    touch ${prefix}.${exo_genome_string}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk -Wversion 2>&1) | sed 's/^.*(GNU Awk) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
