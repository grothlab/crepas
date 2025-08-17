
process MACE_PREPROCESSOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mace:1.2--py27he7e273a_2':
        'biocontainers/mace:1.2_cv1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(chrom_sizes)

    output:
    tuple val(meta), path("*_Forward.wig"), emit: forward_wig
    tuple val(meta), path("*_Reverse.wig"), emit: reverse_wig
    path  "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def treatment  = bam ? "--inputFile ${bam.join(',')}" : ""
    """
    preprocessor.py \\
        ${args} \\
        ${treatment} \\
        --chromSize ${chrom_sizes} \\
        --outPrefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mace_preprocessor: \$(preprocessor.py --version | sed -e "s/mace //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_Forward.wig
    touch ${prefix}_Reverse.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mace_preprocessor: \$(preprocessor.py --version | sed -e "s/mace //g")
    END_VERSIONS
    """
}
