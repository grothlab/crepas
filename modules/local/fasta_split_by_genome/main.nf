process FASTA_SPLIT_BY_GENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data' :
        'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751' }"

    input:
    tuple val(meta), path(fasta)
    val exo_genome_string
    val endo_genome_string

    output:
    tuple val(meta), path ("*.${endo_genome_string}.fasta"), emit: endo_fasta
    tuple val(meta), path ("*.${exo_genome_string}.fasta")  , emit: exo_fasta
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    // Split the number of threads between the two seqkit commands running concurrently in the exo pipe
    def nthreads1 = Math.max(1, (task.cpus / 2) as int)
    def nthreads2 = Math.max(1, (task.cpus - nthreads1) as int)
    """
    seqkit grep \\
        ${args} \\
        --threads ${task.cpus} \\
        -v -r -p "_${exo_genome_string}\$" \\
        ${fasta} \\
        -o ${prefix}.${endo_genome_string}.fasta

    seqkit grep \\
        ${args} \\
        --threads ${nthreads1} \\
        -r -p "_${exo_genome_string}\$" \\
        ${fasta} \\
        | seqkit replace \\
            ${args2} \\
            --threads ${nthreads2} \\
            -p "_${exo_genome_string}\$" -r '' \\
            -o ${prefix}.${exo_genome_string}.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.${endo_genome_string}.fasta
    touch ${prefix}.${exo_genome_string}.fasta
    """
}
