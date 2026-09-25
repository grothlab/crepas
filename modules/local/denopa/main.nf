process DENOPA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://hepingshiming2007/denopa:latest' :
        'docker.io/hepingshiming2007/denopa:latest' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*_ARERs.txt"),           emit: arers
    tuple val(meta), path("*_NFR.txt"),             emit: nfrs
    tuple val(meta), path("*_nucleosomes.txt"),     emit: nucleosomes
    tuple val(meta), path("*_candidates.pkl"),      emit: candidates, optional: true
    tuple val(meta), path("*_frag_len.pkl"),        emit: frag_len,   optional: true
    tuple val(meta), path("*_pileup_signal.hdf"),   emit: pileup,     optional: true
    tuple val(meta), path("*_smooth.hdf"),          emit: smooth,     optional: true
    tuple val("${task.process}"), val('denopa'), eval("/root/miniconda3/bin/python3 -c \"import pkg_resources; print(pkg_resources.get_distribution('deNOPA').version)\""), topic: versions, emit: versions_denopa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args         ?: ""
    def prefix          = task.ext.prefix       ?: "${meta.id}"
    // container's default PATH omits /root/miniconda3/bin, where denopa is installed
    """
    /root/miniconda3/bin/denopa \\
        ${args} \\
        --proc ${task.cpus} \\
        --input ${bam} \\
        --output . \\
        --name ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ARERs.txt
    touch ${prefix}_NFR.txt
    touch ${prefix}_nucleosomes.txt
    """
}
