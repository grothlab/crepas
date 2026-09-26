process TECOUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db7f82ebfe1c9f765a08b89ee98f2a9db9952b8bf6a0508a0305631683d4199c/data':
        'community.wave.seqera.io/library/tetranscripts_pigz:5c9ae6961179bdf5' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(genic_gtf_or_index)
    tuple val(meta3), path(te_gtf_or_index)
    val skip_gz

    output:
    tuple val(meta), path("*.cntTable*"),                        emit: counts
    tuple val("${task.process}"), val('TEcount'), eval("TEcount --version 2>&1 | sed 's/TEcount //g'"), emit: versions_tecount, topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //'"), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args         ?: ""
    def prefix          = task.ext.prefix       ?: "${meta.id}"
    def compress_cmd    = skip_gz ? "" : "pigz -f -p ${task.cpus} ${prefix}.cntTable"

    """
     TEcount \\
        ${args} \\
        --BAM ${bam} \\
        --GTF ${genic_gtf_or_index} \\
        --TE ${te_gtf_or_index} \\
        --project ${prefix} \\
        --outdir ./

    ${compress_cmd}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cmd = skip_gz ? "touch ${prefix}.cntTable" : "touch ${prefix}.cntTable.gz"

    """
    ${cmd}
    """
}
