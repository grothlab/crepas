process BAM_FLAGSTAT_MAPPED {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e4733e27a0c96155d1405ddb1bcb173c68367c0fc396faac8bf7ad1f95436727/data'
        : 'community.wave.seqera.io/library/coreutils_gawk_ripgrep_sed:0be2a714bf896bff'}"

    input:
    tuple val(meta), path(flagstat)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*//'"), emit: versions_awk, topic: versions
    tuple val("${task.process}"), val('ripgrep'), eval("rg --version | head -1 | sed 's/^ripgrep //; s/ (.*//'"), emit: versions_ripgrep, topic: versions
    tuple val("${task.process}"), val('coreutils'), eval("env echo --version | head -1 | sed 's/^echo (GNU coreutils) //'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    TOTAL_READS=\$(rg ${args} '[0-9] mapped \\(' ${flagstat} | awk ${args2} '{print \$1}')
    echo ${args3} \$TOTAL_READS > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.txt
    """
}
