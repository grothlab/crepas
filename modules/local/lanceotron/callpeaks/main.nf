process LANCEOTRON_CALLPEAKS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/lanceotron:1.2.7--pyhdfd78af_0'
        : 'biocontainers/lanceotron:1.2.7--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(bigwig)

    output:
    tuple val(meta), path("*.lanceotron.candidates.tsv"), emit: candidates
    tuple val("${task.process}"), val('lanceotron'), eval('python -c "import importlib.metadata as m; print(m.version(\'lanceotron\'))"'), emit: versions_lanceotron, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export TF_NUM_INTRAOP_THREADS=${task.cpus}
    export TF_NUM_INTEROP_THREADS=1

    lanceotron \\
        callPeaks \\
        ${bigwig} \\
        ${args} \\
        --folder lanceotron_out

    # LanceOtron names its table after the bigWig, up to the first "."; it writes none when it finds no
    # candidate peaks
    table=\$(find lanceotron_out -name '*_L-tron.bed' | head -n 1)
    if [ -n "\$table" ]; then
        mv "\$table" ${prefix}.lanceotron.candidates.tsv
    else
        touch ${prefix}.lanceotron.candidates.tsv
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lanceotron.candidates.tsv
    """
}
