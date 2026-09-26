process HR_REPLISEQ_PLOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b29c6de0ed5b02fcfa32795a17dff4aca14501adafcf995e272722d44a4722c9/data' :
        'community.wave.seqera.io/library/matplotlib-base_numpy_pandas_python:c233ff6bba50fb28' }"

    input:
    tuple val(meta), path(array), path(partition), path(partition_disjoint), path(partition_flanked), path(speeds), path(cluster_rank), path(el_track), path(calls)

    output:
    tuple val(meta), path("*.hr_repliseq.plots.pdf"), emit: plots
    tuple val(meta), path("*_mqc.png")              , emit: mqc_heatmap, optional: true
    tuple val(meta), path("*_mqc.json")             , emit: mqc_sizes  , optional: true
    tuple val("${task.process}"), val('python')    , eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('matplotlib'), eval("python3 -c 'import matplotlib; print(matplotlib.__version__)'"), emit: versions_matplotlib, topic: versions
    tuple val("${task.process}"), val('pandas')    , eval("python3 -c 'import pandas; print(pandas.__version__)'"), emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Drawn above the heatmap when the same sample also has an early/late track; the two
    // designs are independent, so most samples have only one of them.
    def el_arg = el_track ? "--el_track ${el_track}" : ''

    """
    hr_repliseq_plot.py \\
        --array ${array} \\
        --partition ${partition} \\
        --partition_disjoint ${partition_disjoint} \\
        --partition_flanked ${partition_flanked} \\
        --speeds ${speeds} \\
        --cluster_rank ${cluster_rank} \\
        ${el_arg} \\
        --calls ${calls.join(' ')} \\
        --sample_name ${meta.id} \\
        --outdir ./ \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hr_repliseq.plots.pdf
    """
}
