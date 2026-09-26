process HR_REPLISEQ_CALL_FEATURES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e37d089fd2a5e75a8b4baf061893d10ce5744a799786cd752ce74586a8bb24ba/data' :
        'community.wave.seqera.io/library/numpy_pandas_python_scikit-learn:7c7a45cf53383448' }"

    input:
    tuple val(meta), path(array)

    output:
    tuple val(meta), path("*.hr_IZ.bed")          , emit: iz
    tuple val(meta), path("*.hr_TTR.bed")         , emit: ttr
    tuple val(meta), path("*.hr_breakage.overlap.bed")    , emit: breakage
    tuple val(meta), path("*.hr_breakage.disjoint.bed")  , emit: breakage_disjoint
    tuple val(meta), path("*.hr_breakage.flanked.bed")   , emit: breakage_flanked
    tuple val(meta), path("*.hr_termination.bed") , emit: termination
    tuple val(meta), path("*.hr_CTR.bed")         , emit: ctr
    tuple val(meta), path("*.hr_partition.overlap.bed")   , emit: partition
    tuple val(meta), path("*.hr_partition.disjoint.bed") , emit: partition_disjoint
    tuple val(meta), path("*.hr_partition.flanked.bed")  , emit: partition_flanked
    tuple val(meta), path("*.hr_TTR.speed.tsv")          , emit: ttr_speed
    tuple val(meta), path("*.hr_cluster_rank.tsv")       , emit: cluster_rank
    tuple val(meta), path("*.hr_features.qc.txt") , emit: qc
    tuple val(meta), path("hr_repliseq_features_table.tsv"), emit: summary
    tuple val("${task.process}"), val('python')      , eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('numpy')       , eval("python3 -c 'import numpy; print(numpy.__version__)'"), emit: versions_numpy, topic: versions
    tuple val("${task.process}"), val('pandas')      , eval("python3 -c 'import pandas; print(pandas.__version__)'"), emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('scikit-learn'), eval("python3 -c 'import sklearn; print(sklearn.__version__)'"), emit: versions_sklearn, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hr_repliseq_call_features.py \\
        --array ${array} \\
        --sample_name ${meta.id} \\
        --outdir ./ \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hr_IZ.bed
    touch ${prefix}.hr_TTR.bed
    touch ${prefix}.hr_breakage.overlap.bed
    touch ${prefix}.hr_breakage.disjoint.bed
    touch ${prefix}.hr_breakage.flanked.bed
    touch ${prefix}.hr_termination.bed
    touch ${prefix}.hr_CTR.bed
    touch ${prefix}.hr_partition.overlap.bed
    touch ${prefix}.hr_partition.disjoint.bed
    touch ${prefix}.hr_partition.flanked.bed
    touch ${prefix}.hr_TTR.speed.tsv
    touch ${prefix}.hr_cluster_rank.tsv
    touch ${prefix}.hr_features.qc.txt
    touch hr_repliseq_features_table.tsv
    """
}
