process HR_REPLISEQ_MAKE_ARRAY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e37d089fd2a5e75a8b4baf061893d10ce5744a799786cd752ce74586a8bb24ba/data' :
        'community.wave.seqera.io/library/numpy_pandas_python_scikit-learn:7c7a45cf53383448' }"

    input:
    tuple val(meta), path(counts), val(fractions)

    output:
    tuple val(meta), path("*.hr_array.csv")   , emit: array
    tuple val(meta), path("*.hr_array.qc.txt"), emit: qc
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('numpy') , eval("python3 -c 'import numpy; print(numpy.__version__)'"), emit: versions_numpy, topic: versions
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hr_repliseq_make_array.py \\
        --counts ${counts} \\
        --labels ${fractions.join(' ')} \\
        --qc ${prefix}.hr_array.qc.txt \\
        --output ${prefix}.hr_array.csv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hr_array.csv
    touch ${prefix}.hr_array.qc.txt
    """
}
