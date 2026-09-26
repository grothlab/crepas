process REPLISEQ_RTNORMALIZE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7a975bc28eee9e1375f221718022201e184d119495b209b0b133da55d7e3e7c2/data' :
        'community.wave.seqera.io/library/r-argparse_r-base_r-data.table_r-zoo:a9a33411d1a10128' }"

    input:
    tuple val(meta), path(counts), val(phases), val(breps)

    output:
    tuple val(meta), path("*.RT.raw.bedGraph")   , emit: raw
    tuple val(meta), path("*.RT.smooth.bedGraph"), emit: smooth
    tuple val(meta), path("*.RT.raw.covered.bedGraph")   , emit: raw_covered
    tuple val(meta), path("*.RT.smooth.covered.bedGraph"), emit: smooth_covered
    tuple val(meta), path("*.RT_index.raw.bedGraph")   , emit: rt_index_raw   , optional: true
    tuple val(meta), path("*.RT_index.smooth.bedGraph"), emit: rt_index_smooth, optional: true
    tuple val(meta), path("*.coverage.smooth.bedGraph"), emit: replicate_smooth, optional: true
    tuple val(meta), path("*.qc.txt")            , emit: qc
    tuple val(meta), path("rt_summary.tsv")      , emit: summary
    tuple val("${task.process}"), val('r-base'), eval("R --version | head -1 | sed 's/R version //; s/ .*//'"), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    repliseq_rtnormalize.R \\
        --counts ${counts} \\
        --phases ${phases.join(' ')} \\
        --breps ${breps.join(' ')} \\
        --outdir ./ \\
        --prefix ${prefix} \\
        --sample_name ${meta.id} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.RT.raw.bedGraph
    touch ${prefix}.RT.smooth.bedGraph
    touch ${prefix}.RT.raw.covered.bedGraph
    touch ${prefix}.RT.smooth.covered.bedGraph
    touch ${prefix}.coverage.smooth.bedGraph
    touch ${prefix}.qc.txt
    touch rt_summary.tsv
    """
}
