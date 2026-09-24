process REPLISEQ_RT_DOMAINS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f552061e675d909879ccaa1f1a93ba18aea3a92277a0e41c996cf8a3fdb8ce4/data' :
        'community.wave.seqera.io/library/bioconductor-dnacopy_libblas_r-argparse_r-base_pruned:979f7a7b8c7daeaf' }"

    input:
    tuple val(meta), path(bedgraph)

    output:
    tuple val(meta), path("*.RT_domains.bed")            , emit: domains
    tuple val(meta), path("*.RT_domains.early.bed")      , emit: early
    tuple val(meta), path("*.RT_domains.late.bed")       , emit: late
    tuple val(meta), path("*.RT_domains.mid.bed")        , emit: mid, optional: true
    tuple val(meta), path("*.RT_domains.qc.txt")         , emit: qc
    tuple val(meta), path("*.RT_domains.plots.pdf")      , emit: plots, optional: true
    tuple val(meta), path("*_mqc.json")                  , emit: mqc_box, optional: true
    tuple val("${task.process}"), val('r-base'), eval("R --version | head -1 | sed 's/R version //; s/ .*//'"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('bioconductor-dnacopy'), eval("Rscript -e 'cat(as.character(packageVersion(\"DNAcopy\")))'"), emit: versions_dnacopy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    repliseq_rt_domains.R \\
        --bedgraph ${bedgraph} \\
        --sample_name ${meta.id} \\
        --outdir ./ \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.RT_domains.bed
    touch ${prefix}.RT_domains.early.bed
    touch ${prefix}.RT_domains.late.bed
    touch ${prefix}.RT_domains.qc.txt
    """
}
