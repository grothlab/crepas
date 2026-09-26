process REPLISEQ_CLASSIFY_GENES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8a9ad3c9a315b68403d671ecb126c853d599880e51e3a34d22a638eeb5cf9c5/data' :
        'community.wave.seqera.io/library/r-argparse_r-base_r-data.table_r-ggplot2_libblas:4a3be68aad9e7bbe' }"

    input:
    tuple val(meta), path(counts), val(phases), val(breps)

    output:
    tuple val(meta), path("*.gene_RT_class.tsv")               , emit: classes
    tuple val(meta), path("*.gene_RT_class.early.bed")         , emit: early
    tuple val(meta), path("*.gene_RT_class.late.bed")          , emit: late
    tuple val(meta), path("*.gene_RT_class.mid.bed")           , emit: mid, optional: true
    tuple val(meta), path("*.gene_RT_class.unclassified.bed")  , emit: unclassified
    tuple val(meta), path("*.gene_RT_class.qc.txt")            , emit: qc
    tuple val(meta), path("rt_gene_class_counts.tsv")          , emit: summary
    tuple val(meta), path("*.gene_RT_class.plots.pdf")         , emit: plots, optional: true
    tuple val(meta), path("*_mqc.json")                        , emit: mqc_box, optional: true
    tuple val("${task.process}"), val('r-base'), eval("R --version | head -1 | sed 's/R version //; s/ .*//'"), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    repliseq_classify_genes.R \\
        --counts ${counts} \\
        --phases ${phases.join(' ')} \\
        --breps ${breps.join(' ')} \\
        --sample_name ${meta.id} \\
        --outdir ./ \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gene_RT_class.tsv
    touch ${prefix}.gene_RT_class.early.bed
    touch ${prefix}.gene_RT_class.late.bed
    touch ${prefix}.gene_RT_class.unclassified.bed
    touch ${prefix}.gene_RT_class.qc.txt
    touch rt_gene_class_counts.tsv
    """
}
