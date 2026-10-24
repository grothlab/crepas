process TE_COUNTS_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6d/6dd2a477d4dcf63096bf329cf7b6c714ab6e13b56f816f1ec22610fcf49a4aa4/data' :
        'community.wave.seqera.io/library/bioconductor-deseq2_python_r-argparse_r-base_pruned:afff681711e1afd0' }"

    input:
    tuple val(meta), path(counts)
    path pca_header
    path clustering_header

    output:
    tuple val(meta), path("*.plots.pdf")               , emit: pdf
    tuple val(meta), path("*.counts.tsv.gz")           , emit: counts
    tuple val(meta), path("*.normalized_counts.tsv.gz"), emit: normalized_counts
    tuple val(meta), path("*.size_factors.tsv")        , emit: size_factors
    tuple val(meta), path("*.pca.vals.txt")            , emit: pca_txt
    tuple val(meta), path("*.pca.vals_mqc.tsv")        , emit: pca_multiqc
    tuple val(meta), path("*.sample.dists.txt")        , emit: dists_txt
    tuple val(meta), path("*.sample.dists_mqc.tsv")    , emit: dists_multiqc
    tuple val(meta), path("*.dds.rds")                 , emit: dds
    tuple val(meta), path("*.{vst,rlog}.rds")          , emit: transformed
    tuple val(meta), path("*.R_sessionInfo.log")       , emit: log
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^.*R version //; s/ .*\$//'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-deseq2'), eval("Rscript -e 'library(DESeq2); cat(as.character(packageVersion(\"DESeq2\")))'"), topic: versions, emit: versions_deseq2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // MultiQC section IDs must be unique across all TE counting QC runs
    def mqc_id = prefix.replaceAll(/[^A-Za-z0-9_]/, '_')
    """
    te_counts_qc.r \\
        --outprefix ${prefix} \\
        ${args} \\
        ${counts}

    sed -e 's/te_counts_pca/te_counts_pca_${mqc_id}/g' -e 's/TE counts PCA/${prefix}: TE counts PCA/' ${pca_header} > ${prefix}.pca.vals_mqc.tsv
    cat ${prefix}.pca.vals.txt >> ${prefix}.pca.vals_mqc.tsv

    sed -e 's/te_counts_clustering/te_counts_clustering_${mqc_id}/g' -e 's/TE counts sample similarity/${prefix}: TE counts sample similarity/' ${clustering_header} > ${prefix}.sample.dists_mqc.tsv
    cat ${prefix}.sample.dists.txt >> ${prefix}.sample.dists_mqc.tsv
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def transform = args.contains('--rlog') ? 'rlog' : 'vst'
    """
    touch ${prefix}.plots.pdf
    python3 -c "import gzip; [gzip.open(f, 'wb').close() for f in ('${prefix}.counts.tsv.gz', '${prefix}.normalized_counts.tsv.gz')]"
    touch ${prefix}.size_factors.tsv
    touch ${prefix}.pca.vals.txt
    touch ${prefix}.pca.vals_mqc.tsv
    touch ${prefix}.sample.dists.txt
    touch ${prefix}.sample.dists_mqc.tsv
    touch ${prefix}.dds.rds
    touch ${prefix}.${transform}.rds
    touch ${prefix}.R_sessionInfo.log
    """
}
