process DESEQ2_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c37c28e94cb08e049263519febbc9fb2a665247c53487c0bbbe1c78068f84a1/data' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-complexheatmap_bioconductor-deseq2_bioconductor-tximport_pruned:734b06a8a49c9ab7' }"

    input:
    tuple val(meta), path(counts)
    path deseq2_pca_header
    path deseq2_clustering_header

    output:
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*.rds"                , optional:true, emit: rds
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^.*R version //; s/ .*\$//'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-deseq2'), eval("Rscript -e 'library(DESeq2); cat(as.character(packageVersion(\"DESeq2\")))'"), topic: versions, emit: versions_deseq2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    deseq2_qc.r \\
        --count_file $counts \\
        --outdir ./ \\
        --outprefix $prefix \\
        --cores $task.cpus \\
        $args

    sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.pca.txt
    sed -i -e 's/DESeq2 /${meta.id} DESeq2 /g' tmp.pca.txt
    cat tmp.pca.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv

    sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.clustering.txt
    sed -i -e 's/DESeq2 /${meta.id} DESeq2 /g' tmp.clustering.txt
    cat tmp.clustering.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf
    touch ${prefix}.RData
    touch ${prefix}.rds
    touch ${prefix}.pca.vals.txt
    touch ${prefix}.pca.vals_mqc.tsv
    touch ${prefix}.sample.dists.txt
    touch ${prefix}.sample.dists_mqc.tsv
    touch ${prefix}.log
    touch size_factors
    """
}
