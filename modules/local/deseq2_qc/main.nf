process DESEQ2_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-complexheatmap_bioconductor-deseq2_bioconductor-tximport_pruned:95aef8e10f40bd25' :
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
    path "versions.yml"         , emit: versions

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

    sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
    sed -i -e 's/DESeq2 /${meta.id} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv

    sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
    sed -i -e 's/DESeq2 /${meta.id} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
