process MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/7362cb534743c8e26c5aca81c6b04148e43532438a02dbb45778fc69883fa7d0/data' :
        'community.wave.seqera.io/library/coreutils_gawk_r-base:9e0d80fa4a74052b' }"

    input:
    tuple val(meta), path(spp), path(rdata)
    path nsc_header
    path rsc_header
    path correlation_header

    output:
    tuple val(meta), path("*.spp_nsc_mqc.tsv")        , emit: nsc
    tuple val(meta), path("*.spp_rsc_mqc.tsv")        , emit: rsc
    tuple val(meta), path("*.spp_correlation_mqc.tsv"), emit: correlation
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^R version //; s/ .*\$//'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*\$//'"), topic: versions, emit: versions_gawk
    tuple val("${task.process}"), val('coreutils'), eval("cat --version | head -1 | sed 's/^cat (GNU coreutils) //'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp $correlation_header ${prefix}.spp_correlation_mqc.tsv
    Rscript --max-ppsize=500000 -e "load('$rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${prefix}.spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

    awk -v OFS='\t' '{print "${meta.id}", \$9}'  $spp | cat $nsc_header - > ${prefix}.spp_nsc_mqc.tsv
    awk -v OFS='\t' '{print "${meta.id}", \$10}' $spp | cat $rsc_header - > ${prefix}.spp_rsc_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.spp_nsc_mqc.tsv
    touch ${prefix}.spp_rsc_mqc.tsv
    touch ${prefix}.spp_correlation_mqc.tsv
    """
}
