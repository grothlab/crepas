process ROCCO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a0526a7f7adf6fa2e99db65b1da7b13e60dc50568378c5ed02a40f87901725a2/data' :
        'community.wave.seqera.io/library/bedtools_deeptoolsintervals_pybedtools_samtools_pruned:4e193f772646372f' }"

    input:
    tuple val(meta), path(bams_or_bws), path(bamlist), path(bailist)
    tuple val(meta2), path(chrom_sizes)
    tuple val(meta3), path(params_file)
    val effective_genome_size

    output:
    tuple val(meta), path("*.bed"),                             emit: bed
    tuple val(meta), path("*.counts.tsv"),    optional:true,    emit: counts
    tuple val(meta), path("*.narrowPeak"),    optional:true,    emit: narrow_peak
    tuple val(meta), path("*.mps"),           optional:true,    emit: model_mps
    tuple val("${task.process}"), val('rocco'), eval("pip show rocco | sed -n 's/^Version: //p'"), emit: versions_rocco, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args         ?: ""
    def prefix          = task.ext.prefix       ?: "${meta.id}"
    def samples         = bams_or_bws           ? "--input_files ${bams_or_bws.join(' ')}" : ""
    def bamlist_txt     = bamlist               ? "${bamlist.join(',')}" : ""
    def bamlist_arg     = bamlist               ? "--narrowPeak --bamlist_txt ${prefix}.bamlist.txt" : ""
    def sizes_arg       = chrom_sizes           ? "--chrom_sizes_file ${chrom_sizes}" : ""
    def egsize_arg      = effective_genome_size ? "--effective_genome_size ${effective_genome_size}" : ""
    def params_arg      = params_file           ? "--params ${params_file}" : ""
    """
    echo "${bamlist_txt}" | tr ',' '\\n' > ${prefix}.bamlist.txt

    rocco \\
        ${args} \\
        --threads ${task.cpus} \\
        --ecdf_proc ${task.cpus} \\
        --verbose \\
        ${samples} \\
        ${bamlist_arg} \\
        ${sizes_arg} \\
        ${params_arg} \\
        ${egsize_arg} \\
        --outfile ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.counts.tsv
    touch ${prefix}.narrowPeak
    touch ${prefix}.mps
    touch ${prefix}.bamlist.txt
    """
}
