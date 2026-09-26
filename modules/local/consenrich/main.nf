process CONSENRICH {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f4dd21bd3f68dfe71b6d6a8624340b8735b43ef0c63b3a95031b4eb93403790/data'
        : 'community.wave.seqera.io/library/bedtools_deeptools_pybedtools_pybigwig_pruned:01282f183573fac0'}"

    input:
    tuple val(meta), path(treatment_bam), path(treatment_bai), path(control_bam), path(control_bai)
    tuple val(meta2), path(chrom_sizes)
    tuple val(meta3), path(blacklist)
    tuple val(meta4), path(sparsebed)
    tuple val(meta5), path(active_regions)

    output:
    tuple val(meta), path("*consenrich_output.tsv"), emit: results
    tuple val(meta), path("*consenrich_signal_track*.bw"), optional: true, emit: signal_track
    tuple val(meta), path("*consenrich_residuals_track*.bw"), optional: true, emit: residuals_track
    tuple val(meta), path("*consenrich_eratio_track*.bw"), optional: true, emit: eratio_track
    tuple val(meta), path("${prefix}_*.tsv.gz"), optional: true, emit: gain_log
    tuple val(meta), path("consenrich_${prefix}_args.json"), optional: true, emit: args_json
    tuple val(meta), path("*.npz"), optional: true, emit: matrix
    tuple val("${task.process}"), val('consenrich'), eval("pip show consenrich 2>/dev/null | awk '/^Version:/{print \$2}'"), topic: versions, emit: versions_consenrich

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"

    def treatment = treatment_bam ? "--bam_files ${treatment_bam.join(' ')}" : ""
    def control = control_bam ? "--control_files ${control_bam.join(' ')}" : ""
    def blacklist_arg = blacklist ? "--blacklist_file ${blacklist}" : ""
    def sparsebed_arg = sparsebed ? "--sparsebed ${sparsebed}" : ""
    def active_regions_arg = active_regions ? "--active_regions ${active_regions}" : ""
    def no_sparsebed_arg = !sparsebed && !active_regions ? "--no_sparsebed" : ""
    def single_end_arg = meta.single_end ? "--single_end" : ""
    """
    consenrich \\
        ${args} \\
        --threads ${task.cpus} \\
        --n_processes ${task.cpus} \\
        --experiment_id ${prefix} \\
        ${treatment} \\
        ${control} \\
        --sizes_file ${chrom_sizes} \\
        ${blacklist_arg} \\
        ${sparsebed_arg} \\
        ${active_regions_arg} \\
        ${no_sparsebed_arg} \\
        ${single_end_arg}  \\
        --output_file ${prefix}.consenrich_output.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consenrich_output.tsv
    touch ${prefix}.consenrich_signal_track.bw
    touch ${prefix}.consenrich_residuals_track.bw
    touch ${prefix}.consenrich_eratio_track.bw
    touch ${prefix}_consenrich_gain_log.tsv.gz
    touch consenrich_${prefix}_args.json
    touch ${prefix}_consenrich_matrix.npz
    """
}
