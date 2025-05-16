process CONSENRICH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f4dd21bd3f68dfe71b6d6a8624340b8735b43ef0c63b3a95031b4eb93403790/data' :
        'community.wave.seqera.io/library/bedtools_deeptools_pybedtools_pybigwig_pruned:01282f183573fac0' }"

    input:
    tuple val(meta), path(treatment_bam), path(treatment_bai), path(control_bam), path(control_bai)
    path chrom_sizes
    path blacklist
    path sparsebed
    path active_regions

    output:
    tuple val(meta), path("${prefix}.consenrich_output.tsv"),           emit: results
    tuple val(meta), path("*consenrich_signal_track*.bw"),      optional:true, emit: signal_track
    tuple val(meta), path("*consenrich_residuals_track*.bw"),   optional:true, emit: residuals_track
    tuple val(meta), path("*consenrich_eratio_track*.bw"),      optional:true, emit: eratio_track
    tuple val(meta), path("${prefix}_*.tsv.gz"),                        optional:true, emit: gain_log
    tuple val(meta), path("consenrich_${prefix}_args.json"),            optional:true, emit: args_json
    tuple val(meta), path("*.npz"),                                     optional:true, emit: matrix
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args     ?: ""
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    def treatment       = treatment_bam     ? "--bam_files ${treatment_bam.join(' ')}" : ""
    def control         = control_bam       ? "--control_files ${control_bam.join(' ')}" : ""
    def blacklist       = blacklist         ? "--blacklist_file $blacklist" : ""
    def sparsebed       = sparsebed         ? "--sparsebed $sparsebed" : ""
    def active_regions  = active_regions    ? "--active_regions $active_regions" : ""
    def no_sparsebed    = !sparsebed && !active_regions ? "--no_sparsebed" : ""
    def single_end      = meta.single_end   ? "--single_end" : ""
    """
    consenrich \\
        $args \\
        --threads $task.cpus \\
        --n_processes $task.cpus \\
        --experiment_id $prefix \\
        $treatment \\
        $control \\
        --sizes_file $chrom_sizes \\
        $blacklist \\
        $sparsebed \\
        $active_regions \\
        $no_sparsebed \\
        $single_end \\
        --output_file ${prefix}.consenrich_output.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Consenrich: \$(echo \$(Consenrich -h 2>&1) | sed 's/^Consenrich, version //; s/ .*\$//')
    END_VERSIONS
    """
}

// TODO: version parsing
