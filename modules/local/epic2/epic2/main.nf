process EPIC2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f41b3fd2e4992c1a7012a86c5ec7f7c58db2cc262188117eedfb5d078de668f/data'
        : 'community.wave.seqera.io/library/epic2_pyranges_gawk:ae3474e51d3dfc5e'}"

    input:
    tuple val(meta), path(treatment_bam), path(treatment_bai), path(control_bam), path(control_bai)
    tuple val(meta2), path(chrom_sizes)
    val effective_genome_fraction

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.diffusePeak"), emit: peak
    tuple val("${task.process}"), val('epic2'), eval("epic2 --version | sed -e 's/epic2 //g'"), topic: versions, emit: versions_epic2
    tuple val("${task.process}"), val('gawk'), eval("gawk --version 2>&1 | sed 's/^.*GNU Awk //; s/, .*\$//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chromsizes = chrom_sizes ? "--chromsizes ${chrom_sizes}" : ""
    def egf = effective_genome_fraction ? "--effective-genome-fraction ${effective_genome_fraction}" : ""
    def format = meta.single_end ? '' : (args.contains('--guess-bampe') ? '' : '--guess-bampe')
    def treatment = treatment_bam ? "--treatment ${treatment_bam.join(' ')}" : ""
    def control = control_bam ? "--control ${control_bam.join(' ')}" : ""
    """
    epic2 \\
        ${args} \\
        ${format} \\
        ${chromsizes} \\
        ${egf} \\
        ${treatment} \\
        ${control} \\
        --output ${prefix}.bed

    # Convert to HOMER-compatible BED format
    awk \\
        ${args2} \\
        -F'\t' \\
        'NR==1 {next} \\
        { printf "%s\\t%s\\t%s\\tpeak_%d\\t%s\\t%s\\t%s", \$1, \$2, \$3, NR-1, \$5, \$6, \$4; \\
        for (i=7; i<=NF; i++) printf "\\t%s", \$i; \\
        printf "\\n" }' \\
        ${prefix}.bed \\
        > ${prefix}.diffusePeak
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.diffusePeak
    """
}
