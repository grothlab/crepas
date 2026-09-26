process DEEPTOOLS_MULTIBIGWIGSUMMARY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0':
        'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(bigwigs) , val(labels)
    tuple val(meta2), path(blacklist)

    output:
    tuple val(meta), path("*.npz"), emit: matrix
    tuple val(meta), path("*.tsv"), emit: raw_counts, optional: true
    tuple val("${task.process}"), val('deeptools'), eval('multiBigwigSummary --version | sed "s/multiBigwigSummary //g"') , emit: versions_deeptools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "all_bigwig"
    def blacklist_cmd = blacklist ? "--blackListFileName ${blacklist}" : ""
    def outrawcounts_cmd = args.contains('--outRawCounts') ? '' : "--outRawCounts ${prefix}.tsv"
    def label  = labels ? "--labels ${labels.join(' ')}" : ''
    """
    multiBigwigSummary bins \\
        $args \\
        $label \\
        --bwfiles ${bigwigs.join(' ')} \\
        --numberOfProcessors $task.cpus \\
        --outFileName ${prefix}.bigwigSummary.npz \\
        ${outrawcounts_cmd} \\
        $blacklist_cmd
    """
    
    stub:
    def prefix = task.ext.prefix ?: "all_bigwig"
    """
    touch ${prefix}.bigwigSummary.npz
    touch ${prefix}.tsv
    """
}
