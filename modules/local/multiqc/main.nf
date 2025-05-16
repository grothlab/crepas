process MULTIQC {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13a--pyhdfd78af_1':
        'quay.io/biocontainers/multiqc:1.13a--pyhdfd78af_1' }"

    input:
    path workflow_summary
    path multiqc_config
    path mqc_custom_config
    path logo

    path ('fastqc/*')
    path ('trimgalore/fastqc/*')
    path ('trimgalore/*')

    path ('alignment/library/*')
    path ('alignment/library/*')
    path ('alignment/library/*')

    path ('alignment/mergedLibrary/*')
    path ('alignment/mergedLibrary/*')
    path ('alignment/mergedLibrary/*')

    path ('alignment/mergedLibrary/preseq/*')

    path ('alignment/mergedLibrary/umi_dedup/*')
    path ('alignment/mergedLibrary/umi_dedup/*')
    path ('alignment/mergedLibrary/umi_dedup/*')
    path ('alignment/mergedLibrary/umi_dedup/*')

    path('alignment/mergedLibrary/picard_markduplicates/*')
    path('alignment/mergedLibrary/picard_markduplicates/*')
    path('alignment/mergedLibrary/picard_markduplicates/*')
    path('alignment/mergedLibrary/picard_markduplicates/*')

    path ('alignment/mergedLibrary/bam_filter/*')
    path ('alignment/mergedLibrary/bam_filter/*')
    path ('alignment/mergedLibrary/bam_filter/*')

    path ('alignment/mergedLibrary/spikein_split/*')
    path ('alignment/mergedLibrary/spikein_split/*')
    path ('alignment/mergedLibrary/spikein_split/*')

    path ('alignment/mergedLibrary/multimap/*')
    path ('alignment/mergedLibrary/multimap/*')
    path ('alignment/mergedLibrary/multimap/*')

    path ('alignment/mergedLibrary/picard_metrics/*')

    path ('deeptools/*')
    path ('deeptools/*')

    path ('phantompeakqualtools/*')
    path ('phantompeakqualtools/*')
    path ('phantompeakqualtools/*')
    path ('phantompeakqualtools/*')

    path ('macs3/peaks/*')
    path ('macs3/peaks/*')
    path ('macs3/annotation/*')
    path ('macs3/featurecounts/*')

    path ('deseq2/*')
    path ('deseq2/*')

    path ('genrich/peaks/*')
    path ('genrich/peaks/*')
    path ('genrich/annotation/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional: true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc \\
        -f \\
        $args \\
        $custom_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
