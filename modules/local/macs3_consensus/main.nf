/*
 * Consensus peaks across samples, create boolean filtering file, SAF file for featureCounts
 */
process MACS3_CONSENSUS {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6a/6a6c6159c8410c98244fbe18fcd7d46e71cd67ddf5675b6d548eb839c36a3c11/data':
        'community.wave.seqera.io/library/bedtools_biopython_r-optparse_r-upsetr:71c4f71726f54101' }"

    input:
    tuple val(meta), path(peaks)
    val is_narrow_peak

    output:
    tuple val(meta), path("*.bed")          , emit: bed
    tuple val(meta), path("*.saf")          , emit: saf
    tuple val(meta), path("*.pdf")          , emit: pdf
    tuple val(meta), path("*.antibody.txt") , emit: txt
    tuple val(meta), path("*.boolean.txt")  , emit: boolean_txt
    tuple val(meta), path("*.intersect.txt"), emit: intersect_txt
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | head -1 | sed 's/^.*R version //; s/ .*\$//'"), topic: versions, emit: versions_rbase

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline
    def args1        = task.ext.args   ?: ''
    def args2        = task.ext.args2   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def peak_type    = is_narrow_peak  ? 'narrowPeak' : 'broadPeak'
    def mergecols    = is_narrow_peak  ? (2..10).join(',') : (2..9).join(',')
    def collapsecols = is_narrow_peak  ? (['collapse']*9).join(',') : (['collapse']*8).join(',')
    def expandparam  = is_narrow_peak  ? '--is_narrow_peak' : ''
    """
    sort -T '.' -k1,1 -k2,2n ${peaks.collect{peak -> peak.toString()}.sort().join(' ')} \\
        | mergeBed ${args1} -c $mergecols -o $collapsecols > ${prefix}.txt

    macs3_merged_expand.py \\
        ${prefix}.txt \\
        ${peaks.collect{peak -> peak.toString()}.sort().join(',').replaceAll("_peaks.${peak_type}","")} \\
        ${prefix}.boolean.txt \\
        ${args2} \\
        ${expandparam}

    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed

    echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${prefix}.boolean.txt >> ${prefix}.saf

    plot_peak_intersect.r -i ${prefix}.boolean.intersect.txt -o ${prefix}.boolean.intersect.plot.pdf

    echo "${prefix}.bed\t${meta.id}/${prefix}.bed" > ${prefix}.antibody.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.saf
    touch ${prefix}.boolean.intersect.plot.pdf
    touch ${prefix}.antibody.txt
    touch ${prefix}.boolean.txt
    touch ${prefix}.boolean.intersect.txt
    """

}
