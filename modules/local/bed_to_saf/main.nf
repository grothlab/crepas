process BED_TO_SAF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198ab15844726d095c90e454aa9b1cf91b7e3517dd62451791c596e2a2229082/data' :
        'community.wave.seqera.io/library/coreutils_gawk_sed:e167f5dd848b5ae9' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.saf"), emit: saf
    tuple val("${task.process}"), val('gawk'), eval("awk --version | head -1 | sed 's/GNU Awk //; s/,.*//'"), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // BED is 0-based half-open and SAF is 1-based inclusive, so the start shifts by one. Column 4
    // becomes the feature id and column 6 the strand, defaulting to unstranded when absent.
    """
    awk -F'\\t' 'BEGIN { OFS="\\t"; print "GeneID", "Chr", "Start", "End", "Strand" }
        \$0 !~ /^(#|track|browser)/ && NF >= 3 {
            id = (NF >= 4 && \$4 != "") ? \$4 : \$1 "_" \$2 "_" \$3
            strand = (NF >= 6 && (\$6 == "+" || \$6 == "-")) ? \$6 : "."
            print id, \$1, \$2 + 1, \$3, strand
        }' ${bed} > ${prefix}.saf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.saf
    """
}
