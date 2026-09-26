//
// Extract spike-in barcodes from CUT&RUN and CUT&Tag reactions
// More info: https://support.epicypher.com/docs/analyzing-snap-cutana-spike-in-controls-and-results-cut-and-tag
//
process FASTQ_EXTRACT_SPIKEIN_BARCODES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6e009cae040b62b96b2a9750b4d0ac9414d8deb4c87cd80ddbe34b09d8077207/data'
        : 'community.wave.seqera.io/library/coreutils_gawk_pigz_ripgrep_sed:d5a69c5aa4d3b77d'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(barcode_table)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: counts
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion 2>&1 | head -1 | sed 's/^GNU Awk //; s/,.*//'"), emit: versions_awk, topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //'"), emit: versions_pigz, topic: versions
    tuple val("${task.process}"), val('coreutils'), eval("env paste --version | head -1 | sed 's/^paste (GNU coreutils) //'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.spikein_barcode_extract"
    """
    count_barcodes_with_awk() {
        fastq_gz="\$1"
        barcode_tsv="\$2"
        pigz -cd "\$fastq_gz" | awk -F'\\t' -v OFS='\\t' '
            NR==FNR {
                if (FNR==1) next
                n++
                target[n] = \$1
                barcode_id[n] = \$2
                barcode_seq[n] = \$3
                next
            }
            (FNR % 4) == 2 {
                for (i=1; i<=n; i++) {
                    if (index(\$0, barcode_seq[i])) count[i]++
                }
            }
            END {
                for (i=1; i<=n; i++) {
                    print target[i], barcode_id[i], barcode_seq[i], count[i] + 0
                }
            }
        ' "\$barcode_tsv" -
    }

    # Both R1 and R2 are searched because the barcoded side of nucleosome may be ligated in either direction relative to the P5 or P7 adapters.
    count_barcodes_with_awk "${reads[0]}" "${barcode_table}" > r1_counts.tsv
    count_barcodes_with_awk "${reads[1]}" "${barcode_table}" > r2_counts.tsv

    printf "barcode_target\\tbarcode_id\\tbarcode_sequence\\tR1_count\\tR2_count\\n" > ${prefix}.tsv
    paste r1_counts.tsv r2_counts.tsv | awk -F'\\t' -v OFS='\\t' '{ print \$1, \$2, \$3, \$4, \$8 }' >> ${prefix}.tsv

    rm -f r1_counts.tsv r2_counts.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.spikein_barcode_extract"
    """
    touch  ${prefix}.tsv
    """
}
