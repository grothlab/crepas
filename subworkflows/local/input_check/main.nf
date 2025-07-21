//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../../modules/local/samplesheet_check/main'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    seq_center  // string: sequencing center for read group

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it, seq_center) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2, fastq_umi ] ]
def create_fastq_channel(LinkedHashMap row, String seq_center) {
    def meta = [:]
    meta.id                 = row.sample
    meta.single_end         = row.single_end.toBoolean()
    meta.sep_umi_fq         = row.sep_umi_fq.toBoolean()
    meta.okseq_part_file    = row.okseq_part_file
    meta.exp_type           = row.exp_type
    meta.strandedness       = row.strandedness
    meta.antibody           = row.antibody
    meta.control            = row.control
    meta.is_control         = row.is_control.toBoolean()

    def read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id - ~/_T\d+$/}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\'"
    if (seq_center) {
        read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id - ~/_T\d+$/}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\\tCN:${seq_center}\'"
    }
    meta.read_group = read_group

    // add path(s) of the fastq file(s) to the meta map
    // TODO: simplify this (fastq_umi)
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        if (meta.sep_umi_fq) {
            if (!file(row.fastq_umi).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> UMI FastQ file does not exist!\n${row.fastq_umi}"
            }
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_umi) ] ]
        } else {
            fastq_meta = [ meta, [ file(row.fastq_1) ] ]
        }
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        if (meta.sep_umi_fq) {
            if (!file(row.fastq_umi).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> UMI FastQ file does not exist!\n${row.fastq_umi}"
            }
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2), file(row.fastq_umi) ] ]
        } else {
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
        }
    }
    if (meta.okseq_part_file && !file(meta.okseq_part_file).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> OKSeq part file does not exist!\n${meta.okseq_part_file}"
    }
    return fastq_meta
}
