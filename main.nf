#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grothlab/glseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/grothlab/glseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = getGenomeAttribute('fasta')
params.bwa_index        = getGenomeAttribute('bwa')
params.bowtie2_index    = getGenomeAttribute('bowtie2')
params.chromap_index    = getGenomeAttribute('chromap')
params.star_index       = getGenomeAttribute('star')
params.hisat2_index     = getGenomeAttribute('hisat2')
params.gtf              = getGenomeAttribute('gtf')
params.gff              = getGenomeAttribute('gff')
params.gene_bed         = getGenomeAttribute('gene_bed')
params.blacklist        = getGenomeAttribute('blacklist')
params.sparsebed        = getGenomeAttribute('sparsebed')
params.active_regions   = getGenomeAttribute('active_regions')
params.rocco_params     = getGenomeAttribute('rocco_params')
params.splicesites      = getGenomeAttribute('splicesites')
params.initiation_zones = getGenomeAttribute('initiation_zones')
params.macs_gsize       = getMacsGsize(params)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GLSEQ                 } from './workflows/glseq'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_grothlab_glseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_grothlab_glseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow GROTHLAB_GLSEQ {

    main:
    ch_versions = Channel.empty()

    // SUBWORKFLOW: Prepare genome files
    PREPARE_GENOME (
        params.genome,
        params.genomes,
        params.spikein_genome,
        params.aligner,
        params.fasta,
        params.gtf,
        params.gff,
        params.blacklist,
        params.sparsebed,
        params.active_regions,
        params.rocco_params,
        params.gene_bed,
        params.bwa_index,
        params.bowtie2_index,
        params.chromap_index,
        params.star_index,
        params.hisat2_index,
        params.splicesites
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // WORKFLOW: Run grothlab/glseq workflow
    //
    ch_input = Channel.value(file(params.input, checkIfExists: true))

    GLSEQ(
        ch_input,
        ch_versions,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.chrom_sizes_endo,
        PREPARE_GENOME.out.chrom_sizes_exo,
        PREPARE_GENOME.out.scaffolds,
        PREPARE_GENOME.out.filtered_bed,
        PREPARE_GENOME.out.blacklist,
        PREPARE_GENOME.out.sparsebed,
        PREPARE_GENOME.out.active_regions,
        PREPARE_GENOME.out.rocco_params,
        PREPARE_GENOME.out.initiation_zones,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.bowtie2_index,
        PREPARE_GENOME.out.chromap_index,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.splicesites
    )

    emit:
    multiqc_report = GLSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions                // channel: [version1, version2, ...]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    GROTHLAB_GLSEQ ()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        GROTHLAB_GLSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Get macs genome size (macs_gsize)
//
def getMacsGsize(params) {
    def val = null
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey('macs_gsize')) {
            if (params.genomes[ params.genome ][ 'macs_gsize' ].containsKey(params.read_length.toString())) {
                val = params.genomes[ params.genome ][ 'macs_gsize' ][ params.read_length.toString() ]
            }
        }
    }
    return val
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
