#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grothlab/crepas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/grothlab/crepas
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
include { CREPAS                   } from './workflows/crepas'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_grothlab_crepas_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_grothlab_crepas_pipeline'
include { getMacsGsize            } from './subworkflows/local/utils_grothlab_crepas_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_grothlab_crepas_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow GROTHLAB_CREPAS {

    main:
    ch_versions = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        GENOME PARAMETER VALUES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // collect paths from genome attributes file (e.g. iGenomes.config; optional)
    // we cannot overwrite params in the workflow (they stay null as coming from the config file)
    // TODO: simplify, readability
    def fasta               = params.containsKey('fasta') ? params.fasta : (params.refgenie_ignore ? null : getGenomeAttribute('fasta'))
    def bwa_index           = params.containsKey('bwa_index') ? params.bwa_index : (params.refgenie_ignore ? null : getGenomeAttribute('bwa'))
    def bowtie2_index       = params.containsKey('bowtie2_index') ? params.bowtie2_index : (params.refgenie_ignore ? null : getGenomeAttribute('bowtie2_index'))
    def chromap_index       = params.containsKey('chromap_index') ? params.chromap_index : (params.refgenie_ignore ? null : getGenomeAttribute('chromap'))
    def star_index          = params.containsKey('star_index') ? params.star_index : (params.refgenie_ignore ? null : getGenomeAttribute('star'))
    def hisat2_index        = params.containsKey('hisat2_index') ? params.hisat2_index : (params.refgenie_ignore ? null : getGenomeAttribute('hisat2'))
    def gtf                 = params.containsKey('gtf') ? params.gtf : (params.refgenie_ignore ? null : getGenomeAttribute('gtf'))
    def gff                 = params.containsKey('gff') ? params.gff : (params.refgenie_ignore ? null : getGenomeAttribute('gff'))
    def gene_bed            = params.containsKey('gene_bed') ? params.gene_bed : (params.refgenie_ignore ? null : getGenomeAttribute('gene_bed'))
    def blacklist           = params.containsKey('blacklist') ? params.blacklist : (params.refgenie_ignore ? null : getGenomeAttribute('blacklist'))
    def sparsebed           = params.containsKey('sparsebed') ? params.sparsebed : (params.refgenie_ignore ? null : getGenomeAttribute('sparsebed'))
    def active_regions      = params.containsKey('active_regions') ? params.active_regions : (params.refgenie_ignore ? null : getGenomeAttribute('active_regions'))
    def rocco_params        = params.containsKey('rocco_params') ? params.rocco_params : (params.refgenie_ignore ? null : getGenomeAttribute('rocco_params'))
    def splicesites         = params.containsKey('splicesites') ? params.splicesites : (params.refgenie_ignore ? null : getGenomeAttribute('splicesites'))
    def okseq_rfd_file      = params.containsKey('okseq_rfd_file') ? params.okseq_rfd_file : (params.refgenie_ignore ? null : getGenomeAttribute('okseq_rfd_file'))
    def initiation_zones    = params.containsKey('initiation_zones') ? params.initiation_zones : (params.refgenie_ignore ? null : getGenomeAttribute('initiation_zones'))
    def tecount_gene_index  = params.containsKey('tecount_gene_index') ? params.tecount_gene_index : (params.refgenie_ignore ? null : getGenomeAttribute('tecount_gene_index'))
    def telocal_gene_index  = params.containsKey('telocal_gene_index') ? params.telocal_gene_index : (params.refgenie_ignore ? null : getGenomeAttribute('telocal_gene_index'))
    def te_gtf              = params.containsKey('te_gtf') ? params.te_gtf : (params.refgenie_ignore ? null : getGenomeAttribute('te_gtf'))
    def tecount_te_index    = params.containsKey('tecount_te_index') ? params.tecount_te_index : (params.refgenie_ignore ? null : getGenomeAttribute('tecount_te_index'))
    def telocal_te_index    = params.containsKey('telocal_te_index') ? params.telocal_te_index : (params.refgenie_ignore ? null : getGenomeAttribute('telocal_te_index'))
    def macs_gsize          = params.containsKey('macs_gsize') ? params.macs_gsize : getMacsGsize(params)


    //
    // SUBWORKFLOW: Prepare reference genome files
    //
    PREPARE_GENOME (
        params.genome,
        params.genomes,
        params.spikein_genome,
        params.aligner,
        fasta,
        gtf,
        gff,
        blacklist,
        params.read_length,
        macs_gsize,
        sparsebed,
        active_regions,
        rocco_params,
        params.skip_gtf_index,
        gene_bed,
        bwa_index,
        bowtie2_index,
        chromap_index,
        star_index,
        hisat2_index,
        splicesites,
        okseq_rfd_file,
        initiation_zones,
        params.skip_te_counting,
        params.skip_telocal,
        tecount_gene_index,
        telocal_gene_index,
        te_gtf,
        tecount_te_index,
        telocal_te_index
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // WORKFLOW: Run grothlab/crepas workflow
    //
    ch_samplesheet = channel.value(file(params.input, checkIfExists: true))
    CREPAS (
        ch_samplesheet,
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
        PREPARE_GENOME.out.okseq_rfd_file,
        PREPARE_GENOME.out.initiation_zones,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.bowtie2_index,
        PREPARE_GENOME.out.chromap_index,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.splicesites
    )

    emit:
    multiqc_report = CREPAS.out.multiqc_report // channel: /path/to/multiqc_report.html
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
        params.validate_params,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    GROTHLAB_CREPAS ()

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
        GROTHLAB_CREPAS.out.multiqc_report
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
