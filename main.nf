#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grothlab/crepas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/grothlab/crepas
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CREPAS                  } from './workflows/crepas'
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
    def fasta                 = params.fasta ?: (params.refgenie_ignore ? null : getGenomeAttribute('fasta'))
    def spikein_fasta         = params.spikein_fasta ?: (params.refgenie_ignore ? null : getGenomeAttribute('spikein_fasta'))
    def hybrid_fasta          = params.hybrid_fasta ?: (params.refgenie_ignore ? null : getGenomeAttribute('hybrid_fasta'))
    def spikein_barcode_table = params.spikein_barcode_table ?: (params.refgenie_ignore ? null : getGenomeAttribute('spikein_barcode_table'))
    def bwa_index             = params.bwa_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('bwa'))
    def bwamem2_index         = params.bwamem2_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('bwamem2'))
    def minibwa_index         = params.minibwa_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('minibwa'))
    def bowtie_index          = params.bowtie_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('bowtie'))
    def bowtie2_index         = params.bowtie2_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('bowtie2'))
    def chromap_index         = params.chromap_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('chromap'))
    def star_index            = params.star_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('star'))
    def hisat2_index          = params.hisat2_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('hisat2'))
    def minimap2_index        = params.minimap2_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('minimap2'))
    def strobealign_index     = params.strobealign_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('strobealign_'))
    def gtf                   = params.gtf ?: (params.refgenie_ignore ? null : getGenomeAttribute('gtf'))
    def gff                   = params.gff ?: (params.refgenie_ignore ? null : getGenomeAttribute('gff'))
    def gene_bed              = params.gene_bed ?: (params.refgenie_ignore ? null : getGenomeAttribute('gene_bed'))
    def blacklist             = params.blacklist ?: (params.refgenie_ignore ? null : getGenomeAttribute('blacklist'))
    def sparsebed             = params.sparsebed ?: (params.refgenie_ignore ? null : getGenomeAttribute('sparsebed'))
    def active_regions        = params.active_regions ?: (params.refgenie_ignore ? null : getGenomeAttribute('active_regions'))
    def rocco_params          = params.rocco_params ?: (params.refgenie_ignore ? null : getGenomeAttribute('rocco_params'))
    def splicesites           = params.splicesites ?: (params.refgenie_ignore ? null : getGenomeAttribute('splicesites'))
    def okseq_rfd_file        = params.okseq_rfd_file ?: (params.refgenie_ignore ? null : getGenomeAttribute('okseq_rfd_file'))
    def initiation_zones      = params.initiation_zones ?: (params.refgenie_ignore ? null : getGenomeAttribute('initiation_zones'))
    def te_counting_gene_gtf  = params.te_counting_gene_gtf ?: (params.refgenie_ignore ? null : getGenomeAttribute('te_counting_gene_gtf'))
    def tecount_gene_index    = params.tecount_gene_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('tecount_gene_index'))
    def telocal_gene_index    = params.telocal_gene_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('telocal_gene_index'))
    def te_gtf                = params.te_gtf ?: (params.refgenie_ignore ? null : getGenomeAttribute('te_gtf'))
    def tecount_te_index      = params.tecount_te_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('tecount_te_index'))
    def telocal_te_index      = params.telocal_te_index ?: (params.refgenie_ignore ? null : getGenomeAttribute('telocal_te_index'))
    def macs_gsize            = params.macs_gsize ?: getMacsGsize(params)


    //
    // SUBWORKFLOW: Prepare reference genome files
    //
    PREPARE_GENOME (
        params.genome,
        params.spikein_genome,
        params.aligner,
        fasta,
        spikein_fasta,
        hybrid_fasta,
        gtf,
        gff,
        blacklist,
        spikein_barcode_table,
        params.read_length,
        macs_gsize,
        sparsebed,
        active_regions,
        rocco_params,
        params.skip_gtf_index,
        gene_bed,
        bwa_index,
        bwamem2_index,
        minibwa_index,
        bowtie_index,
        bowtie2_index,
        chromap_index,
        star_index,
        hisat2_index,
        minimap2_index,
        strobealign_index,
        splicesites,
        okseq_rfd_file,
        initiation_zones,
        params.skip_te_counting,
        params.skip_telocal,
        te_counting_gene_gtf,
        tecount_gene_index,
        telocal_gene_index,
        te_gtf,
        tecount_te_index,
        telocal_te_index
    )

    //
    // WORKFLOW: Run grothlab/crepas workflow
    //
    CREPAS (
        ch_versions,
        PREPARE_GENOME.out.fasta_to_align,
        PREPARE_GENOME.out.fai_to_align,
        PREPARE_GENOME.out.chromsizes_to_align,
        PREPARE_GENOME.out.endo_chromsizes,
        PREPARE_GENOME.out.spikein_chromsizes,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.effective_gsize,
        PREPARE_GENOME.out.effective_gfraction,
        PREPARE_GENOME.out.whitelist,
        PREPARE_GENOME.out.blacklist,
        PREPARE_GENOME.out.spikein_barcode_table,
        PREPARE_GENOME.out.sparsebed,
        PREPARE_GENOME.out.active_regions,
        PREPARE_GENOME.out.rocco_params,
        PREPARE_GENOME.out.okseq_rfd_file,
        PREPARE_GENOME.out.initiation_zones,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.bwamem2_index,
        PREPARE_GENOME.out.minibwa_index,
        PREPARE_GENOME.out.bowtie_index,
        PREPARE_GENOME.out.bowtie2_index,
        PREPARE_GENOME.out.chromap_index,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.minimap2_index,
        PREPARE_GENOME.out.strobealign_index,
        PREPARE_GENOME.out.splicesites,
        PREPARE_GENOME.out.tecount_gene_index,
        PREPARE_GENOME.out.telocal_gene_index,
        PREPARE_GENOME.out.tecount_te_index,
        PREPARE_GENOME.out.telocal_te_index
    )
    emit:
    multiqc_report = CREPAS.out.multiqc_report // channel: /path/to/multiqc_report.html
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
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
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
        GROTHLAB_CREPAS.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
