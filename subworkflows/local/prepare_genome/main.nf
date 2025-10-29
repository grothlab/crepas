//
// Uncompress and prepare reference genome files
//

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_SPARSEBED
    GUNZIP as GUNZIP_ACTIVE_REGIONS
    GUNZIP as GUNZIP_ROCCO_PARAMS
    GUNZIP as GUNZIP_BLACKLIST
    GUNZIP as GUNZIP_INITIATION_ZONES
    GUNZIP as GUNZIP_OKSEQ_RFD_FILE
    GUNZIP as GUNZIP_SPLICESITES
    GUNZIP as GUNZIP_TECOUNT_GENE_INDEX
    GUNZIP as GUNZIP_TELOCAL_GENE_INDEX
    GUNZIP as GUNZIP_TE_GTF
    GUNZIP as GUNZIP_TECOUNT_TE_INDEX
    GUNZIP as GUNZIP_TELOCAL_TE_INDEX
    } from '../../../modules/nf-core/gunzip/main'

include {
    UNTAR as UNTAR_BWA_INDEX
    UNTAR as UNTAR_BOWTIE2_INDEX
    UNTAR as UNTAR_STAR_INDEX
    UNTAR as UNTAR_CHROMAP_INDEX
    UNTAR as UNTAR_HISAT2_INDEX
    } from '../../../modules/nf-core/untar/main'

include { GFFREAD              } from '../../../modules/nf-core/gffread/main'
include { CUSTOM_GETCHROMSIZES } from '../../../modules/nf-core/custom/getchromsizes/main'
include { BWA_INDEX            } from '../../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD        } from '../../../modules/nf-core/bowtie2/build/main'
include { CHROMAP_INDEX        } from '../../../modules/nf-core/chromap/index/main'
include { STAR_GENOMEGENERATE      } from '../../../modules/nf-core/star/genomegenerate/main'
include { HISAT2_BUILD       } from '../../../modules/nf-core/hisat2/build/main'
include { HISAT2_EXTRACTSPLICESITES } from '../../../modules/nf-core/hisat2/extractsplicesites/main'
include { KHMER_UNIQUEKMERS        } from '../../../modules/nf-core/khmer/uniquekmers/main'

include { GTF2BED                  } from '../../../modules/local/gtf2bed/main'
include { GENOME_WHITELIST_REGIONS } from '../../../modules/local/genome_whitelist_regions/main'
include { CHROM_SIZES_SPIKEIN_SPLIT  } from '../../../modules/local/chrom_sizes_spikein_split/main'

include {
    TETRANSCRIPTS_INDEXER as TETRANSCRIPTS_INDEXER_GENE
    TETRANSCRIPTS_INDEXER as TETRANSCRIPTS_INDEXER_TE
    } from '../../../modules/local/tetranscripts/indexer/main'

include {
    TELOCAL_INDEXER as TELOCAL_INDEXER_GENE
    TELOCAL_INDEXER as TELOCAL_INDEXER_TE
    } from '../../../modules/local/telocal/indexer/main'

workflow PREPARE_GENOME {
    take:
    genome             //    string: genome name
    spikein_genome     //    string: spikein genome name
    prepare_tool_index //    string  : tool to prepare index for
    fasta              //    path: path to genome fasta file
    gtf                //    file: /path/to/genome.gtf
    gff                //    file: /path/to/genome.gff
    blacklist          //    file: /path/to/blacklist.bed
    read_length        //    integer: read length for khmer
    macs_gsize         //    string: genome size for MACS2
    sparsebed          //    file: /path/to/sparsebed.bed
    active_regions     //    file: /path/to/active_regions.bed
    rocco_params       //    file: /path/to/rocco_params.yml
    gene_bed           //    file: /path/to/gene.bed
    bwa_index          //    file: /path/to/bwa/index/
    bowtie2_index      //    file: /path/to/bowtie2/index/
    chromap_index      //    file: /path/to/chromap/index/
    star_index         //    file: /path/to/star/index/
    hisat2_index       //    file: /path/to/hisat2/index/
    splicesites        //    file: /path/to/splicesites.txt
    okseq_rfd_file     //    file: /path/to/okseq_rfd_file.bed
    initiation_zones   //    file: /path/to/initiation_zones.bed
    skip_te_counting   //    boolean: skip TE counting
    skip_telocal    //    boolean: skip TElocal indexing
    tecount_gene_index //    file: /path/to/tecount_gene_index.Ind
    telocal_gene_index //    file: /path/to/telocal_gene_index.Ind
    te_gtf     //    file: /path/to/te_gtf.gtf
    tecount_te_index   //    file: /path/to/tecount_te_index.Ind
    telocal_te_index   //    file: /path/to/telocal_te_index.locInd


    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [id:'fasta'], file(fasta, checkIfExists: true) ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value([ [ id:'fasta' ], file(fasta, checkIfExists: true) ])
    }
    
    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    ch_gtf = Channel.empty()
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [id:'gtf'], file(gtf, checkIfExists: true) ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value( [ [id:'gtf'], file(gtf, checkIfExists: true) ] )
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [id:'gff'], file(gff, checkIfExists: true) ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.value( [ [id:'gff'], file(gff, checkIfExists: true) ] )
        }
        ch_gtf      = GFFREAD ( ch_gff, ch_fasta.map{ it[1] } ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    ch_sparsebed = Channel.empty()
    if (sparsebed) {
        if (sparsebed.endsWith('.gz')) {
            ch_sparsebed = GUNZIP_SPARSEBED ( [ [id:'sparsebed'], file(sparsebed, checkIfExists: true) ] ).gunzip
            ch_versions  = ch_versions.mix(GUNZIP_SPARSEBED.out.versions)
        } else {
            ch_sparsebed = Channel.value( [ [id:'sparsebed'], file(sparsebed, checkIfExists: true) ] )
        }
    }

    ch_active_regions = Channel.empty()
    if (active_regions) {
        if (active_regions.endsWith('.gz')) {
            ch_active_regions = GUNZIP_ACTIVE_REGIONS ( [ [id:'active_regions'], file(active_regions, checkIfExists: true) ] ).gunzip
            ch_versions  = ch_versions.mix(GUNZIP_ACTIVE_REGIONS.out.versions)
        } else {
            ch_active_regions = Channel.value( [ [id:'active_regions'], file(active_regions, checkIfExists: true) ] )
        }
    }

    ch_rocco_params = Channel.empty()
    if (rocco_params) {
        if (rocco_params.endsWith('.gz')) {
            ch_rocco_params = GUNZIP_ROCCO_PARAMS ( [ [id:'rocco_params'], file(rocco_params, checkIfExists: true) ] ).gunzip
            ch_versions     = ch_versions.mix(GUNZIP_ROCCO_PARAMS.out.versions)
        } else {
            ch_rocco_params = Channel.value( [ [id:'rocco_params'], file(rocco_params, checkIfExists: true) ] )
        }
    } 

    // Create dummy file 
    // https://github.com/nf-core/sarek/blob/a7679b9b5c178351b1e96a3ffe7ee81ddf9aad06/main.nf#L201
    //ch_dummy_file = file("$baseDir/assets/dummy_file.txt", checkIfExists: true)

    //ch_blacklist = Channel.value( [ [id:'blacklist'], ch_dummy_file ] )
    // Uncompress blacklist file if required
    ch_blacklist = Channel.empty()
    if (blacklist) {
        if (blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST ( [ [id:'blacklist'], file(blacklist, checkIfExists: true) ] ).gunzip
            ch_versions  = ch_versions.mix(GUNZIP_BLACKLIST.out.versions)
        } else {
            ch_blacklist = Channel.value( [ [id:'blacklist'], file(blacklist, checkIfExists: true) ] )
        }
    }

    ch_okseq_rfd_file = Channel.empty().first() // .first() ensures it is a value channel
    if (okseq_rfd_file) {
        if (okseq_rfd_file.endsWith('.gz')) {
            ch_okseq_rfd_file = GUNZIP_OKSEQ_RFD_FILE ( [ [id:'okseq_rfd_file'], file(okseq_rfd_file, checkIfExists: true) ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_OKSEQ_RFD_FILE.out.versions)
        } else {
            ch_okseq_rfd_file = Channel.value( [ [id:'okseq_rfd_file'], file(okseq_rfd_file, checkIfExists: true) ] )
        }
    }

    ch_initiation_zones = Channel.empty().first() // .first() ensures it is a value channel
    if (initiation_zones) {
        if (initiation_zones.endsWith('.gz')) {
            ch_initiation_zones = GUNZIP_INITIATION_ZONES ( [ [id:'initiation_zones'], file(initiation_zones, checkIfExists: true) ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_INITIATION_ZONES.out.versions)
        } else {
            ch_initiation_zones = Channel.value( [ [id:'initiation_zones'], file(initiation_zones, checkIfExists: true) ] )
        }
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (!gene_bed) {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    } else {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [id:'gene_bed'], file(gene_bed, checkIfExists: true) ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.value( [ [id:'gene_bed'], file(gene_bed, checkIfExists: true) ] )
        }
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta )
    ch_chrom_sizes_endo = CUSTOM_GETCHROMSIZES.out.sizes
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Create endogenous genome chromosome sizes file
    //
    ch_chrom_sizes_exo = Channel.empty()
    if (spikein_genome) {
        CHROM_SIZES_SPIKEIN_SPLIT ( ch_chrom_sizes_endo, spikein_genome, genome )
        ch_chrom_sizes_endo = CHROM_SIZES_SPIKEIN_SPLIT.out.endo_sizes.map { [ it[0] + [ genome: genome ], it[1] ] }
        ch_chrom_sizes_exo = CHROM_SIZES_SPIKEIN_SPLIT.out.exo_sizes.map { [ it[0] + [ genome: spikein_genome ], it[1] ] }
        ch_versions        = ch_versions.mix(CHROM_SIZES_SPIKEIN_SPLIT.out.versions)
    }

    //
    // MODULE: Calculate genome size with khmer
    //

    // TODO: genome size is calculated with khmer even when not needed (no chipseq samples)
    // this is could be a workaround (https://github.com/nextflow-io/nextflow/discussions/5102#discussioncomment-9939140)
    ch_effective_gsize = Channel.empty()
    if (!macs_gsize) {
        KHMER_UNIQUEKMERS (
            ch_fasta,
            read_length
        )
        ch_effective_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it[1].text.trim() }
    }

    // Create a channel with the effective genome fraction
    ch_chrom_sizes_endo
        .map { meta, bed ->
            bed.splitCsv(header: false, sep: '\t')
        }
        .flatMap { bed ->
            bed.collect { chr, size ->
                [size.toLong()]
            }
        }
        .sum()
        .combine(ch_effective_gsize)
        .map { size, egs ->
            egs.toDouble() / size.toDouble()
        }
        .set { ch_effective_gfraction }

    // TODO: Print to file for debuggin
    ch_effective_gfraction
        .map { egf ->
            "${egf}"
        }
        .collectFile(name: 'ch_effective_gfraction.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/PREPARE_GENOME")



    //
    // Prepare genome intervals for filtering by removing regions in blacklist file
    //
    ch_whitelist = Channel.empty()
    GENOME_WHITELIST_REGIONS (
        ch_chrom_sizes_endo,
        ch_blacklist//.ifEmpty([[:], []])
    )
    ch_whitelist = GENOME_WHITELIST_REGIONS.out.bed
    ch_versions = ch_versions.mix(GENOME_WHITELIST_REGIONS.out.versions)


    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index = Channel.empty()
    if (prepare_tool_index == 'bwa') {
        if (bwa_index) {
            if (bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [id:'bwa_index'], file(bwa_index, checkIfExists: true) ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                ch_bwa_index = Channel.value( [ [id:'bwa_index'], file(bwa_index, checkIfExists: true) ] )
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta ).index
            ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    //
    // Uncompress Bowtie2 index or generate from scratch if required
    //
    ch_bowtie2_index = Channel.empty()
    if (prepare_tool_index == 'bowtie2') {
        if (bowtie2_index) {
            if (bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( [ [id:'bowtie2_index'], file(bowtie2_index, checkIfExists: true) ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
            } else {
                ch_bowtie2_index = Channel.value( [ [id:'bowtie2_index'], file(bowtie2_index, checkIfExists: true) ] )
            }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta ).index
            ch_versions      = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }

    //
    // Uncompress CHROMAP index or generate from scratch if required
    //
    ch_chromap_index = Channel.empty()
    if (prepare_tool_index == 'chromap') {
        if (chromap_index) {
            if (chromap_index.endsWith('.tar.gz')) {
                ch_chromap_index = UNTAR_CHROMAP_INDEX ( [ [id:'chromap_index'], file(chromap_index, checkIfExists: true) ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_CHROMAP_INDEX.out.versions)
            } else {
                ch_chromap_index = Channel.value( [ [id:'chromap_index'], file(chromap_index, checkIfExists: true) ] )
            }
        } else {
            ch_chromap_index = CHROMAP_INDEX ( ch_fasta ).index
            ch_versions  = ch_versions.mix(CHROMAP_INDEX.out.versions)
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if (prepare_tool_index == 'star') {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [id:'star_index'], file(star_index, checkIfExists: true) ] ).untar
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value( [ [id:'star_index'], file(star_index, checkIfExists: true) ] )
            }
        } else {
            ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
            ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()
    if (prepare_tool_index == 'hisat2') {
        if (!splicesites) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            if (splicesites.endsWith('.gz')) {
                ch_splicesites = GUNZIP_SPLICESITES ( [ [id:'splicesites'], file(splicesites, checkIfExists: true) ] ).gunzip
                ch_versions    = ch_versions.mix(GUNZIP_SPLICESITES.out.versions)
            } else {
                ch_splicesites = Channel.value( [ [id:'splicesites'], file(splicesites, checkIfExists: true) ] )
            }
        }
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [id:'hisat2_index'], file(hisat2_index, checkIfExists: true) ] ).untar
                ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            } else {
                ch_hisat2_index = Channel.value( [ [id:'hisat2_index'], file(hisat2_index, checkIfExists: true) ] )
            }
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_fasta, ch_gtf, ch_splicesites ).index
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    //
    // Uncompress TEcount or TElocal indices or generate from scratch if required
    //
    ch_tecount_gene_index = Channel.empty()
    ch_telocal_gene_index = Channel.empty()
    ch_tecount_te_index = Channel.empty()
    ch_telocal_te_index = Channel.empty()
    if (!skip_te_counting) {
        if (tecount_gene_index) {
            if (tecount_gene_index.endsWith('.gz')) {
                ch_tecount_gene_index = GUNZIP_TECOUNT_GENE_INDEX ( [ [id:'tecount_gene_index'], file(tecount_gene_index, checkIfExists: true) ] ).gunzip
                ch_versions = ch_versions.mix(GUNZIP_TECOUNT_GENE_INDEX.out.versions)
            } else {
                ch_tecount_gene_index = Channel.value( [ [id:'tecount_gene_index'], file(tecount_gene_index, checkIfExists: true) ] )
            }
        } else {
            ch_tecount_gene_index = TETRANSCRIPTS_INDEXER_GENE ( ch_gtf, 'gene' ).index
            ch_versions = ch_versions.mix(TETRANSCRIPTS_INDEXER_GENE.out.versions)
        }
        if (tecount_te_index) {
            if (tecount_te_index.endsWith('.gz')) {
                ch_tecount_te_index = GUNZIP_TECOUNT_TE_INDEX ( [ [id:'tecount_te_index'], file(tecount_te_index, checkIfExists: true) ] ).gunzip
                ch_versions = ch_versions.mix(GUNZIP_TECOUNT_TE_INDEX.out.versions)
            } else {
                ch_tecount_te_index = Channel.value( [ [id:'tecount_te_index'], file(tecount_te_index, checkIfExists: true) ] )
            }
        } else {
            if (te_gtf.endsWith('.gz')) {
                ch_te_gtf = GUNZIP_TE_GTF ( [ [id:'te_gtf'], file(te_gtf, checkIfExists: true) ] ).gunzip
                ch_versions = ch_versions.mix(GUNZIP_TE_GTF.out.versions)
            } else {
                ch_te_gtf = Channel.value( [ [id:'te_gtf'], file(te_gtf, checkIfExists: true) ] )
            }
            ch_tecount_te_index = TETRANSCRIPTS_INDEXER_TE ( ch_te_gtf, 'te' ).index
            ch_versions = ch_versions.mix(TETRANSCRIPTS_INDEXER_TE.out.versions)
        }

        if (!skip_telocal) {
            if (telocal_gene_index) {
                if (telocal_gene_index.endsWith('.gz')) {
                    ch_telocal_gene_index = GUNZIP_TELOCAL_GENE_INDEX ( [ [id:'telocal_gene_index'], file(telocal_gene_index, checkIfExists: true) ] ).gunzip
                    ch_versions = ch_versions.mix(GUNZIP_TELOCAL_GENE_INDEX.out.versions)
                } else {
                    ch_telocal_gene_index = Channel.value( [ [id:'telocal_gene_index'], file(telocal_gene_index, checkIfExists: true) ] )
                }
            } else {
                ch_telocal_gene_index = TELOCAL_INDEXER_GENE ( ch_gtf, 'gene' ).index
                ch_versions = ch_versions.mix(TELOCAL_INDEXER_GENE.out.versions)
            }
            if (telocal_te_index) {
                if (telocal_te_index.endsWith('.gz')) {
                    ch_telocal_te_index = GUNZIP_TELOCAL_TE_INDEX ( [ [id:'telocal_te_index'], file(telocal_te_index, checkIfExists: true) ] ).gunzip
                    ch_versions = ch_versions.mix(GUNZIP_TELOCAL_TE_INDEX.out.versions)
                } else {
                    ch_telocal_te_index = Channel.value( [ [id:'telocal_te_index'], file(telocal_te_index, checkIfExists: true) ] )
                }
            } else if (te_gtf) {
                ch_telocal_te_index = TELOCAL_INDEXER_TE ( ch_te_gtf, 'TE' ).index
                ch_versions = ch_versions.mix(TELOCAL_INDEXER_TE.out.versions)
            }
        }
    }

        
    emit:
    fasta                  = ch_fasta                  //    channel: [ val(meta), [ genome.fasta ]]
    fai                    = ch_fai                    //    channel: [ val(meta), [ genome.fai ]]
    gtf                    = ch_gtf                    //    channel: [ val(meta), [ genome.gtf ]]
    gene_bed               = ch_gene_bed               //    channel: [ val(meta), [ gene.bed ]]
    chrom_sizes_endo       = ch_chrom_sizes_endo       //    channel: [ val(meta), [ genome_endo.sizes ]]
    chrom_sizes_exo        = ch_chrom_sizes_exo        //    channel: [ val(meta), [ genome_exo.sizes ]]
    effective_gsize        = ch_effective_gsize        //    channel: [ val(meta), [ effective_genome_size.txt ]]
    effective_gfraction    = ch_effective_gfraction
    whitelist              = ch_whitelist              //    channel: [ val(meta), [ *.include_regions.bed ]]
    blacklist              = ch_blacklist              //    channel: [  blacklist.bed ]
    sparsebed              = ch_sparsebed              //    channel: [ val(meta), [ sparsebed.bed ]]
    active_regions         = ch_active_regions         //    channel: [ val(meta), [ active_regions.bed ]]
    rocco_params           = ch_rocco_params           //    channel: [ val(meta), [ rocco_params.yml ]]
    okseq_rfd_file         = ch_okseq_rfd_file         //    channel: [ val(meta), [ okseq_rfd_file.bed ]]
    initiation_zones       = ch_initiation_zones       //    channel: [ val(meta), [ initiation_zones.bed ]]
    bwa_index              = ch_bwa_index              //    path: bwa/index/
    bowtie2_index          = ch_bowtie2_index          //    channel: [ val(meta), [ bowtie2/index/ ]]
    chromap_index          = ch_chromap_index          //    channel: [ val(meta), [ chromap/index/ ]]
    star_index             = ch_star_index             //    channel: [ val(meta), [ star/index/ ]]
    hisat2_index           = ch_hisat2_index           //    channel: [ val(meta), [ hisat2/index/ ]]
    splicesites            = ch_splicesites            //    channel: [ val(meta), [ splicesites.txt ]]
    tecount_gene_index     = ch_tecount_gene_index     //    channel: [ val(meta), [ tecount_gene_index.Ind ]]
    telocal_gene_index     = ch_telocal_gene_index     //    channel: [ val(meta), [ telocal_gene_index.Ind ]]
    tecount_te_index       = ch_tecount_te_index       //    channel: [ val(meta), [ tecount_te_index.Ind ]]
    telocal_te_index       = ch_telocal_te_index       //    channel: [ val(meta), [ telocal_te_index.locInd ]]
    versions               = ch_versions                //    channel: [ versions.yml ]
}
