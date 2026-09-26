//
// Uncompress and prepare reference genome files
//

include {
    GUNZIP as GUNZIP_ENDO_FASTA
    GUNZIP as GUNZIP_SPIKEIN_FASTA
    GUNZIP as GUNZIP_HYBRID_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_SPARSEBED
    GUNZIP as GUNZIP_ACTIVE_REGIONS
    GUNZIP as GUNZIP_ROCCO_PARAMS
    GUNZIP as GUNZIP_BLACKLIST
    GUNZIP as GUNZIP_SPIKEIN_BARCODE_TABLE
    GUNZIP as GUNZIP_INITIATION_ZONES
    GUNZIP as GUNZIP_OKSEQ_RFD_FILE
    GUNZIP as GUNZIP_SPLICESITES
    GUNZIP as GUNZIP_TE_COUNTING_GENE_GTF
    GUNZIP as GUNZIP_TECOUNT_GENE_INDEX
    GUNZIP as GUNZIP_TELOCAL_GENE_INDEX
    GUNZIP as GUNZIP_TE_GTF
    GUNZIP as GUNZIP_TECOUNT_TE_INDEX
    GUNZIP as GUNZIP_TELOCAL_TE_INDEX
    GUNZIP as GUNZIP_STROBEALIGN_INDEX
    } from '../../../modules/nf-core/gunzip/main'

include {
    UNTAR as UNTAR_BWA_INDEX
    UNTAR as UNTAR_BWAMEM2_INDEX
    UNTAR as UNTAR_MINIBWA_INDEX
    UNTAR as UNTAR_BOWTIE_INDEX
    UNTAR as UNTAR_BOWTIE2_INDEX
    UNTAR as UNTAR_STAR_INDEX
    UNTAR as UNTAR_CHROMAP_INDEX
    UNTAR as UNTAR_HISAT2_INDEX
    UNTAR as UNTAR_MINIMAP2_INDEX
    } from '../../../modules/nf-core/untar/main'

include { BUILD_HYBRID_FASTA            } from '../../../modules/local/build_hybrid_fasta/main'
include { GFFREAD                       } from '../../../modules/nf-core/gffread/main'
include {
    SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_ENDO
    SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SPIKEIN
    SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_HYBRID
} from '../../../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX                     } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX                 } from '../../../modules/nf-core/bwamem2/index/main'
include { MINIBWA_INDEX                 } from '../../../modules/nf-core/minibwa/index/main'
include { BOWTIE_BUILD                  } from '../../../modules/nf-core/bowtie/build/main'
include { BOWTIE2_BUILD                 } from '../../../modules/nf-core/bowtie2/build/main'
include { CHROMAP_INDEX                 } from '../../../modules/nf-core/chromap/index/main'
include { STAR_GENOMEGENERATE           } from '../../../modules/nf-core/star/genomegenerate/main'
include { HISAT2_BUILD                  } from '../../../modules/nf-core/hisat2/build/main'
include { HISAT2_EXTRACTSPLICESITES     } from '../../../modules/nf-core/hisat2/extractsplicesites/main'
include { MINIMAP2_INDEX                } from '../../../modules/nf-core/minimap2/index/main'
include { KHMER_UNIQUEKMERS             } from '../../../modules/nf-core/khmer/uniquekmers/main'
include { RFD_TO_IZ                     } from '../../../modules/local/rfd_to_iz/main'

include { GFF3SORT                      } from '../../../modules/local/gff3sort/main'
include { HTSLIB_BGZIPTABIX             } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { EAUTILS_GTF2BED               } from '../../../modules/nf-core/ea-utils/gtf2bed/main'
include { GENOME_WHITELIST_REGIONS      } from '../../../modules/local/genome_whitelist_regions/main'
include { FASTA_SPLIT_BY_GENOME         } from '../../../modules/local/fasta_split_by_genome/main'

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
    spikein_fasta      //    path: path to spike-in genome fasta file
    hybrid_fasta       //    path: path to hybrid fasta file
    gtf                //    file: /path/to/genome.gtf
    gff                //    file: /path/to/genome.gff
    blacklist          //    file: /path/to/blacklist.bed
    spikein_barcode_table //    file: /path/to/spikein_barcode_table.tsv
    read_length        //    integer: read length for khmer
    macs_gsize         //    string: genome size for MACS2
    sparsebed          //    file: /path/to/sparsebed.bed
    active_regions     //    file: /path/to/active_regions.bed
    rocco_params       //    file: /path/to/rocco_params.yml
    skip_gtf_index      //    boolean: skip GTF indexing
    gene_bed           //    file: /path/to/gene.bed
    bwa_index          //    file: /path/to/bwa/index/
    bwamem2_index      //    file: /path/to/bwamem2/index/
    minibwa_index      //    file: /path/to/minibwa/index/
    bowtie_index        //    file: /path/to/bowtie/index/
    bowtie2_index      //    file: /path/to/bowtie2/index/
    chromap_index      //    file: /path/to/chromap/index/
    star_index         //    file: /path/to/star/index/
    hisat2_index       //    file: /path/to/hisat2/index/
    minimap2_index     //    file: /path/to/minimap2/index/
    strobealign_index  //    file: /path/to/strobealign_index
    splicesites        //    file: /path/to/splicesites.txt
    okseq_rfd_file     //    file: /path/to/okseq_rfd_file.bed
    initiation_zones   //    file: /path/to/initiation_zones.bed
    skip_te_counting   //    boolean: skip TE counting
    skip_telocal    //    boolean: skip TElocal indexing
    te_counting_gene_gtf //    file: /path/to/te_counting_gene_gtf.gtf
    tecount_gene_index //    file: /path/to/tecount_gene_index.Ind
    telocal_gene_index //    file: /path/to/telocal_gene_index.Ind
    te_gtf     //    file: /path/to/te_gtf.gtf
    tecount_te_index   //    file: /path/to/tecount_te_index.Ind
    telocal_te_index   //    file: /path/to/telocal_te_index.locInd


    main:

    //
    // Uncompress genome endogenous fasta file if required
    //
    ch_fasta_to_align = channel.empty()
    ch_endo_fasta = channel.empty()
    if (fasta) {
        if (fasta.endsWith('.gz')) {
            ch_endo_fasta    = GUNZIP_ENDO_FASTA ( [ [id:'fasta'], file(fasta, checkIfExists: true) ] ).gunzip
        } else {
            ch_endo_fasta    = channel.value([ [ id:'fasta' ], file(fasta, checkIfExists: true) ])
        }
        ch_fasta_to_align = ch_endo_fasta
    }

    //
    // Uncompress or generate hybrid genome fasta file if required
    //
    ch_hybrid_fasta = channel.empty()
    ch_spikein_fasta = channel.empty()
    if (hybrid_fasta) {
        if (hybrid_fasta.endsWith('.gz')) {
            ch_hybrid_fasta    = GUNZIP_HYBRID_FASTA ( [ [id:'hybrid_fasta'], file(hybrid_fasta, checkIfExists: true) ] ).gunzip
        } else {
            ch_hybrid_fasta = channel.value([ [ id:'hybrid_fasta' ], file(hybrid_fasta, checkIfExists: true) ])
        }
        ch_fasta_to_align = ch_hybrid_fasta

        //
        // MODULE: Split hybrid fasta file into endogenous and spike-in fasta files
        //
        FASTA_SPLIT_BY_GENOME ( ch_hybrid_fasta, spikein_genome, genome )
        ch_endo_fasta = FASTA_SPLIT_BY_GENOME.out.endo_fasta.map { it -> [ it[0] + [ genome: genome ], it[1] ] }
        ch_spikein_fasta = FASTA_SPLIT_BY_GENOME.out.exo_fasta.map { it -> [ it[0] + [ genome: spikein_genome ], it[1] ] }

    } else if (spikein_fasta) {
        if (spikein_fasta.endsWith('.gz')) {
            ch_spikein_fasta    = GUNZIP_SPIKEIN_FASTA ( [ [id:'spikein_fasta'], file(spikein_fasta, checkIfExists: true) ] ).gunzip
        } else {
            ch_spikein_fasta = channel.value([ [ id:'spikein_fasta' ], file(spikein_fasta, checkIfExists: true) ])
        }

        //
        // MODULE: Create hybrid fasta file with endogenous and spike-in chromosomes
        //
        BUILD_HYBRID_FASTA (
            ch_endo_fasta.map { meta, fa -> [ meta, fa, genome ] },
            ch_spikein_fasta.map { meta, fa -> [ meta, fa, spikein_genome ] }
        )
        ch_hybrid_fasta = BUILD_HYBRID_FASTA.out.fasta
        ch_fasta_to_align = ch_hybrid_fasta

    }

    // Create channel: [ val(meta), fasta, fai ]
    // For FAIDX, since we do not have a fai yet
    ch_endo_fasta_fai = ch_endo_fasta.combine(channel.value([[]])).first()
    ch_spikein_fasta_fai = ch_spikein_fasta.combine(channel.value([[]])).first()
    ch_hybrid_fasta_fai = ch_hybrid_fasta.combine(channel.value([[]])).first()

    //
    // MODULE: Create fai and chromosome sizes files
    //
    SAMTOOLS_FAIDX_ENDO ( ch_endo_fasta_fai, true )
    SAMTOOLS_FAIDX_SPIKEIN ( ch_spikein_fasta_fai, true )
    SAMTOOLS_FAIDX_HYBRID ( ch_hybrid_fasta_fai, true )

    ch_endo_chromsizes = SAMTOOLS_FAIDX_ENDO.out.sizes
    ch_endo_fai = SAMTOOLS_FAIDX_ENDO.out.fai
    ch_spikein_chromsizes = SAMTOOLS_FAIDX_SPIKEIN.out.sizes
    ch_spikein_fai = SAMTOOLS_FAIDX_SPIKEIN.out.fai
    ch_hybrid_chromsizes = SAMTOOLS_FAIDX_HYBRID.out.sizes
    ch_hybrid_fai = SAMTOOLS_FAIDX_HYBRID.out.fai


    if (hybrid_fasta || spikein_fasta) {
        ch_fai_to_align = ch_hybrid_fai
        ch_chromsizes_to_align = ch_hybrid_chromsizes
    } else if (fasta) {
        ch_fai_to_align = ch_endo_fai
        ch_chromsizes_to_align = ch_endo_chromsizes
    }


    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    ch_gtf = channel.empty()
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            gtf = file(gtf, checkIfExists: true)
            ch_gtf      = GUNZIP_GTF ( [ [id:"${gtf.getBaseName(2)}"], gtf ] ).gunzip
        } else {
            gtf = file(gtf, checkIfExists: true)
            ch_gtf = channel.value( [ [id:"${gtf.getBaseName(1)}"], gtf ] )
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            gff = file(gff, checkIfExists: true)
            ch_gff      = GUNZIP_GFF ( [ [id:"${gff.getBaseName(2)}"], gff ] ).gunzip
        } else {
            gff = file(gff, checkIfExists: true)
            ch_gff = channel.value( [ [id:"${gff.getBaseName(1)}"], gff ] )
        }
        ch_gtf      = GFFREAD ( ch_gff, ch_endo_fasta.map{ it -> it[1] } ).gtf
    }

    if (!skip_gtf_index) {

        //
        // MODULE: Sort GTF file with gff3sort
        //
        GFF3SORT ( ch_gtf )
        ch_gtf = GFF3SORT.out.sorted

        //
        // MODULE: Compress and index sorted GTF file with bgzip/tabix
        //
        HTSLIB_BGZIPTABIX (
            ch_gtf.map { meta, file -> [ meta, file, [], [] ] },
            'compress',
            true,
            'gtf'
        )

    }

    if (spikein_barcode_table) {
        if (spikein_barcode_table.endsWith('.gz')) {
            ch_spikein_barcode_table = GUNZIP_SPIKEIN_BARCODE_TABLE ( [ [id:'spikein_barcode_table'], file(spikein_barcode_table, checkIfExists: true) ] ).gunzip
        } else {
            ch_spikein_barcode_table = channel.value( [ [id:'spikein_barcode_table'], file(spikein_barcode_table, checkIfExists: true) ] )
        }
    } else {
        ch_spikein_barcode_table = channel.value( [ [id:'spikein_barcode_table'], file("${projectDir}/assets/barcodes/sc_spikein_barcodes.tsv", checkIfExists: true) ] )
    }

    ch_sparsebed = channel.empty().first()
    if (sparsebed) {
        if (sparsebed.endsWith('.gz')) {
            ch_sparsebed = GUNZIP_SPARSEBED ( [ [id:'sparsebed'], file(sparsebed, checkIfExists: true) ] ).gunzip
        } else {
            ch_sparsebed = channel.value( [ [id:'sparsebed'], file(sparsebed, checkIfExists: true) ] )
        }
    }

    ch_active_regions = channel.empty().first()
    if (active_regions) {
        if (active_regions.endsWith('.gz')) {
            ch_active_regions = GUNZIP_ACTIVE_REGIONS ( [ [id:'active_regions'], file(active_regions, checkIfExists: true) ] ).gunzip
        } else {
            ch_active_regions = channel.value( [ [id:'active_regions'], file(active_regions, checkIfExists: true) ] )
        }
    }

    ch_rocco_params = channel.empty()
    if (rocco_params) {
        if (rocco_params.endsWith('.gz')) {
            ch_rocco_params = GUNZIP_ROCCO_PARAMS ( [ [id:'rocco_params'], file(rocco_params, checkIfExists: true) ] ).gunzip
        } else {
            ch_rocco_params = channel.value( [ [id:'rocco_params'], file(rocco_params, checkIfExists: true) ] )
        }
    }

    // Create dummy file
    // https://github.com/nf-core/sarek/blob/a7679b9b5c178351b1e96a3ffe7ee81ddf9aad06/main.nf#L201
    //ch_dummy_file = file("$baseDir/assets/dummy_file.txt", checkIfExists: true)

    //ch_blacklist = channel.value( [ [id:'blacklist'], ch_dummy_file ] )
    // Uncompress blacklist file if required
    ch_blacklist = channel.empty().first()
    if (blacklist) {
        if (blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST ( [ [id:'blacklist'], file(blacklist, checkIfExists: true) ] ).gunzip
        } else {
            ch_blacklist = channel.value( [ [id:'blacklist'], file(blacklist, checkIfExists: true) ] )
        }
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (!gene_bed) {
        ch_gene_bed = EAUTILS_GTF2BED ( ch_gtf ).bed
    } else {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [id:'gene_bed'], file(gene_bed, checkIfExists: true) ] ).gunzip
        } else {
            ch_gene_bed = channel.value( [ [id:'gene_bed'], file(gene_bed, checkIfExists: true) ] )
        }
    }


    ch_okseq_rfd_file = channel.empty().first()
    if (okseq_rfd_file) {
        if (okseq_rfd_file.endsWith('.gz')) {
            ch_okseq_rfd_file = GUNZIP_OKSEQ_RFD_FILE ( [ [id:'okseq_rfd_file'], file(okseq_rfd_file, checkIfExists: true) ] ).gunzip
        } else {
            ch_okseq_rfd_file = channel.value( [ [id:'okseq_rfd_file'], file(okseq_rfd_file, checkIfExists: true) ] )
        }
    }

    ch_initiation_zones = channel.empty().first()
    if (initiation_zones) {
        if (initiation_zones.endsWith('.gz')) {
            ch_initiation_zones = GUNZIP_INITIATION_ZONES ( [ [id:'initiation_zones'], file(initiation_zones, checkIfExists: true) ] ).gunzip
        } else {
            ch_initiation_zones = channel.value( [ [id:'initiation_zones'], file(initiation_zones, checkIfExists: true) ] )
        }
    } else if (okseq_rfd_file) {
        //
        // MODULE: Process OK-seq RFD file to get initiation zones
        //
        RFD_TO_IZ (
            ch_okseq_rfd_file,
            ch_blacklist.ifEmpty([[:], []]),
            ch_endo_chromsizes
        )
        ch_initiation_zones = RFD_TO_IZ.out.iz_bed
    }


    //
    // MODULE: Calculate genome size with khmer
    //

    // TODO: genome size is calculated with khmer even when not needed (no chipseq samples)
    // this is could be a workaround (https://github.com/nextflow-io/nextflow/discussions/5102#discussioncomment-9939140)
    ch_effective_gsize = channel.empty()
    if (!macs_gsize) {
        KHMER_UNIQUEKMERS (
            ch_endo_fasta,
            read_length
        )
        ch_effective_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it -> it[1].text.trim() }
    }

    // Create a channel with the effective genome fraction
    ch_endo_chromsizes
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
        .first()
        .map { size, egs ->
            // khmer's stub just touches an empty kmers.txt, so egs.toDouble() would throw "empty String" under -stub
            workflow.stubRun ? 1.0 : egs.toDouble() / size.toDouble()
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
    ch_whitelist = channel.empty()
    GENOME_WHITELIST_REGIONS (
        ch_endo_chromsizes,
        ch_blacklist//.ifEmpty([[:], []])
    )
    ch_whitelist = GENOME_WHITELIST_REGIONS.out.bed


    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index = channel.empty()
    if (prepare_tool_index == 'bwa') {
        if (bwa_index) {
            if (bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [id:'bwa_index'], file(bwa_index, checkIfExists: true) ] ).untar
            } else {
                ch_bwa_index = channel.value( [ [id:'bwa_index'], file(bwa_index, checkIfExists: true) ] )
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta_to_align ).index
        }
    }

    //
    // Uncompress BWAMEM2 index or generate from scratch if required
    //
    ch_bwamem2_index = channel.empty()
    if (prepare_tool_index == 'bwamem2') {
        if (bwamem2_index) {
            if (bwamem2_index.endsWith('.tar.gz')) {
                ch_bwamem2_index = UNTAR_BWAMEM2_INDEX ( [ [id:'bwamem2_index'], file(bwamem2_index, checkIfExists: true) ] ).untar
            } else {
                ch_bwamem2_index = channel.value( [ [id:'bwamem2_index'], file(bwamem2_index, checkIfExists: true) ] )
            }
        } else {
            ch_bwamem2_index = BWAMEM2_INDEX ( ch_fasta_to_align ).index
        }
    }

    //
    // Uncompress minibwa index or generate from scratch if required
    //
    ch_minibwa_index = channel.empty()
    if (prepare_tool_index == 'minibwa') {
        if (minibwa_index) {
            if (minibwa_index.endsWith('.tar.gz')) {
                ch_minibwa_index = UNTAR_MINIBWA_INDEX ( [ [id:'minibwa_index'], file(minibwa_index, checkIfExists: true) ] ).untar
            } else {
                ch_minibwa_index = channel.value( [ [id:'minibwa_index'], file(minibwa_index, checkIfExists: true) ] )
            }
        } else {
            ch_minibwa_index = MINIBWA_INDEX ( ch_fasta_to_align ).index
        }
    }

    //
    // Uncompress Bowtie index or generate from scratch if required
    //
    ch_bowtie_index = channel.empty()
    if (prepare_tool_index == 'bowtie') {
        if (bowtie_index) {
            if (bowtie_index.endsWith('.tar.gz')) {
                ch_bowtie_index = UNTAR_BOWTIE_INDEX ( [ [id:'bowtie_index'], file(bowtie_index, checkIfExists: true) ] ).untar
            } else {
                ch_bowtie_index = channel.value( [ [id:'bowtie_index'], file(bowtie_index, checkIfExists: true) ] )
            }
        } else {
            ch_bowtie_index = BOWTIE_BUILD ( ch_fasta_to_align ).index
        }
    }

    //
    // Uncompress Bowtie2 index or generate from scratch if required
    //
    ch_bowtie2_index = channel.empty()
    if (prepare_tool_index == 'bowtie2') {
        if (bowtie2_index) {
            if (bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( [ [id:'bowtie2_index'], file(bowtie2_index, checkIfExists: true) ] ).untar
            } else {
                ch_bowtie2_index = channel.value( [ [id:'bowtie2_index'], file(bowtie2_index, checkIfExists: true) ] )
            }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta_to_align ).index
        }
    }

    //
    // Uncompress CHROMAP index or generate from scratch if required
    //
    ch_chromap_index = channel.empty()
    if (prepare_tool_index == 'chromap') {
        if (chromap_index) {
            if (chromap_index.endsWith('.tar.gz')) {
                ch_chromap_index = UNTAR_CHROMAP_INDEX ( [ [id:'chromap_index'], file(chromap_index, checkIfExists: true) ] ).untar
            } else {
                ch_chromap_index = channel.value( [ [id:'chromap_index'], file(chromap_index, checkIfExists: true) ] )
            }
        } else {
            ch_chromap_index = CHROMAP_INDEX ( ch_fasta_to_align ).index
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = channel.empty()
    if (prepare_tool_index == 'star') {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [id:'star_index'], file(star_index, checkIfExists: true) ] ).untar
            } else {
                ch_star_index = channel.value( [ [id:'star_index'], file(star_index, checkIfExists: true) ] )
            }
        } else {
            ch_star_index = STAR_GENOMEGENERATE ( ch_fasta_to_align, ch_gtf ).index
        }
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = channel.empty()
    ch_hisat2_index = channel.empty()
    if (prepare_tool_index == 'hisat2') {
        if (!splicesites) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf ).txt
        } else {
            if (splicesites.endsWith('.gz')) {
                ch_splicesites = GUNZIP_SPLICESITES ( [ [id:'splicesites'], file(splicesites, checkIfExists: true) ] ).gunzip
            } else {
                ch_splicesites = channel.value( [ [id:'splicesites'], file(splicesites, checkIfExists: true) ] )
            }
        }
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [id:'hisat2_index'], file(hisat2_index, checkIfExists: true) ] ).untar
            } else {
                ch_hisat2_index = channel.value( [ [id:'hisat2_index'], file(hisat2_index, checkIfExists: true) ] )
            }
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_fasta_to_align, ch_gtf, ch_splicesites ).index
        }
    }

    //
    // Uncompress Minimap2 index or generate from scratch if required
    //
    ch_minimap2_index = channel.empty()
    if (prepare_tool_index == 'minimap2') {
        if (minimap2_index) {
            if (minimap2_index.endsWith('.tar.gz')) {
                ch_minimap2_index = UNTAR_MINIMAP2_INDEX ( [ [id:'minimap2_index'], file(minimap2_index, checkIfExists: true) ] ).untar
            } else {
                ch_minimap2_index = channel.value( [ [id:'minimap2_index'], file(minimap2_index, checkIfExists: true) ] )
            }
        } else {
            ch_minimap2_index = MINIMAP2_INDEX ( ch_fasta_to_align ).index
        }
    }

    //
    // Uncompress Strobealign index or generate from scratch if required
    //
    ch_strobealign_index = channel.empty().first()
    if (prepare_tool_index == 'strobealign') {
        if (strobealign_index) {
            if (strobealign_index.endsWith('.gz')) {
                ch_strobealign_index = GUNZIP_STROBEALIGN_INDEX ( [ [id:'strobealign_index'], file(strobealign_index, checkIfExists: true) ] ).gunzip
            } else {
                ch_strobealign_index = channel.value( [ [id:'strobealign_index'], file(strobealign_index, checkIfExists: true) ] )
            }
        } else {
            ch_strobealign_index = channel.value([[:], []])
        }
    }



    //
    // Uncompress gene GTF for TE counting or use existing gene GTF
    // With this, it is possible to provide a different GTF specifically for TE counting
    //
    ch_te_counting_gene_gtf = ch_gtf
    if (te_counting_gene_gtf) {
        if (te_counting_gene_gtf.endsWith('.gz')) {
            ch_te_counting_gene_gtf = GUNZIP_TE_COUNTING_GENE_GTF ( [ [id:'te_counting_gene_gtf'], file(te_counting_gene_gtf, checkIfExists: true) ] ).gunzip
        } else {
            ch_te_counting_gene_gtf = channel.value( [ [id:'te_counting_gene_gtf'], file(te_counting_gene_gtf, checkIfExists: true) ] )
        }
    }

    //
    // Uncompress TEcount or TElocal indices or generate from scratch if required
    //
    ch_tecount_gene_index = channel.empty()
    ch_telocal_gene_index = channel.empty()
    ch_tecount_te_index = channel.empty()
    ch_telocal_te_index = channel.empty()
    if (!skip_te_counting) {
        if (tecount_gene_index) {
            if (tecount_gene_index.endsWith('.gz')) {
                ch_tecount_gene_index = GUNZIP_TECOUNT_GENE_INDEX ( [ [id:'tecount_gene_index'], file(tecount_gene_index, checkIfExists: true) ] ).gunzip
            } else {
                ch_tecount_gene_index = channel.value( [ [id:'tecount_gene_index'], file(tecount_gene_index, checkIfExists: true) ] )
            }
        } else {
            ch_tecount_gene_index = TETRANSCRIPTS_INDEXER_GENE ( ch_te_counting_gene_gtf, 'gene' ).index
        }
        if (tecount_te_index) {
            if (tecount_te_index.endsWith('.gz')) {
                ch_tecount_te_index = GUNZIP_TECOUNT_TE_INDEX ( [ [id:'tecount_te_index'], file(tecount_te_index, checkIfExists: true) ] ).gunzip
            } else {
                ch_tecount_te_index = channel.value( [ [id:'tecount_te_index'], file(tecount_te_index, checkIfExists: true) ] )
            }
        } else {
            if (te_gtf.endsWith('.gz')) {
                ch_te_gtf = GUNZIP_TE_GTF ( [ [id:'te_gtf'], file(te_gtf, checkIfExists: true) ] ).gunzip
            } else {
                ch_te_gtf = channel.value( [ [id:'te_gtf'], file(te_gtf, checkIfExists: true) ] )
            }
            ch_tecount_te_index = TETRANSCRIPTS_INDEXER_TE ( ch_te_gtf, 'te' ).index
        }

        if (!skip_telocal) {
            if (telocal_gene_index) {
                if (telocal_gene_index.endsWith('.gz')) {
                    ch_telocal_gene_index = GUNZIP_TELOCAL_GENE_INDEX ( [ [id:'telocal_gene_index'], file(telocal_gene_index, checkIfExists: true) ] ).gunzip
                } else {
                    ch_telocal_gene_index = channel.value( [ [id:'telocal_gene_index'], file(telocal_gene_index, checkIfExists: true) ] )
                }
            } else {
                ch_telocal_gene_index = TELOCAL_INDEXER_GENE ( ch_te_counting_gene_gtf, 'gene' ).index
            }
            if (telocal_te_index) {
                if (telocal_te_index.endsWith('.gz')) {
                    ch_telocal_te_index = GUNZIP_TELOCAL_TE_INDEX ( [ [id:'telocal_te_index'], file(telocal_te_index, checkIfExists: true) ] ).gunzip
                } else {
                    ch_telocal_te_index = channel.value( [ [id:'telocal_te_index'], file(telocal_te_index, checkIfExists: true) ] )
                }
            } else if (te_gtf) {
                ch_telocal_te_index = TELOCAL_INDEXER_TE ( ch_te_gtf, 'TE' ).index
            }
        }
    }


    emit:
    fasta_to_align         = ch_fasta_to_align
    fai_to_align           = ch_fai_to_align
    chromsizes_to_align    = ch_chromsizes_to_align
    endo_fasta             = ch_endo_fasta                  //    channel: [ val(meta), [ genome.fasta ]]
    endo_fai               = ch_endo_fai                    //    channel: [ val(meta), [ genome.fai ]]
    spikein_fasta          = ch_spikein_fasta          //    channel: [ val(meta), [ spikein_genome.fasta ]]
    spikein_fai            = ch_spikein_fai                    //    channel: [ val(meta), [ spikein_genome.fai ]]
    hybrid_fasta           = ch_hybrid_fasta          //    channel: [ val(meta), [ hybrid_genome.fasta ]]
    hybrid_fai             = ch_hybrid_fai            //    channel: [ val(meta), [ hybrid_genome.fai ]]
    gtf                    = ch_gtf                    //    channel: [ val(meta), [ genome.gtf ]]
    gene_bed               = ch_gene_bed               //    channel: [ val(meta), [ gene.bed ]]
    endo_chromsizes        = ch_endo_chromsizes       //    channel: [ val(meta), [ genome.sizes ]]
    spikein_chromsizes     = ch_spikein_chromsizes        //    channel: [ val(meta), [ spikein_genome.sizes ]]
    hybrid_chromsizes      = ch_hybrid_chromsizes   //    channel: [ val(meta), [ hybrid_genome.sizes ]]
    effective_gsize        = ch_effective_gsize        //    channel: [ val(meta), [ effective_genome_size.txt ]]
    effective_gfraction    = ch_effective_gfraction
    whitelist              = ch_whitelist              //    channel: [ val(meta), [ *.include_regions.bed ]]
    blacklist              = ch_blacklist              //    channel: [  blacklist.bed ]
    spikein_barcode_table  = ch_spikein_barcode_table  //    channel: [ val(meta), [ spikein_barcode_table.tsv ]]
    sparsebed              = ch_sparsebed              //    channel: [ val(meta), [ sparsebed.bed ]]
    active_regions         = ch_active_regions         //    channel: [ val(meta), [ active_regions.bed ]]
    rocco_params           = ch_rocco_params           //    channel: [ val(meta), [ rocco_params.yml ]]
    okseq_rfd_file         = ch_okseq_rfd_file         //    channel: [ val(meta), [ okseq_rfd_file.bed ]]
    initiation_zones       = ch_initiation_zones       //    channel: [ val(meta), [ initiation_zones.bed ]]
    bwa_index              = ch_bwa_index              //    path: bwa/index/
    bwamem2_index          = ch_bwamem2_index          //    channel: [ val(meta), [ bwamem2/index/ ]]
    minibwa_index          = ch_minibwa_index          //    channel: [ val(meta), [ minibwa/index/ ]]
    bowtie_index           = ch_bowtie_index            //    channel: [ val(meta), [ bowtie/index/ ]]
    bowtie2_index          = ch_bowtie2_index          //    channel: [ val(meta), [ bowtie2/index/ ]]
    chromap_index          = ch_chromap_index          //    channel: [ val(meta), [ chromap/index/ ]]
    star_index             = ch_star_index             //    channel: [ val(meta), [ star/index/ ]]
    hisat2_index           = ch_hisat2_index           //    channel: [ val(meta), [ hisat2/index/ ]]
    minimap2_index         = ch_minimap2_index         //    channel: [ val(meta), [ minimap2/index/ ]]
    strobealign_index      = ch_strobealign_index         //    channel: [ val(meta), [ strobealign/index/ ]]
    splicesites            = ch_splicesites            //    channel: [ val(meta), [ splicesites.txt ]]
    te_counting_gene_gtf   = ch_te_counting_gene_gtf   //    channel: [ val(meta), [ te_counting_gene_gtf.gtf ]]
    tecount_gene_index     = ch_tecount_gene_index     //    channel: [ val(meta), [ tecount_gene_index.Ind ]]
    telocal_gene_index     = ch_telocal_gene_index     //    channel: [ val(meta), [ telocal_gene_index.Ind ]]
    tecount_te_index       = ch_tecount_te_index       //    channel: [ val(meta), [ tecount_te_index.Ind ]]
    telocal_te_index       = ch_telocal_te_index       //    channel: [ val(meta), [ telocal_te_index.locInd ]]
}
