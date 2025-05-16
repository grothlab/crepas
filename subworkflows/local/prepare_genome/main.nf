//
// Uncompress and prepare reference genome files
//

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_SPARSEBED
    GUNZIP as GUNZIP_BLACKLIST } from '../../../modules/nf-core/gunzip/main'

include {
    UNTAR as UNTAR_BWA_INDEX
    UNTAR as UNTAR_BOWTIE2_INDEX
    UNTAR as UNTAR_STAR_INDEX
    UNTAR as UNTAR_CHROMAP_INDEX
    UNTAR as UNTAR_HISAT2_INDEX  } from '../../../modules/nf-core/untar/main'

include { GFFREAD              } from '../../../modules/nf-core/gffread/main'
include { CUSTOM_GETCHROMSIZES } from '../../../modules/nf-core/custom/getchromsizes/main'
include { BWA_INDEX            } from '../../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD        } from '../../../modules/nf-core/bowtie2/build/main'
include { CHROMAP_INDEX        } from '../../../modules/nf-core/chromap/index/main'
include { STAR_GENOMEGENERATE      } from '../../../modules/nf-core/star/genomegenerate/main'
include { HISAT2_BUILD       } from '../../../modules/nf-core/hisat2/build/main'
include { HISAT2_EXTRACTSPLICESITES } from '../../../modules/nf-core/hisat2/extractsplicesites/main'

include { GTF2BED                  } from '../../../modules/local/gtf2bed/main'
include { GENOME_BLACKLIST_REGIONS } from '../../../modules/local/genome_blacklist_regions/main'
include { CHROM_SIZES_SPIKEIN_SPLIT  } from '../../../modules/local/chrom_sizes_spikein_split/main'

workflow PREPARE_GENOME {
    take:
    genome             //    string: genome name
    genomes            //    map: genome attributes
    spikein_genome     //    string: spikein genome name
    prepare_tool_index //    string  : tool to prepare index for
    fasta              //    path: path to genome fasta file
    gtf                //    file: /path/to/genome.gtf
    gff                //    file: /path/to/genome.gff
    blacklist          //    file: /path/to/blacklist.bed
    sparsebed          //    file: /path/to/sparsebed.bed
    gene_bed           //    file: /path/to/gene.bed
    bwa_index          //    file: /path/to/bwa/index/
    bowtie2_index      //    file: /path/to/bowtie2/index/
    chromap_index      //    file: /path/to/chromap/index/
    star_index         //    file: /path/to/star/index/
    hisat2_index       //    file: /path/to/hisat2/index/
    splicesites        //    file: /path/to/splicesites.txt


    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [id:'fasta'], params.fasta ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.of([ [ id:'fasta' ], file(params.fasta) ])
    }

    // ch_fasta_exo = Channel.empty()
    // if (params.spikein_fasta && params.spikein_fasta.endsWith('.gz')) {
    //     ch_fasta_exo = GUNZIP_FASTA ( [ [id:'spikein_fasta'], params.spikein_fasta ] ).gunzip
    //     ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    // } else if (params.spikein_fasta) {
    //     ch_fasta_exo = Channel.of([ [id:'spikein_fasta'], file(params.spikein_fasta) ])
    // }


    // Make fasta file available if reference saved or IGV is run
    // if (params.save_reference || !params.skip_igv) {
    //     file("${params.outdir}/genome/").mkdirs()
    //     // copy fasta file (second element of tuple) to output directory
    //     ch_fasta.map{ it[1] }.collect{ it.copyTo("${params.outdir}/genome/") }
    //     //ch_fasta_exo.map{ it[1] }.collect{ it.copyTo("${params.outdir}/genome/") }
    // }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [id:'gtf'], params.gtf ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.of( [ [id:'gtf'], file(params.gtf) ] )
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [id:'gff'], params.gff ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.of( [ [id:'gff'], file(params.gff) ] )
        }
        ch_gtf      = GFFREAD ( ch_gff, ch_fasta.map{ it[1] } ).gtf.map{ [ [id:'gtf'], it[1] ] }
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    ch_sparsebed = Channel.empty()
    if (params.sparsebed) {
        if (params.sparsebed.endsWith('.gz')) {
            ch_sparsebed = GUNZIP_SPARSEBED ( [ [id:'sparsebed'], params.sparsebed ] ).gunzip
            ch_versions  = ch_versions.mix(GUNZIP_SPARSEBED.out.versions)
        } else {
            ch_sparsebed = Channel.of( [ [id:'sparsebed'], file(params.sparsebed) ] )
        }
    }

    // Create dummy file 
    // https://github.com/nf-core/sarek/blob/a7679b9b5c178351b1e96a3ffe7ee81ddf9aad06/main.nf#L201
    ch_dummy_file = file("$baseDir/assets/dummy_file.txt", checkIfExists: true)

    ch_blacklist = Channel.of( [ [id:'blacklist'], ch_dummy_file ] )
    // Uncompress blacklist file if required
    if (params.blacklist) {
        if (params.blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST ( [ [id:'blacklist'], params.blacklist ] ).gunzip
            ch_versions  = ch_versions.mix(GUNZIP_BLACKLIST.out.versions)
        } else {
            ch_blacklist = Channel.of( [ [id:'blacklist'], file(params.blacklist) ] )
        }
    }

    ch_initiation_zones = Channel.of( [ [id:'initiation_zones'], ch_dummy_file ] )
    if (params.initiation_zones) {
        ch_initiation_zones = Channel.of( [ [id:'initiation_zones'], file(params.initiation_zones) ] )
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    // If --gtf is supplied along with --genome
    // Make gene bed from supplied --gtf instead of using iGenomes one automatically
    def make_bed = false
    if (!params.gene_bed) {
        make_bed = true
    } else if (params.genome && params.gtf) {
        if (params.genomes[ params.genome ].gtf != params.gtf) {
            make_bed = true
        }
    }

    if (make_bed) {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed.map{ [ [id:'gene_bed'], it[1] ] }
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    } else {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [id:'gene_bed'], params.gene_bed ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.of( [ [id:'gene_bed'], file(params.ch_gene_bed) ] )
        }
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta )
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Create endogenous genome chromosome sizes file
    //
    ch_chrom_sizes_endo = ch_chrom_sizes
    ch_chrom_sizes_exo = Channel.empty()
    if (spikein_genome) {
        CHROM_SIZES_SPIKEIN_SPLIT ( ch_chrom_sizes, spikein_genome, genome )
        ch_chrom_sizes_endo = CHROM_SIZES_SPIKEIN_SPLIT.out.endo_sizes.map { [ it[0] + [ genome: genome ], it[1] ] }
        ch_chrom_sizes_exo = CHROM_SIZES_SPIKEIN_SPLIT.out.exo_sizes.map { [ it[0] + [ genome: spikein_genome ], it[1] ] }
        ch_versions        = ch_versions.mix(CHROM_SIZES_SPIKEIN_SPLIT.out.versions)
    }

    // get list of scaffolds in CHROM_SIZES_SPIKEIN_SPLIT.out.sizes
    // this are the ones that contain a dot in the first column of the tab-separated file
    ch_scaffolds = Channel.empty()
    ch_chrom_sizes
        .map{
            meta, bed ->
                bed.splitCsv(header: false, sep: '\t').findAll{ it[0].contains('.') }
        }
        .map { it[0] } // Extract the scaffold IDs
        .set { ch_scaffolds }

    // TODO: remove channel output to file for debugging
    // ch_scaffolds
    //     .map {
    //         scaffolds ->
    //             "${scaffolds}"
    //     }
    //     .collectFile( name: 'ch_scaffolds.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    //
    // Prepare genome intervals for filtering by removing regions in blacklist file
    //
    ch_genome_filtered_bed = Channel.empty()

    GENOME_BLACKLIST_REGIONS (
        ch_chrom_sizes,
        // if second element of tuple is empty, use [] as input for GENOME_BLACKLIST_REGIONS
        ch_blacklist.map{ it[1] }.ifEmpty([])
    )
    ch_genome_filtered_bed = GENOME_BLACKLIST_REGIONS.out.bed
    ch_versions = ch_versions.mix(GENOME_BLACKLIST_REGIONS.out.versions)


    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index = Channel.empty()
    if (prepare_tool_index == 'bwa') {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [:], params.bwa_index ] ).untar.map{ it[1] }
                ch_versions  = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                ch_bwa_index = file(params.bwa_index)
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
        if (params.bowtie2_index) {
            if (params.bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( [ [:], params.bowtie2_index ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
            } else {
                ch_bowtie2_index = Channel.of( [ [:], file(params.bowtie2_index) ] )
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
        if (params.chromap_index) {
            if (params.chromap_index.endsWith('.tar.gz')) {
                ch_chromap_index = UNTAR_CHROMAP_INDEX ( [ [:], params.chromap_index ] ).untar
                ch_versions  = ch_versions.mix(UNTAR.out.versions)
            } else {
                ch_chromap_index = Channel.of( [ [:], file(params.chromap_index) ] )
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
        if (params.star_index) {
            if (params.star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], params.star_index ] ).untar
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.of( [ [:], file(params.star_index) ] )
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
            ch_splicesites = Channel.value(file(splicesites))
        }
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [:], hisat2_index ] ).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            } else {
                ch_hisat2_index = Channel.value(file(hisat2_index))
            }
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    } 

    emit:
    fasta         = ch_fasta                  //    channel: [ val(meta), [ genome.fasta ]]
    fai           = ch_fai                    //    channel: [ val(meta), [ genome.fai ]]
    gtf           = ch_gtf                    //    channel: [ val(meta), [ genome.gtf ]]
    gene_bed      = ch_gene_bed               //    channel: [ val(meta), [ gene.bed ]]
    chrom_sizes   = ch_chrom_sizes            //    channel: [ val(meta), [ genome.sizes ]]
    chrom_sizes_endo = ch_chrom_sizes_endo //    channel: [ val(meta), [ genome_endo.sizes ]]
    chrom_sizes_exo = ch_chrom_sizes_exo //    channel: [ val(meta), [ genome_exo.sizes ]]
    scaffolds  = ch_scaffolds              //    channel: [ scaffolds ]
    filtered_bed  = ch_genome_filtered_bed    //    channel: [ val(meta), [ *.include_regions.bed ]]
    blacklist     = ch_blacklist              //    channel: [  blacklist.bed ]
    sparsebed     = ch_sparsebed              //    channel: [ val(meta), [ sparsebed.bed ]]
    initiation_zones = ch_initiation_zones    //    channel: [ val(meta), [ initiation_zones.bed ]]
    bwa_index     = ch_bwa_index              //    path: bwa/index/
    bowtie2_index = ch_bowtie2_index          //    channel: [ val(meta), [ bowtie2/index/ ]]
    chromap_index = ch_chromap_index          //    channel: [ val(meta), [ chromap/index/ ]]
    star_index    = ch_star_index             //    channel: [ val(meta), [ star/index/ ]]
    hisat2_index  = ch_hisat2_index           //    channel: [ val(meta), [ hisat2/index/ ]]
    splicesites   = ch_splicesites            //    channel: [ val(meta), [ splicesites.txt ]]
    versions      = ch_versions.ifEmpty(null) //    channel: [ versions.yml ]
}
