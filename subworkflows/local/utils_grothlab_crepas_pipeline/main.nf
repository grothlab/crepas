//
// Subworkflow with functionality specific to the grothlab/crepas pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../../subworkflows/nf-core/utils_nfschema_plugin'
include { UTILS_NFCORE_PIPELINE     } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../../subworkflows/nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion        } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { logColours                } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap          } from 'plugin/nf-schema'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    emit:
    versions = ch_versions

}


/*
========================================================================================
    SUBWORKFLOW FOR SAMPLESHEET CHECKING
========================================================================================
*/

// List of cases to handle

// 1. Value in `input_control` column must be present in the `sample` column
// 2. strandedness must be specified for SCAR-seq and OK-seq samples and must not be specified for other exp_types
// 3. Antibody must be specified for ChIP-seq and ChIP-exo samples, but must not be specified if input_control is set or if exp_type is not ChIP-seq or ChIP-exo
// 4. TODO: check that technical replicate is not duplicated within a biological replicate

workflow INPUT_CHECK {
    take:
    ch_fastq        // channel: [ val(meta), [ fastq_1, fastq_2, fastq_umi ] ]
    ch_seq_center  // string: sequencing center for read group

    main:


    // TODO: print for debugging
    ch_fastq.map { meta, fastqs -> "${meta}\t${fastqs}" }
        .collectFile(name: 'ch_fastq_1.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/INPUT_CHECK")


    // Check if within each biological replicate all technical replicates are set
    ch_fastq
        .map { meta, fastqs -> [ meta.id, meta.brep, meta, fastqs ] }
        .groupTuple(by: [0,1])
        .map { id, brep, metas, fastq_lists ->
            def new_metas = metas
            if (metas.any { it -> it.trep }) {
                //println "There is a technical replicate set for sample group: ${sample_group}"
                if (!metas.every { it -> it.trep }) {
                    error(
                        """
                        ERROR: If any technical replicate within a biological replicate is assigned an ID, then all the technical replicates within that biological replicate must have an ID. 

                        Check biological replicate '${brep}' of sample ${id} in the samplesheet.

                        """.stripIndent()
                    )
                }
            } else {
                //println "No technical replicates set for sample group: ${sample_group}"
                // Assign a technical replicate number
                def trep_counter = 1
                new_metas = metas.collect { meta ->
                    def meta_clone = meta.clone()
                    meta_clone.trep = trep_counter
                    // TODO: is_input_control has not been defined here yet
                    meta_clone.input_trep = !meta.is_input_control ? meta_clone.trep : []
                    trep_counter += 1
                    return meta_clone
                }
            }
            // Check that all meta.trep are unique
            def trep_ids = new_metas.collect { it -> it.trep }
            if (trep_ids.size() != trep_ids.unique().size()) {
                error(
                    """
                    ERROR: Technical replicate IDs must be unique within each biological replicate.

                    Check biological replicate '${brep}' of sample ${id} in the samplesheet.

                    """.stripIndent()
                )
            }
            
            return [ id, brep, new_metas, fastq_lists ]
        }
        .transpose()
        .map { id, brep, meta, fastqs -> [ meta, fastqs ] }
        .set { ch_fastq }

    // TODO: print for debugging
    ch_fastq.map { meta, fastqs -> "${meta}\t${fastqs}" }
        .collectFile(name: 'ch_fastq_2.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/INPUT_CHECK")

    // Count technical replicates per biological replicate to avoid .groupTuple() bottlenecks downstream
    // See: https://nextflow-io.github.io/nf-schema/latest/samplesheets/examples/#combining-a-channel
    ch_fastq
        .map { meta, fastqs -> 
            def id_brep = "${meta.id}_${meta.brep}"
            [ id_brep ]
        }
        .flatten() // avoid e.g., "[id_brep]:1"
        .reduce([:]) { trep_count, id_brep ->
            trep_count[id_brep] = (trep_count[id_brep] ?: 0) + 1
            trep_count
        }
        .combine(ch_fastq)
        .map { trep_count, meta, fastqs -> 
            def meta_clone = meta.clone()
            def id_brep = "${meta.id}_${meta.brep}"
            meta_clone.trep_count = trep_count[id_brep]
            [ meta_clone, fastqs ]
        }
        .set { ch_fastq }

    // TODO: print for debugging
    ch_fastq.map { meta, fastqs -> "${meta}\t${fastqs}" }
        .collectFile(name: 'ch_fastq_3.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/INPUT_CHECK")


    ch_fastq
        //.combine(seq_center)
        .map { meta, fastqs -> //, seq_center ->
            def meta_clone = meta.clone()
            meta_clone.id = "${meta.exp_type}_${meta.id}_bRep_${meta.brep}_tRep_${meta.trep}"
            meta_clone.original_id = meta.id
            if (meta.input_control) {
                meta_clone.input_control = "${meta.exp_type}_${meta.input_control}_bRep_${meta.input_brep}_tRep_${meta.input_trep}"
            }

            def read_group = "\'@RG\\tID:${meta_clone.id}\\tSM:${meta_clone.id - ~/_tRep_.*$/}\\tPL:ILLUMINA\\tLB:${meta_clone.id}\\tPU:1\'"
            if (ch_seq_center) {
                read_group = "\'@RG\\tID:${meta_clone.id}\\tSM:${meta_clone.id - ~/_tRep_.*$/}\\tPL:ILLUMINA\\tLB:${meta_clone.id}\\tPU:1\\tCN:${ch_seq_center}\'"
            }
            meta_clone.read_group = read_group
            [ meta_clone, fastqs ]
        }
        .set { ch_fastq } 

    

    // TODO: print for debugging
    ch_fastq.map { meta, fastqs -> "${meta}\t${fastqs}" }
        .collectFile(name: 'ch_fastq_4.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/INPUT_CHECK")

    //
    // PARSING EXTRA META FIELDS
    //
    // get list of input control samples
    ch_fastq
        .map { meta, fastqs -> meta.input_control }
        .unique()
        .collect()
        .map { it -> [it] } // unflatten
        .set { ch_ipcontrol_list }

    // TODO: print for debugging
    ch_ipcontrol_list
        .map { ipcontrol -> "${ipcontrol}" }
        .collectFile(name: 'ch_ipcontrol_list.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/INPUT_CHECK")


    // Add meta.is_input_control to ch_fastq
    ch_fastq
        .combine(ch_ipcontrol_list.ifEmpty([[]]))// .ifEmpty([[]]) is to avoid empty ch_fastq if there are no input controls
        .map { meta, fastqs, ipcontrol_list ->
            def meta_clone = meta.clone()
            meta_clone.is_input_control = ipcontrol_list.contains(meta.id)
            [ meta_clone, fastqs ]
        }
        .set { ch_fastq }

    // TODO: print for debugging
    ch_fastq.map { meta, fastqs -> "${meta}\t${fastqs}" }
        .collectFile(name: 'ch_fastq_5.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/INPUT_CHECK")
    
    // Create list of samples
    ch_fastq
        .map { meta, fastqs -> meta.id }
        .unique()
        .collect()
        .map { it -> [it] } // unflatten
        .set { ch_samples }

    // TODO: print for debugging
    ch_samples
        .map { sample_list -> "${sample_list}" }
        .collectFile(name: 'ch_samples.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/INPUT_CHECK")

    ch_ipcontrol_list
        .flatten()
        .combine(ch_samples)
        .map { ip_control, sample_list ->
            def ip_control_exists = sample_list.contains(ip_control)
            if (!ip_control_exists) {
                error(
                """
                ERROR:  The `input_control` sample '${ip_control}' does not exist in your sample list:

                ${sample_list}

                Please check your samplesheet and ensure all `input_control` values refer to an existing ID in the `sample` column.

                If you specified an `input_control_biological_replicate` or `input_control_technical_replicate`, also ensure that there is an input control with those values in the `biological_replicate` and `technical_replicate` columns.

                """.stripIndent()
                )
            }
            return [ ip_control ]
        }
        .collect()
        .map { it -> [it] } // unflatten
        .set { ch_ipcontrols }

    // TODO: print for debugging
    ch_ipcontrols
        .map { ip_control_list -> "${ip_control_list}" }
        .collectFile(name: 'ch_ipcontrols.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/INPUT_CHECK")
    
    // filter ch_fastq to only include samples whose meta.input_control is in ch_ipcontrols 
    ch_fastq
        .combine(ch_ipcontrols.ifEmpty([[]]))
        .filter { meta, fastqs, ipcontrol_list ->
            ipcontrol_list.contains(meta.input_control) || !meta.input_control
        }
        .map { meta, fastqs, ipcontrol_list ->
            // Input controls checks
            if (meta.is_input_control && meta.input_control) {
                error("ERROR: Input control samples must not have an `input_control` value set. Check sample: ${meta.id}")
            }
            // SCAR-seq and OK-seq checks
            if (['SCAR-seq', 'OK-seq'].contains(meta.exp_type)) {
                if (!meta.strandedness) {
                    error("ERROR: `strandedness` must be specified for SCAR-seq and OK-seq samples. Check sample: ${meta.id}")
                }
                if (meta.exp_type == 'SCAR-seq') {
                    if (!(params.containsKey('initiation_zones') || params.containsKey('okseq_rfd_file'))) {
                        if (params.refgenie_ignore && params.igenomes_ignore) {
                            error("ERROR: a SCAR-seq sample has been inputted, but neither `--initiation_zones` nor `--okseq_rfd_file` have been provided, and reference genomes are being ignored (`--refgenie_ignore true` and `--igenomes_ignore true`). You should provide either an initiation zones file or an OK-seq RFD file.")
                        } else if (!getGenomeAttribute('initiation_zones') && !getGenomeAttribute('okseq_rfd_file')) {
                            error("ERROR: a SCAR-seq sample has been inputted, but neither `--initiation_zones` nor `--okseq_rfd_file` have been found among reference genomes (iGenomes or Refgenie). You should provide either an initiation zones file or an OK-seq RFD file.")
                        }
                    error("ERROR: a SCAR-seq sample has been inputted, but neither `--initiation_zones` nor `--okseq_rfd_file` have been provided. You should provide either an initiation zones file or an OK-seq RFD file.")
                    }
                }
            } else if (meta.strandedness) {
                error("ERROR: `strandedness` must not be specified for samples other than SCAR-seq and OK-seq. Check sample: ${meta.id}")
            }
            // Antibody checks
            if (['ChIP-seq', 'ChIP-exo', 'ChOR-seq', 'SCAR-seq', 'CUTandTag', 'CUTandRUN', 'TIP-seq'].contains(meta.exp_type)) {
                if (!meta.antibody && !meta.is_input_control) {
                    log.warn("`antibody` should be specified for non-input ChIP-seq, ChIP-exo, ChOR-seq, SCAR-seq, CUTandTag, CUTandRUN, and TIP-seq samples. Check sample: ${meta.id}. Ignore this warning if the sample is an input control but you are not actually using it as the input control of another sample.")
                }
                if (meta.antibody && meta.is_input_control) {
                    error("ERROR: `antibody` must not be specified for input control samples. Check sample: ${meta.id}")
                }
            } else {
                if (meta.antibody) {
                    error("ERROR: `antibody` must not be specified for samples other than ChIP-seq, ChIP-exo, ChOR-seq, SCAR-seq, CUT&Tag, CUT&RUN, and TIP-seq. Check sample: ${meta.id}")
                }
            }
            return [ meta, fastqs ]
        }
        .set { ch_fastq }


    // Parse UMI meta fields
    ch_fastq
        .map { meta, fastqs ->
            if (!meta.with_umi) {
                if (meta.extract_umi) {
                    error("ERROR: `extract_umi` is set to true, but `with_umi` is set to false. Check sample: ${meta.id}")
                }
                if (meta.sep_umi_fq) {
                    if (!meta.extract_umi) {
                        log.warn("ERROR: `fastq_umi` has been provided in the samplesheet, but `extract_umi` is set to false. UMIs will not be transfered for sample: ${meta.id}")
                    }
                    log.warn("WARNING: `fastq_umi` has been provided in the samplesheet, but `with_umi` is set to false. UMIs will not be used for deduplication. Check sample: ${meta.id}")
                }
            }
            return [ meta, fastqs ]
        }
        .set { ch_fastq }

    // TODO: print for debugging
    ch_fastq.map { meta, fastqs -> "${meta}\t${fastqs}" }
        .collectFile(name: 'ch_fastq_6.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/INPUT_CHECK")


    emit:
    fastq = ch_fastq                                    // channel: [ val(meta), [ reads ] ]
    versions = channel.empty() // channel: [ versions.yml ]
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_report.toList(),
            )
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        } else {
            return null
        }
    } else {
        return null
    }
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

//
// Check and validate pipeline parameters
//
def validateInputParameters() {

    // commenting-out this is only a temporary fix to avoid
    // the pipeline from failing due to missing genome in igenomes
    //genomeExistsError()

    if (!params.fasta) {
        error("Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file.")
    }

    if (!params.gtf && !params.containsKey('gff')) {
        error("No GTF or GFF3 annotation specified! The pipeline requires at least one of these files.")
    }

    if (params.gtf && params.containsKey('gff')) {
        gtfGffWarn(log)
    }

    if (!params.skip_flTbl) {
        if (!params.containsKey('blacklist')) {
            if (params.refgenie_ignore && params.igenomes_ignore) {
                error("Blacklist filtering is enabled (`--skip_flTbl false`), a blacklist file (`--blacklist`) has not been provided, and reference genomes are being ignored (`--refgenie_ignore true` and `--igenomes_ignore true`). You should set the pipeline to skip blacklist filtering (`--skip_flTbl`) or provide a blacklist.")
            } else if (!getGenomeAttribute('blacklist')) {
                error("Blacklist filtering is enabled (`--skip_flTbl false`) but no valid blacklist file has been found among reference genomes (iGenomes or Refgenie). You should set the pipeline to skip blacklist filtering (`--skip_flTbl`) or provide a blacklist.")
            }
        error("Blacklist filtering is enabled (`--skip_flTbl false`) but no valid blacklist file has been provided. You should set the pipeline to skip blacklist filtering (`--skip_flTbl`) or provide a blacklist.")
        }
    }

    if (!params.skip_te_counting) {
        if (!params.containsKey('tecount_te_index') && !params.containsKey('te_gtf')) {
            if (params.refgenie_ignore && params.igenomes_ignore) {
                error("TE counting is enabled (`--skip_te_counting false`), a TEcount TE index file (`--tecount_te_index`) has not been provided, and reference genomes are being ignored (`--refgenie_ignore true` and `--igenomes_ignore true`). You should set the pipeline to skip TE counting (`--skip_te_counting`) or provide a TEcount TE index.")
            } else if (!getGenomeAttribute('tecount_te_index') && !getGenomeAttribute('te_gtf')) {
                error("TE counting is enabled (`--skip_te_counting false`) but no valid TEcount TE index file has been found among reference genomes (iGenomes or Refgenie). You should set the pipeline to skip TE counting (`--skip_te_counting`) or provide a TEcount TE index.")
            }
        }
        if (!params.skip_telocal) {
            if (!params.containsKey('telocal_te_index') && !params.containsKey('te_gtf')) {
                if (params.refgenie_ignore && params.igenomes_ignore) {
                    error("TElocal counting is enabled (`--skip_telocal false`), a TElocal TE index file (`--telocal_te_index`) has not been provided, and reference genomes are being ignored (`--refgenie_ignore true` and `--igenomes_ignore true`). You should set the pipeline to skip TE local counting (`--skip_telocal`) or provide a TElocal TE index.")
                } else if (!getGenomeAttribute('telocal_te_index') && !getGenomeAttribute('te_gtf')) {
                    error("TElocal counting is enabled (`--skip_telocal false`) but no valid TElocal TE index file has been found among reference genomes (iGenomes or Refgenie). You should set the pipeline to skip TE local counting (`--skip_telocal`) or provide a TElocal TE index.")
                }
            }
        }
    }

    if (!params.containsKey('macs_gsize')) {
        macsGsizeWarn(log)
    }

    if (params.umi_dedup_tool == 'umicollapse' && !['directional', 'adjacency', 'cluster'].contains(params.umi_grouping_method)) {
        error("The `-umi_grouping_method` parameter must be set to 'directional', 'adjacency' or 'cluster' when using the 'umicollapse' UMI deduplication tool.")
    }

    if (params.hardtrim3_length && params.hardtrim5_length) {
        error("Both `--hardtrim3_length` and `--hardtrim5_length` parameters have been provided. Please provide only one.")
    }

    if (!params.read_length && !params.macs_gsize) {
        error ("Both `--read_length` and `--macs_gsize` not specified! Please specify either to infer MACS3 genome size for peak calling.")
    }

    if (params.trim_q_cutoff && params.trim_nextseq) {
        error("Both `--trim_q_cutoff` and `--trim_nextseq` parameters have been provided (or one is set by default and the other was provided on top). Please provide only one.")
    }

    // if --read_length and either hardtrim3_length or hardtrim5_length are provided, then check if they are equal
    if (params.read_length && (params.hardtrim3_length || params.hardtrim5_length)) {
        if (params.hardtrim3_length && params.hardtrim3_length != params.read_length) {
            error("The `--read_length` and `--hardtrim3_length` parameters must be equal.")
        }
        if (params.hardtrim5_length && params.hardtrim5_length != params.read_length) {
            error("The `--read_length` and `--hardtrim5_length` parameters must be equal.")
        }
    }

    if (params.map_n_multimappers) {
        if (!['chromap', 'bowtie2', 'hisat2', 'star', 'bowtie', 'strobealign', 'minimap2'].contains(params.aligner)) {
            error("The `--map_n_multimappers` parameter requires the aligner to be set to 'chromap', 'bowtie2', 'hisat2', 'star', 'bowtie', 'strobealign', or 'minimap2'.")
        }
    }

    if (params.multimap_allocation_method == 'chromap' && params.aligner != 'chromap') {
        error("Allocating multimapping reads with 'chromap' requires the aligner to be set to 'chromap'.")
    }
    
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()
    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Print a warning if both GTF and GFF have been provided
//
def gtfGffWarn(log) {
    log.warn "=============================================================================\n" +
        "  Both '--gtf' and '--gff' parameters have been provided.\n" +
        "  Using GTF file as priority.\n" +
        "==================================================================================="
}

//
// Print a warning if macs_gsize parameter has not been provided
//
def macsGsizeWarn(log) {
    log.warn "=============================================================================\n" +
        "  --macs_gsize parameter has not been provided.\n" +
        "  It will be auto-calculated by 'khmer unique-kmers.py' using the '--read_length' parameter.\n" +
        "  Explicitly provide '--macs_gsize macs3_genome_size' to change this behaviour.\n" +
        "==================================================================================="
}