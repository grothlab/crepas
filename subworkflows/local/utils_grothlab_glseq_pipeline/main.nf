//
// Subworkflow with functionality specific to the grothlab/glseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../../subworkflows/nf-core/utils_nfvalidation_plugin'
include { UTILS_NFCORE_PIPELINE     } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../../subworkflows/nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion        } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { logColours                } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap          } from 'plugin/nf-validation'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

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
    pre_help_text = grothLabGlseqLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --genome GRCh37 --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
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
    // commenting-out this is only a temporary fix to avoid
    // the pipeline from failing due to missing genome in igenomes
    //validateInputParameters()

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
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
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
// Check and validate pipeline parameters
//
def validateInputParameters() {

    genomeExistsError()

    if (!params.fasta) {
        error("Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file.")
    }

    if (!params.gtf && !params.gff) {
        error("No GTF or GFF3 annotation specified! The pipeline requires at least one of these files.")
    }

    if (params.gtf && params.gff) {
        gtfGffWarn(log)
    }

    if (!params.macs_gsize) {
        macsGsizeWarn(log)
    }

    if (params.umi_dedup_tool == 'umicollapse' && !['directional', 'adjacency', 'cluster'].contains(params.umi_grouping_method)) {
        error("The '--umi_grouping_method' parameter must be set to 'directional', 'adjacency' or 'cluster' when using the 'umicollapse' UMI deduplication tool.")
    }

    if (params.hardtrim3_length && params.hardtrim5_length) {
        error("Both '--hardtrim3_length' and '--hardtrim5_length' parameters have been provided. Please provide only one.")
    }

    if (!params.read_length && !params.macs_gsize) {
        error ("Both '--read_length' and '--macs_gsize' not specified! Please specify either to infer MACS3 genome size for peak calling.")
    }

    if (params.trim_q_cutoff && params.trim_nextseq) {
        error("Both '--trim_q_cutoff' and '--trim_nextseq' parameters have been provided (or one is set by default and the other was provided on top). Please provide only one.")
    }

    // if --read_length and either hardtrim3_length or hardtrim5_length are provided, then check if they are equal
    if (params.read_length && (params.hardtrim3_length || params.hardtrim5_length)) {
        if (params.hardtrim3_length && params.hardtrim3_length != params.read_length) {
            error("The '--read_length' and '--hardtrim3_length' parameters must be equal.")
        }
        if (params.hardtrim5_length && params.hardtrim5_length != params.read_length) {
            error("The '--read_length' and '--hardtrim5_length' parameters must be equal.")
        }
    }

    if (params.map_n_multimappers) {
        if (!['chromap', 'bowtie2', 'hisat2'].contains(params.aligner)) {
            error("The '--map_n_multimappers' parameter requires the aligner to be set to 'chromap', 'bowtie2', or 'hisat2'.")
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
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
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


def grothLabGlseqLogo(monochrome_logs = true) {
    Map colors = logColours(monochrome_logs)
    String.format(
        """\n
        ${dashedLine(monochrome_logs)}
${colors.white}█████████████████████████████████████████████████████████████████████████████${colors.reset}
${colors.white}█▌/////////////////////////////////////////////////////////////////////////▐█${colors.reset}
${colors.white}█▌/////////////////////////////////////////////////////////////////////////▐█${colors.reset}
${colors.white}█▌////██████╗/██████╗//██████╗/████████╗██╗//██╗██╗//////█████╗/██████╗////▐█${colors.reset}
${colors.white}█▌///██╔════╝/██╔══██╗██╔═══██╗╚══██╔══╝██║//██║██║/////██╔══██╗██╔══██╗///▐█${colors.reset}
${colors.white}█▌///██║//███╗██████╔╝██║///██║///██║///███████║██║/////███████║██████╔╝///▐█${colors.reset}
${colors.white}█▌///██║///██║██╔══██╗██║///██║///██║///██╔══██║██║/////██╔══██║██╔══██╗///▐█${colors.reset}
${colors.white}█▌///╚██████╔╝██║//██║╚██████╔╝///██║///██║//██║███████╗██║//██║██████╔╝///▐█${colors.reset}
${colors.white}█▌////╚═════╝/╚═╝//╚═╝/╚═════╝////╚═╝///╚═╝//╚═╝╚══════╝╚═╝//╚═╝╚═════╝////▐█${colors.reset}
${colors.white}█▌/////////////////////////////////////////////////////////////////////////▐█${colors.reset}
${colors.white}█▌////////////////██████╗/██╗/////███████╗███████╗/██████╗/////////////////▐█${colors.reset}
${colors.white}█▌///////////////██╔════╝/██║/////██╔════╝██╔════╝██╔═══██╗////////////////▐█${colors.reset}
${colors.white}█▌///////////////██║//███╗██║/////███████╗█████╗//██║///██║////////////////▐█${colors.reset}
${colors.white}█▌///////////////██║///██║██║/////╚════██║██╔══╝//██║▄▄/██║////////////////▐█${colors.reset}
${colors.white}█▌///////////////╚██████╔╝███████╗███████║███████╗╚██████╔╝////////////////▐█${colors.reset}
${colors.white}█▌////////////////╚═════╝/╚══════╝╚══════╝╚══════╝/╚══▀▀═╝/////////////////▐█${colors.reset}
${colors.white}█▌/////////////////////////////////////////////////////////////////////////▐█${colors.reset}
${colors.white}█▌/////////////////////////////////////////////////////////////////////////▐█${colors.reset}
${colors.white}█████████████████████████████████████████████████████████████████████████████${colors.reset}

${colors.purple}  ${workflow.manifest.name} ${getWorkflowVersion()}${colors.reset}
        ${dashedLine(monochrome_logs)}
        """.stripIndent()
    )
}

//
// Return dashed line
//
def dashedLine(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    return "-${colors.dim}----------------------------------------------------${colors.reset}-"
}

//
// Citation string for pipeline
//
def workflowCitation() {
    def temp_doi_ref = ""
    String[] manifest_doi = workflow.manifest.doi.tokenize(",")
    // Using a loop to handle multiple DOIs
    // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
    // Removing ` ` since the manifest.doi is a string and not a proper list
    for (String doi_ref: manifest_doi) temp_doi_ref += "  https://doi.org/${doi_ref.replace('https://doi.org/', '').replace(' ', '')}\n"
    return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
        "* The pipeline\n" +
        temp_doi_ref + "\n" +
        "* The nf-core framework\n" +
        "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
        "* Software dependencies\n" +
        "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
}