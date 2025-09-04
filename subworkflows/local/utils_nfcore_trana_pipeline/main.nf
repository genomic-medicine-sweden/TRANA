//
// Subworkflow with functionality specific to the genomic-medicine-sweden/TRANA pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENERATE_INPUT             } from '../../../modules/local/generate_input/main.nf'
include { MERGE_BARCODES             } from '../../../modules/local/merge_barcodes/main.nf'
include { MERGE_BARCODES_SAMPLESHEET } from '../../../modules/local/merge_barcodes_samplesheet/main.nf'

include { UTILS_NFSCHEMA_PLUGIN      } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap           } from 'plugin/nf-schema'
include { samplesheetToList          } from 'plugin/nf-schema'
include { completionEmail            } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary          } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification             } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE      } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE    } from '../../nf-core/utils_nextflow_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version                 // boolean: Display version and exit
    validate_params         // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs         // boolean: Do not use coloured log outputs
    nextflow_cli_args       // array: List of positional nextflow CLI args
    outdir                  // string: The output directory where the results will be saved
    samplesheet             // string: Path to input samplesheet
    barcodes_samplesheet    // string: Path to input barcodes samplesheet
    merge_fastq_pass        // string: Path to directory ONT fastq_pass dir

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

    // Check input path parameters to see if they exist
    def checkPathParamList = !merge_fastq_pass ? [params.input] : []
    checkPathParamList += [params.multiqc_config, params.fasta]
    checkPathParamList.findAll { it }.each { path ->
        if (!file(path).exists()) {
            exit 1, "ERROR: File path does not exist: ${path}"
        }
    }

    //
    // MODULE: Concatenate input read files
    //
    if ( merge_fastq_pass && !barcodes_samplesheet ) {
        MERGE_BARCODES(merge_fastq_pass)
        ch_versions = ch_versions.mix(MERGE_BARCODES.out.versions)
        GENERATE_INPUT(MERGE_BARCODES.out.fastq_dir_merged).sample_sheet_merged.set{ ch_samplesheet_path }
        ch_versions = ch_versions.mix(GENERATE_INPUT.out.versions)

    } else if ( merge_fastq_pass && barcodes_samplesheet ) {
        MERGE_BARCODES_SAMPLESHEET(barcodes_samplesheet, merge_fastq_pass)
        ch_versions = ch_versions.mix(MERGE_BARCODES_SAMPLESHEET.out.versions)
        GENERATE_INPUT(MERGE_BARCODES_SAMPLESHEET.out.fastq_dir_merged).sample_sheet_merged.set{ ch_samplesheet_path }
        ch_versions = ch_versions.mix(GENERATE_INPUT.out.versions)
    } else if ( !merge_fastq_pass && !barcodes_samplesheet && samplesheet ) {
        ch_samplesheet_path = channel.value(samplesheet)
    } else {
        error "Invalid input. Please specify either '--input' or '--merge_fastq_pass' (and '--barcodes_samplesheet' if available)."
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_samplesheet_path
        .flatMap { samplesheet_path ->
            samplesheetToList(samplesheet_path, "${projectDir}/assets/schema_input.json")
        }
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { validateInputSamplesheet( it ) }
        .set { ch_reads }

    emit:
    reads       = ch_reads            // channel: [ val(meta), [ reads ] ]
    samplesheet = ch_samplesheet_path // channel: [ val(meta), [ samplesheet ] ]
    versions    = ch_versions         // channel: [ versions.yml ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    workflow.onComplete {
        //
        // Completion email and summary
        //
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report)
        }
        completionSummary(monochrome_logs)
        //
        // Instant message notification
        //
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs[0] ]
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
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
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

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {

    println "ROW: ${row}"

    // Create meta map
    def meta = [:]
    meta.id                 = row.sample
    meta.sequencing_run     = row.sequencing_run
    meta.single_end         = row.single_end.toBoolean()
    meta.fw_primer          = row.FW_primer
    meta.rv_primer          = row.RV_primer

    // Add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

