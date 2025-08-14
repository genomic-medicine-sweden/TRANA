/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                       } from 'plugin/nf-schema'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline/main.nf'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline/main.nf'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_taco_pipeline/main.nf'
include { GENERATE_MASTER_HTML                   } from '../modules/local/generate_master_html/main.nf'
include { EMU_ABUNDANCE                          } from '../modules/local/emu/abundance/main.nf'
include { KRONA_KTIMPORTTAXONOMY                 } from '../modules/nf-core/krona/ktimporttaxonomy/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main.nf'
include { FASTQC                                 } from '../modules/nf-core/fastqc/main.nf'
include { CUTADAPT                               } from '../modules/nf-core/cutadapt/main.nf'
include { NANOPLOT as NANOPLOT_UNPROCESSED_READS } from '../modules/nf-core/nanoplot/main.nf'
include { NANOPLOT as NANOPLOT_PROCESSED_READS   } from '../modules/nf-core/nanoplot/main.nf'
include { PORECHOP_ABI                           } from '../modules/nf-core/porechop/abi/main.nf'
include { SEQTK_SAMPLE                           } from '../modules/nf-core/seqtk/sample/main.nf'
include { FILTLONG                               } from '../modules/nf-core/filtlong/main.nf'
include { EMU_COMBINE_OUTPUTS                    } from '../modules/local/emu/combine_outputs/main.nf'
include { COMBINE_REPORTS                        } from '../modules/local/phyloseq/main.nf'
include { PHYLOSEQ_OBJECT                        } from '../modules/local/phyloseq/main.nf'
include { ASSIGNMENT_HEATMAP                     } from '../modules/local/assignment_heatmap/main.nf'
include { CTRL_COMPARISON                        } from '../modules/local/ctrl_comparison/main.nf'
include { TRANSLATE_TAXIDS                       } from '../modules/local/translate_taxids/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TACO {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input
    ch_reads        // channel: reads from PIPELINE_INITIALISATION

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // MODULE: Run FastQC
    //
    FASTQC (ch_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    if (params.seqtype == "map-ont") {

        //
        // MODULE: Run NANOPLOT_UNPROCESSED_READS
        //
        NANOPLOT_UNPROCESSED_READS(ch_reads)
        ch_versions = ch_versions.mix(NANOPLOT_UNPROCESSED_READS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            NANOPLOT_UNPROCESSED_READS.out.txt.collect{ it[1] }
        )

        //
        // Optional adapter trimming
        //
        if (params.adapter_trimming) {
            //
            // MODULE: Run Porechop to trim ONT adapters
            //
            PORECHOP_ABI(ch_reads, [])

            PORECHOP_ABI.out.reads.map {
                meta, reads -> [meta + [single_end: 1], reads]
            }.set{ ch_optionally_trimmed_reads }

            ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(
                PORECHOP_ABI.out.log.collect{ it[1] }
            )
        } else {
            ch_reads.set{ ch_optionally_trimmed_reads }
        }

        //
        // Optional quality filtering
        //
        if (params.quality_filtering) {
            //
            // MODULE: Run filtlong to filter on read length
            //
            FILTLONG(ch_optionally_trimmed_reads.map {
                meta, reads -> [meta, [], reads]
            }).reads.set { ch_processed_reads }

            ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(
                FILTLONG.out.log.collect{ it[1] }
            )
        } else {
            ch_optionally_trimmed_reads.set{ ch_processed_reads }
        }

        //
        // MODULE: run NANOPLOT_PROCESSED_READS
        //
        NANOPLOT_PROCESSED_READS(ch_processed_reads)
        ch_versions = ch_versions.mix(NANOPLOT_PROCESSED_READS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            NANOPLOT_PROCESSED_READS.out.txt.collect{ it[1] }
        )

    } else if (params.seqtype == "sr") {

        //
        // MODULE: Run cutadapt to trim short-reads adapters
        //
        if (!params.skip_cutadapt) {
            CUTADAPT(ch_reads)
            CUTADAPT.out.reads.set { ch_processed_reads }
            ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(
                CUTADAPT.out.log.collect{ it[1] }
            )
        } else {
            ch_reads.set { ch_processed_reads }
        }
    } else {
        error "Invalid seqtype. Please specify either 'map-ont' or 'sr'."
    }

    if ( params.sample_size ) {
        //
        // MODULE: Downsample reads
        SEQTK_SAMPLE(
            ch_processed_reads.map { meta, reads -> tuple(meta, reads, params.sample_size) }
            ).reads.set { ch_processed_optionally_sampled_reads }
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
    } else {
        ch_processed_reads.set{ ch_processed_optionally_sampled_reads }
    }

    //
    // MODULE: run EMU abundance calculation
    //
    EMU_ABUNDANCE(ch_processed_optionally_sampled_reads)
    ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions)
    if (params.run_krona) {
        //
        // MODULE: Run krona plot
        //
        KRONA_KTIMPORTTAXONOMY(
            EMU_ABUNDANCE.out.report,
            file(params.krona_taxonomy_tab, checkExists: true)
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY.out.versions.first())
    }

    //MODULE: run translate_taxids
    if (params.make_heatmap) {
        TRANSLATE_TAXIDS(EMU_ABUNDANCE.out.assignment_report)
        ch_versions = ch_versions.mix(TRANSLATE_TAXIDS.out.versions)

        //Module: run assignment_heatmap
        ASSIGNMENT_HEATMAP(TRANSLATE_TAXIDS.out.assignment_translated_report)
        ch_versions = ch_versions.mix(ASSIGNMENT_HEATMAP.out.versions)
    }
    // MODULE: run emu combine-outputs
    ch_emu_combine_input_files = Channel.empty()
    // Collect all reports into a single list containing the paths
    ch_emu_combine_input_files = EMU_ABUNDANCE.out.report
        .map { it[1] }  // Extract only the file path from the tuple (meta, path)
        .collect()
        .set { collected_files }
    EMU_COMBINE_OUTPUTS(collected_files)
    ch_versions = ch_versions.mix(EMU_COMBINE_OUTPUTS.out.versions)

    // MODULE: run ctrl_comparison
    if (params.ctrl_1) {
        CTRL_COMPARISON(EMU_COMBINE_OUTPUTS.out.combined_report,EMU_COMBINE_OUTPUTS.out.combined_counts_report)
        ch_versions = ch_versions.mix(CTRL_COMPARISON.out.versions)
    }
    //
    // MODULE: Generate PHYLOSEQ object
    if (params.phyloseq)  {
        ch_tax_file = Channel.fromPath("$projectDir/assets/databases/emu_database/taxonomy.tsv", checkIfExists: true)
        report_ch = EMU_ABUNDANCE.out.report
        all_reports_ch = report_ch
            .map { meta, path -> path }
            .collect()
        COMBINE_REPORTS (
            all_reports_ch
        )

        ch_tax_file = Channel.fromPath(
            "$projectDir/assets/databases/emu_database/taxonomy.tsv",
            checkIfExists: true
        )
        PHYLOSEQ_OBJECT (
            COMBINE_REPORTS.out.combinedreport,
            ch_tax_file
        )

        ch_versions = ch_versions.mix(PHYLOSEQ_OBJECT.out.versions)
    }

    // collect tool versions.
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT REPORTS & MultiQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    //
    // MODULE: MultiQC
    //

    ch_multiqc_config           = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config    = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo             = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_workflow_summary         = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_methods_description      = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    multiqc_report = MULTIQC.out.report.toList()

    //
    // MODULE: Generate master.html from reports
    //

    ch_reads
        .map {
            meta, reads -> meta
        }
        .first()
        .set {ch_meta}
    GENERATE_MASTER_HTML(ch_meta, ch_samplesheet)
    ch_versions = ch_versions.mix(GENERATE_MASTER_HTML.out.versions)

    emit:
    master_html             = GENERATE_MASTER_HTML.out.html      // channel: [ path(master.html) ]
    versions                = ch_versions                        // channel: [ path(versions.yml) ]
    nanostats_unprocessed   = (params.seqtype == "map-ont") ? NANOPLOT_UNPROCESSED_READS.out.txt : Channel.empty()  // channel: [ path(master.html) ]
    nanostats_processed     = (params.seqtype == "map-ont") ? NANOPLOT_PROCESSED_READS.out.txt   : Channel.empty()  // channel: [ path(master.html) ]
    multiqc_report          = multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
