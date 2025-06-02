//
// This file holds several functions specific to the workflow/cdnaseq.nf in the nf-core/cdnaseq pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowCdnaseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)

        if (!params.fasta) {
            Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMap(workflow, params) {
        def summary_log = ""
        def summary = [:]
        if (workflow.revision) summary['revision'] = workflow.revision
        summary['runName']      = workflow.runName
        if (workflow.containerEngine) summary['containerEngine'] = workflow.containerEngine
        if (workflow.container) summary['container'] = workflow.container
        summary['nextflow'] = [:]
        summary['nextflow']['version'] = workflow.nextflow.version
        summary['nextflow']['build'] = workflow.nextflow.build
        summary['nextflow']['timestamp'] = workflow.nextflow.timestamp
        summary['parameters'] = [:]
        summary['parameters']['input_samplesheet'] = params.input_samplesheet
        summary['parameters']['ref_dir'] = params.ref_dir
        summary['parameters']['output_dir'] = params.output_dir
        summary['parameters']['genome_build'] = params.genome_build
        summary['parameters']['perform_mut_ref_alignment'] = params.perform_mut_ref_alignment
        summary['parameters']['skip_fastqc'] = params.skip_fastqc
        summary['parameters']['skip_trimming'] = params.skip_trimming
        summary['parameters']['skip_variant_calling'] = params.skip_variant_calling
        summary['parameters']['skip_splicing_analysis'] = params.skip_splicing_analysis
        summary['parameters']['skip_quantification'] = params.skip_quantification
        summary['parameters']['skip_multiqc'] = params.skip_multiqc
        if (params.targets_bed) summary['parameters']['targets_bed'] = params.targets_bed
        if (params.adapter_fasta) summary['parameters']['adapter_fasta'] = params.adapter_fasta
        summary['parameters']['max_memory'] = params.max_memory
        summary['parameters']['max_cpus'] = params.max_cpus
        summary['parameters']['max_time'] = params.max_time
        
        return groovy.json.JsonBuilder(summary).toPrettyString()
    }

    //
    // Generate methods description for MultiQC
    //
    public static String methodsDescriptionText(workflow, mqc_config, params) {
        // TODO nf-core: Optional to add a citation of the pipeline here
        def meta = [:]
        meta.workflow = workflow
        meta.params   = params
        meta.mqc_config = mqc_config

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methodsDescriptionTemplate()).make(meta)

        return description_html.toString()
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

    //
    // Generate methods description template
    //
    private static String methodsDescriptionTemplate() {
        return """
        <h4>Methods</h4>
        <h3>Data processing</h3>
        <p>Raw sequencing reads were processed using the cdnaseq-nf pipeline (version ${workflow.manifest.version}). 
        Quality control of raw sequencing data was performed using FastQC (v${params.versions?.fastqc ?: 'N/A'})
        ${!params.skip_trimming ? "and adapter trimming was conducted using BBDuk from BBTools suite" : ""}.
        Reads were aligned to the ${params.genome_build ?: 'reference'} genome using STAR aligner 
        ${params.perform_mut_ref_alignment ? "including alignment to patient-specific mutated references generated from provided VCF files" : ""}.
        ${!params.skip_star_second_pass ? "A two-pass alignment strategy was employed where splice junctions identified in the first pass were used to guide the second pass alignment." : "A single-pass alignment strategy was employed for fast processing."}
        ${!params.skip_variant_calling ? "Variant calling was performed using GATK HaplotypeCaller in RNA-seq mode." : ""}
        ${!params.skip_splicing_analysis ? "Splicing analysis was conducted using Regtools for junction extraction and annotation, and Portcullis for comprehensive splice site analysis." : ""}
        ${!params.skip_quantification ? "Gene expression quantification was performed using featureCounts." : ""}
        All quality control metrics were summarized using MultiQC (v${params.versions?.multiqc ?: 'N/A'}).
        </p>
        """
    }
}
