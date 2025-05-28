/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS/FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WorkflowCdnaseq } from '../lib/WorkflowCdnaseq'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters
if (params.input_samplesheet) { 
    ch_input_samplesheet = file(params.input_samplesheet, checkIfExists: true) 
} else { 
    exit 1, 'Input samplesheet not specified!' 
}

// Check reference directory
if (params.ref_dir) {
    if (!file(params.ref_dir).exists()) {
        exit 1, "Reference directory does not exist: ${params.ref_dir}"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MultiQC config file
ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_dummy_file              = Channel.fromPath("$projectDir/assets/adapters.fa", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local modules
//
include { FASTQC_RAW                  } from '../modules/local/fastqc'
include { FASTQC_TRIM                 } from '../modules/local/fastqc'
include { BBDUK_TRIM                  } from '../modules/local/bbduk_trim'
include { CREATE_MUT_REF_FASTA        } from '../modules/local/gatk_tools'
include { INDEX_FASTA                 } from '../modules/local/gatk_tools'
include { CREATE_SEQ_DICT             } from '../modules/local/gatk_tools'
include { SPLITNCIGARREADS            } from '../modules/local/gatk_tools'
include { GATK_HAPLOTYPECALLER_RNA    } from '../modules/local/gatk_tools'
include { STAR_INDEX                  } from '../modules/local/star_index'
include { STAR_ALIGN_P1_REF           } from '../modules/local/star_align'
include { STAR_ALIGN_P1_MUT           } from '../modules/local/star_align'
include { AGGREGATE_SJ_P1             } from '../modules/local/star_align'
include { STAR_ALIGN_P2               } from '../modules/local/star_align'
include { MARK_DUPLICATES_PICARD      } from '../modules/local/bam_processing'
include { SAMTOOLS_INDEX              } from '../modules/local/bam_processing'
include { SAMTOOLS_SORT               } from '../modules/local/bam_processing'
include { QUALIMAP_RNASEQ             } from '../modules/local/bam_processing'
include { REGTOOLS_JUNCTIONS_EXTRACT  } from '../modules/local/splicing_tools'
include { REGTOOLS_JUNCTIONS_ANNOTATE } from '../modules/local/splicing_tools'
include { PORTCULLIS                  } from '../modules/local/splicing_tools'
include { FEATURECOUNTS               } from '../modules/local/quantification'
include { MULTIQC                     } from '../modules/local/multiqc'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/local/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { TRIMGALORE                  } from '../modules/nf-core/trimgalore/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CDNASEQ {

    ch_versions = Channel.empty()

    //
    // Parse samplesheet
    //
    Channel
        .fromPath(params.input_samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [:]
            meta.id = row.sample_id
            meta.patient_id = row.patient_id
            meta.panel = row.panel
            meta.strandedness = row.strandedness ?: 'unstranded'
            meta.single_end = false  // Assuming paired-end by default

            def reads = []
            if (row.fastq_r1) reads.add(file(row.fastq_r1, checkIfExists: true))
            if (row.fastq_r2) reads.add(file(row.fastq_r2, checkIfExists: true))
            
            def vcf = row.vcf && row.vcf != '' ? file(row.vcf, checkIfExists: true) : null
            
            return [meta, reads, vcf]
        }
        .set { ch_input }

    // Set up reference files
    def reference_fasta = file("${params.ref_dir}/*.primary_assembly.genome.fa", checkIfExists: true)[0]
    def reference_gtf = file("${params.ref_dir}/*.annotation.gtf", checkIfExists: true)[0]
    def reference_fai = file("${params.ref_dir}/*.primary_assembly.genome.fa.fai", checkIfExists: true)[0]
    def reference_dict = file("${params.ref_dir}/*.primary_assembly.genome.dict", checkIfExists: true)[0]

    //
    // MODULE: Run FastQC on raw reads
    //
    if (!params.skip_fastqc) {
        FASTQC_RAW (
            ch_input.map { meta, reads, vcf -> [meta, reads] }
        )
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    //
    // MODULE: Trim adapters with BBDuk
    //
    if (!params.skip_trimming) {
        BBDUK_TRIM (
            ch_input.map { meta, reads, vcf -> [meta, reads] },
            file(params.adapter_fasta, checkIfExists: true)
        )
        ch_trimmed_reads = BBDUK_TRIM.out.reads
        ch_versions = ch_versions.mix(BBDUK_TRIM.out.versions.first())

        // Run FastQC on trimmed reads
        if (!params.skip_fastqc) {
            FASTQC_TRIM (
                ch_trimmed_reads
            )
        }
    } else {
        ch_trimmed_reads = ch_input.map { meta, reads, vcf -> [meta, reads] }
    }

    //
    // Create mutated references for patients with VCFs
    //
    ch_input
        .filter { meta, reads, vcf -> vcf != null }
        .map { meta, reads, vcf -> 
            def vcf_index = file(vcf.toString() + ".tbi", checkIfExists: false)
            return [meta, vcf, vcf_index]
        }
        .groupTuple(by: [0])
        .map { meta, vcfs, vcf_indices -> [meta, vcfs[0], vcf_indices[0]] }
        .set { ch_vcf_for_mutref }

    if (params.perform_mut_ref_alignment) {
        CREATE_MUT_REF_FASTA (
            ch_vcf_for_mutref,
            reference_fasta,
            reference_fai,
            reference_dict
        )
        ch_versions = ch_versions.mix(CREATE_MUT_REF_FASTA.out.versions.first())
    }

    //
    // STAR Index Creation
    //
    STAR_INDEX (
        reference_fasta,
        reference_gtf
    )
    ch_star_index = STAR_INDEX.out.index
    ch_versions = ch_versions.mix(STAR_INDEX.out.versions.first())

    //
    // STAR First Pass Alignment - Reference
    //
    STAR_ALIGN_P1_REF (
        ch_trimmed_reads,
        ch_star_index,
        reference_gtf
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_P1_REF.out.versions.first())

    //
    // STAR First Pass Alignment - Mutated Reference (if VCF provided)
    //
    if (params.perform_mut_ref_alignment) {
        // Create STAR index for mutated reference
        CREATE_MUT_REF_FASTA.out.fasta
            .combine(Channel.value(reference_gtf))
            .set { ch_mut_ref_for_star }

        STAR_INDEX.alias('STAR_INDEX_MUT') (
            CREATE_MUT_REF_FASTA.out.fasta,
            reference_gtf
        )
        ch_star_index_mut = STAR_INDEX.alias('STAR_INDEX_MUT').out.index

        // Align to mutated reference - only for samples with VCF
        ch_trimmed_reads
            .join(ch_input.map { meta, reads, vcf -> [meta.id, vcf] }, by: [0])
            .filter { meta, reads, vcf -> vcf != null }
            .map { meta, reads, vcf -> [meta, reads] }
            .combine(ch_star_index_mut)
            .set { ch_reads_mut_index }

        STAR_ALIGN_P1_MUT (
            ch_reads_mut_index.map { meta, reads, index -> [meta, reads] },
            ch_reads_mut_index.map { meta, reads, index -> index },
            reference_gtf
        )
    }

    //
    // Aggregate Splice Junctions from First Pass
    //
    ch_sj_p1_ref = STAR_ALIGN_P1_REF.out.sj_tab
    
    if (params.perform_mut_ref_alignment) {
        ch_sj_p1_mut = STAR_ALIGN_P1_MUT.out.sj_tab
        ch_sj_p1_all = ch_sj_p1_ref.mix(ch_sj_p1_mut)
    } else {
        ch_sj_p1_all = ch_sj_p1_ref
    }

    AGGREGATE_SJ_P1 (
        ch_sj_p1_all.map { meta, sj -> sj }.collect()
    )
    ch_filtered_sj = AGGREGATE_SJ_P1.out.filtered_sj
    ch_versions = ch_versions.mix(AGGREGATE_SJ_P1.out.versions.first())

    //
    // STAR Second Pass Alignment
    //
    STAR_ALIGN_P2 (
        ch_trimmed_reads,
        ch_star_index,
        reference_gtf,
        ch_filtered_sj
    )
    ch_bam_sorted = STAR_ALIGN_P2.out.bam_sorted
    ch_versions = ch_versions.mix(STAR_ALIGN_P2.out.versions.first())

    //
    // BAM Processing Pipeline
    //
    
    // Mark duplicates
    MARK_DUPLICATES_PICARD (
        ch_bam_sorted
    )
    ch_bam_markdup = MARK_DUPLICATES_PICARD.out.bam
    ch_versions = ch_versions.mix(MARK_DUPLICATES_PICARD.out.versions.first())

    // Index BAM files
    SAMTOOLS_INDEX (
        ch_bam_markdup
    )
    ch_bam_indexed = ch_bam_markdup.join(SAMTOOLS_INDEX.out.bai)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // Variant Calling Pipeline
    //
    if (!params.skip_variant_calling) {
        // Split N CIGAR reads
        SPLITNCIGARREADS (
            ch_bam_indexed.map { meta, bam, bai -> [meta, bam, bai] },
            reference_fasta,
            reference_fai,
            reference_dict
        )
        ch_bam_split = SPLITNCIGARREADS.out.bam
        ch_versions = ch_versions.mix(SPLITNCIGARREADS.out.versions.first())

        // GATK HaplotypeCaller
        GATK_HAPLOTYPECALLER_RNA (
            ch_bam_split,
            reference_fasta,
            reference_fai,
            reference_dict,
            params.targets_bed ? file(params.targets_bed, checkIfExists: true) : []
        )
        ch_versions = ch_versions.mix(GATK_HAPLOTYPECALLER_RNA.out.versions.first())
    }

    //
    // Splicing Analysis
    //
    if (!params.skip_splicing_analysis) {
        // Extract junctions with regtools
        REGTOOLS_JUNCTIONS_EXTRACT (
            ch_bam_indexed.map { meta, bam, bai -> [meta, bam, bai] }
        )
        ch_versions = ch_versions.mix(REGTOOLS_JUNCTIONS_EXTRACT.out.versions.first())

        // Annotate junctions with regtools
        REGTOOLS_JUNCTIONS_ANNOTATE (
            REGTOOLS_JUNCTIONS_EXTRACT.out.junctions,
            reference_fasta,
            reference_gtf
        )
        ch_versions = ch_versions.mix(REGTOOLS_JUNCTIONS_ANNOTATE.out.versions.first())

        // Run Portcullis for comprehensive splice junction analysis
        PORTCULLIS (
            ch_bam_indexed.map { meta, bam, bai -> [meta, bam] },
            reference_fasta
        )
        ch_versions = ch_versions.mix(PORTCULLIS.out.versions.first())
    }

    //
    // Quantification
    //
    if (!params.skip_quantification) {
        FEATURECOUNTS (
            ch_bam_indexed.map { meta, bam, bai -> [meta, bam] },
            reference_gtf
        )
        ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions.first())
    }

    //
    // Collect software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MultiQC Report
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowCdnaseq.paramsSummaryMap(workflow, params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowCdnaseq.methodsDescriptionText(workflow, ch_dummy_file, params)
        ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        
        if (!params.skip_fastqc) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
            if (!params.skip_trimming) {
                ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIM.out.zip.collect{it[1]}.ifEmpty([]))
            }
        }
        
        if (!params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(BBDUK_TRIM.out.log.collect{it[1]}.ifEmpty([]))
        }
        
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_P2.out.log_final.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(MARK_DUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]))
        
        if (!params.skip_quantification) {
            ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]))
        }

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)
    }

    // Print completion message
    log.info """
    ========================================
    cdnaseq-nf pipeline execution completed
    ========================================
    Output directory: ${params.output_dir}
    """.stripIndent()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}
