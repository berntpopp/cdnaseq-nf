/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS/FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WorkflowCdnaseq class from lib/ is auto-loaded by Nextflow

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
include { STAR_INDEX as STAR_INDEX_STD} from '../modules/local/star_index'
include { STAR_INDEX as STAR_INDEX_MUT} from '../modules/local/star_index'
include { STAR_INDEX as STAR_INDEX_P2 } from '../modules/local/star_index'
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
    IMPORT ADDITIONAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local modules 
//
// Remove unused aliases that can cause confusion
// Previously: include { FASTQC_RAW as FASTQC        } from '../modules/local/fastqc'
// Previously: include { BBDUK_TRIM as TRIMGALORE    } from '../modules/local/bbduk_trim'

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
    ch_empty_sj_file = Channel.empty()
    
    STAR_INDEX_STD (
        reference_fasta,
        reference_gtf,
        params.sjdb_overhang,
        "standard_ref_index",
        ch_empty_sj_file
    )
    ch_star_index_std = STAR_INDEX_STD.out.index
    ch_versions = ch_versions.mix(STAR_INDEX_STD.out.versions.first())

    //
    // STAR First Pass Alignment - Reference
    //
    // 1. Create a VCF info channel keyed by sample ID
    ch_input 
        .map { meta_orig, _, vcf_file_orig -> // reads_orig are not needed here directly
            def vcf_idx_orig = vcf_file_orig ? file(vcf_file_orig.toString() + ".tbi", checkIfExists: false) : null
            [ meta_orig.id, meta_orig, vcf_file_orig, vcf_idx_orig ] // [id, meta, vcf, vcf_index]
        }
        .set { ch_vcf_info_by_id }

    // 2. Prepare ch_trimmed_reads for joining (assuming ch_trimmed_reads is [meta_obj, reads_list])
    ch_trimmed_reads
        .map { meta_trim, reads_trim_list -> 
            [ meta_trim.id, meta_trim, reads_trim_list ] // [id, meta, reads_list]
        }
        .set { ch_trimmed_reads_by_id }

    // 3. Join trimmed reads with VCF info
    ch_trimmed_reads_by_id
        .join( ch_vcf_info_by_id, by: [0] ) // Join on meta.id
        .map { id, meta_trim_obj, reads_trim_list, _, vcf_file_final, vcf_idx_final -> // _ is ignored meta from vcf_info
            // Ensure vcf_file_final and vcf_idx_final are passed as null if they were originally null
            [ meta_trim_obj, reads_trim_list, vcf_file_final, vcf_idx_final ] 
        }
        .set { ch_reads_and_vcf_for_alignment } // Emits: [meta, reads_list, vcf_path_or_null, vcf_index_path_or_null]

    // --- STAR_ALIGN_P1_REF ---
    STAR_ALIGN_P1_REF (
        ch_reads_and_vcf_for_alignment.map { meta, reads_list, vcf_file, vcf_idx_file -> [meta, reads_list] }, // meta, reads
        ch_star_index_std,
        ch_reads_and_vcf_for_alignment.map { meta, reads_list, vcf_file, vcf_idx_file -> vcf_file ?: [] },       // vcf
        ch_reads_and_vcf_for_alignment.map { meta, reads_list, vcf_file, vcf_idx_file -> vcf_idx_file ?: [] }    // vcf_index
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

        STAR_INDEX_MUT (
            CREATE_MUT_REF_FASTA.out.fasta,
            reference_gtf,
            params.sjdb_overhang,
            CREATE_MUT_REF_FASTA.out.fasta.map { meta, fasta -> "${meta.patient_id}_mut_index" },
            ch_empty_sj_file
        )
        ch_star_indices_mut = STAR_INDEX_MUT.out.index

        // --- STAR_ALIGN_P1_MUT ---
        // ch_star_indices_mut is assumed to be [meta_patient_obj_from_mutref, mut_star_index_path]
        // We need to join based on patient_id
        
        ch_reads_and_vcf_for_alignment // [meta_rna, reads_list, vcf_file, vcf_idx_file]
            .map { meta_rna, reads_list, vcf, vcf_idx -> 
                [meta_rna.patient_id, meta_rna, reads_list, vcf, vcf_idx] // Key by patient_id
            }
            .join( ch_star_indices_mut.map { meta_mut_idx, mut_idx_path -> [meta_mut_idx.patient_id, mut_idx_path] } , by: [0] ) // Join on patient_id
            .map { patient_id, meta_rna_obj, reads_list, vcf_file, vcf_idx_file, mut_idx_path ->
                 // Structure for STAR_ALIGN_P1_MUT: [meta, reads, index, vcf_optional, vcf_index_optional]
                [ meta_rna_obj, reads_list, mut_idx_path, vcf_file, vcf_idx_file ]
            }
            .set { ch_for_star_align_p1_mut_input }

        STAR_ALIGN_P1_MUT (
            ch_for_star_align_p1_mut_input.map { meta, reads_list, mut_idx_path, vcf_file, vcf_idx_file -> [meta, reads_list] }, // meta, reads
            ch_for_star_align_p1_mut_input.map { meta, reads_list, mut_idx_path, vcf_file, vcf_idx_file -> mut_idx_path },       // index (mutated)
            ch_for_star_align_p1_mut_input.map { meta, reads_list, mut_idx_path, vcf_file, vcf_idx_file -> vcf_file ?: [] },     // vcf
            ch_for_star_align_p1_mut_input.map { meta, reads_list, mut_idx_path, vcf_file, vcf_idx_file -> vcf_idx_file ?: [] }  // vcf_index
        )
    }

    //
    // Aggregate Splice Junctions from First Pass
    //
    ch_sj_p1_ref = STAR_ALIGN_P1_REF.out.sj
    
    if (params.perform_mut_ref_alignment) {
        ch_sj_p1_mut = STAR_ALIGN_P1_MUT.out.sj
        ch_sj_p1_all = ch_sj_p1_ref.mix(ch_sj_p1_mut)
    } else {
        ch_sj_p1_all = ch_sj_p1_ref
    }

    // Create unified BAM channel for downstream processing
    if (params.skip_star_second_pass) {
        // Use first-pass alignments as final
        // Prioritize mutated reference alignment if available, otherwise use reference alignment
        if (params.perform_mut_ref_alignment) {
            ch_p1_ref_bams = STAR_ALIGN_P1_REF.out.bam
            ch_p1_mut_bams = STAR_ALIGN_P1_MUT.out.bam
            
            // Create a channel with all P1 BAMs, tagged by alignment type
            ch_p1_ref_tagged = ch_p1_ref_bams.map { meta, bam -> 
                [meta.patient_id ?: meta.id, meta, bam, 'ref'] 
            }
            ch_p1_mut_tagged = ch_p1_mut_bams.map { meta, bam -> 
                [meta.patient_id ?: meta.id, meta, bam, 'mut'] 
            }
            
            // Combine and select best alignment per sample
            ch_all_p1_bams = ch_p1_ref_tagged.mix(ch_p1_mut_tagged)
                .groupTuple(by: 0) // Group by patient_id/sample_id
                .map { patient_id, metas, bams, types ->
                    // Prefer mutated alignment if available, otherwise use reference
                    def mut_idx = types.findIndexOf { it == 'mut' }
                    def selected_idx = mut_idx >= 0 ? mut_idx : 0
                    [metas[selected_idx], bams[selected_idx]]
                }
            ch_final_aligned_bams_for_processing = ch_all_p1_bams
        } else {
            // Only reference alignment available
            ch_final_aligned_bams_for_processing = STAR_ALIGN_P1_REF.out.bam
        }
    } else {
        // Traditional two-pass approach
        AGGREGATE_SJ_P1 (
            ch_sj_p1_all.map { meta, sj -> sj }.collect()
        )
        ch_filtered_sj = AGGREGATE_SJ_P1.out.filtered_sj
        ch_versions = ch_versions.mix(AGGREGATE_SJ_P1.out.versions.first())

        //
        // STAR Index for P2 with aggregated junctions
        //
        STAR_INDEX_P2 (
            reference_fasta,
            reference_gtf,
            params.sjdb_overhang,
            "p2_index_with_junctions",
            ch_filtered_sj
        )
        ch_star_index_p2 = STAR_INDEX_P2.out.index

        //
        // STAR Second Pass Alignment
        //
        STAR_ALIGN_P2 (
            ch_reads_and_vcf_for_alignment.map { meta, reads_list, vcf_file, vcf_idx_file -> [meta, reads_list] }, // meta, reads
            ch_star_index_p2,
            reference_gtf,
            ch_filtered_sj,
            ch_reads_and_vcf_for_alignment.map { meta, reads_list, vcf_file, vcf_idx_file -> vcf_file ?: [] },     // vcf_optional
            ch_reads_and_vcf_for_alignment.map { meta, reads_list, vcf_file, vcf_idx_file -> vcf_idx_file ?: [] }  // vcf_index_optional
        )
        ch_final_aligned_bams_for_processing = STAR_ALIGN_P2.out.bam_sorted
        ch_versions = ch_versions.mix(STAR_ALIGN_P2.out.versions.first())
    }

    //
    // BAM Processing Pipeline
    //
    
    // Mark duplicates
    MARK_DUPLICATES_PICARD (
        ch_final_aligned_bams_for_processing
    )
    ch_bam_markdup = MARK_DUPLICATES_PICARD.out.bam
    ch_versions = ch_versions.mix(MARK_DUPLICATES_PICARD.out.versions.first())

    // Index BAM files
    SAMTOOLS_INDEX (
        ch_bam_markdup.map { meta, bam, bai -> [meta, bam] } // Extract just meta and bam 
    )
    ch_bam_indexed = SAMTOOLS_INDEX.out.bam // This already contains [meta, bam, bai]
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
        
        // Include STAR alignment logs based on skip_star_second_pass setting
        if (params.skip_star_second_pass) {
            ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_P1_REF.out.log_final.collect{it[1]}.ifEmpty([]))
            if (params.perform_mut_ref_alignment) {
                ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_P1_MUT.out.log_final.collect{it[1]}.ifEmpty([]))
            }
        } else {
            ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_P2.out.log_final.collect{it[1]}.ifEmpty([]))
        }
        
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
    // Define summary params if needed
    def summary_params = WorkflowCdnaseq.paramsSummaryMap(workflow, params)
    
    if (params.email || params.email_on_fail) {
        // Comment out the NfcoreTemplate.email call as it may not be available
        // NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
        log.info "Pipeline completed. Email notification was requested but NfcoreTemplate may not be available."
    }
    
    // Use a try-catch block to handle potential errors with the summary method
    try {
        // NfcoreTemplate.summary(workflow, params, log)
        log.info "Pipeline execution summary not available (NfcoreTemplate may not be present)"
    } catch (Exception e) {
        log.warn "Could not generate execution summary: ${e.message}"
    }
    
    if (params.hook_url) {
        // Comment out the notification call as it may not be available
        // NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
        log.info "Pipeline completed. IM notification was requested but NfcoreTemplate may not be available."
    }
}
