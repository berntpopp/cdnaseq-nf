# Output

## Introduction

This document describes the output produced by the cdnaseq-nf pipeline. Most of the plots are taken from the MultiQC report, which summarizes results at the end of the pipeline.

## Pipeline Overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read quality control
- [Trimming](#trimming) - Adapter and quality trimming
- [STAR Alignment](#star-alignment) - Two-pass RNA-seq alignment
- [BAM Processing](#bam-processing) - BAM file processing and indexing
- [Variant Calling](#variant-calling) - RNA-seq variant detection
- [Splicing Analysis](#splicing-analysis) - Splice junction analysis
- [Quantification](#quantification) - Gene expression quantification
- [MultiQC](#multiqc) - Aggregate report

## Output Directory Structure

```
results/
├── fastqc/
│   ├── raw/                    # FastQC results for raw reads
│   └── trimmed/                # FastQC results for trimmed reads
├── trimming/
│   ├── bbduk/                  # BBDuk trimming results
│   └── logs/                   # Trimming log files
├── alignment/
│   ├── pass1_ref/              # First-pass reference alignment
│   ├── pass1_mut/              # First-pass mutated reference alignment
│   ├── splice_junctions/       # Aggregated splice junctions
│   ├── pass2/                  # Second-pass alignment results
│   └── star_index/             # STAR index files
├── bam_processing/
│   ├── markdup/                # Duplicate-marked BAM files
│   ├── split_reads/            # Split N CIGAR reads
│   └── qc/                     # BAM quality control
├── variant_calling/
│   └── gatk/                   # GATK HaplotypeCaller results
├── splicing_analysis/
│   ├── regtools/               # Regtools junction analysis
│   └── portcullis/             # Portcullis splice site analysis
├── quantification/
│   └── featurecounts/          # Gene expression counts
├── multiqc/
│   └── multiqc_report.html     # Comprehensive QC report
└── pipeline_info/              # Pipeline execution information
```

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) generates quality control reports for sequencing data.

<details markdown="1">
<summary>Output files</summary>

- `fastqc/raw/`
  - `*_fastqc.html`: FastQC report containing quality metrics
  - `*_fastqc.zip`: Zip file containing the FastQC report and data
- `fastqc/trimmed/`
  - `*_fastqc.html`: FastQC report for trimmed reads
  - `*_fastqc.zip`: Zip file containing the FastQC report and data

</details>

FastQC provides information about:
- Basic statistics (e.g. number of reads, sequence length)
- Per base sequence quality
- Per tile sequence quality
- Per sequence quality scores
- Per base sequence content
- Per sequence GC content
- Per base N content
- Sequence length distribution
- Sequence duplication levels
- Overrepresented sequences
- Adapter content

## Trimming

[BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) is used for adapter trimming and quality filtering.

<details markdown="1">
<summary>Output files</summary>

- `trimming/bbduk/`
  - `*.trimmed.fastq.gz`: Trimmed FastQ files
  - `*.trimming_report.txt`: Trimming statistics
- `trimming/logs/`
  - `*.bbduk.log`: Detailed trimming logs

</details>

BBDuk removes:
- Adapter sequences
- Low quality bases
- Short reads after trimming

## STAR Alignment

[STAR](https://github.com/alexdobin/STAR) performs splice-aware alignment of RNA-seq reads using a two-pass approach.

<details markdown="1">
<summary>Output files</summary>

- `alignment/pass1_ref/`
  - `*.p1.Aligned.sortedByCoord.out.bam`: First-pass reference alignment
  - `*.p1.SJ.out.tab`: Splice junctions from first pass
  - `*.p1.Log.final.out`: Alignment summary statistics
- `alignment/pass1_mut/` (if VCF provided)
  - `*.p1mut.Aligned.sortedByCoord.out.bam`: First-pass mutated reference alignment
  - `*.p1mut.SJ.out.tab`: Splice junctions from mutated reference
- `alignment/splice_junctions/`
  - `filtered_junctions.tab`: Aggregated and filtered splice junctions
- `alignment/pass2/`
  - `*.p2.Aligned.sortedByCoord.out.bam`: Final alignment with junction guidance
  - `*.p2.SJ.out.tab`: Final splice junctions
  - `*.p2.Log.final.out`: Final alignment statistics

</details>

STAR alignment provides:
- Uniquely mapped reads percentage
- Multi-mapping reads percentage
- Unmapped reads percentage
- Splice junction discovery
- Alignment statistics

## BAM Processing

BAM files are processed using [Picard Tools](https://broadinstitute.github.io/picard/) and [Samtools](http://www.htslib.org/).

<details markdown="1">
<summary>Output files</summary>

- `bam_processing/markdup/`
  - `*.markdup.bam`: Duplicate-marked BAM files
  - `*.markdup.metrics.txt`: Duplicate marking statistics
- `bam_processing/split_reads/`
  - `*.split.bam`: BAM files with split N CIGAR reads
- `bam_processing/qc/`
  - `*.qualimap/`: Qualimap RNA-seq QC results

</details>

BAM processing includes:
- Duplicate marking and removal
- Split N CIGAR reads (for variant calling)
- Quality control metrics
- Indexing for downstream analysis

## Variant Calling

[GATK HaplotypeCaller](https://gatk.broadinstitute.org/) performs variant calling optimized for RNA-seq data.

<details markdown="1">
<summary>Output files</summary>

- `variant_calling/gatk/`
  - `*.vcf.gz`: Called variants in VCF format
  - `*.vcf.gz.tbi`: VCF index files
  - `*.variant_calling_metrics.txt`: Variant calling statistics

</details>

Variant calling features:
- RNA-seq specific filtering
- Support for target region analysis
- Quality score recalibration
- Allele-specific expression analysis

## Splicing Analysis

### Regtools

[Regtools](https://regtools.readthedocs.io/) extracts and annotates splice junctions.

<details markdown="1">
<summary>Output files</summary>

- `splicing_analysis/regtools/`
  - `*.junctions.bed`: Extracted splice junctions
  - `*.annotated.junctions.bed`: Annotated splice junctions
  - `*.junction_stats.txt`: Junction statistics

</details>

### Portcullis

[Portcullis](https://portcullis.readthedocs.io/) performs comprehensive splice site analysis.

<details markdown="1">
<summary>Output files</summary>

- `splicing_analysis/portcullis/`
  - `*.junctions.tab`: Portcullis junction calls
  - `*.junctions.bed`: BED format junctions
  - `*.portcullis.log`: Analysis log

</details>

Splicing analysis provides:
- Novel splice junction discovery
- Junction annotation with known databases
- Splice site strength prediction
- Alternative splicing detection

## Quantification

[featureCounts](http://subread.sourceforge.net/) quantifies gene expression levels.

<details markdown="1">
<summary>Output files</summary>

- `quantification/featurecounts/`
  - `*.featureCounts.txt`: Gene-level count matrix
  - `*.featureCounts.txt.summary`: Counting statistics
  - `*.featureCounts.txt.jcounts`: Junction counts (if requested)

</details>

Quantification features:
- Gene-level expression counts
- Multi-mapping read handling
- Strand-specific counting
- Feature-level statistics

## MultiQC

[MultiQC](http://multiqc.info) aggregates results from multiple tools into a single report.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: Comprehensive quality control report
  - `multiqc_data/`: Directory containing parsed data from tools

</details>

The MultiQC report includes:
- Sample quality metrics
- Alignment statistics
- Duplication rates
- Gene expression overview
- Splice junction summary
- Tool version information

## Pipeline Information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report.html`: Nextflow execution report
  - `execution_timeline.html`: Timeline of pipeline execution
  - `execution_trace.txt`: Resource usage trace
  - `pipeline_dag.svg`: Pipeline directed acyclic graph
  - `software_versions.yml`: Software versions used

</details>

## Quality Control Recommendations

### Read Quality
- **Good**: >90% of bases with Q30+
- **Acceptable**: >80% of bases with Q30+
- **Poor**: <80% of bases with Q30+

### Alignment Quality
- **Good**: >85% uniquely mapped reads
- **Acceptable**: >70% uniquely mapped reads
- **Poor**: <70% uniquely mapped reads

### Duplication Rate
- **Good**: <20% duplicates
- **Acceptable**: <30% duplicates
- **High**: >30% duplicates (may indicate library prep issues)

### Expression Detection
- **Good**: >15,000 genes with >1 count
- **Acceptable**: >10,000 genes with >1 count
- **Poor**: <10,000 genes with >1 count

## Troubleshooting

### Low Alignment Rates
- Check for adapter contamination
- Verify reference genome build
- Assess read quality

### High Duplication Rates
- May indicate low input material
- Check library construction protocol
- Consider PCR optimization

### Low Gene Detection
- Check strandedness settings
- Verify GTF annotation compatibility
- Assess sequencing depth

For additional support, please refer to the [usage documentation](usage.md) or create an issue on the project repository.
