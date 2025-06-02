# Usage

## Introduction

This document describes how to run the cdnaseq-nf pipeline for RNA-seq analysis with splicing aberration detection.

## Samplesheet Input

The pipeline requires a CSV samplesheet with information about the samples to be processed.

### Full Samplesheet

The samplesheet must contain the following columns:

| Column         | Description                                               |
| -------------- | --------------------------------------------------------- |
| `sample_id`    | Custom sample name (unique identifier)                   |
| `patient_id`   | Patient identifier for grouping samples                  |
| `panel`        | Panel name for sample grouping                           |
| `fastq_r1`     | Full path to FastQ file for read 1                       |
| `fastq_r2`     | Full path to FastQ file for read 2 (empty for SE)        |
| `vcf`          | Full path to VCF file (optional, for mutated refs)       |
| `strandedness` | Sample strand-specificity (unstranded/forward/reverse)   |

### Example Samplesheet

```csv
sample_id,patient_id,panel,fastq_r1,fastq_r2,vcf,strandedness
SAMPLE1_T,PATIENT1,PANEL_A,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz,/data/patient1.vcf.gz,unstranded
SAMPLE1_N,PATIENT1,PANEL_A,/data/sample1_normal_R1.fastq.gz,/data/sample1_normal_R2.fastq.gz,,unstranded
SAMPLE2_T,PATIENT2,PANEL_A,/data/sample2_R1.fastq.gz,/data/sample2_R2.fastq.gz,/data/patient2.vcf.gz,unstranded
```

## Reference Preparation

Before running the pipeline, you need to prepare reference files:

### Using the Reference Preparation Script

```bash
# Download hg38 references and create indices
python bin/prepare_references.py \\
    --genome_build hg38 \\
    --output_dir /path/to/references \\
    --gencode_version_hg38 45 \\
    --star_threads 8

# For hg19
python bin/prepare_references.py \\
    --genome_build hg19 \\
    --output_dir /path/to/references \\
    --gencode_version_hg19 19 \\
    --star_threads 8

# Skip STAR index building if desired
python bin/prepare_references.py \\
    --genome_build hg38 \\
    --output_dir /path/to/references \\
    --skip_star_index
```

### Manual Reference Preparation

If you prefer to prepare references manually:

```bash
# Create reference directory
mkdir -p /path/to/references
cd /path/to/references

# Download FASTA and GTF files (example for GRCh38)
wget ftp://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz

# Decompress
gunzip *.gz

# Index FASTA file
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Create sequence dictionary
gatk CreateSequenceDictionary \\
    -R Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
    -O Homo_sapiens.GRCh38.dna.primary_assembly.dict
```

## Running the Pipeline

### Basic Command

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

### Profile Selection

#### Conda Environment
```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

#### SLURM Cluster
```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile slurm \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

#### Docker
```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile docker \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

#### Singularity
```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile singularity \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

> **Note**: Docker and Singularity profiles use pre-built container images with default fallback containers (`nfcore/base:2.1`). Individual processes have specific containers defined in their modules. You can customize container settings in `conf/docker.config` or `conf/singularity.config`.

### Minimal First-Pass Alignment
For faster analysis focusing only on QC, trimming, and first-pass alignments:

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda,minimal_alignment \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

This profile:

- Skips the second-pass STAR alignment (`skip_star_second_pass = true`)
- Performs FastQC and trimming
- Generates first-pass STAR alignments (P1-Ref and conditionally P1-Mut if VCFs provided)
- Skips most downstream analyses (variant calling, detailed splicing analysis, quantification)
- Includes a final MultiQC report
- Uses reduced resource allocation suitable for smaller systems

> **Note**: Docker and Singularity profiles use pre-built container images with default fallback containers (`nfcore/base:2.1`). Individual processes have specific containers defined in their modules. You can customize container settings in `conf/docker.config` or `conf/singularity.config`.

### Advanced Usage

#### Skip Specific Steps

```bash
# Skip adapter trimming
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --skip_trimming

# Skip second-pass STAR alignment (use first-pass only)
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --skip_star_second_pass

# Skip variant calling
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --skip_variant_calling

# Skip splicing analysis
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --skip_splicing_analysis
```

#### Custom Resource Limits

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --max_memory 256.GB \\
    --max_cpus 32 \\
    --max_time 480.h
```

#### Variant Calling with Target Regions

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --targets_bed /path/to/targets.bed
```

#### Output Publishing Options

By default, the pipeline copies output files to the results directory. You can change this behavior using the `publish_dir_mode` parameter:

```bash
# Use symbolic links instead of copying (faster, saves disk space)
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --publish_dir_mode symlink

# Use hard links (fast, saves space, but files must be on same filesystem)
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --publish_dir_mode link
```

**Available options:**
- `copy` (default): Copy files to output directory
- `symlink`: Create symbolic links (saves disk space)
- `link`: Create hard links (saves space, same filesystem only)  
- `move`: Move files to output directory
- `rellink`: Create relative symbolic links
- `copyNoFollow`: Copy files without following symbolic links

#### STAR Indexing RAM Control

For systems with limited memory, you can control RAM usage during in-workflow STAR indexing (mutated references and P2 index):

```bash
# Example with 30GB RAM limit and reduced index size
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --star_index_limit_genome_generate_ram 30000000000 \\
    --star_index_genome_sa_sparse_d 2 \\
    --star_index_genome_chr_bin_nbits 16
```

**RAM Control Parameters:**
- `--star_index_limit_genome_generate_ram`: Direct RAM limit in bytes (e.g., 30000000000 for 30GB)
- `--star_index_genome_sa_sparse_d`: Higher values (2-3) reduce index size and RAM requirements
- `--star_index_genome_chr_bin_nbits`: Lower values (15-16) reduce RAM for human genome

These parameters are particularly useful when the pipeline generates indices for mutated references, which are full genome indices and can be memory-intensive.

**Note:** These parameters control STAR indexing *within* the pipeline workflow. They are different from the RAM control options available in the standalone `bin/prepare_references.py` script, which is used for initial reference preparation.

**Alternatively**, you can use the provided low-memory configuration file:

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    -c conf/low_memory.config \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

### Resume Failed Runs

If a pipeline run fails, you can resume from where it left off:

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    -resume
```

## Configuration

### Custom Configuration

Create a custom configuration file (`custom.config`):

```groovy
// Custom resource configuration
process {
    withLabel: process_high {
        cpus = 32
        memory = 128.GB
        time = 24.h
    }
    
    withName: STAR_ALIGN_P2 {
        cpus = 16
        memory = 64.GB
        time = 12.h
    }
}

// Custom parameters
params {
    max_memory = '256.GB'
    max_cpus = 32
    max_time = '240.h'
}
```

Use the custom configuration:

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    -c custom.config \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

## Test Run

Run the pipeline with test data:

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile test,conda \\
    --output_dir ./test_results
```

## Troubleshooting

### Common Issues

1. **Out of Memory Errors**: Increase memory allocation for specific processes
2. **STAR Index Issues**: Ensure STAR index was built with sufficient memory
3. **VCF File Issues**: Ensure VCF files are properly indexed (`.tbi` files)
4. **Permission Issues**: Check file permissions and paths

### Debug Mode

Run with debug information:

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    -with-trace \\
    -with-report \\
    -with-timeline
```

### Resource Monitoring

Monitor resource usage:

```bash
nextflow run berntpopp/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    -with-trace trace.txt
```

## Parameter Reference

For a complete list of parameters, see the [README.md](../README.md#pipeline-parameters) or run:

```bash
nextflow run berntpopp/cdnaseq-nf --help
```
