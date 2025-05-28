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
nextflow run laborberlin/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

### Profile Selection

#### Conda Environment
```bash
nextflow run laborberlin/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

#### SLURM Cluster
```bash
nextflow run laborberlin/cdnaseq-nf \\
    -profile slurm \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

#### Docker
```bash
nextflow run laborberlin/cdnaseq-nf \\
    -profile docker \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

### Advanced Usage

#### Skip Specific Steps

```bash
# Skip adapter trimming
nextflow run laborberlin/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --skip_trimming

# Skip variant calling
nextflow run laborberlin/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --skip_variant_calling

# Skip splicing analysis
nextflow run laborberlin/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --skip_splicing_analysis
```

#### Custom Resource Limits

```bash
nextflow run laborberlin/cdnaseq-nf \\
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
nextflow run laborberlin/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    --targets_bed /path/to/targets.bed
```

### Resume Failed Runs

If a pipeline run fails, you can resume from where it left off:

```bash
nextflow run laborberlin/cdnaseq-nf \\
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
nextflow run laborberlin/cdnaseq-nf \\
    -profile conda \\
    -c custom.config \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results
```

## Test Run

Run the pipeline with test data:

```bash
nextflow run laborberlin/cdnaseq-nf \\
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
nextflow run laborberlin/cdnaseq-nf \\
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
nextflow run laborberlin/cdnaseq-nf \\
    -profile conda \\
    --input_samplesheet samplesheet.csv \\
    --ref_dir /path/to/references \\
    --output_dir ./results \\
    -with-trace trace.txt
```

## Parameter Reference

For a complete list of parameters, see the [README.md](../README.md#pipeline-parameters) or run:

```bash
nextflow run laborberlin/cdnaseq-nf --help
```
