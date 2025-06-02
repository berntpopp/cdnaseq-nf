# cdnaseq-nf

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://apptainer.org/)

## Introduction

**cdnaseq-nf** is a comprehensive Nextflow DSL2 pipeline for RNA sequencing analysis with a focus on splicing aberrations using a "double alignment with mutated reference" strategy. The pipeline is designed for detection of splice-affecting variants and quantification of aberrant splicing events in RNA-seq data.

The pipeline performs:

- Quality control and adapter trimming
- Two-pass STAR alignment with splice junction discovery
- Optional alignment to patient-specific mutated references (VCF-based)
- BAM processing and duplicate marking
- Variant calling using GATK HaplotypeCaller (RNA-seq mode)
- Comprehensive splicing analysis using Regtools and Portcullis
- Gene expression quantification
- Multi-sample quality control reporting

## Pipeline Summary

1. **Quality Control** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. **Adapter Trimming** ([`BBDuk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/))
3. **Reference Preparation** (Optional mutated reference generation from VCF)
4. **STAR Indexing** ([`STAR`](https://github.com/alexdobin/STAR))
5. **First-Pass Alignment** (Reference and optionally mutated reference)
6. **Splice Junction Aggregation** (Custom filtering and merging)
7. **Second-Pass Alignment** (Junction-guided alignment)
8. **BAM Processing** ([`Picard`](https://broadinstitute.github.io/picard/), [`Samtools`](http://www.htslib.org/))
9. **Variant Calling** ([`GATK HaplotypeCaller`](https://gatk.broadinstitute.org/))
10. **Splicing Analysis** ([`Regtools`](https://regtools.readthedocs.io/), [`Portcullis`](https://portcullis.readthedocs.io/))
11. **Quantification** ([`featureCounts`](http://subread.sourceforge.net/))
12. **Quality Control Report** ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Prepare your samplesheet (see [usage docs](#samplesheet) for details):

   ```csv
   sample_id,patient_id,panel,fastq_r1,fastq_r2,vcf,strandedness
   SAMPLE1,PATIENT1,PANEL_A,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,/path/to/patient1.vcf.gz,unstranded
   SAMPLE2,PATIENT2,PANEL_A,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,,unstranded
   ```

4. Prepare reference files:

   ```bash
   # Download and prepare reference files for hg38
   python bin/prepare_references.py \\
       --genome_build hg38 \\
       --output_dir /path/to/references \\
       --gencode_version_hg38 45 \\
       --star_threads 8

   # Or for hg19
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

5. Run the pipeline:

   ```bash
   nextflow run laborberlin/cdnaseq-nf \\
       -profile conda \\
       --input_samplesheet samplesheet.csv \\
       --ref_dir /path/to/references \\
       --output_dir ./results
   ```

## Documentation

The cdnaseq-nf pipeline comes with documentation about the pipeline: [usage](docs/usage.md) and [output](docs/output.md).

## Samplesheet

The input samplesheet must be a CSV file with the following required columns:

| Column         | Description                                               |
| -------------- | --------------------------------------------------------- |
| `sample_id`    | Custom sample name                                        |
| `patient_id`   | Patient identifier for grouping samples                  |
| `panel`        | Panel name for sample grouping                           |
| `fastq_r1`     | Full path to FastQ file for read 1                       |
| `fastq_r2`     | Full path to FastQ file for read 2 (leave empty for SE)  |
| `vcf`          | Full path to VCF file (optional, for mutated references) |
| `strandedness` | Sample strand-specificity                                |

An example samplesheet is provided in `assets/samplesheet_example.csv`.

## Pipeline Parameters

### Input/Output Options

| Parameter            | Description                           | Type     | Default                   | Required |
| -------------------- | ------------------------------------- | -------- | ------------------------- | -------- |
| `input_samplesheet`  | Path to input samplesheet            | `string` |                           | ✓        |
| `ref_dir`            | Directory containing reference files  | `string` |                           | ✓        |
| `output_dir`         | Path for output files                 | `string` | `./results`               |          |

### Reference Options

| Parameter                      | Description                                    | Type      | Default   | Required |
| ------------------------------ | ---------------------------------------------- | --------- | --------- | -------- |
| `genome_build`                 | Genome build (GRCh37, GRCh38, GRCm38, etc.)   | `string`  | `GRCh38`  |          |
| `perform_mut_ref_alignment`    | Perform alignment to mutated references        | `boolean` | `true`    |          |

### STAR Indexing RAM Control (In-Workflow)

For STAR indexing steps performed *within* the pipeline (e.g., for mutated references or the P2 index), you can control RAM usage with the following parameters:

| Parameter                               | Description                                           | Type      | Default | Required |
| --------------------------------------- | ----------------------------------------------------- | --------- | ------- | -------- |
| `star_index_limit_genome_generate_ram`  | Limits RAM used by STAR for index generation (bytes) | `integer` | `null`  |          |
| `star_index_genome_sa_sparse_d`         | STAR genomeSAsparseD parameter (higher = less RAM)   | `integer` | `null`  |          |
| `star_index_genome_chr_bin_nbits`       | STAR genomeChrBinNbits parameter (lower = less RAM)  | `integer` | `null`  |          |

**Examples:**
- `--star_index_limit_genome_generate_ram 30000000000` (30GB RAM limit)
- `--star_index_genome_sa_sparse_d 2` (Reduce index size and RAM requirements)
- `--star_index_genome_chr_bin_nbits 16` (For human genome with limited RAM)

These parameters are particularly useful when generating indices for full mutated genomes on systems with limited memory.

**Pre-configured option:** Use the provided low-memory configuration file with `-c conf/low_memory.config` for systems with limited RAM.

**Pre-configured option:** Use the provided low-memory configuration file with `-c conf/low_memory.config` for systems with limited RAM.

**Pre-configured option:** Use the provided low-memory configuration file with `--c conf/low_memory.config` for systems with limited RAM.

### Analysis Options

| Parameter                     | Description                              | Type      | Default | Required |
| ----------------------------- | ---------------------------------------- | --------- | ------- | -------- |
| `skip_fastqc`                 | Skip FastQC quality control              | `boolean` | `false` |          |
| `skip_trimming`               | Skip adapter trimming                    | `boolean` | `false` |          |
| `skip_variant_calling`        | Skip GATK variant calling                | `boolean` | `false` |          |
| `skip_splicing_analysis`      | Skip splicing analysis                   | `boolean` | `false` |          |
| `skip_quantification`         | Skip gene expression quantification     | `boolean` | `false` |          |
| `skip_multiqc`                | Skip MultiQC report generation           | `boolean` | `false` |          |

### Additional Options

| Parameter        | Description                           | Type     | Default                            | Required |
| ---------------- | ------------------------------------- | -------- | ---------------------------------- | -------- |
| `targets_bed`    | BED file with target regions for variant calling | `string` |                                    |          |
| `adapter_fasta`  | FASTA file with adapter sequences     | `string` | `$projectDir/assets/adapters.fa`   |          |

### Resource Limits

| Parameter    | Description              | Type     | Default  | Required |
| ------------ | ------------------------ | -------- | -------- | -------- |
| `max_memory` | Maximum memory per job   | `string` | `128.GB` |          |
| `max_cpus`   | Maximum CPUs per job     | `int`    | `16`     |          |
| `max_time`   | Maximum time per job     | `string` | `240.h`  |          |

## Profiles

The pipeline supports multiple execution profiles:

- `conda`: Use Conda for software management
- `mamba`: Use Mamba (faster Conda) for software management
- `docker`: Use Docker containers with default fallback images
- `singularity`: Use Singularity containers with default fallback images
- `slurm`: Configure for SLURM scheduler
- `test`: Run with test data

### Container Profiles

The `docker` and `singularity` profiles are now fully functional with:

- **Default container images**: Uses `nfcore/base:2.1` as fallback for processes without specific containers
- **Per-process containers**: Individual modules specify appropriate biocontainer images
- **Customizable configurations**: Modify `conf/docker.config` or `conf/singularity.config` for advanced settings

Example usage:
```bash
# Docker
nextflow run laborberlin/cdnaseq-nf -profile docker --input_samplesheet samplesheet.csv --ref_dir /path/to/refs

# Singularity  
nextflow run laborberlin/cdnaseq-nf -profile singularity --input_samplesheet samplesheet.csv --ref_dir /path/to/refs
```

## Output

The pipeline generates various output files organized in the following structure:

```text
results/
├── fastqc/                     # FastQC reports
├── trimming/                   # Trimmed reads and logs
├── alignment/
│   ├── pass1_ref/              # First-pass reference alignment
│   ├── pass1_mut/              # First-pass mutated reference alignment
│   ├── splice_junctions/       # Aggregated splice junctions
│   └── pass2/                  # Second-pass alignment results
├── bam_processing/             # Processed BAM files
├── variant_calling/            # GATK variant calls
├── splicing_analysis/          # Regtools and Portcullis results
├── quantification/             # featureCounts results
├── multiqc/                    # MultiQC report
└── pipeline_info/              # Pipeline execution info
```

For detailed output descriptions, see [docs/output.md](docs/output.md).

## Citations

If you use cdnaseq-nf for your analysis, please cite:

- **Nextflow**: Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820.

- **STAR**: Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21.

- **GATK**: McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303.

- **MultiQC**: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8.

## Credits

This pipeline was developed by the Labor Berlin team for comprehensive RNA-seq splicing analysis.

## Support

For support, please:

1. Check the [documentation](docs/)
2. Search existing [issues](https://github.com/laborberlin/cdnaseq-nf/issues)
3. Create a new [issue](https://github.com/laborberlin/cdnaseq-nf/issues/new) with detailed information

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
