#!/usr/bin/env python3

"""
prepare_references.py

Python script for downloading and preparing reference genomes and annotations
for the cdnaseq-nf pipeline.

This script automates the download of FASTA files, GTF annotations from GENCODE,
and optionally builds STAR indices for specified genome builds.
"""

import argparse
import os
import sys
import subprocess
import urllib.request
import gzip
import shutil
from pathlib import Path

def run_command(cmd, description=""):
    """Execute a shell command and handle errors."""
    print(f"Running: {description}")
    print(f"Command: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(f"Success: {description}")
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error: {description} failed")
        print(f"Exit code: {e.returncode}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        sys.exit(1)

def download_file(url, output_path, description=""):
    """Download a file from URL."""
    print(f"Downloading: {description}")
    print(f"URL: {url}")
    print(f"Output: {output_path}")
    
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"Downloaded successfully: {output_path}")
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        sys.exit(1)

def decompress_gz(gz_path, output_path):
    """Decompress a .gz file."""
    print(f"Decompressing: {gz_path} -> {output_path}")
    
    try:
        with gzip.open(gz_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Decompressed successfully: {output_path}")
    except Exception as e:
        print(f"Error decompressing {gz_path}: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Download and prepare reference genomes for cdnaseq-nf pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download hg38 references with GENCODE v45
  python prepare_references.py --genome_build hg38 --output_dir /path/to/refs --gencode_version_hg38 45

  # Download hg19 references with GENCODE v19
  python prepare_references.py --genome_build hg19 --output_dir /path/to/refs --gencode_version_hg19 19

  # Download without building STAR index
  python prepare_references.py --genome_build hg38 --output_dir /path/to/refs --skip_star_index
        """
    )

    parser.add_argument(
        "--genome_build",
        choices=["hg19", "hg38"],
        required=True,
        help="Genome build to download (hg19 or hg38)"
    )
    
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for reference files"
    )
    
    parser.add_argument(
        "--gencode_version_hg19",
        type=int,
        default=19,
        help="GENCODE version for hg19 (default: 19)"
    )
    
    parser.add_argument(
        "--gencode_version_hg38",
        type=int,
        default=45,
        help="GENCODE version for hg38 (default: 45)"
    )
    
    parser.add_argument(
        "--skip_star_index",
        action="store_true",
        help="Skip building STAR index"
    )
    
    parser.add_argument(
        "--star_threads",
        type=int,
        default=8,
        help="Number of threads for STAR index building (default: 8)"
    )
    
    parser.add_argument(
        "--star_sjdb_overhang",
        type=int,
        default=149,
        help="STAR sjdbOverhang parameter (default: 149, for 150bp reads)"
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set URLs and file names based on genome build
    if args.genome_build == "hg38":
        gencode_version = args.gencode_version_hg38
        fasta_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"
        gtf_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version}/gencode.v{gencode_version}.annotation.gtf.gz"
        genome_name = "GRCh38"
    else:  # hg19
        gencode_version = args.gencode_version_hg19
        fasta_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.primary_assembly.genome.fa.gz"
        gtf_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version}/gencode.v{gencode_version}.annotation.gtf.gz"
        genome_name = "GRCh37"

    # Define file paths
    fasta_gz = output_dir / f"{genome_name}.primary_assembly.genome.fa.gz"
    fasta_fa = output_dir / f"{genome_name}.primary_assembly.genome.fa"
    gtf_gz = output_dir / f"gencode.v{gencode_version}.annotation.gtf.gz"
    gtf_file = output_dir / f"gencode.v{gencode_version}.annotation.gtf"
    fasta_fai = output_dir / f"{genome_name}.primary_assembly.genome.fa.fai"
    fasta_dict = output_dir / f"{genome_name}.primary_assembly.genome.dict"

    print(f"Preparing {args.genome_build} references with GENCODE v{gencode_version}")
    print(f"Output directory: {output_dir}")

    # Download FASTA
    if not fasta_gz.exists():
        download_file(fasta_url, fasta_gz, f"{genome_name} FASTA")
    else:
        print(f"FASTA already exists: {fasta_gz}")

    # Download GTF
    if not gtf_gz.exists():
        download_file(gtf_url, gtf_gz, f"GENCODE v{gencode_version} GTF")
    else:
        print(f"GTF already exists: {gtf_gz}")

    # Decompress FASTA
    if not fasta_fa.exists():
        decompress_gz(fasta_gz, fasta_fa)
    else:
        print(f"Decompressed FASTA already exists: {fasta_fa}")

    # Decompress GTF
    if not gtf_file.exists():
        decompress_gz(gtf_gz, gtf_file)
    else:
        print(f"Decompressed GTF already exists: {gtf_file}")

    # Index FASTA with samtools
    if not fasta_fai.exists():
        run_command(f"samtools faidx {fasta_fa}", "Indexing FASTA with samtools")
    else:
        print(f"FASTA index already exists: {fasta_fai}")

    # Create sequence dictionary with GATK/Picard
    if not fasta_dict.exists():
        run_command(
            f"gatk CreateSequenceDictionary -R {fasta_fa} -O {fasta_dict}",
            "Creating sequence dictionary"
        )
    else:
        print(f"Sequence dictionary already exists: {fasta_dict}")

    # Build STAR index if requested
    if not args.skip_star_index:
        star_index_dir = output_dir / f"STAR_index_{genome_name}_gencode.v{gencode_version}"
        
        if not star_index_dir.exists() or not (star_index_dir / "SA").exists():
            star_index_dir.mkdir(exist_ok=True)
            
            star_cmd = f"""STAR --runMode genomeGenerate \\
                --genomeDir {star_index_dir} \\
                --genomeFastaFiles {fasta_fa} \\
                --sjdbGTFfile {gtf_file} \\
                --sjdbOverhang {args.star_sjdb_overhang} \\
                --runThreadN {args.star_threads}"""
            
            run_command(star_cmd, "Building STAR index")
        else:
            print(f"STAR index already exists: {star_index_dir}")
    else:
        print("Skipping STAR index building")

    print("\n" + "="*60)
    print("Reference preparation completed successfully!")
    print("="*60)
    print(f"Reference files location: {output_dir}")
    print(f"FASTA: {fasta_fa}")
    print(f"GTF: {gtf_file}")
    print(f"FASTA index: {fasta_fai}")
    print(f"Sequence dictionary: {fasta_dict}")
    
    if not args.skip_star_index:
        star_index_dir = output_dir / f"STAR_index_{genome_name}_gencode.v{gencode_version}"
        print(f"STAR index: {star_index_dir}")

if __name__ == "__main__":
    main()
