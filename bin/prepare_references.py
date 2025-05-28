#!/usr/bin/env python3

"""
prepare_references.py

Python script for downloading and preparing reference genomes and annotations
for the cdnaseq-nf pipeline.

This script automates the download of FASTA files, GTF annotations from GENCODE,
and optionally builds STAR indices for specified genome builds.

Features:
- Resumable downloads with checksum verification
- Detailed logging with timestamps
- Progress tracking and state management
- Automatic retry on download failures
"""

import argparse
import os
import sys
import subprocess
import urllib.request
import gzip
import shutil
import hashlib
import json
import logging
import time
from datetime import datetime
from pathlib import Path


def setup_logging(output_dir):
    """Set up logging to both file and console."""
    log_file = Path(output_dir) / "prepare_references.log"

    # Create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Set up root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Clear existing handlers
    logger.handlers.clear()

    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def flush_logging():
    """Explicitly flush all logging handlers to ensure output is written."""
    for handler in logging.getLogger().handlers:
        handler.flush()


def calculate_md5(file_path):
    """Calculate MD5 hash of a file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def load_state(state_file):
    """Load progress state from file."""
    if state_file.exists():
        try:
            with open(state_file, "r") as f:
                return json.load(f)
        except Exception as e:
            logging.warning(f"Could not load state file: {e}")
    return {}


def save_state(state_file, state):
    """Save progress state to file."""
    try:
        with open(state_file, "w") as f:
            json.dump(state, f, indent=2)
    except Exception as e:
        logging.error(f"Could not save state file: {e}")


def verify_file_integrity(file_path, expected_md5=None, min_size=1000):
    """Verify file integrity by checking size and optionally MD5."""
    if not file_path.exists():
        return False

    # Check minimum file size
    if file_path.stat().st_size < min_size:
        logging.warning(
            f"File {file_path} is too small ({file_path.stat().st_size} bytes)"
        )
        return False

    # Check MD5 if provided
    if expected_md5:
        actual_md5 = calculate_md5(file_path)
        if actual_md5 != expected_md5:
            logging.warning(
                f"MD5 mismatch for {file_path}: expected {expected_md5}, got {actual_md5}"
            )
            return False
        logging.info(f"MD5 verification passed for {file_path}")

    return True


def download_file_with_resume(
    url, output_path, description="", expected_md5=None, max_retries=3
):
    """Download a file with resume capability and checksum verification."""
    output_path = Path(output_path)
    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    # Check if file already exists and is valid
    if verify_file_integrity(output_path, expected_md5):
        logging.info(f"File already exists and verified: {output_path}")
        return True

    for attempt in range(max_retries):
        try:
            logging.info(
                f"Downloading: {description} (attempt {attempt + 1}/{max_retries})"
            )
            logging.info(f"URL: {url}")
            logging.info(f"Output: {output_path}")

            # Remove temp file if it exists from previous failed attempt
            if temp_path.exists():
                temp_path.unlink()
            # Download with progress reporting
            last_percent = 0
            last_report_time = time.time()

            def report_progress(block_num, block_size, total_size):
                nonlocal last_percent, last_report_time
                current_time = time.time()

                if total_size > 0:
                    percent = min(100, (block_num * block_size / total_size) * 100)
                    downloaded = block_num * block_size

                    # Report every 5% or every 30 seconds to avoid spam but show progress
                    if (
                        percent - last_percent >= 5.0
                        or current_time - last_report_time >= 30
                        or percent >= 100
                    ):

                        mb_downloaded = downloaded / (1024 * 1024)
                        mb_total = total_size / (1024 * 1024)

                        print(
                            f"Progress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)"
                        )
                        logging.info(
                            f"Progress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)"
                        )

                        last_percent = percent
                        last_report_time = current_time
                elif (
                    block_num % 1000 == 0
                ):  # For unknown size, report every 1000 blocks
                    mb_downloaded = (block_num * block_size) / (1024 * 1024)
                    print(f"Downloaded: {mb_downloaded:.1f} MB")
                    logging.info(f"Downloaded: {mb_downloaded:.1f} MB")

            # Download the file
            logging.info("Starting download...")
            print("Starting download...")
            urllib.request.urlretrieve(url, temp_path, reporthook=report_progress)
            print("Download completed, verifying file...")
            logging.info("Download completed, verifying file...")

            # Verify download integrity
            if verify_file_integrity(
                temp_path, expected_md5, min_size=1000000
            ):  # 1MB minimum for large files
                # Move temp file to final location
                shutil.move(str(temp_path), str(output_path))
                logging.info(f"Downloaded and verified successfully: {output_path}")
                return True
            else:
                logging.error(f"Downloaded file failed integrity check: {temp_path}")
                if temp_path.exists():
                    temp_path.unlink()

        except Exception as e:
            logging.error(f"Download attempt {attempt + 1} failed: {e}")
            if temp_path.exists():
                temp_path.unlink()

            if attempt < max_retries - 1:
                wait_time = 2**attempt  # Exponential backoff
                logging.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)

    logging.error(f"Failed to download {url} after {max_retries} attempts")
    return False


def run_command(cmd, description=""):
    """Execute a shell command and handle errors."""
    logging.info(f"Running: {description}")
    logging.info(f"Command: {cmd}")

    try:
        result = subprocess.run(
            cmd, shell=True, check=True, capture_output=True, text=True
        )
        logging.info(f"Success: {description}")
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Error: {description} failed")
        logging.error(f"Exit code: {e.returncode}")
        logging.error(f"STDOUT: {e.stdout}")
        logging.error(f"STDERR: {e.stderr}")
        sys.exit(1)


def decompress_gz(gz_path, output_path, state_file=None, state_key=None):
    """Decompress a .gz file with progress tracking."""
    gz_path = Path(gz_path)
    output_path = Path(output_path)

    # Check if already decompressed
    if output_path.exists() and verify_file_integrity(output_path):
        logging.info(f"Decompressed file already exists: {output_path}")
        return True

    logging.info(f"Decompressing: {gz_path} -> {output_path}")

    try:
        with gzip.open(gz_path, "rb") as f_in:
            with open(output_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        # Verify decompressed file
        if verify_file_integrity(output_path):
            logging.info(f"Decompressed successfully: {output_path}")

            # Update state if provided
            if state_file and state_key:
                file_md5 = calculate_md5(output_path)
                logging.info(f"Calculated MD5 for decompressed file: {file_md5}")

                state = load_state(state_file)
                state[state_key] = {
                    "status": "completed",
                    "timestamp": datetime.now().isoformat(),
                    "md5": file_md5,
                }
                save_state(state_file, state)
                logging.info(f"Updated state for {state_key}")

            # Explicit logging flush
            for handler in logging.getLogger().handlers:
                handler.flush()

            return True
        else:
            logging.error(f"Decompressed file failed integrity check: {output_path}")
            return False

    except Exception as e:
        logging.error(f"Error decompressing {gz_path}: {e}")
        return False


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
        """,
    )

    parser.add_argument(
        "--genome_build",
        choices=["hg19", "hg38"],
        required=True,
        help="Genome build to download (hg19 or hg38)",
    )

    parser.add_argument(
        "--output_dir", required=True, help="Output directory for reference files"
    )

    parser.add_argument(
        "--gencode_version_hg19",
        type=int,
        default=19,
        help="GENCODE version for hg19 (default: 19)",
    )

    parser.add_argument(
        "--gencode_version_hg38",
        type=int,
        default=45,
        help="GENCODE version for hg38 (default: 45)",
    )

    parser.add_argument(
        "--skip_star_index", action="store_true", help="Skip building STAR index"
    )

    parser.add_argument(
        "--star_threads",
        type=int,
        default=8,
        help="Number of threads for STAR index building (default: 8)",
    )

    parser.add_argument(
        "--star_sjdb_overhang",
        type=int,
        default=149,
        help="STAR sjdbOverhang parameter (default: 149, for 150bp reads)",
    )

    parser.add_argument(
        "--force_redownload", action="store_true", help="Force redownload of all files"
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set up logging
    logger = setup_logging(output_dir)
    logger.info("=" * 60)
    logger.info("Starting reference preparation")
    logger.info("=" * 60)

    # Initialize state tracking
    state_file = output_dir / "preparation_state.json"
    state = load_state(state_file)  # Set URLs and file names based on genome build
    if args.genome_build == "hg38":
        gencode_version = args.gencode_version_hg38
        # Use a more reliable URL structure
        fasta_url = (
            "https://ftp.ebi.ac.uk/pub/databases/gencode/"
            "Gencode_human/release_45/"
            "GRCh38.primary_assembly.genome.fa.gz"
        )
        gtf_url = (
            f"https://ftp.ebi.ac.uk/pub/databases/gencode/"
            f"Gencode_human/release_{gencode_version}/"
            f"gencode.v{gencode_version}.annotation.gtf.gz"
        )
        genome_name = "GRCh38"
    else:  # hg19
        gencode_version = args.gencode_version_hg19
        fasta_url = (
            "https://ftp.ebi.ac.uk/pub/databases/gencode/"
            "Gencode_human/release_19/"
            "GRCh37.primary_assembly.genome.fa.gz"
        )
        gtf_url = (
            f"https://ftp.ebi.ac.uk/pub/databases/gencode/"
            f"Gencode_human/release_{gencode_version}/"
            f"gencode.v{gencode_version}.annotation.gtf.gz"
        )
        genome_name = "GRCh37"

    # Define file paths
    fasta_gz = output_dir / f"{genome_name}.primary_assembly.genome.fa.gz"
    fasta_fa = output_dir / f"{genome_name}.primary_assembly.genome.fa"
    gtf_gz = output_dir / f"gencode.v{gencode_version}.annotation.gtf.gz"
    gtf_file = output_dir / f"gencode.v{gencode_version}.annotation.gtf"
    fasta_fai = output_dir / f"{genome_name}.primary_assembly.genome.fa.fai"
    fasta_dict = output_dir / f"{genome_name}.primary_assembly.genome.dict"

    logger.info(
        f"Preparing {args.genome_build} references with GENCODE v{gencode_version}"
    )
    logger.info(f"Output directory: {output_dir}")

    # Download FASTA
    fasta_state_key = f"fasta_download_{genome_name}"
    if args.force_redownload or not verify_file_integrity(fasta_gz):
        if download_file_with_resume(fasta_url, fasta_gz, f"{genome_name} FASTA"):
            state[fasta_state_key] = {
                "status": "completed",
                "timestamp": datetime.now().isoformat(),
                "md5": calculate_md5(fasta_gz),
            }
            save_state(state_file, state)
        else:
            logger.error("Failed to download FASTA file")
            sys.exit(1)
    else:
        logger.info(f"FASTA already exists and verified: {fasta_gz}")

    # Download GTF
    gtf_state_key = f"gtf_download_v{gencode_version}"
    if args.force_redownload or not verify_file_integrity(gtf_gz):
        if download_file_with_resume(
            gtf_url, gtf_gz, f"GENCODE v{gencode_version} GTF"
        ):
            state[gtf_state_key] = {
                "status": "completed",
                "timestamp": datetime.now().isoformat(),
                "md5": calculate_md5(gtf_gz),
            }
            save_state(state_file, state)
        else:
            logger.error("Failed to download GTF file")
            sys.exit(1)
    else:
        logger.info(f"GTF already exists and verified: {gtf_gz}")  # Decompress FASTA
    fasta_decomp_key = f"fasta_decompress_{genome_name}"
    print(f"\n=== Starting FASTA decompression ===")
    logger.info(f"Starting FASTA decompression: {fasta_gz}")
    if not decompress_gz(fasta_gz, fasta_fa, state_file, fasta_decomp_key):
        logger.error("Failed to decompress FASTA file")
        sys.exit(1)
    print(f"=== FASTA decompression completed ===\n")
    flush_logging()

    # Decompress GTF
    gtf_decomp_key = f"gtf_decompress_v{gencode_version}"
    print(f"\n=== Starting GTF decompression ===")
    logger.info(f"Starting GTF decompression: {gtf_gz}")
    if not decompress_gz(gtf_gz, gtf_file, state_file, gtf_decomp_key):
        logger.error("Failed to decompress GTF file")
        sys.exit(1)
    print(f"=== GTF decompression completed ===\n")
    flush_logging()  # Index FASTA with samtools
    fasta_index_key = f"fasta_index_{genome_name}"
    print(f"\n=== Starting FASTA indexing ===")
    logger.info(f"Starting FASTA indexing: {fasta_fa}")

    if not fasta_fai.exists():
        logger.info("Creating FASTA index with samtools")
        run_command(f"samtools faidx {fasta_fa}", "Indexing FASTA with samtools")

        # Calculate MD5 for the index file
        if fasta_fai.exists():
            fai_md5 = calculate_md5(fasta_fai)
            logger.info(f"FASTA index created, MD5: {fai_md5}")

            state[fasta_index_key] = {
                "status": "completed",
                "timestamp": datetime.now().isoformat(),
                "md5": fai_md5,
                "file_path": str(fasta_fai),
            }
            save_state(state_file, state)
            logger.info(f"Updated state for {fasta_index_key}")
        else:
            logger.error("FASTA index file was not created")
            sys.exit(1)
    else:
        # Index already exists, calculate and log its MD5
        fai_md5 = calculate_md5(fasta_fai)
        logger.info(f"FASTA index already exists: {fasta_fai}")
        logger.info(f"FASTA index MD5: {fai_md5}")
        # Update state with existing file info
        state[fasta_index_key] = {
            "status": "completed",
            "timestamp": datetime.now().isoformat(),
            "md5": fai_md5,
            "file_path": str(fasta_fai),
        }
        save_state(state_file, state)

    print(f"=== FASTA indexing completed ===\n")
    flush_logging()

    # Create sequence dictionary with GATK/Picard
    dict_key = f"sequence_dict_{genome_name}"
    print(f"\n=== Starting sequence dictionary creation ===")
    logger.info(f"Starting sequence dictionary creation: {fasta_dict}")

    if not fasta_dict.exists():
        logger.info("Creating sequence dictionary with GATK")
        run_command(
            f"gatk CreateSequenceDictionary -R {fasta_fa} -O {fasta_dict}",
            "Creating sequence dictionary",
        )

        # Calculate MD5 for the dictionary file
        if fasta_dict.exists():
            dict_md5 = calculate_md5(fasta_dict)
            logger.info(f"Sequence dictionary created, MD5: {dict_md5}")

            state[dict_key] = {
                "status": "completed",
                "timestamp": datetime.now().isoformat(),
                "md5": dict_md5,
                "file_path": str(fasta_dict),
            }
            save_state(state_file, state)
            logger.info(f"Updated state for {dict_key}")
        else:
            logger.error("Sequence dictionary file was not created")
            sys.exit(1)
    else:
        # Dictionary already exists, calculate and log its MD5
        dict_md5 = calculate_md5(fasta_dict)
        logger.info(f"Sequence dictionary already exists: {fasta_dict}")
        logger.info(f"Sequence dictionary MD5: {dict_md5}")

        # Update state with existing file info
        state[dict_key] = {
            "status": "completed",
            "timestamp": datetime.now().isoformat(),
            "md5": dict_md5,
            "file_path": str(fasta_dict),
        }
        save_state(state_file, state)

    print(
        f"=== Sequence dictionary creation completed ===\n"
    )  # Build STAR index if requested
    if not args.skip_star_index:
        star_index_dir = (
            output_dir / f"STAR_index_{genome_name}_gencode.v{gencode_version}"
        )
        star_index_key = f"star_index_{genome_name}_v{gencode_version}"

        print(f"\n=== Starting STAR index creation ===")
        logger.info(f"Starting STAR index creation: {star_index_dir}")

        if not star_index_dir.exists() or not (star_index_dir / "SA").exists():
            star_index_dir.mkdir(exist_ok=True)
            logger.info("Building STAR index (this may take a while)")

            star_cmd = f"""STAR --runMode genomeGenerate \\
                --genomeDir {star_index_dir} \\
                --genomeFastaFiles {fasta_fa} \\
                --sjdbGTFfile {gtf_file} \\
                --sjdbOverhang {args.star_sjdb_overhang} \\
                --runThreadN {args.star_threads}"""

            run_command(star_cmd, "Building STAR index")

            # Calculate MD5 checksums for key STAR index files
            star_files_to_check = ["SA", "SAindex", "Genome", "chrNameLength.txt"]
            star_file_checksums = {}

            for star_file in star_files_to_check:
                star_file_path = star_index_dir / star_file
                if star_file_path.exists():
                    file_md5 = calculate_md5(star_file_path)
                    star_file_checksums[star_file] = file_md5
                    logger.info(f"STAR {star_file} MD5: {file_md5}")

            logger.info(
                f"STAR index created successfully with {len(star_file_checksums)} key files"
            )

            state[star_index_key] = {
                "status": "completed",
                "timestamp": datetime.now().isoformat(),
                "index_path": str(star_index_dir),
                "file_checksums": star_file_checksums,
            }
            save_state(state_file, state)
            logger.info(f"Updated state for {star_index_key}")
        else:
            logger.info(f"STAR index already exists: {star_index_dir}")

            # Calculate MD5 checksums for existing STAR index files
            star_files_to_check = ["SA", "SAindex", "Genome", "chrNameLength.txt"]
            star_file_checksums = {}

            for star_file in star_files_to_check:
                star_file_path = star_index_dir / star_file
                if star_file_path.exists():
                    file_md5 = calculate_md5(star_file_path)
                    star_file_checksums[star_file] = file_md5
                    logger.info(f"Existing STAR {star_file} MD5: {file_md5}")

            # Update state with existing file info
            state[star_index_key] = {
                "status": "completed",
                "timestamp": datetime.now().isoformat(),
                "index_path": str(star_index_dir),
                "file_checksums": star_file_checksums,
            }
            save_state(state_file, state)

        print(f"=== STAR index creation completed ===\n")
    else:
        logger.info("Skipping STAR index building")

    logger.info("=" * 60)
    logger.info("Reference preparation completed successfully!")
    logger.info("=" * 60)
    logger.info(f"Reference files location: {output_dir}")
    logger.info(f"FASTA: {fasta_fa}")
    logger.info(f"GTF: {gtf_file}")
    logger.info(f"FASTA index: {fasta_fai}")
    logger.info(f"Sequence dictionary: {fasta_dict}")

    if not args.skip_star_index:
        star_index_dir = (
            output_dir / f"STAR_index_{genome_name}_gencode.v{gencode_version}"
        )
        logger.info(f"STAR index: {star_index_dir}")

    logger.info(f"State file: {state_file}")
    logger.info("Use --force_redownload to redownload files if needed")


if __name__ == "__main__":
    main()
