#!/usr/bin/env python3

"""
Pipeline validation script for cdnaseq-nf
Checks that all required files are present and properly configured
"""

import os
import sys
import yaml
import json
from pathlib import Path


def check_file_exists(file_path, description=""):
    """Check if a file exists and report result"""
    if os.path.exists(file_path):
        print(f"✓ {description}: {file_path}")
        return True
    else:
        print(f"✗ {description}: {file_path} - NOT FOUND")
        return False


def check_directory_structure():
    """Check the expected directory structure"""
    print("Checking directory structure...")

    required_dirs = [
        "assets",
        "bin",
        "conf",
        "docs",
        "lib",
        "modules/local",
        "subworkflows/local",
        "workflows",
        ".github/workflows",
    ]

    all_exist = True
    for dir_path in required_dirs:
        if os.path.exists(dir_path):
            print(f"✓ Directory: {dir_path}")
        else:
            print(f"✗ Directory: {dir_path} - NOT FOUND")
            all_exist = False

    return all_exist


def check_core_files():
    """Check core pipeline files"""
    print("\nChecking core pipeline files...")

    core_files = [
        ("main.nf", "Main pipeline script"),
        ("nextflow.config", "Main configuration"),
        ("README.md", "Documentation"),
        ("LICENSE", "License file"),
    ]

    all_exist = True
    for file_path, description in core_files:
        if not check_file_exists(file_path, description):
            all_exist = False

    return all_exist


def check_config_files():
    """Check configuration files"""
    print("\nChecking configuration files...")

    config_files = [
        ("conf/base.config", "Base configuration"),
        ("conf/conda.config", "Conda configuration"),
        ("conf/modules.config", "Modules configuration"),
        ("conf/slurm.config", "SLURM configuration"),
        ("conf/test.config", "Test configuration"),
    ]

    all_exist = True
    for file_path, description in config_files:
        if not check_file_exists(file_path, description):
            all_exist = False

    return all_exist


def check_module_files():
    """Check module files"""
    print("\nChecking module files...")

    module_files = [
        ("modules/local/fastqc.nf", "FastQC module"),
        ("modules/local/bbduk_trim.nf", "BBDuk trimming module"),
        ("modules/local/gatk_tools.nf", "GATK tools module"),
        ("modules/local/star_index.nf", "STAR index module"),
        ("modules/local/star_align.nf", "STAR alignment module"),
        ("modules/local/bam_processing.nf", "BAM processing module"),
        ("modules/local/splicing_tools.nf", "Splicing tools module"),
        ("modules/local/quantification.nf", "Quantification module"),
        ("modules/local/multiqc.nf", "MultiQC module"),
    ]

    all_exist = True
    for file_path, description in module_files:
        if not check_file_exists(file_path, description):
            all_exist = False

    return all_exist


def check_workflow_files():
    """Check workflow files"""
    print("\nChecking workflow files...")

    workflow_files = [
        ("workflows/cdnaseq.nf", "Main workflow"),
        ("lib/WorkflowCdnaseq.groovy", "Workflow helper"),
    ]

    all_exist = True
    for file_path, description in workflow_files:
        if not check_file_exists(file_path, description):
            all_exist = False

    return all_exist


def check_asset_files():
    """Check asset files"""
    print("\nChecking asset files...")

    asset_files = [
        ("assets/samplesheet_example.csv", "Example samplesheet"),
        ("assets/samplesheet_test.csv", "Test samplesheet"),
        ("assets/adapters.fa", "Adapter sequences"),
        ("assets/multiqc_config.yml", "MultiQC configuration"),
    ]

    all_exist = True
    for file_path, description in asset_files:
        if not check_file_exists(file_path, description):
            all_exist = False

    return all_exist


def check_script_files():
    """Check script files"""
    print("\nChecking script files...")

    script_files = [
        ("bin/prepare_references.py", "Reference preparation script"),
        ("bin/filter_splice_junctions.py", "Junction filtering script"),
    ]

    all_exist = True
    for file_path, description in script_files:
        if not check_file_exists(file_path, description):
            all_exist = False

    return all_exist


def check_documentation():
    """Check documentation files"""
    print("\nChecking documentation...")

    doc_files = [
        ("docs/usage.md", "Usage documentation"),
        ("docs/output.md", "Output documentation"),
    ]

    all_exist = True
    for file_path, description in doc_files:
        if not check_file_exists(file_path, description):
            all_exist = False

    return all_exist


def check_github_workflows():
    """Check GitHub workflow files"""
    print("\nChecking GitHub workflows...")

    workflow_files = [
        (".github/workflows/ci.yml", "CI workflow"),
    ]

    all_exist = True
    for file_path, description in workflow_files:
        if not check_file_exists(file_path, description):
            all_exist = False

    return all_exist


def validate_nextflow_config():
    """Basic validation of nextflow.config syntax"""
    print("\nValidating nextflow.config syntax...")

    try:
        with open("nextflow.config", "r") as f:
            content = f.read()

        # Basic checks for required sections
        required_sections = ["params", "process", "profiles"]
        missing_sections = []

        for section in required_sections:
            if section not in content:
                missing_sections.append(section)

        if missing_sections:
            print(f"✗ Missing sections in nextflow.config: {missing_sections}")
            return False
        else:
            print("✓ nextflow.config appears valid")
            return True

    except Exception as e:
        print(f"✗ Error reading nextflow.config: {e}")
        return False


def validate_multiqc_config():
    """Validate MultiQC configuration"""
    print("\nValidating MultiQC configuration...")

    try:
        with open("assets/multiqc_config.yml", "r") as f:
            config = yaml.safe_load(f)

        required_keys = ["title", "module_order"]
        missing_keys = []

        for key in required_keys:
            if key not in config:
                missing_keys.append(key)

        if missing_keys:
            print(f"✗ Missing keys in MultiQC config: {missing_keys}")
            return False
        else:
            print("✓ MultiQC configuration appears valid")
            return True

    except Exception as e:
        print(f"✗ Error validating MultiQC config: {e}")
        return False


def main():
    """Main validation function"""
    print("cdnaseq-nf Pipeline Validation")
    print("=" * 40)

    # Change to script directory
    script_dir = Path(__file__).parent.absolute()
    os.chdir(script_dir)

    checks = [
        check_directory_structure(),
        check_core_files(),
        check_config_files(),
        check_module_files(),
        check_workflow_files(),
        check_asset_files(),
        check_script_files(),
        check_documentation(),
        check_github_workflows(),
        validate_nextflow_config(),
        validate_multiqc_config(),
    ]

    print("\n" + "=" * 40)

    if all(checks):
        print("✓ All validation checks passed!")
        print("\nPipeline appears to be properly configured.")
        print("\nNext steps:")
        print("1. Test the pipeline with: nextflow run . -profile test,conda")
        print("2. Prepare reference files using: python bin/prepare_references.py")
        print("3. Run with your data using the instructions in docs/usage.md")
        return 0
    else:
        print("✗ Some validation checks failed!")
        print("\nPlease fix the missing files/issues before running the pipeline.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
