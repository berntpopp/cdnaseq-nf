#!/usr/bin/env python

# This script takes the output from various tools' versions and arranges it
# into a YAML file for MultiQC.

import os
import yaml
import re

regexes = {
    "cdnaseq": ["v_pipeline.txt", r"(\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
    "FastQC": ["v_fastqc.txt", r"FastQC v(\S+)"],
    "MultiQC": ["v_multiqc.txt", r"multiqc, version (\S+)"],
    "BBMAP": ["v_bbmap.txt", r"BBMap version (\S+)"],
    "STAR": ["v_star.txt", r"(\S+)"],
    "GATK": ["v_gatk.txt", r"(\S+)"],
    "Picard": ["v_picard.txt", r"(\S+)"],
    "Samtools": ["v_samtools.txt", r"samtools (\S+)"],
    "Regtools": ["v_regtools.txt", r"(\S+)"],
    "Portcullis": ["v_portcullis.txt", r"(\S+)"],
    "featureCounts": ["v_featurecounts.txt", r"featureCounts v(\S+)"],
    "QualiMap": ["v_qualimap.txt", r"QualiMap v(\S+)"],
}

results = {"cdnaseq-nf": {"pipeline": "1.0.0", "tools": {}}}

# Search each file using its regex
for tool, value in regexes.items():
    path, regex = value
    if os.path.isfile(path):
        with open(path) as x:
            versions = x.read().strip()
            match = re.search(regex, versions)
            if match:
                results["cdnaseq-nf"]["tools"][tool] = match.group(1)

# Dump to YAML
print(
    yaml.dump(results, default_flow_style=False)
)
