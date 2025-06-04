process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc=1.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.15--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0' }"

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    
    """
    multiqc \\
        --force \\
        $args \\
        $config \\
        $extra_config \\
        $logo \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "bioconda::multiqc=1.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.15--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0' }"

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def _args = task.ext.args ?: ''
    """
    #!/usr/bin/env python

    # This script takes the output from various tools' versions and arranges it
    # into a YAML file for MultiQC.

    import os
    import yaml
    import re
    import glob

    # Collect versions from all version.yml files
    versions_dict = {}
    for version_file in glob.glob("$versions"):
        with open(version_file) as f:
            versions_dict.update(yaml.safe_load(f))

    # Create a simplified version for software_versions.yml and MultiQC
    software_versions = {}
    for tool, version in versions_dict.items():
        # Clean up the tool name
        tool_clean = tool.split(':')[-1].strip()  # Remove any process prefix
        # Handle the case where version is a dict
        if isinstance(version, dict):
            for subtool, subversion in version.items():
                software_versions[subtool] = subversion
        else:
            software_versions[tool_clean] = version

    # Write to YAML files
    with open("software_versions.yml", "w") as f:
        yaml.dump(software_versions, f, default_flow_style=False)

    # Create MultiQC YAML with a more MultiQC-friendly structure
    mqc_versions = {'section_name': 'Software Versions'}
    software_mqc = {}
    for tool, version in software_versions.items():
        software_mqc[tool] = {'software': version}
    mqc_versions['software_versions'] = software_mqc

    with open("software_versions_mqc.yml", "w") as f:
        yaml.dump(mqc_versions, f, default_flow_style=False)

    # Output versions file for this process
    import platform
    with open("versions.yml", "w") as f:
        yaml.dump({"${task.process}": {"python": platform.python_version()}}, f, default_flow_style=False)
    """
}
