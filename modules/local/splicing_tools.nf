process REGTOOLS_JUNCTIONS_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::regtools=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/regtools:1.0.0--h7d7f7ad_2' :
        'quay.io/biocontainers/regtools:1.0.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference_fasta)
    path(reference_fai)

    output:
    tuple val(meta), path("*.junc"), emit: junctions
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    regtools junctions extract \\
        $args \\
        -s 0 \\
        -o ${prefix}.junc \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(regtools --version 2>&1 | grep -o 'Version.*' | cut -f2 -d' ')
    END_VERSIONS
    """
}

process REGTOOLS_JUNCTIONS_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::regtools=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/regtools:1.0.0--h7d7f7ad_2' :
        'quay.io/biocontainers/regtools:1.0.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(junctions)
    path(reference_fasta)
    path(reference_fai)
    path(gtf)

    output:
    tuple val(meta), path("*.annotated.junc"), emit: annotated_junctions
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    regtools junctions annotate \\
        $args \\
        -o ${prefix}.annotated.junc \\
        $junctions \\
        $reference_fasta \\
        $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(regtools --version 2>&1 | grep -o 'Version.*' | cut -f2 -d' ')
    END_VERSIONS
    """
}

process PORTCULLIS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::portcullis=1.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/portcullis:1.2.2--h1b792b2_1' :
        'quay.io/biocontainers/portcullis:1.2.2--h1b792b2_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference_fasta)
    path(reference_fai)
    path(gtf)

    output:
    tuple val(meta), path("portcullis_out/*.pass.junctions.tab"), emit: junctions
    tuple val(meta), path("portcullis_out/*.pass.junctions.bed"), emit: bed
    tuple val(meta), path("portcullis_out")                     , emit: results
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Prepare BAM list
    echo $bam > bam_list.txt

    # Full Portcullis pipeline
    portcullis full \\
        $args \\
        --threads $task.cpus \\
        --output portcullis_out \\
        --reference $reference_fasta \\
        --annotation $gtf \\
        bam_list.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        portcullis: \$(portcullis --version 2>&1 | grep -o 'portcullis [0-9.]*' | cut -f2 -d' ')
    END_VERSIONS
    """
}
