process MARK_DUPLICATES_PICARD {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.0--hdfd78af_0' :
        'quay.io/biocontainers/picard:3.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.marked.bam"), path("*.marked.bai"), emit: bam
    tuple val(meta), path("*.MarkDuplicates.metrics.txt")       , emit: metrics
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 'ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    mkdir -p tmp

    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args \\
        INPUT=$bam \\
        OUTPUT=${prefix}.marked.bam \\
        METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt

    # Index the marked BAM file
    picard \\
        -Xmx${avail_mem}M \\
        BuildBamIndex \\
        INPUT=${prefix}.marked.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}

process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.9--h10a08f8_12' :
        'biocontainers/samtools:1.9--h10a08f8_12' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("*.bai"), emit: bam
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    samtools index $args $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.9--h10a08f8_12' :
        'biocontainers/samtools:1.9--h10a08f8_12' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    samtools sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.sorted.bam \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process QUALIMAP_RNASEQ {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::qualimap=2.2.2d"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--1' :
        'quay.io/biocontainers/qualimap:2.2.2d--1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(gtf)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = 'non-strand-specific'
    def memory = task.memory ? "-Xmx${task.memory.toGiga()}G" : ''
    
    if (meta.strandedness == 'forward') {
        strandedness = 'strand-specific-forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = 'strand-specific-reverse'
    }

    """
    qualimap \\
        $memory \\
        rnaseq \\
        $args \\
        -bam $bam \\
        -gtf $gtf \\
        -p $strandedness \\
        -outdir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap //; s/Built.*\$//')
    END_VERSIONS
    """
}
