process FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::subread=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0' :
        'quay.io/biocontainers/subread:2.0.1--hed695b0_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(gtf)

    output:
    tuple val(meta), path("*.featureCounts.txt")        , emit: counts
    tuple val(meta), path("*.featureCounts.txt.summary"), emit: summary
    path  "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Map strandedness to featureCounts parameter
    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    
    def paired_end = meta.single_end ? '' : '-p'
    
    """
    featureCounts \\
        $args \\
        -a $gtf \\
        -o ${prefix}.featureCounts.txt \\
        -T $task.cpus \\
        -s $strandedness \\
        $paired_end \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -o 'featureCounts v[0-9.]*' | sed 's/featureCounts v//')
    END_VERSIONS
    """
}
