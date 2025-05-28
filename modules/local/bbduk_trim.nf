process BBDUK_TRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0' :
        'biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(reads)
    path(adapters)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")           , emit: log
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 minlen=36'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? 
        "in=${reads[0]}" : 
        "in1=${reads[0]} in2=${reads[1]}"
    def output_command = meta.single_end ? 
        "out=${prefix}.trim.fastq.gz" : 
        "out1=${prefix}.trim.R1.fastq.gz out2=${prefix}.trim.R2.fastq.gz"

    """
    bbduk.sh \\
        $reads_command \\
        $output_command \\
        ref=$adapters \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g \\
        2> ${prefix}.bbduk.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh 2>&1 | grep -o 'BBMap version [0-9.]*' | cut -d' ' -f3)
    END_VERSIONS
    """
}
