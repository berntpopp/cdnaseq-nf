process STAR_INDEX {
    tag "$index_name"
    label 'process_high'

    conda "bioconda::star=2.7.10a"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.10a--h9ee0642_0' :
        'quay.io/biocontainers/star:2.7.10a--h9ee0642_0' }"

    input:
    path(fasta)
    path(gtf)
    val(sjdb_overhang)
    val(index_name)
    path(sjdb_file)

    output:
    path("star_index_${index_name}"), emit: index
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sjdb_cmd = sjdb_file ? "--sjdbFileChrStartEnd $sjdb_file" : ""
    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ""
    
    """
    mkdir star_index_${index_name}

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index_${index_name} \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf \\
        --sjdbOverhang $sjdb_overhang \\
        --runThreadN $task.cpus \\
        $sjdb_cmd \\
        $memory \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
