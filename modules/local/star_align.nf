process STAR_ALIGN_P1_REF {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.10a bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
    tuple val(meta), path(reads)
    path(index)
    path(vcf)
    path(vcf_index)

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("*.SJ.out.tab")                   , emit: sj
    tuple val(meta), path("*.Log.final.out")                , emit: log_final
    tuple val(meta), path("*.Log.out")                      , emit: log_out
    tuple val(meta), path("*.Log.progress.out")             , emit: log_progress
    path  "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? 
        "--readFilesIn ${reads[0]}" : 
        "--readFilesIn ${reads[0]} ${reads[1]}"
    def vcf_cmd = vcf ? "--varVCFfile $vcf" : ""
    def compression_cmd = reads[0].toString().endsWith('.gz') ? '--readFilesCommand zcat' : ''
    def memory = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ""
    
    """
    STAR \\
        --genomeDir $index \\
        $reads_command \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}.p1ref. \\
        $compression_cmd \\
        $vcf_cmd \\
        $memory \\
        $args

    # Index the BAM file
    samtools index ${prefix}.p1ref.Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process STAR_ALIGN_P1_MUT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.10a bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
    tuple val(meta), path(reads)
    path(index)
    path(vcf)
    path(vcf_index)

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("*.SJ.out.tab")                   , emit: sj
    tuple val(meta), path("*.Log.final.out")                , emit: log_final
    tuple val(meta), path("*.Log.out")                      , emit: log_out
    tuple val(meta), path("*.Log.progress.out")             , emit: log_progress
    path  "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? 
        "--readFilesIn ${reads[0]}" : 
        "--readFilesIn ${reads[0]} ${reads[1]}"
    def vcf_cmd = vcf ? "--varVCFfile $vcf" : ""
    def compression_cmd = reads[0].toString().endsWith('.gz') ? '--readFilesCommand zcat' : ''
    def memory = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ""
    
    """
    STAR \\
        --genomeDir $index \\
        $reads_command \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}.p1mut. \\
        $compression_cmd \\
        $vcf_cmd \\
        $memory \\
        $args

    # Index the BAM file
    samtools index ${prefix}.p1mut.Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process AGGREGATE_SJ_P1 {
    tag "aggregate_sj"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path(sj_files)

    output:
    path("filtered_junctions.tab"), emit: filtered_sj
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--min_unique_reads 3 --min_samples 1'
    
    """
    python $projectDir/bin/filter_splice_junctions.py \\
        --input_files ${sj_files.join(' ')} \\
        --output filtered_junctions.tab \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

process STAR_ALIGN_P2 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.10a bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
    tuple val(meta), path(reads)
    path(index)
    path(gtf)
    path(sj_filtered)
    path(vcf_optional)
    path(vcf_index_optional)

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: bam_sorted
    tuple val(meta), path("*.SJ.out.tab")                   , emit: sj_tab
    tuple val(meta), path("*.Log.final.out")                , emit: log_final
    tuple val(meta), path("*.Log.out")                      , emit: log_out
    tuple val(meta), path("*.Log.progress.out")             , emit: log_progress
    path  "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? 
        "--readFilesIn ${reads[0]}" : 
        "--readFilesIn ${reads[0]} ${reads[1]}"
    def compression_cmd = reads[0].toString().endsWith('.gz') ? '--readFilesCommand zcat' : ''
    def memory = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ""
    def sj_cmd = sj_filtered ? "--sjdbFileChrStartEnd $sj_filtered" : ""
    def vcf_cmd = vcf_optional ? "--varVCFfile $vcf_optional" : ""
    
    """
    STAR \\
        --genomeDir $index \\
        $reads_command \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}.p2. \\
        --sjdbGTFfile $gtf \\
        $sj_cmd \\
        $compression_cmd \\
        $vcf_cmd \\
        $memory \\
        $args

    # Index the BAM file
    samtools index ${prefix}.p2.Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
