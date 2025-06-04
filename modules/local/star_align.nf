process STAR_ALIGN_P1_REF {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.10a bioconda::samtools=1.9"
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
    
    // Custom temporary directory handling for FIFO file compatibility with unique subdirectory per task
    def effective_star_temp_dir_path = ""
    def temp_dir_cmd = ""
    // Check if params.star_temp_dir is not null, not empty string, and not the string "null"
    if (params.star_temp_dir && params.star_temp_dir.toString().trim() != "" && params.star_temp_dir.toString() != "null") {
        effective_star_temp_dir_path = "${params.star_temp_dir}/${task.process.replaceAll(':','_')}_${meta.id}_${task.attempt}"
        temp_dir_cmd = "--outTmpDir ${effective_star_temp_dir_path}"
    }
    
    """
    # Prepare temporary directory path if specified
    if [[ -n "${effective_star_temp_dir_path}" ]]; then
        # Ensure parent directory exists
        mkdir -p "\$(dirname "${effective_star_temp_dir_path}")"
        
        # Remove the specific temp directory if it already exists from a previous run
        if [[ -d "${effective_star_temp_dir_path}" ]]; then
            rm -rf "${effective_star_temp_dir_path}"
            echo "Removed existing temporary directory: ${effective_star_temp_dir_path}"
        fi
        
        echo "Using custom STAR temporary directory: ${effective_star_temp_dir_path}"
    fi

    STAR \\
        --genomeDir $index \\
        $reads_command \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}.p1ref. \\
        $compression_cmd \\
        $vcf_cmd \\
        $memory \\
        $temp_dir_cmd \\
        $args

    # Index the BAM file
    samtools index ${prefix}.p1ref.Aligned.sortedByCoord.out.bam

    # Clean up the specific temp dir after STAR finishes, if it was created
    if [[ -n "${effective_star_temp_dir_path}" ]]; then
        rm -rf "${effective_star_temp_dir_path}" || echo "Warning: Could not remove temp dir ${effective_star_temp_dir_path}"
    fi

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

    conda "bioconda::star=2.7.10a bioconda::samtools=1.9"
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
    
    // Custom temporary directory handling for FIFO file compatibility with unique subdirectory per task
    def effective_star_temp_dir_path = ""
    def temp_dir_cmd = ""
    // Check if params.star_temp_dir is not null, not empty string, and not the string "null"
    if (params.star_temp_dir && params.star_temp_dir.toString().trim() != "" && params.star_temp_dir.toString() != "null") {
        effective_star_temp_dir_path = "${params.star_temp_dir}/${task.process.replaceAll(':','_')}_${meta.id}_${task.attempt}"
        temp_dir_cmd = "--outTmpDir ${effective_star_temp_dir_path}"
    }
    
    """
    # Prepare temporary directory path if specified
    if [[ -n "${effective_star_temp_dir_path}" ]]; then
        # Ensure parent directory exists
        mkdir -p "\$(dirname "${effective_star_temp_dir_path}")"
        
        # Remove the specific temp directory if it already exists from a previous run
        if [[ -d "${effective_star_temp_dir_path}" ]]; then
            rm -rf "${effective_star_temp_dir_path}"
            echo "Removed existing temporary directory: ${effective_star_temp_dir_path}"
        fi
        
        echo "Using custom STAR temporary directory: ${effective_star_temp_dir_path}"
    fi

    STAR \\
        --genomeDir $index \\
        $reads_command \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}.p1mut. \\
        $compression_cmd \\
        $vcf_cmd \\
        $memory \\
        $temp_dir_cmd \\
        $args

    # Index the BAM file
    samtools index ${prefix}.p1mut.Aligned.sortedByCoord.out.bam

    # Clean up the specific temp dir after STAR finishes, if it was created
    if [[ -n "${effective_star_temp_dir_path}" ]]; then
        rm -rf "${effective_star_temp_dir_path}" || echo "Warning: Could not remove temp dir ${effective_star_temp_dir_path}"
    fi

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
    path(filter_script)

    output:
    path("filtered_junctions.tab"), emit: filtered_sj
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--min_unique_reads 3 --min_samples 1'
    
    """
    python ${filter_script} \\
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

    conda "bioconda::star=2.7.10a bioconda::samtools=1.18"
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
    
    // Custom temporary directory handling for FIFO file compatibility with unique subdirectory per task
    def effective_star_temp_dir_path = ""
    def temp_dir_cmd = ""
    // Check if params.star_temp_dir is not null, not empty string, and not the string "null"
    if (params.star_temp_dir && params.star_temp_dir.toString().trim() != "" && params.star_temp_dir.toString() != "null") {
        effective_star_temp_dir_path = "${params.star_temp_dir}/${task.process.replaceAll(':','_')}_${meta.id}_${task.attempt}"
        temp_dir_cmd = "--outTmpDir ${effective_star_temp_dir_path}"
    }
    
    """
    # Prepare temporary directory path if specified
    if [[ -n "${effective_star_temp_dir_path}" ]]; then
        # Ensure parent directory exists
        mkdir -p "\$(dirname "${effective_star_temp_dir_path}")"
        
        # Remove the specific temp directory if it already exists from a previous run
        if [[ -d "${effective_star_temp_dir_path}" ]]; then
            rm -rf "${effective_star_temp_dir_path}"
            echo "Removed existing temporary directory: ${effective_star_temp_dir_path}"
        fi
        
        echo "Using custom STAR temporary directory: ${effective_star_temp_dir_path}"
    fi

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
        $temp_dir_cmd \\
        $args

    # Index the BAM file
    samtools index ${prefix}.p2.Aligned.sortedByCoord.out.bam

    # Clean up the specific temp dir after STAR finishes, if it was created
    if [[ -n "${effective_star_temp_dir_path}" ]]; then
        rm -rf "${effective_star_temp_dir_path}" || echo "Warning: Could not remove temp dir ${effective_star_temp_dir_path}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
