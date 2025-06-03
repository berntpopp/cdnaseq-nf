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
    
    // RAM limiting logic
    def limit_ram_cmd = ""
    if (params.star_index_limit_genome_generate_ram) {
        limit_ram_cmd = "--limitGenomeGenerateRAM ${params.star_index_limit_genome_generate_ram}"
        log.info "[STAR_INDEX] Using --limitGenomeGenerateRAM ${params.star_index_limit_genome_generate_ram}"
    } else if (task.memory) {
        // Fallback to task.memory if specific param not set, ensuring it's a positive value
        def calculated_ram = task.memory.toBytes() - 100000000 // Subtract 100MB buffer
        if (calculated_ram > 0) {
            limit_ram_cmd = "--limitGenomeGenerateRAM ${calculated_ram}"
            log.info "[STAR_INDEX] Using task.memory based --limitGenomeGenerateRAM ${calculated_ram}"
        } else {
            log.info "[STAR_INDEX] task.memory too low for fallback RAM limit, STAR will use its default."
        }    } else {
        log.info "[STAR_INDEX] No RAM limit specified, STAR will use its default."
    }
    
    def sparse_sa_cmd = ""
    if (params.star_index_genome_sa_sparse_d != null) {
        sparse_sa_cmd = "--genomeSAsparseD ${params.star_index_genome_sa_sparse_d}"
        log.info "[STAR_INDEX] Using --genomeSAsparseD ${params.star_index_genome_sa_sparse_d}"
    }
    
    def chr_bin_nbits_cmd = ""
    if (params.star_index_genome_chr_bin_nbits != null) {
        chr_bin_nbits_cmd = "--genomeChrBinNbits ${params.star_index_genome_chr_bin_nbits}"
        log.info "[STAR_INDEX] Using --genomeChrBinNbits ${params.star_index_genome_chr_bin_nbits}"
    }
    
    // Custom temporary directory handling for FIFO file compatibility with unique subdirectory per task
    def effective_star_temp_dir_path = ""
    def temp_dir_cmd = ""
    // Check if params.star_temp_dir is not null, not empty string, and not the string "null"
    if (params.star_temp_dir && params.star_temp_dir.toString().trim() != "" && params.star_temp_dir.toString() != "null") {
        effective_star_temp_dir_path = "${params.star_temp_dir}/${task.process.replaceAll(':','_')}_${index_name}_${task.attempt}"
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
    
    mkdir star_index_${index_name}

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index_${index_name} \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf \\
        --sjdbOverhang $sjdb_overhang \\
        --runThreadN $task.cpus \\
        $sjdb_cmd \\
        $limit_ram_cmd \\
        $sparse_sa_cmd \\
        $chr_bin_nbits_cmd \\
        $temp_dir_cmd \\
        $args

    # Clean up the specific temp dir after STAR finishes, if it was created
    if [[ -n "${effective_star_temp_dir_path}" ]]; then
        rm -rf "${effective_star_temp_dir_path}" || echo "Warning: Could not remove temp dir ${effective_star_temp_dir_path}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
