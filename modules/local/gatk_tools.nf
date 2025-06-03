process CREATE_MUT_REF_FASTA {
    tag "$meta.patient_id"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'broadinstitute/gatk:4.4.0.0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(reference_fasta)
    path(reference_fai)
    path(reference_dict)

    output:
    tuple val(meta), path("*.mutated.fa")    , emit: fasta
    tuple val(meta), path("*.mutated.fa.fai"), emit: fai
    tuple val(meta), path("*.mutated.dict")  , emit: dict
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    
    """
    gatk FastaAlternateReferenceMaker \\
        -R $reference_fasta \\
        -O ${prefix}.mutated.fa \\
        -V $vcf \\
        $args

    # Clean FASTA headers to remove problematic characters
    sed -i 's/:/_/g' ${prefix}.mutated.fa

    # Index the mutated reference
    samtools faidx ${prefix}.mutated.fa

    # Create sequence dictionary
    gatk CreateSequenceDictionary \\
        -R ${prefix}.mutated.fa \\
        -O ${prefix}.mutated.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process INDEX_FASTA {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::samtools=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    path(fasta)

    output:
    path("*.fai")        , emit: fai
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools faidx $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process CREATE_SEQ_DICT {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'broadinstitute/gatk:4.4.0.0' }"

    input:
    path(fasta)

    output:
    path("*.dict")       , emit: dict
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def dict_name = fasta.baseName + ".dict"
    
    """
    gatk CreateSequenceDictionary \\
        -R $fasta \\
        -O $dict_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process SPLITNCIGARREADS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'broadinstitute/gatk:4.4.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference_fasta)
    path(reference_fai)
    path(reference_dict)

    output:
    tuple val(meta), path("*.splitn.bam"), path("*.splitn.bam.bai"), emit: bam
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--rf ReassignOneMappingQuality --RMQF 255 --RMQT 60 --U ALLOW_N_CIGAR_READS'
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    gatk SplitNCigarReads \\
        -R $reference_fasta \\
        -I $bam \\
        -O ${prefix}.splitn.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process GATK_HAPLOTYPECALLER_RNA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'broadinstitute/gatk:4.4.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference_fasta)
    path(reference_fai)
    path(reference_dict)
    path(target_bed)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_cmd = target_bed ? "-L $target_bed" : ""
    
    """
    gatk HaplotypeCaller \\
        -R $reference_fasta \\
        -I $bam \\
        -O ${prefix}.vcf.gz \\
        $intervals_cmd \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
