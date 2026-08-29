process FASTA_QC {
    tag "$meta.id"
    label 'process_single'

    // Set the conda environment or container to use for this process
    conda "bioconda::seqkit=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--he881be0_0' :
        'quay.io/biocontainers/seqkit:2.13.0--he881be0_0' }"

    input:
    tuple val(meta), path(fasta)
    val split_mode  // 'size' (split_value = sequences/chunk) or 'parts' (split_value = chunk count)
    val split_value

    output:
    tuple val(meta), path("*_clean.fasta")      , emit: fasta
    tuple val(meta), path("*_chunks/*.fasta")   , emit: chunks
    path "versions.yml"                         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fasta_qc.sh ${fasta} ${prefix}_clean.fasta ${split_mode} ${split_value}

    # nf-core modules: Add version information to output
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | awk '{print \$2}')
    END_VERSIONS
    """
}