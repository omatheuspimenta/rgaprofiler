process DEEPCOIL2 {
    tag "$meta.id"
    label 'process_medium' // Change this label to 'process_medium' to reflect the resource requirements of the process

    // Set container to use for this process
    container 'deepcoil:nextflow' 
    // container 'ghcr.io/omatheuspimenta/deepcoil2:2.0.2'

    input:
    // Load the input fasta file and metadata. One chuck of fasta file is processed at a time.
    // Nextflow will handle the parallelization of the chunks.
    tuple val(meta), path(fasta)

    output:
    // DeepCoil2 will output a directory with the predictions for each input fasta file. 
    // The output will be a tuple containing the metadata and the path to the results directory.
    tuple val(meta), path("results/*"), emit: predictions
    path "versions.yml"               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results/

    # Call deepcoil with the input fasta file and specify the output path
    deepcoil -i ${fasta} -out_path results/

    # Get the version of deepcoil and write it to a versions.yml file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepcoil: \$(pip show deepcoil | grep Version | awk '{print \$2}')
    END_VERSIONS
    """
}