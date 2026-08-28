process DEEPCOIL2 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu' // DeepCoil2 supports GPU execution; see params.use_gpu / docs/software-setup.md

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
    tuple val(meta), path("results")  , emit: results_dir // whole directory, for tools (rga_classify) that expect a --deepcoil DIR rather than a file list
    path "versions.yml"               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // deepcoil defaults to CPU: --gpu is a real, separate opt-in flag on top of the
    // container just having device access (--gpus all alone doesn't make deepcoil use it).
    def gpu_flag = task.ext.use_gpu ? '--gpu' : "-n_cpu ${task.cpus}"
    """
    mkdir -p results/

    # Call deepcoil with the input fasta file and specify the output path
    deepcoil -i ${fasta} -out_path results/ ${gpu_flag}

    # Get the version of deepcoil and write it to a versions.yml file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepcoil: \$(pip show deepcoil | grep Version | awk '{print \$2}')
    END_VERSIONS
    """
}