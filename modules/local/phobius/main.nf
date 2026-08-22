process PHOBIUS {
    tag "$meta.id"
    label 'process_single'

    // Set the docker container to use for this process.
    container 'quay.io/phobius:local' 

    input:
    tuple val(meta), path(fasta)

    output:
    // Phobius will output a TSV file with the predictions for each input fasta file.
    tuple val(meta), path("*_phobius.tsv"), emit: predictions
    path "versions.yml"                   , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run Phobius with the input fasta file and specify the output path
    phobius.pl -short ${fasta} > ${prefix}_phobius.tsv

    # Get the version of Phobius and write it to a versions.yml file hardcoded to 1.01 since Phobius does not provide a version command.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phobius: 1.01
    END_VERSIONS
    """
}