process INTERPROSCAN {
    tag "$meta.id"
    label 'process_high' // InterProScan is resource-intensive, so we label it as 'process_high' to allocate more resources.

    container 'ghcr.io/omatheuspimenta/interproscan:5.78-109.0'
    // container 'quay.io/interproscan_base:local' // local dev build

    input:
    tuple val(meta), path(fasta)
    path ipr_dir // The path to the InterProScan directory, which contains the interproscan.sh script and necessary databases.

    output:
    tuple val(meta), path("*_interpro.tsv"), emit: tsv
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Since we inject the path, the interproscan.sh script will be in the root of ipr_dir
    # The -dp flag disables online searches (Pre-calculated match lookup), essential to prevent crashes in clusters
    ${ipr_dir}/interproscan.sh \\
        -i ${fasta} \\
        -f TSV \\
        -o ${prefix}_interpro.tsv \\
        -cpu ${task.cpus} \\
        -dp \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(${ipr_dir}/interproscan.sh --version | grep 'InterProScan version' | sed 's/.*version //')
    END_VERSIONS
    """
}
