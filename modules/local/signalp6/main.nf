process SIGNALP6 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu' // SignalP6 supports GPU execution; see params.use_gpu / docs/software-setup.md.
                         // Unlike deeploc2, SignalP6 has no --device flag: it auto-detects GPU from
                         // whether the *weight files themselves* were GPU-converted (signalp6_convert_models),
                         // which softwares/SignalP6/'s default install is not -- so this only takes effect
                         // once a GPU-converted model_dir is supplied. See docs/software-setup.md.

    // Set container to use for this process
    container 'signalp6:baseline'
    // container 'ghcr.io/omatheuspimenta/signalp6:6.0h'

    input:
    tuple val(meta), path(fasta)
    path signalp6_models // softwares_dir/SignalP6/signalp-6-package/models -- DTU model weights (license-gated, not baked into the image)

    output:
    tuple val(meta), path("results/*"), emit: predictions
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Nextflow's docker profile runs containers as the host UID, which has no
    # passwd entry in the image and thus no \$HOME -- point it at the task work
    # dir (always writable) so matplotlib doesn't warn about an unwritable
    # ~/.config/matplotlib.
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD

    # This repo's softwares/SignalP6/ only carries the 'slow-sequential' mode's
    # weights (the only mode whose files are actually installed -- 'fast' and
    # 'slow' are separate downloads, see docs/software-setup.md), hence the
    # --mode default here. --organism other and --format none match how the
    # committed R570 ground truth was generated (models/README.md, signalp
    # --organism default).
    signalp6 \\
        --fastafile ${fasta} \\
        --output_dir results/ \\
        --model_dir ${signalp6_models} \\
        --organism other \\
        --mode slow-sequential \\
        --format none \\
        ${args}
    mv results/prediction_results.txt results/${prefix}_signalp6_predictions.txt
    mv results/output.gff3 results/${prefix}_signalp6.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp6: \$(signalp6 --version | awk '{print \$NF}')
    END_VERSIONS
    """
}
