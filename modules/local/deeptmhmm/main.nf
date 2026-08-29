process DEEPTMHMM {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu' // DeepTMHMM auto-detects GPU via torch.cuda.is_available() -- no CLI flag,
                         // see params.use_gpu / docs/software-setup.md.

    // Set container to use for this process
    container 'ghcr.io/omatheuspimenta/deeptmhmm:1.0'
    // container 'deeptmhmm:baseline' // local dev build

    input:
    tuple val(meta), path(fasta)
    path deeptmhmm_weights // softwares_dir/DeepTMHMM/DeepTMHMM-Academic-License-v1.0 -- DTU/BioLib model weights
                            // (CC BY-NC-SA 4.0 academic license, not baked into the image)

    output:
    tuple val(meta), path("results/*"), emit: predictions
    path "versions.yml"                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Nextflow's docker profile runs containers as the host UID, which has no
    # passwd entry in the image and thus no \$HOME -- point it at the task work
    # dir (always writable) so matplotlib doesn't warn about an unwritable
    # ~/.config/matplotlib.
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD
    task_dir=\$PWD

    # predict.py resolves its 8 weight files (5 CV checkpoints + 3 ESM1b files)
    # by bare relative filename in the CWD, with no --model-dir-style override
    # (unlike deeploc2/signalp6) -- so it must be run from /opt/deeptmhmm/ with
    # the real weights symlinked in there first.
    for f in deeptmhmm_cv_0.model deeptmhmm_cv_1.model deeptmhmm_cv_2.model \\
             deeptmhmm_cv_3.model deeptmhmm_cv_4.model esm_model_args.pt \\
             esm_model_alphabet.pt esm_model_state_dict.pt; do
        ln -sfn "\$(realpath "\${task_dir}/${deeptmhmm_weights}/\${f}")" /opt/deeptmhmm/\${f}
    done

    # predict.py also creates --output-dir itself and errors if it already exists.
    cd /opt/deeptmhmm
    python predict.py --fasta "\${task_dir}/${fasta}" --output-dir "\${task_dir}/results"
    cd "\${task_dir}"

    mv results/TMRs.gff3 results/${prefix}_deeptmhmm.gff3
    mv results/predicted_topologies.3line results/${prefix}_predicted_topologies.3line

    # Hardcoded: predict.py has no --version flag (only --fasta/--output-dir),
    # and DeepTMHMM is vendored as source, not pip installed. 1.0 matches the
    # DTU/BioLib "DeepTMHMM 1.0 - Academic Version" release this was vendored from.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptmhmm: 1.0
    END_VERSIONS
    """
}
