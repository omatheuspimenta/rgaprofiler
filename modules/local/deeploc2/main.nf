process DEEPLOC2 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu' // DeepLoc2 supports GPU execution; see params.use_gpu / docs/software-setup.md

    // Set container to use for this process
    container 'ghcr.io/omatheuspimenta/deeploc2:1.0.0'
    // container 'deeploc2:baseline' // local dev build

    input:
    tuple val(meta), path(fasta)
    path deeploc2_models // softwares_dir/DeepLoc2/DeepLoc2/models -- DTU model checkpoints (license-gated, not baked into the image)
    path torch_cache      // softwares_dir/DeepLoc2/torch_cache -- ESM1b base encoder weights (fair-esm, not baked into the image)

    output:
    tuple val(meta), path("results/*.csv"), emit: predictions
    tuple val(meta), path("results")      , emit: results_dir // whole directory, for DEEPLOC2_MERGE to fan back in across chunks
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def device = task.ext.use_gpu ? 'cuda' : 'cpu'
    """
    mkdir -p results/

    # Nextflow's docker profile runs containers as the host UID, which has no
    # passwd entry in the image and thus no \$HOME -- point it at the task work
    # dir (always writable) so matplotlib/huggingface don't warn about an
    # unwritable ~/.cache.
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD

    # DeepLoc2's package_data model checkpoints and torch hub's checkpoint cache
    # are both resolved from fixed, code-relative locations (pkg_resources /
    # \$TORCH_HOME) rather than a CLI flag -- symlink the staged inputs there
    # instead of baking the license-gated/oversized weights into the image.
    ln -sfn \$(realpath ${deeploc2_models}) /opt/deeploc2/DeepLoc2/models
    export TORCH_HOME=\$(realpath ${torch_cache})

    deeploc2 -f ${fasta} -o results/ -m Fast -d ${device} ${args}
    mv results/results_*.csv results/${prefix}_deeploc2.csv

    # Hardcoded: DeepLoc2 is vendored as source (see docker/deeploc2/), not pip
    # installed, so there's no package metadata to query at runtime. 1.0.0 is
    # the version pinned in softwares/DeepLoc2/env.yml and DeepLoc2.egg-info.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeploc2: 1.0.0
    END_VERSIONS
    """
}
