process DEEPLOC2_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Merges per-chunk DeepLoc2 CSVs back into one CSV per sample: a fixed one-line
    // header followed by one row per protein -- keep the first chunk's header,
    // concatenate every chunk's data rows. Reuses deeploc2's own image rather than
    // building a dedicated one for this text-file operation.
    container 'ghcr.io/omatheuspimenta/deeploc2:1.0.0'
    // container 'deeploc2:baseline' // local dev build

    input:
    // stageAs with a wildcard avoids a name collision: every chunk's DEEPLOC2 task
    // independently names its output directory "results" (same directory name, since
    // every chunk shares the same sample meta), so without this they'd all try to
    // stage under the identical directory name in this task's work dir.
    tuple val(meta), path(dirs, stageAs: 'chunk_?')

    output:
    tuple val(meta), path("results/*.csv"), emit: predictions
    tuple val(meta), path("results")      , emit: results_dir

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results/

    first_csv=\$(ls chunk_*/*_deeploc2.csv | head -n1)
    head -n1 "\$first_csv" > results/${prefix}_deeploc2.csv
    for f in chunk_*/*_deeploc2.csv; do
        tail -n +2 "\$f" >> results/${prefix}_deeploc2.csv
    done
    """
}
