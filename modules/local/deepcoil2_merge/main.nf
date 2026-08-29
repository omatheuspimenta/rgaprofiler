process DEEPCOIL2_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Just collects already-computed per-protein .out files, one per input sequence
    // (DeepCoil2's own output shape -- see modules/local/deepcoil2), into a single
    // directory -- no name collisions are possible since every protein appears in
    // exactly one chunk. Reuses deepcoil2's own image rather than building a
    // dedicated one for `cp`.
    container 'ghcr.io/omatheuspimenta/deepcoil2:2.0.2'
    // container 'deepcoil:nextflow' // local dev build

    input:
    // stageAs with a wildcard avoids a name collision: every chunk's DEEPCOIL2 task
    // independently names its output directory "results" (same directory name, since
    // every chunk shares the same sample meta), so without this they'd all try to
    // stage under the identical directory name in this task's work dir.
    tuple val(meta), path(dirs, stageAs: 'chunk_?')

    output:
    tuple val(meta), path("results/*"), emit: predictions
    tuple val(meta), path("results")  , emit: results_dir

    script:
    """
    mkdir -p results/
    cp chunk_*/*.out results/
    """
}
