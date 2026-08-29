process SIGNALP6_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Merges per-chunk SignalP6 results back into one result per sample, respecting
    // each file's own record structure -- never an arbitrary line/byte split:
    //   - *_signalp6_predictions.txt: a fixed 2-line header followed by one row per
    //     protein -- keep the first chunk's header, concatenate every chunk's data rows.
    //   - *_signalp6.gff3 / region_output.gff3: a single '##gff-version 3' header
    //     followed by one block per protein -- same header-once + concatenate approach.
    //   - processed_entries.fasta: plain FASTA concatenation is lossless here since
    //     every chunk contributes disjoint, complete records (each protein was scored
    //     in exactly one chunk).
    //   - output.json: a single JSON object per chunk (not a per-protein list), so it
    //     has no lossless line-level merge -- every chunk's copy is kept, individually,
    //     rather than silently dropped or naively concatenated into invalid JSON.
    // Reuses signalp6's own image rather than building a dedicated one for these
    // text-file operations.
    container 'ghcr.io/omatheuspimenta/signalp6:6.0h'
    // container 'signalp6:baseline' // local dev build

    input:
    // stageAs with a wildcard avoids a name collision: every chunk's SIGNALP6 task
    // independently names its output directory "results" (same directory name, since
    // every chunk shares the same sample meta), so without this they'd all try to
    // stage under the identical directory name in this task's work dir.
    tuple val(meta), path(dirs, stageAs: 'chunk_?')

    output:
    tuple val(meta), path("results/*"), emit: predictions
    tuple val(meta), path("results")  , emit: results_dir

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results/

    first_predictions=\$(ls chunk_*/*_signalp6_predictions.txt | head -n1)
    head -n2 "\$first_predictions" > results/${prefix}_signalp6_predictions.txt
    for f in chunk_*/*_signalp6_predictions.txt; do
        tail -n +3 "\$f" >> results/${prefix}_signalp6_predictions.txt
    done

    first_gff3=\$(ls chunk_*/*_signalp6.gff3 | head -n1)
    head -n1 "\$first_gff3" > results/${prefix}_signalp6.gff3
    for f in chunk_*/*_signalp6.gff3; do
        tail -n +2 "\$f" >> results/${prefix}_signalp6.gff3
    done

    first_region=\$(ls chunk_*/region_output.gff3 | head -n1)
    head -n1 "\$first_region" > results/region_output.gff3
    for f in chunk_*/region_output.gff3; do
        tail -n +2 "\$f" >> results/region_output.gff3
    done

    cat chunk_*/processed_entries.fasta > results/processed_entries.fasta

    i=1
    for f in chunk_*/output.json; do
        cp "\$f" results/chunk_\${i}_output.json
        i=\$((i+1))
    done
    """
}
