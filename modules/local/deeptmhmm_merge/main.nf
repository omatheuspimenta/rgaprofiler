process DEEPTMHMM_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Merges per-chunk DeepTMHMM results back into one result per sample, respecting
    // each file's own record structure -- never an arbitrary line/byte split:
    //   - *_deeptmhmm.gff3: a single '##gff-version 3' header followed by one
    //     '# ... // '-delimited block per protein -- keep the first chunk's header,
    //     concatenate every chunk's per-protein blocks.
    //   - *_predicted_topologies.3line: no shared header, just repeated
    //     header+seq+topology triplets -- plain concatenation is already lossless.
    //   - embeddings/ and probabilities/: per-protein intermediate files named by a
    //     content hash -- collected into one directory each (a hash collision, if it
    //     ever happened, would only occur for byte-identical content, so overwriting
    //     is harmless).
    //   - deeptmhmm_results.md: a run summary, not a per-protein record -- one copy per
    //     chunk is kept (informational only, not consumed downstream) rather than
    //     naively concatenated or dropped.
    // Reuses deeptmhmm's own image rather than building a dedicated one for these
    // text-file operations.
    container 'ghcr.io/omatheuspimenta/deeptmhmm:1.0'
    // container 'deeptmhmm:baseline' // local dev build

    input:
    // stageAs with a wildcard avoids a name collision: every chunk's DEEPTMHMM task
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
    mkdir -p results/embeddings results/probabilities results/summaries

    first_gff3=\$(ls chunk_*/*_deeptmhmm.gff3 | head -n1)
    head -n1 "\$first_gff3" > results/${prefix}_deeptmhmm.gff3
    for f in chunk_*/*_deeptmhmm.gff3; do
        tail -n +2 "\$f" >> results/${prefix}_deeptmhmm.gff3
    done

    cat chunk_*/*_predicted_topologies.3line > results/${prefix}_predicted_topologies.3line

    # embeddings/ always has one file per protein in practice; probabilities/ can be
    # genuinely empty (observed with the real DeepTMHMM image/weights this pipeline
    # uses) -- '[ -e ... ] &&' skips a directory with nothing to copy instead of
    # erroring on an unmatched glob (harmless under `set -e`: bash only aborts on the
    # failure of the *last* command in an && list, not an earlier one).
    for f in chunk_*/embeddings/*; do [ -e "\$f" ] && cp "\$f" results/embeddings/; done
    for f in chunk_*/probabilities/*; do [ -e "\$f" ] && cp "\$f" results/probabilities/; done

    i=1
    for f in chunk_*/deeptmhmm_results.md; do
        cp "\$f" results/summaries/chunk_\${i}_deeptmhmm_results.md
        i=\$((i+1))
    done
    """
}
