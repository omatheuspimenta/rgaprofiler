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
    //     content hash -- collected into one directory each. The same name DOES occur in
    //     several chunks whenever identical sequences land in different chunks (179 times
    //     across 8 chunks of 2400 real R570 proteins, always byte-identical content); the
    //     last chunk's copy wins, as with a plain per-file copy loop.
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

    # glob array instead of 'ls | head -n1': head closing the pipe early makes ls die
    # of SIGPIPE (exit 141), which 'pipefail' then turns into a task failure.
    gff3_files=(chunk_*/*_deeptmhmm.gff3)
    first_gff3="\${gff3_files[0]}"
    head -n1 "\$first_gff3" > results/${prefix}_deeptmhmm.gff3
    for f in chunk_*/*_deeptmhmm.gff3; do
        tail -n +2 "\$f" >> results/${prefix}_deeptmhmm.gff3
    done

    # Every multi-file command below goes through a bash array + builtin printf |
    # xargs, never a glob handed straight to an external command: with enough chunks
    # (or one per-protein file per sequence, as in embeddings/) that argument list
    # passes the kernel's ARG_MAX ("Argument list too long", exit 126). Arrays and
    # printf never exec; xargs splits the list into as many calls as fit, in order.
    topology_files=(chunk_*/*_predicted_topologies.3line)
    printf '%s\\0' "\${topology_files[@]}" | xargs -0 cat -- > results/${prefix}_predicted_topologies.3line

    # embeddings/ always has one file per protein in practice; probabilities/ can be
    # genuinely empty (observed with the real DeepTMHMM image/weights this pipeline
    # uses) -- nullglob turns an unmatched glob into an empty array, which is then
    # skipped instead of handing cp a literal, nonexistent pattern.
    shopt -s nullglob
    embedding_files=(chunk_*/embeddings/*)
    probability_files=(chunk_*/probabilities/*)
    shopt -u nullglob
    # One cp call refuses to write the same destination name twice ("will not
    # overwrite just-created"), so keep only the last occurrence of each file name, in
    # glob order -- the copy a per-file cp loop would have left in place.
    last_by_name() { local -A last=(); local f; for f in "\$@"; do last[\${f##*/}]=\$f; done; printf '%s\\0' "\${last[@]}"; }
    if (( \${#embedding_files[@]} )); then
        last_by_name "\${embedding_files[@]}" | xargs -0 cp -t results/embeddings/ --
    fi
    if (( \${#probability_files[@]} )); then
        last_by_name "\${probability_files[@]}" | xargs -0 cp -t results/probabilities/ --
    fi

    # Also picks up summaries/*.md -- this process's own renamed copies -- since a sample
    # with too many chunks for one task to stage is merged in two levels (see
    # splitForMerge in workflows/rgaprofiler.nf), feeding merged results back in here.
    shopt -s nullglob
    summary_files=(chunk_*/deeptmhmm_results.md chunk_*/summaries/*.md)
    shopt -u nullglob
    i=1
    for f in "\${summary_files[@]}"; do
        cp "\$f" results/summaries/chunk_\${i}_deeptmhmm_results.md
        i=\$((i+1))
    done
    """
}
