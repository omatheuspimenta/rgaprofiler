process PHOBIUS_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Merges per-chunk Phobius results back into the exact file a single unchunked
    // Phobius run over the whole clean FASTA writes: one 'SEQENCE ID ...' header, then
    // one row per protein in the FASTA's own record order. Plain concatenation would NOT
    // give that order: FASTA_QC's seqkit split2 deals records out to the chunks
    // round-robin, so rows are re-emitted in the order of the clean FASTA (passed in
    // alongside) instead. Phobius scores every sequence independently, so the rows
    // themselves don't depend on chunking -- verified byte-identical (same md5, 299,732
    // lines) against the unchunked run on the full R570 proteome with --num_blocks 1000.
    // Reuses phobius's own image rather than building a dedicated one for mawk.
    container 'ghcr.io/omatheuspimenta/phobius:1.01'
    // container 'quay.io/phobius:local' // local dev build

    input:
    // stageAs with a wildcard avoids a name collision: every chunk's PHOBIUS task
    // independently names its output "<meta.id>_phobius.tsv" (same prefix, since every
    // chunk shares the same sample meta).
    tuple val(meta), path(tsvs, stageAs: 'chunk_?/*'), path(fasta)

    output:
    tuple val(meta), path("*_phobius.tsv"), emit: predictions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // ext.partial: set for PHOBIUS_MERGE_BATCH (conf/modules.config), which only ever
    // sees one batch of chunks (see splitForMerge in workflows/rgaprofiler.nf) -- so the
    // "every FASTA record has a row" check is only enforced by the final merge.
    def require_all = task.ext.partial ? 0 : 1
    """
    # The chunk TSV list comes from a bash array (never exec'd) through the builtin
    # printf into a file list read by awk itself, so no single command line grows with
    # the chunk count (see modules/local/deepcoil2_merge for the ARG_MAX failure mode).
    tsv_files=(chunk_*/*)
    printf '%s\\n' "\${tsv_files[@]}" > tsv_list.txt

    awk -v require_all=${require_all} '
        # 1st file: the clean FASTA -> record order, keyed by the ID Phobius reports
        # (the header up to the first whitespace).
        FILENAME == ARGV[1] { if (/^>/) order[++n] = substr(\$1, 2); next }
        # 2nd file: the list of chunk TSVs -> read every row of each.
        {
            file = \$0
            first = 1
            while ((getline line < file) > 0) {
                if (first) { first = 0; if (!have_header) { header = line; have_header = 1 }; continue }
                split(line, f, " ")
                if (f[1] in row) { print "ERROR: protein " f[1] " reported by more than one chunk" > "/dev/stderr"; failed = 1; exit 1 }
                row[f[1]] = line
            }
            close(file)
        }
        END {
            # awk still runs END after an exit in a rule above -- stop right away then.
            if (failed) exit 1
            if (!have_header) { print "ERROR: no Phobius rows to merge" > "/dev/stderr"; exit 1 }
            print header
            for (i = 1; i <= n; i++) {
                if (order[i] in row) { print row[order[i]]; printed++ }
                else if (require_all) { print "ERROR: no Phobius row for FASTA record " order[i] > "/dev/stderr"; exit 1 }
            }
            nrows = 0; for (k in row) nrows++
            if (printed != nrows) { print "ERROR: " (nrows - printed) " Phobius row(s) match no FASTA record" > "/dev/stderr"; exit 1 }
        }
    ' ${fasta} tsv_list.txt > ${prefix}_phobius.tsv
    rm -f tsv_list.txt
    """
}
