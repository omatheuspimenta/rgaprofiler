process INTERPROSCAN_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Just concatenates already-computed TSVs (InterProScan's TSV format has no header
    // row, confirmed against a real run -- plain `cat` is a correct, exact merge here).
    // Reuses interproscan's own image rather than building a dedicated one for `cat`.
    container 'ghcr.io/omatheuspimenta/interproscan:5.78-109.0'
    // container 'quay.io/interproscan_base:local' // local dev build

    input:
    // stageAs with a wildcard avoids a name collision: every chunk's INTERPROSCAN task
    // independently names its output "<meta.id>_interpro.tsv" (same prefix, since every
    // chunk shares the same sample meta), so without this they'd all try to stage under
    // the identical filename in this task's work dir.
    tuple val(meta), path(tsvs, stageAs: 'chunk_?/*')

    output:
    tuple val(meta), path("*_interpro.tsv"), emit: tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Not 'cat \${tsvs}': Nextflow expands that into one argument per chunk on a single
    # cat call, which passes the kernel's ARG_MAX ("Argument list too long", exit 126)
    # once there are enough chunks. A glob into a bash array and the builtin printf
    # never exec; xargs splits the list into as many cat calls as fit, in order.
    tsv_files=(chunk_*/*)
    printf '%s\\0' "\${tsv_files[@]}" | xargs -0 cat -- > ${prefix}_interpro.tsv
    """
}
