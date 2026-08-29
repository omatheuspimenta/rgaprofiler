process INTERPROSCAN_MERGE {
    tag "$meta.id"
    label 'process_single'

    // Just concatenates already-computed TSVs (InterProScan's TSV format has no header
    // row, confirmed against a real run -- plain `cat` is a correct, exact merge here).
    // Reuses interproscan's own image rather than building a dedicated one for `cat`.
    container 'quay.io/interproscan_base:local'

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
    cat ${tsvs} > ${prefix}_interpro.tsv
    """
}
