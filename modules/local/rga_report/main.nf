process RGA_REPORT {
    tag "$meta.id"
    label 'process_single'

    // Reuses rga_classify's image rather than building a new one: this is a small,
    // pipeline-authored script (bin/rga_report.py, not vendored from anywhere) that only
    // needs pandas + PyYAML to read rga_classify's own harmonised outputs, and that image
    // already has both in its uv-managed venv (see docker/rga_classify/Dockerfile).
    container 'rga_classify:baseline'

    input:
    tuple val(meta), path(rga_predictions_tsv), path(rga_summary_tsv), path(versions_yml)

    output:
    // Internal dir deliberately not named 'summary_report' -- the withName override in
    // conf/modules.config already publishes this process to .../summary_report/, and
    // publishDir preserves the emitted path's own subdirectory, so matching names would
    // double up into .../summary_report/summary_report/report.html.
    tuple val(meta), path("report_out/*"), emit: report
    path "versions.yml"                  , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\$PWD

    rga_report.py \\
        --predictions ${rga_predictions_tsv} \\
        --summary ${rga_summary_tsv} \\
        --versions ${versions_yml} \\
        --sample-name ${prefix} \\
        --outdir report_out/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(uv run --project /opt/rga_classify python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
