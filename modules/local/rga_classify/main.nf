process RGA_CLASSIFY {
    tag "$meta.id"
    label 'process_medium' // CPU-only; the reference R570 run is ~72s / 2.84GB peak RSS at --workers 4 (ARCHITECTURE.md §14)

    // Set container to use for this process
    container 'ghcr.io/omatheuspimenta/rga_classify:0.0.2'
    // container 'rga_classify:baseline' // local dev build

    input:
    tuple val(meta), path(interproscan_tsv), path(phobius_tsv), path(deeptmhmm_gff3), path(signalp_txt), path(deeploc_csv), path(deepcoil_dir, stageAs: 'deepcoil_data'), val(run_settings)

    output:
    tuple val(meta), path("rga_out/*"), emit: results
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Pipeline-level settings that change some tools' numbers (chunking, CPU vs GPU --
    // see docs/usage.md#chunking-and-reproducibility) but that rga_classify itself never
    // sees. Kept on one line: a multi-line value would break the script's indentation.
    def settings_json = groovy.json.JsonOutput.toJson(run_settings ?: [:])
    """
    # Own output dir deliberately not named 'results' -- DEEPCOIL2's directory-shaped input
    # (staged here as deepcoil_data/, see the stageAs above) is itself named 'results' at
    # its source (modules/local/deepcoil2), which would otherwise collide.
    mkdir -p rga_out/

    # rga_classify is vendored (docker/rga_classify/, from omatheuspimenta/rgapredictor)
    # and run exactly as documented upstream -- explicit --<tool> paths (not
    # --input-dir glob discovery) point it straight at each of this pipeline's own tool
    # outputs, since all six already match the formats it expects (see PLAN.md Stage 6).
    uv run --project /opt/rga_classify python /opt/rga_classify/code/rgas/rgas_prediction.py \\
        --interproscan ${interproscan_tsv} \\
        --phobius ${phobius_tsv} \\
        --deeptmhmm ${deeptmhmm_gff3} \\
        --signalp ${signalp_txt} \\
        --deeploc ${deeploc_csv} \\
        --deepcoil deepcoil_data \\
        --outdir rga_out/ \\
        --organism-name ${prefix} \\
        --workers ${task.cpus} \\
        ${args}

    # Fold those settings into rga_classify's own run_metadata.json, so one file
    # records everything needed to reproduce (or explain differences between) runs.
    cat <<-'END_SETTINGS' > pipeline_settings.json
    ${settings_json}
    END_SETTINGS
    uv run --project /opt/rga_classify python - <<-'END_PY'
    import json
    path = "rga_out/run_metadata.json"
    with open(path) as f:
        metadata = json.load(f)
    with open("pipeline_settings.json") as f:
        metadata["pipeline_settings"] = json.load(f)
    with open(path, "w") as f:
        f.write(json.dumps(metadata, indent=2, sort_keys=True, default=str))
    END_PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rga_classify: \$(uv run --project /opt/rga_classify python /opt/rga_classify/code/rgas/rgas_prediction.py --version | awk '{print \$NF}')
    END_VERSIONS
    """
}
