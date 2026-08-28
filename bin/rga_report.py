#!/opt/rga_classify/.venv/bin/python3
"""Render this pipeline's lightweight custom summary report (Stage 7).

Shebang points directly at rga_classify's uv-managed venv (this script's only
container, see modules/local/rga_report/main.nf) rather than `/usr/bin/env
python3`, since Nextflow's automatic bin/ PATH-staging makes this script
callable by bare name from any process, but does not put that venv's own bin/
on PATH -- a generic env-python3 shebang would hit the venv-less system
interpreter and fail to import pandas/yaml.

Per the maintainer's decision, this pipeline does not use MultiQC (most of its
tools have no native MultiQC module anyway and would need custom parsers
regardless). Instead this reads rga_classify's own harmonised outputs
(rga_predictions.tsv / rga_summary_counts.tsv -- already a single source of
truth combining all six tools' evidence per protein, see docs/rga/README.md
Sec7.1 vendored at docker/rga_classify/src/) plus the pipeline's collated
software_versions.yml, and writes one self-contained HTML page: no external
CSS/JS/CDN, matching the house style already set by rga_classify's own
report.html.

Pipeline-authored code, not vendored from anywhere.
"""

from __future__ import annotations

import argparse
import html
from pathlib import Path

import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True, help="rga_predictions.tsv")
    parser.add_argument("--summary", type=Path, required=True, help="rga_summary_counts.tsv")
    parser.add_argument("--versions", type=Path, required=True, help="collated software_versions.yml")
    parser.add_argument("--sample-name", required=True, help="sample/organism label")
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def load_versions(path: Path) -> list[tuple[str, str, str]]:
    """Flatten the collated versions.yml into (process, tool, version) rows."""
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    rows = []
    for process, tools in raw.items():
        if not isinstance(tools, dict):
            continue
        for tool, version in tools.items():
            rows.append((process, str(tool), str(version)))
    return sorted(rows)


#: Harmonised features that only ever come from InterProScan's own analyses (Pfam/SMART/
#: Gene3D/etc. domain models, and Coils via its 'Coil' accession). NOT 'CC' as a whole --
#: a CC call can also come from DeepCoil2 alone, so 'CC' in features_found doesn't by
#: itself mean InterProScan contributed anything; cc_coils/cc_rx_domain say that precisely.
_INTERPROSCAN_ONLY_FEATURES = ("NB-ARC", "TIR", "RPW8", "LRR", "STTK", "LysM")


def tool_contributions(predictions: pd.DataFrame) -> list[tuple[str, int, int]]:
    """One row per underlying annotation tool: (tool, n_with_evidence, n_total).

    Deliberately reads the specific per-tool boolean/count columns rga_classify already
    computed (sp_phobius, cc_deepcoil, ...) rather than the harmonised `features_found`
    string, which mixes evidence from every tool behind each of the nine controlled-
    vocabulary features (e.g. 'CC' alone doesn't say whether DeepCoil2, InterProScan
    Coils, or the Rx domain model -- or several of them -- called it).
    """
    n = len(predictions)

    def count(mask: pd.Series) -> int:
        return int(mask.sum())

    has_domain_feature = predictions["features_found"].apply(
        lambda cell: cell != "NA"
        and any(feat in cell.split(";") for feat in _INTERPROSCAN_ONLY_FEATURES)
    )
    interproscan_evidence = has_domain_feature | predictions["cc_coils"] | predictions["cc_rx_domain"]

    return [
        ("InterProScan (domain hit, incl. Coils/Rx domain)", count(interproscan_evidence), n),
        (
            "Phobius (signal peptide or TM)",
            count((predictions["sp_phobius"]) | (predictions["n_tm_phobius"].astype(int) > 0)),
            n,
        ),
        ("DeepTMHMM (TM helix)", count(predictions["n_tm_deeptmhmm"].astype(int) > 0), n),
        ("SignalP 6.0 (signal peptide)", count(predictions["sp_signalp"]), n),
        ("DeepCoil2 (coiled coil)", count(predictions["cc_deepcoil"]), n),
        (
            "DeepLoc 2.0 (localisation predicted)",
            count(predictions["predicted_localization"] != "NA"),
            n,
        ),
    ]


def render_html(
    sample_name: str,
    predictions: pd.DataFrame,
    summary: pd.DataFrame,
    versions: list[tuple[str, str, str]],
) -> str:
    n_total = len(predictions)
    n_rga = int((predictions["is_rga"] == True).sum())  # noqa: E712 (pandas bool column)
    pct_rga = f"{100 * n_rga / n_total:.1f}" if n_total else "0.0"

    family_rows = summary[summary["level"] == "family"].sort_values(
        "n_proteins", ascending=False
    )
    subclass_rows = summary[summary["level"] == "subclass"].sort_values(
        "n_proteins", ascending=False
    )

    def esc(x) -> str:
        return html.escape(str(x))

    def table(headers: list[str], rows: list[list]) -> str:
        thead = "".join(f"<th>{esc(h)}</th>" for h in headers)
        tbody = "".join(
            "<tr>" + "".join(f"<td>{esc(c)}</td>" for c in row) + "</tr>" for row in rows
        )
        return f"<table><thead><tr>{thead}</tr></thead><tbody>{tbody}</tbody></table>"

    family_table = table(
        ["RGA family", "proteins", "% of proteome"],
        [
            [r.rga_family, r.n_proteins, f"{r.percent_of_proteome:.2f}"]
            for r in family_rows.itertuples()
        ],
    )
    subclass_table = table(
        ["RGA subclass", "proteins", "% of proteome"],
        [
            [r.rga_subclass, r.n_proteins, f"{r.percent_of_proteome:.2f}"]
            for r in subclass_rows.itertuples()
            if r.rga_family != "Non-RGA"
        ],
    )
    tool_table = table(
        ["Tool", "proteins with evidence", "% of proteome"],
        [
            [tool, n_hit, f"{100 * n_hit / n_total:.1f}" if n_total else "0.0"]
            for tool, n_hit, n_total in tool_contributions(predictions)
        ],
    )
    versions_table = table(
        ["Process", "Tool", "Version"], [[p, t, v] for p, t, v in versions]
    )

    links = "".join(
        f'<li><a href="{href}">{label}</a></li>'
        for label, href in [
            ("Full RGA classification report (rga_classify)", "../rga/rga_out/report.html"),
            ("Per-protein predictions (rga_predictions.tsv)", "../rga/rga_out/rga_predictions.tsv"),
            ("InterProScan domain hits", "../interproscan/"),
            ("Phobius", "../phobius/"),
            ("DeepTMHMM", "../deeptmhmm/"),
            ("SignalP 6.0", "../signalp6/"),
            ("DeepLoc 2.0", "../deeploc2/"),
            ("DeepCoil2", "../deepcoil2/"),
        ]
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>rgaprofiler summary -- {esc(sample_name)}</title>
<style>
  body {{ font-family: -apple-system, Helvetica, Arial, sans-serif; max-width: 960px;
         margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; background: #fff; }}
  h1 {{ border-bottom: 3px solid #2c5f2d; padding-bottom: .3rem; }}
  h2 {{ margin-top: 2rem; color: #2c5f2d; }}
  .stat {{ display: inline-block; margin: .5rem 2rem .5rem 0; }}
  .stat .n {{ font-size: 2rem; font-weight: bold; display: block; }}
  .stat .label {{ color: #555; font-size: .9rem; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; }}
  th, td {{ border: 1px solid #ddd; padding: .4rem .6rem; text-align: left; }}
  th {{ background: #2c5f2d; color: #fff; }}
  tr:nth-child(even) {{ background: #f6f6f6; }}
  code {{ background: #f0f0f0; padding: .1rem .3rem; border-radius: 3px; }}
</style>
</head>
<body>
<h1>rgaprofiler summary report</h1>
<p>Sample / organism: <code>{esc(sample_name)}</code></p>

<div class="stat"><span class="n">{n_total}</span><span class="label">proteins analysed</span></div>
<div class="stat"><span class="n">{n_rga}</span><span class="label">RGA candidates ({pct_rga}%)</span></div>

<h2>RGA families</h2>
{family_table}

<h2>RGA subclasses (excluding Non-RGA)</h2>
{subclass_table}

<h2>Per-tool contribution</h2>
<p>Proteins for which each underlying tool reported evidence feeding into the RGA call
(see <code>rga_predictions.tsv</code> for the per-protein detail).</p>
{tool_table}

<h2>Detailed outputs</h2>
<ul>
{links}
</ul>

<h2>Software versions</h2>
{versions_table}

</body>
</html>
"""


def main() -> int:
    args = parse_args()
    predictions = pd.read_csv(args.predictions, sep="\t", keep_default_na=False)
    summary = pd.read_csv(args.summary, sep="\t", keep_default_na=False)
    versions = load_versions(args.versions)

    args.outdir.mkdir(parents=True, exist_ok=True)
    report_path = args.outdir / "report.html"
    report_path.write_text(
        render_html(args.sample_name, predictions, summary, versions), encoding="utf-8"
    )
    print(f"wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
