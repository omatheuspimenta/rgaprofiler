"""Summary statistics and report generation (HTML + Markdown).

Every number shown in ``report.html``, ``report.md``, ``rga_summary_counts.tsv``
and ``run_metadata.json`` comes from the single :class:`Summary` object built
here. Nothing is ever recomputed independently for a particular output, which
is what makes the cross-output consistency assertion meaningful.
"""

from __future__ import annotations

import html
import logging
import shlex
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import pandas as pd

LOGGER = logging.getLogger(__name__)

#: Bibliography rendered at the end of both reports.
REFERENCES: tuple[tuple[str, str], ...] = (
    (
        "Rody HVS, Bombardelli RGH, Creste S, Camargo LEA, Van Sluys M-A, "
        "Monteiro-Vitorello CB (2019). Genome survey of resistance gene analogs in "
        "sugarcane: genomic features and differential expression of the innate immune "
        "system from a smut-resistant genotype. BMC Genomics 20:809.",
        "10.1186/s12864-019-6207-y",
    ),
    (
        "Li P, Quan X, Jia G, Xiao J, Cloutier S, You FM (2016). RGAugury: a pipeline "
        "for genome-wide prediction of resistance gene analogs (RGAs) in plants. "
        "BMC Genomics 17:852.",
        "10.1186/s12864-016-3197-x",
    ),
    (
        "Sekhwal MK, Li P, Lam I, Wang X, Cloutier S, You FM (2015). Disease resistance "
        "gene analogs (RGAs) in plants. Int J Mol Sci 16:19248-19290.",
        "10.3390/ijms160819248",
    ),
    (
        "Kourelis J, Sakai T, Adachi H, Kamoun S (2021). RefPlantNLR is a comprehensive "
        "collection of experimentally validated plant disease resistance proteins from "
        "the NLR family. PLoS Biology 19(10):e3001124.",
        "10.1371/journal.pbio.3001124",
    ),
    (
        "Smith M, Jones JT, Hein I (2025). Resistify: a novel NLR classifier that reveals "
        "Helitron-associated NLR expansion in Solanaceae. "
        "Bioinform Biol Insights 19:11779322241308944.",
        "10.1177/11779322241308944",
    ),
    (
        "Shiu S-H, Bleecker AB (2003). Expansion of the receptor-like kinase/Pelle gene "
        "family and receptor-like proteins in Arabidopsis. Plant Physiol 132:530-543.",
        "10.1104/pp.103.021964",
    ),
    (
        "Jones JDG, Dangl JL (2006). The plant immune system. Nature 444:323-329.",
        "10.1038/nature05286",
    ),
    (
        "Jones P et al. (2014). InterProScan 5: genome-scale protein function "
        "classification. Bioinformatics 30:1236-1240.",
        "10.1093/bioinformatics/btu031",
    ),
    (
        "Blum M et al. (2025). InterPro: the protein sequence classification resource "
        "in 2025. Nucleic Acids Res 53:D444-D456.",
        "10.1093/nar/gkae1082",
    ),
    (
        "Paysan-Lafosse T et al. (2025). The Pfam protein families database: embracing "
        "AI/ML. Nucleic Acids Res 53:D523-D534.",
        "10.1093/nar/gkae997",
    ),
    (
        "Kall L, Krogh A, Sonnhammer ELL (2004). A combined transmembrane topology and "
        "signal peptide prediction method. J Mol Biol 338:1027-1036.",
        "10.1016/j.jmb.2004.03.016",
    ),
    (
        "Hallgren J et al. (2022). DeepTMHMM predicts alpha and beta transmembrane "
        "proteins using deep neural networks. bioRxiv.",
        "10.1101/2022.04.08.487609",
    ),
    (
        "Teufel F et al. (2022). SignalP 6.0 predicts all five types of signal peptides "
        "using protein language models. Nat Biotechnol 40:1023-1025.",
        "10.1038/s41587-021-01156-3",
    ),
    (
        "Thumuluri V et al. (2022). DeepLoc 2.0: multi-label subcellular localization "
        "prediction using protein language models. Nucleic Acids Res 50:W228-W234.",
        "10.1093/nar/gkac278",
    ),
    (
        "Ludwiczak J, Winski A, Szczepaniak K, Alva V, Dunin-Horkawicz S (2019). "
        "DeepCoil - a fast and accurate prediction of coiled-coil domains in protein "
        "sequences. Bioinformatics 35(16):2790-2795.",
        "10.1093/bioinformatics/bty1062",
    ),
    (
        "Lupas A, Van Dyke M, Stock J (1991). Predicting coiled coils from protein "
        "sequences. Science 252:1162-1164.",
        "10.1126/science.252.5009.1162",
    ),
    (
        "Simm D, Hatje K, Waack S, Kollmar M (2021). Critical assessment of "
        "coiled-coil predictions based on protein structure data. Scientific "
        "Reports 11:12439.",
        "10.1038/s41598-021-91886-w",
    ),
)


@dataclass
class Summary:
    """Every count reported by the run, computed exactly once.

    Attributes
    ----------
    n_proteins : int
        Number of proteins in the output table.
    n_rga : int
        Number of proteins with ``is_rga = True``.
    family_counts, subclass_counts : pandas.DataFrame
        Counts and percentages per family / per family+subclass.
    subclass_confidence : pandas.DataFrame
        Subclass x confidence cross-table. The headline count of a class and how
        much of it is trustworthy are different questions, and for TM-CC in
        particular they have very different answers.
    confidence_counts : pandas.DataFrame
        Counts per confidence level among RGAs.
    architecture_counts : pandas.DataFrame
        Most frequent domain architectures among RGAs.
    warning_counts : pandas.DataFrame
        Frequency of each distinct warning.
    """

    n_proteins: int
    n_rga: int
    family_counts: pd.DataFrame
    subclass_counts: pd.DataFrame
    confidence_counts: pd.DataFrame
    architecture_counts: pd.DataFrame
    warning_counts: pd.DataFrame
    subclass_confidence: pd.DataFrame = field(default_factory=pd.DataFrame)
    extras: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        """Return the summary as JSON-serialisable counts for ``run_metadata.json``."""
        return {
            "n_proteins": self.n_proteins,
            "n_rga": self.n_rga,
            "by_family": dict(
                zip(self.family_counts["rga_family"], self.family_counts["n_proteins"])
            ),
            "by_subclass": dict(
                zip(
                    self.subclass_counts["rga_subclass"],
                    self.subclass_counts["n_proteins"],
                )
            ),
            "by_confidence": dict(
                zip(
                    self.confidence_counts["confidence"],
                    self.confidence_counts["n_proteins"],
                )
            ),
            **{key: value for key, value in self.extras.items()},
        }


def summarize(predictions: pd.DataFrame, top_n: int = 20) -> Summary:
    """Compute every count the run reports.

    Parameters
    ----------
    predictions : pandas.DataFrame
        The per-protein prediction table.
    top_n : int, default 20
        Number of domain architectures kept in the architecture table.

    Returns
    -------
    Summary
        Count tables shared by all outputs.
    """
    total = len(predictions)
    rgas = predictions[predictions["is_rga"]]
    family, subclass = _class_counts(predictions, total)

    confidence = (
        rgas.groupby("confidence", dropna=False)
        .size()
        .rename("n_proteins")
        .reset_index()
        .sort_values("n_proteins", ascending=False, kind="stable")
    )
    architecture = (
        rgas.groupby("domain_architecture", dropna=False)
        .size()
        .rename("n_proteins")
        .reset_index()
        .sort_values(
            ["n_proteins", "domain_architecture"],
            ascending=[False, True],
            kind="stable",
        )
        .head(top_n)
    )
    subclass_confidence = _subclass_confidence(rgas)
    warnings = _warning_counts(predictions)
    LOGGER.info("Summary: %d proteins, %d RGAs", total, len(rgas))
    return Summary(
        n_proteins=total,
        n_rga=int(len(rgas)),
        family_counts=family.reset_index(drop=True),
        subclass_counts=subclass.reset_index(drop=True),
        confidence_counts=confidence.reset_index(drop=True),
        subclass_confidence=subclass_confidence,
        architecture_counts=architecture.reset_index(drop=True),
        warning_counts=warnings,
    )


def _subclass_confidence(rgas: pd.DataFrame) -> pd.DataFrame:
    """Cross-tabulate subclass against confidence, largest class first.

    The count of a class and the trustworthiness of that count are separate
    facts, and reporting only the first invites a reader to quote a number the
    grading already qualified.
    """
    if rgas.empty:
        return pd.DataFrame(
            columns=["rga_subclass", "high", "medium", "low", "n_proteins"]
        )
    table = (
        pd.crosstab(rgas["rga_subclass"], rgas["confidence"])
        .reindex(columns=["high", "medium", "low"], fill_value=0)
        .reset_index()
    )
    table["n_proteins"] = table[["high", "medium", "low"]].sum(axis=1)
    return (
        table.sort_values(
            ["n_proteins", "rga_subclass"], ascending=[False, True], kind="stable"
        )
        .reset_index(drop=True)
        .astype({"high": int, "medium": int, "low": int, "n_proteins": int})
    )


def _class_counts(
    predictions: pd.DataFrame, total: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Count proteins per family and per family x subclass, with percentages."""
    family = (
        predictions.groupby("rga_family", dropna=False)
        .size()
        .rename("n_proteins")
        .reset_index()
        .sort_values(
            ["n_proteins", "rga_family"], ascending=[False, True], kind="stable"
        )
    )
    family["percent_of_proteome"] = (100.0 * family["n_proteins"] / total).round(4)

    subclass = (
        predictions.groupby(["rga_family", "rga_subclass"], dropna=False)
        .size()
        .rename("n_proteins")
        .reset_index()
        .sort_values(
            ["rga_family", "n_proteins", "rga_subclass"],
            ascending=[True, False, True],
            kind="stable",
        )
    )
    subclass["percent_of_proteome"] = (100.0 * subclass["n_proteins"] / total).round(4)
    return family, subclass


def _warning_counts(predictions: pd.DataFrame) -> pd.DataFrame:
    """Count how often each distinct warning was raised."""
    tally: dict[str, int] = {}
    for value in predictions["warnings"].dropna():
        for warning in str(value).split(";"):
            warning = warning.strip()
            if warning:
                tally[warning] = tally.get(warning, 0) + 1
    frame = pd.DataFrame(
        {"warning": list(tally), "n_proteins": [tally[k] for k in tally]}
    )
    if frame.empty:
        return pd.DataFrame(
            {"warning": pd.Series(dtype=str), "n_proteins": pd.Series(dtype=int)}
        )
    return frame.sort_values(
        ["n_proteins", "warning"], ascending=[False, True], kind="stable"
    ).reset_index(drop=True)


def summary_table(summary: Summary) -> pd.DataFrame:
    """Build ``rga_summary_counts.tsv`` from the shared :class:`Summary`."""
    family = summary.family_counts.assign(level="family", rga_subclass="ALL")
    subclass = summary.subclass_counts.assign(level="subclass")
    columns = [
        "level",
        "rga_family",
        "rga_subclass",
        "n_proteins",
        "percent_of_proteome",
    ]
    return pd.concat([family[columns], subclass[columns]], ignore_index=True)


# ---------------------------------------------------------------------------
# rendering helpers
# ---------------------------------------------------------------------------


def _md_table(frame: pd.DataFrame, missing: str = "NA") -> str:
    """Render a DataFrame as a GitHub-flavoured Markdown table."""
    if frame.empty:
        return "_(no rows)_\n"
    header = "| " + " | ".join(str(c) for c in frame.columns) + " |"
    divider = "| " + " | ".join("---" for _ in frame.columns) + " |"
    body = [
        "| " + " | ".join(missing if pd.isna(v) else str(v) for v in row) + " |"
        for row in frame.itertuples(index=False)
    ]
    return "\n".join([header, divider, *body]) + "\n"


def _html_table(frame: pd.DataFrame, missing: str = "NA") -> str:
    """Render a DataFrame as an HTML table with escaped cell contents."""
    if frame.empty:
        return "<p class='muted'>(no rows)</p>"
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in frame.columns)
    rows = []
    for row in frame.itertuples(index=False):
        cells = "".join(
            f"<td>{missing if pd.isna(v) else html.escape(str(v))}</td>" for v in row
        )
        rows.append(f"<tr>{cells}</tr>")
    return (
        "<div class='scroll'><table><thead><tr>"
        + head
        + "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></div>"
    )


def _svg_bars(labels: Sequence[str], values: Sequence[int], width: int = 640) -> str:
    """Render a horizontal bar chart as inline SVG (no plotting library needed)."""
    if not labels:
        return "<p class='muted'>(nothing to plot)</p>"
    row_height, label_width, pad = 26, 150, 8
    height = row_height * len(labels) + pad
    top = max(values) or 1
    bars = []
    for index, (label, value) in enumerate(zip(labels, values)):
        y = index * row_height + pad
        bar = max(1, int((width - label_width - 90) * value / top))
        bars.append(
            f"<text x='0' y='{y + 14}' class='lab'>{html.escape(str(label))}</text>"
            f"<rect x='{label_width}' y='{y + 3}' width='{bar}' height='16' rx='3'/>"
            f"<text x='{label_width + bar + 6}' y='{y + 16}' class='val'>{value:,}</text>"
        )
    return (
        f"<svg viewBox='0 0 {width} {height}' width='100%' height='{height}' "
        f"role='img' class='chart'>" + "".join(bars) + "</svg>"
    )


_CSS = """
/* A deliberately light-only palette. The report is a printable record that gets
   screenshotted into slides and pasted into theses, so it renders the same for
   everyone rather than following the reader's system theme. `color-scheme:light`
   stops the browser from auto-darkening form controls and scrollbars around it. */
:root { color-scheme: light;
        --fg:#1b1f23; --muted:#6a737d; --line:#e1e4e8; --accent:#2f6f4f; --bg:#ffffff;
        --card:#f6f8fa; --stripe:#fbfcfd; --warn-bg:#fff8e6; --warn-line:#e6c65c; }
* { box-sizing: border-box; }
body { margin:0; padding:0 0 4rem; font-family: -apple-system, BlinkMacSystemFont,
       "Segoe UI", Roboto, Helvetica, Arial, sans-serif; color:var(--fg);
       background:var(--bg); line-height:1.55; }
.wrap { max-width: 1040px; margin: 0 auto; padding: 0 1.25rem; }
header { background:var(--accent); color:#fff; padding:2rem 0 1.5rem; margin-bottom:2rem; }
header h1 { margin:0 0 .3rem; font-size:1.7rem; }
header p { margin:0; opacity:.9; font-size:.95rem; }
h2 { margin-top:2.5rem; padding-bottom:.3rem; border-bottom:2px solid var(--line);
     font-size:1.25rem; }
h3 { margin-top:1.6rem; font-size:1.02rem; }
table { border-collapse:collapse; width:100%; font-size:.87rem; }
th, td { border:1px solid var(--line); padding:.35rem .5rem; text-align:left;
         vertical-align:top; }
th { background:var(--card); font-weight:600; position:sticky; top:0; }
tbody tr:nth-child(even) { background:var(--stripe); }
.scroll { overflow-x:auto; }
.muted { color:var(--muted); }
.card { background:var(--card); border:1px solid var(--line); border-radius:6px;
        padding:.85rem 1rem; margin:1rem 0; }
.kpis { display:flex; flex-wrap:wrap; gap:.75rem; margin:1rem 0; }
.kpi { flex:1 1 150px; background:var(--card); border:1px solid var(--line);
       border-radius:6px; padding:.7rem .9rem; }
.kpi .n { font-size:1.5rem; font-weight:700; display:block; }
.kpi .t { font-size:.8rem; color:var(--muted); }
.chart rect { fill:var(--accent); }
.chart .lab { font-size:12px; fill:var(--fg); }
.chart .val { font-size:12px; fill:var(--muted); }
code, .mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
              font-size:.85em; }
a { color:var(--accent); }
nav.toc { background:var(--card); border:1px solid var(--line); border-radius:6px;
          padding:.8rem 1rem; margin:1.5rem 0; }
nav.toc ol { margin:.4rem 0 0; padding-left:1.4rem; columns:2; column-gap:2rem; }
nav.toc li { margin:.15rem 0; font-size:.9rem; }
pre.cmd { background:var(--card); border:1px solid var(--line); border-radius:6px;
          padding:.8rem 1rem; overflow-x:auto; font-size:.82rem; line-height:1.5;
          margin:.6rem 0 1rem; }
pre.cmd .c { color:var(--muted); }
nav.toc strong { font-size:.85rem; text-transform:uppercase; letter-spacing:.04em;
                 color:var(--muted); }
.note { background:var(--warn-bg); border:1px solid var(--warn-line);
        border-left-width:4px; border-radius:4px; padding:.8rem 1rem; margin:1rem 0; }
.note p { margin:.4rem 0; }
.note p:first-child { margin-top:0; }
.note p:last-child { margin-bottom:0; }
.scroll { max-height:32rem; overflow:auto; border:1px solid var(--line);
          border-radius:4px; }
.scroll table { border:0; }
h2 a.anchor { text-decoration:none; color:var(--muted); font-weight:400;
              opacity:0; padding-left:.4rem; }
h2:hover a.anchor { opacity:1; }
.grade-high { color:#1a7f37; font-weight:600; }
.grade-low  { color:#9a6700; }
@media print {
  header { background:none; color:var(--fg); border-bottom:2px solid var(--line); }
  nav.toc, .scroll { break-inside:avoid; max-height:none; overflow:visible; }
  h2 { break-after:avoid; }
  table { font-size:.75rem; }
  a { text-decoration:none; color:var(--fg); }
}
ol.refs li { margin-bottom:.45rem; font-size:.88rem; }
"""


_METHODS = (
    "Six independent annotation tools are harmonised into a single controlled "
    "vocabulary of protein features (NB-ARC, TIR, RPW8, CC, LRR, kinase, LysM, "
    "transmembrane helix, signal peptide). Protein domains come from InterProScan "
    "and are matched by accession only, never by description text, because "
    "descriptions change between releases and match unrelated entries. Overlapping "
    "hits reported by several signature databases for the same region are merged "
    "before anything is counted, so one LRR seen by Pfam, SMART and Gene3D counts "
    "once. Transmembrane helices are taken from Phobius and DeepTMHMM, signal "
    "peptides from SignalP 6.0 and Phobius, and coiled coils from three channels: a "
    "domain-level profile HMM plus the DeepCoil2 and InterProScan Coils predictors. "
    "A helix predicted inside the signal "
    "peptide is discarded, because signal peptides are routinely mistaken for "
    "transmembrane helices. Each protein is then passed through an ordered list of "
    "mutually exclusive rules and receives the first class that fits, together with "
    "a written justification citing the exact signatures behind the call. Subcellular "
    "localisation from DeepLoc 2.0 never decides a class; it only raises or lowers "
    "the reported confidence and flags inconsistencies."
)

_CC_NOTE = (
    "The coiled coil is the least reliable feature in every published RGA pipeline, "
    "and it is the one that decides CNL against NL. Three channels are used here, "
    "and they are not of equal weight. The leading one is a curated profile HMM for "
    "a named domain (the Rx N-terminal domain, PF18052 / IPR041118), which carries "
    "the same kind of evidence as the NB-ARC model every NLR call already rests on. "
    "The other two, DeepCoil2 and InterProScan Coils, are biophysical propensity "
    "predictors: neither publishes a recommended score cut-off, and Simm et al. "
    "(2021), benchmarking coiled-coil predictors against the whole PDB, found a "
    "30-fold spread in how many coiled coils they call and agreement with structure "
    "close to random. They are kept because they cover proteins no domain model "
    "reaches, and a call resting on them alone is graded down rather than hidden. "
    "The tables below show how much the channels disagree and how far the subclass "
    "counts move with the policy."
)


def render_markdown(ctx: Mapping[str, Any], summary: Summary) -> str:
    """Render the run report as Markdown.

    Parameters
    ----------
    ctx : mapping
        Run context assembled by the orchestrator.
    summary : Summary
        The shared count object; no count is recomputed here.

    Returns
    -------
    str
        The complete Markdown document.
    """
    lines = (
        _md_preamble(ctx, summary) + _md_counts(summary) + _md_diagnostics(ctx, summary)
    )
    lines.extend(
        f"{index}. {text} doi:{doi}"
        for index, (text, doi) in enumerate(REFERENCES, start=1)
    )
    lines.append("")
    return "\n".join(lines)


def _md_preamble(ctx: Mapping[str, Any], summary: Summary) -> list[str]:
    """Sections 1-4 of the Markdown report: what, how, metadata, rules."""
    return [
        f"# RGA prediction report -- {ctx['organism']}",
        "",
        f"*Generated {ctx['timestamp']} by `rgas_prediction.py` v{ctx['version']} "
        f"(config v{ctx['config_version']}).*",
        "",
        "## 1. What this report shows",
        "",
        f"{summary.n_proteins:,} proteins were examined and {summary.n_rga:,} "
        f"({100.0 * summary.n_rga / max(summary.n_proteins, 1):.2f}%) carry at least one "
        "feature associated with plant immune receptors. These are *candidates* "
        "identified from protein domains and topology: they are not experimentally "
        "validated resistance genes.",
        "",
        "## 2. How the call was made",
        "",
        _METHODS,
        "",
        "## 3. Run metadata",
        "",
        "### Reproduce this run",
        "",
        "The exact command that produced this report, quoted as it was invoked, so it "
        "can be pasted back into a shell from the repository root:",
        "",
        "```bash",
        _command_lines(ctx["command"]),
        "```",
        "",
        "### Settings",
        "",
        f"- Output directory: `{ctx['outdir']}`",
        f"- Consensus policies: TM `{ctx['options']['policies']['tm']}`, "
        f"SP `{ctx['options']['policies']['sp']}`, CC `{ctx['options']['policies']['cc']}`",
        f"- Coiled-coil calling: threshold {ctx['options']['cc_threshold']}, "
        f"minimum length {ctx['options']['cc_min_length']} residues, "
        f"maximum gap {ctx['options']['cc_max_gap']} residues",
        f"- Minimum LRR copies: {ctx['options']['min_lrr_copies']}",
        "",
        "### Input files",
        "",
        _md_table(ctx["inputs"]),
        "### Evidence channels",
        "",
        _md_table(ctx["availability"]),
        "## 4. Rules applied",
        "",
        _md_table(ctx["rules"]),
    ]


def _md_counts(summary: Summary) -> list[str]:
    """Section 5 of the Markdown report: the count tables."""
    return [
        "## 5. Counts",
        "",
        "### By family",
        "",
        _md_table(summary.family_counts),
        "### By subclass",
        "",
        _md_table(summary.subclass_counts),
        "### Confidence of RGA calls",
        "",
        _md_table(summary.confidence_counts),
        "",
        "### Confidence by subclass",
        "",
        _md_table(summary.subclass_confidence),
        "### Most frequent domain architectures among RGAs",
        "",
        _md_table(summary.architecture_counts),
    ]


def _md_diagnostics(ctx: Mapping[str, Any], summary: Summary) -> list[str]:
    """Sections 6-10 of the Markdown report: CC evidence, IDs, warnings, top RGAs."""
    return [
        "## 6. Coiled-coil evidence",
        "",
        _CC_NOTE,
        "",
        "### DeepCoil2 versus InterProScan Coils (whole proteome)",
        "",
        _md_table(ctx["cc_contingency_table"]),
        "### Subclass counts under each `--cc-policy`",
        "",
        _md_table(ctx["cc_policy_counts"]),
        "### Sensitivity to the segment-calling parameters",
        "",
        _md_table(ctx["cc_sensitivity"]),
        "## 7. Identifier reconciliation",
        "",
        _md_table(ctx["id_report"]),
        "## 8. Warnings",
        "",
        _md_table(summary.warning_counts),
        f"## 9. Top {len(ctx['top_rgas'])} RGA candidates",
        "",
        _md_table(ctx["top_rgas"]),
        "## 10. References",
        "",
    ]


#: HTML section anchors and titles, in page order. One source for the headings
#: and the table of contents, so the two can never drift apart.
_SECTIONS: tuple[tuple[str, str], ...] = (
    ("what", "What this report shows"),
    ("how", "How the call was made"),
    ("metadata", "Run metadata"),
    ("rules", "Rules applied"),
    ("counts", "Counts"),
    ("confidence", "How much of each class is trustworthy"),
    ("cc", "Coiled-coil evidence"),
    ("ids", "Identifier reconciliation"),
    ("warnings", "Warnings"),
    ("top", "Top RGA candidates"),
    ("refs", "References"),
)


def _section(anchor: str) -> str:
    """Open one numbered ``<h2>`` with a stable anchor and a self link."""
    number = [a for a, _ in _SECTIONS].index(anchor) + 1
    title = dict(_SECTIONS)[anchor]
    return (
        f"<h2 id='{anchor}'>{number}. {html.escape(title)}"
        f"<a class='anchor' href='#{anchor}' aria-label='link to this section'>#</a></h2>"
    )


def _toc() -> str:
    """Render the table of contents from :data:`_SECTIONS`."""
    items = "".join(
        f"<li><a href='#{anchor}'>{html.escape(title)}</a></li>"
        for anchor, title in _SECTIONS
    )
    return f"<nav class='toc'><strong>Contents</strong><ol>{items}</ol></nav>"


def render_html(ctx: Mapping[str, Any], summary: Summary) -> str:
    """Render the run report as a self-contained HTML page.

    The page inlines its own CSS and draws its charts with inline SVG: it loads
    no external stylesheet, script, font or image, so it can be opened from a
    file system or emailed as a single attachment.
    """
    body = (
        _html_preamble(ctx, summary)
        + _html_counts(summary)
        + _html_diagnostics(ctx, summary)
    )
    return (
        "<!DOCTYPE html>\n<html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        f"<title>RGA prediction &mdash; {html.escape(ctx['organism'])}</title>"
        f"<style>{_CSS}</style></head><body>" + "".join(body) + "</body></html>"
    )


def _html_preamble(ctx: Mapping[str, Any], summary: Summary) -> list[str]:
    """Header, headline numbers, methodology and run metadata."""
    percent = 100.0 * summary.n_rga / max(summary.n_proteins, 1)
    return [
        "<header><div class='wrap'>",
        f"<h1>Resistance Gene Analog prediction &mdash; {html.escape(ctx['organism'])}</h1>",
        f"<p>Generated {html.escape(ctx['timestamp'])} by rgas_prediction.py "
        f"v{html.escape(ctx['version'])} (config v{html.escape(str(ctx['config_version']))})</p>",
        "</div></header><div class='wrap'>",
        "<div class='kpis'>",
        _kpi(f"{summary.n_proteins:,}", "proteins examined"),
        _kpi(f"{summary.n_rga:,}", "RGA candidates"),
        _kpi(f"{percent:.2f}%", "of the proteome"),
        _kpi(str(len(ctx["rules"])), "rules applied"),
        "</div>",
        _toc(),
        _section("what"),
        "<p>Every protein in the proteome was checked for the protein domains and "
        "membrane topology that characterise plant immune receptors. A protein listed "
        "here as an RGA is a <em>candidate</em> identified from sequence features "
        "alone; it is not an experimentally validated resistance gene.</p>",
        "<div class='note'>"
        "<p><strong>Before quoting a number from this report.</strong></p>"
        "<p>Counts are per <em>protein</em>. In a polyploid genome the same gene "
        "appears several times, so treat them as upper bounds and use "
        "<code>rga_predictions_by_locus.tsv</code> for a locus-level view.</p>"
        "<p>The classes that depend on a coiled coil &mdash; <code>CNL</code>, "
        "<code>CN</code> and <code>TM-CC</code> &mdash; move with the "
        "<code>--cc-policy</code> setting, so state the policy alongside the number. "
        f"This run used <code>{html.escape(ctx['options']['policies']['cc'])}</code>; "
        "section 7 shows what the other settings would have given.</p>"
        "<p><code>TM-CC</code> is a screening bucket rather than a class: it is "
        "defined by the two least specific features, and section 6 shows how much of "
        "it the grading already flags as weak.</p>"
        "</div>",
        _section("how"),
        f"<p>{html.escape(_METHODS)}</p>",
        _section("metadata"),
        "<h3>Reproduce this run</h3>",
        "<p>The exact command that produced this report. It is quoted as it was "
        "actually invoked, so it can be pasted back into a shell from the repository "
        "root; the input checksums in <code>run_metadata.json</code> say which files it "
        "read.</p>",
        _command_block(ctx["command"]),
        "<h3>Settings</h3>",
        "<div class='card mono'>",
        f"<div>output directory: {html.escape(str(ctx['outdir']))}</div>",
        f"<div>policies: TM={html.escape(ctx['options']['policies']['tm'])}, "
        f"SP={html.escape(ctx['options']['policies']['sp'])}, "
        f"CC={html.escape(ctx['options']['policies']['cc'])}</div>",
        f"<div>coiled coil: threshold {ctx['options']['cc_threshold']}, "
        f"min length {ctx['options']['cc_min_length']}, "
        f"max gap {ctx['options']['cc_max_gap']}</div>",
        f"<div>min LRR copies: {ctx['options']['min_lrr_copies']}</div>",
        "</div>",
        "<h3>Input files</h3>",
        _html_table(ctx["inputs"]),
        "<h3>Evidence channels</h3>",
        _html_table(ctx["availability"]),
        _section("rules"),
        _html_table(ctx["rules"]),
    ]


def _html_counts(summary: Summary) -> list[str]:
    """The count tables and their inline-SVG bar charts."""
    families = summary.family_counts
    subclasses = summary.subclass_counts[
        summary.subclass_counts["rga_family"] != "Non-RGA"
    ]
    return [
        _section("counts"),
        "<h3>By family</h3>",
        _svg_bars(list(families["rga_family"]), list(families["n_proteins"])),
        _html_table(families),
        "<h3>By subclass (RGA families only)</h3>",
        _svg_bars(list(subclasses["rga_subclass"]), list(subclasses["n_proteins"])),
        _html_table(summary.subclass_counts),
        "<h3>Most frequent domain architectures among RGAs</h3>",
        _html_table(summary.architecture_counts),
        _section("confidence"),
        "<p>Every call starts at <code>high</code> and is demoted once per triggered "
        "caveat. A large class that is mostly <code>low</code> is a screening result, "
        "not a finding; <code>confidence_demotions</code> in "
        "<code>rga_predictions.tsv</code> names the caveats that fired for each "
        "protein.</p>",
        "<h3>Overall</h3>",
        _svg_bars(
            list(summary.confidence_counts["confidence"]),
            list(summary.confidence_counts["n_proteins"]),
        ),
        "<h3>By subclass</h3>",
        _html_table(summary.subclass_confidence),
    ]


def _html_diagnostics(ctx: Mapping[str, Any], summary: Summary) -> list[str]:
    """Coiled-coil diagnostics, identifier reconciliation, warnings, references."""
    return [
        _section("cc"),
        f"<p>{html.escape(_CC_NOTE)}</p>",
        "<h3>Channel agreement (whole proteome)</h3>",
        _html_table(ctx["cc_contingency_table"]),
        "<h3>Subclass counts under each --cc-policy</h3>",
        _html_table(ctx["cc_policy_counts"]),
        "<h3>Sensitivity to the segment-calling parameters</h3>",
        _html_table(ctx["cc_sensitivity"]),
        _section("ids"),
        _html_table(ctx["id_report"]),
        _section("warnings"),
        _html_table(summary.warning_counts),
        _section("top"),
        f"<p class='muted'>The {len(ctx['top_rgas'])} highest-confidence, "
        "longest RGA candidates. The complete table is "
        "<code>rga_predictions.tsv</code>.</p>",
        _html_table(ctx["top_rgas"]),
        _section("refs"),
        "<ol class='refs'>",
        *(
            f"<li>{html.escape(text)} "
            f"<a href='https://doi.org/{html.escape(doi)}'>doi:{html.escape(doi)}</a></li>"
            for text, doi in REFERENCES
        ),
        "</ol>",
        "</div>",
    ]


def _command_lines(command: str) -> str:
    """Wrap a command one flag per line, for a fenced code block."""
    parts = shlex.split(command)
    lines: list[str] = []
    current = "uv run python"
    for part in parts:
        if part.startswith("--"):
            lines.append(current)
            current = "    " + part
        else:
            current = f"{current} {shlex.quote(part)}"
    lines.append(current)
    return " \\\n".join(lines)


def _command_block(command: str) -> str:
    """Render a command as a copy-pasteable block, wrapped at each flag.

    One flag per line with a trailing backslash: long enough to be unreadable on
    one line, and this is the form a reader is meant to paste back.
    """
    return f"<pre class='cmd'>{html.escape(_command_lines(command))}</pre>"


def _kpi(number: str, text: str) -> str:
    """Render one headline number for the top of the HTML report."""
    return (
        f"<div class='kpi'><span class='n'>{html.escape(number)}</span>"
        f"<span class='t'>{html.escape(text)}</span></div>"
    )
