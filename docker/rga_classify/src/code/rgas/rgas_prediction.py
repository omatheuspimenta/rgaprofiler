"""Genome-wide prediction and classification of Resistance Gene Analogs (RGAs).

Command-line entry point. The pipeline reads pre-computed protein-annotation
outputs (InterProScan, Phobius, DeepTMHMM, SignalP 6.0, DeepLoc 2.0, DeepCoil2),
harmonises them into a single feature vocabulary, applies an ordered list of
mutually exclusive classification rules, and writes machine-readable tables, a
human-readable report and a reproducibility record.

Nothing in this file is specific to sugarcane or to any particular tool version:
every accession, threshold and rule lives in ``code/rgas/config/rga_config.yaml``.

Examples
--------
Run on the bundled R570 data with the default configuration::

    uv run python code/rgas_prediction.py \\
        --input-dir data/rgas --organism-name "Saccharum R570"
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from contextlib import contextmanager
from random import choice
from rich_argparse import RichHelpFormatter

from rich.console import Console
from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
import platform
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from rga import (
    __version__,
    evidence as evidence_mod,
    parsers,
    report as report_mod,
    rules,
)
from rga.config import Config, load_config, normalize_id
from rga.progress import ProgressCallback, null_progress

#: Tools the pipeline can read; only ``interproscan`` is mandatory.
TOOLS: tuple[str, ...] = (
    "interproscan",
    "phobius",
    "deeptmhmm",
    "signalp",
    "deeploc",
    "deepcoil",
)

# ---------------------------------------------------------------------------
# console and logging
# ---------------------------------------------------------------------------

#: One console shared by the log handler and the progress bar. Sharing it is
#: what lets rich print log lines *above* a live progress bar instead of
#: scribbling over it.
CONSOLE = Console()

#: ``markup`` is deliberately off: log messages carry bracketed Python lists of
#: accessions (``['PF00931', ...]``) that rich would otherwise parse as style
#: tags and silently delete.
_RICH_HANDLER_KWARGS: dict[str, Any] = {
    "console": CONSOLE,
    "show_time": False,
    "show_path": False,
    "markup": False,
    "rich_tracebacks": True,
}

logging.basicConfig(
    format="%(message)s",
    level="INFO",
    handlers=[RichHandler(**_RICH_HANDLER_KWARGS)],
)

LOGGER = logging.getLogger(__name__)


def _make_progress() -> Progress:
    """Build the progress bar used for the long parsing stages.

    Returns
    -------
    rich.progress.Progress
        A transient progress bar bound to :data:`CONSOLE`. The spinner glyph is
        chosen at random purely for cosmetics; it has no effect on the run and
        is not part of the reproducibility record.
    """
    spinners = [
        "aesthetic",
        "shark",
        "dots",
        "line",
        "bouncingBall",
        "moon",
        "earth",
        "monkey",
        "runner",
        "pong",
        "weather",
        "clock",
    ]
    return Progress(
        SpinnerColumn(choice(spinners)),
        TaskProgressColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
        console=CONSOLE,
        transient=True,
    )


@contextmanager
def tracked(description: str, total: int | None = None):
    """Run a stage under a transient progress bar.

    Yields the :class:`~rga.progress.ProgressCallback` that the stage should
    report to. This is the *only* place where ``rich`` meets the pipeline
    stages: everything in ``code/rga/`` reports through the plain callback and
    stays importable without a console.

    Parameters
    ----------
    description : str
        Text shown beside the bar.
    total : int, optional
        Units of work, when the caller already knows it. A stage that discovers
        its own total later announces it through the callback.

    Yields
    ------
    ProgressCallback
        Pass this as the stage's ``on_progress`` argument.

    Examples
    --------
    >>> with tracked("classifying", total=len(records)) as report:  # doctest: +SKIP
    ...     predictions = rules.classify_proteome(cfg, ev, options, report)
    """
    with _make_progress() as progress:
        task = progress.add_task(description, total=total)

        def report(advance: int = 0, total: int | None = None) -> None:
            """Forward a stage's progress report to the rich task."""
            if total is not None:
                progress.update(task, total=total)
            if advance:
                progress.advance(task, advance)

        yield report


# CLI
def build_parser() -> argparse.ArgumentParser:
    """Construct the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="rgas_prediction.py",
        description=(
            "[bold cyan]Genome-wide prediction and classification "
            "of Resistance Gene Analogs[/bold cyan]\n\n"
            "[dim]Combines predictions from InterProScan, Phobius, "
            "DeepTMHMM, SignalP, DeepLoc, and DeepCoil2.[/dim]"
        ),
        formatter_class=RichHelpFormatter,
        epilog=(
            "[bold]Examples:[/bold]\n"
            "  [cyan]rgas_prediction.py[/cyan] "
            "[green]--input-dir[/green] results "
            "[green]--outdir[/green] predictions\n"
            "  [cyan]rgas_prediction.py[/cyan] "
            "[green]--input-dir[/green] results "
            "[green]--rga-only[/green]\n"
        ),
    )

    # Input / output
    io_group = parser.add_argument_group(
        "Input and output"
    )

    io_group.add_argument(
        "--input-dir",
        type=Path,
        metavar="DIR",
        help="directory containing the tool outputs",
    )
    io_group.add_argument(
        "--outdir",
        type=Path,
        metavar="DIR",
        help="directory for generated output files",
    )
    io_group.add_argument(
        "--config",
        type=Path,
        metavar="FILE",
        default=Path(__file__).resolve().parent
        / "config"
        / "rga_config.yaml",
        help="pipeline configuration file",
    )
    io_group.add_argument(
        "--organism-name",
        metavar="NAME",
        default="unnamed organism",
        help="organism label used in reports",
    )

    # Per-tool inputs
    tool_group = parser.add_argument_group(
        "Tool inputs",
        description=(
            "[dim]Explicit paths override the corresponding files "
            "found under --input-dir.[/dim]"
        ),
    )

    tool_group.add_argument(
        "--interproscan",
        type=Path,
        metavar="TSV",
        help="[bold]InterProScan[/bold] TSV [yellow](required)[/yellow]",
    )
    tool_group.add_argument(
        "--phobius",
        type=Path,
        metavar="FILE",
        help="[bold]Phobius[/bold] short-format output",
    )
    tool_group.add_argument(
        "--deeptmhmm",
        type=Path,
        metavar="GFF3",
        help="[bold]DeepTMHMM[/bold] TMRs.gff3",
    )
    tool_group.add_argument(
        "--signalp",
        type=Path,
        metavar="FILE",
        help="[bold]SignalP 6.0[/bold] prediction_results.txt",
    )
    tool_group.add_argument(
        "--deeploc",
        type=Path,
        metavar="CSV",
        help="[bold]DeepLoc 2.0[/bold] results CSV",
    )
    tool_group.add_argument(
        "--deepcoil",
        type=Path,
        metavar="PATH",
        help="[bold]DeepCoil2[/bold] directory or .tar.xz archive",
    )

    
    # Consensus
    consensus_group = parser.add_argument_group(
        "Consensus policies and thresholds"
    )

    consensus_group.add_argument(
        "--tm-policy",
        metavar="POLICY",
        choices=["union", "intersection", "deeptmhmm", "phobius"],
        help="transmembrane prediction policy",
    )
    consensus_group.add_argument(
        "--sp-policy",
        metavar="POLICY",
        choices=["union", "intersection", "signalp", "phobius"],
        help="signal peptide prediction policy",
    )
    consensus_group.add_argument(
        "--cc-policy",
        metavar="POLICY",
        choices=["rx_domain", "deepcoil", "coils", "union", "intersection"],
        help="coiled-coil prediction policy",
    )
    consensus_group.add_argument(
        "--cc-threshold",
        type=float,
        metavar="SCORE",
        help="DeepCoil2 plateau-score cut-off",
    )
    consensus_group.add_argument(
        "--cc-min-length",
        type=int,
        metavar="RESIDUES",
        help="minimum coiled-coil segment length",
    )
    consensus_group.add_argument(
        "--cc-max-gap",
        type=int,
        metavar="RESIDUES",
        help="maximum gap merged between coiled-coil segments",
    )
    consensus_group.add_argument(
        "--min-lrr-copies",
        type=int,
        metavar="N",
        help="minimum number of merged LRR hits",
    )

    
    # Output filtering
    output_group = parser.add_argument_group(
        "Output filtering"
    )

    output = output_group.add_mutually_exclusive_group()

    output.add_argument(
        "--keep-non-rga",
        dest="keep_non_rga",
        action="store_true",
        default=True,
        help="write all proteins to rga_predictions.tsv [dim](default)[/dim]",
    )
    output.add_argument(
        "--rga-only",
        dest="keep_non_rga",
        action="store_false",
        help="write only predicted RGA candidates",
    )

    
    # Performance / diagnostics
    runtime_group = parser.add_argument_group(
        "Runtime and diagnostics"
    )

    runtime_group.add_argument(
        "--workers",
        type=int,
        metavar="N",
        default=4,
        help="number of processes used to read DeepCoil2",
    )
    runtime_group.add_argument(
        "--refresh-deepcoil-cache",
        action="store_true",
        help="re-parse DeepCoil2 even when a cached segment table exists",
    )
    runtime_group.add_argument(
        "--log-level",
        metavar="LEVEL",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="logging verbosity",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def setup_logging(level: str, logfile: Path) -> None:
    """Send logs to the console (rich) and to ``<outdir>/logs/run.log`` (plain).

    Parameters
    ----------
    level : str
        Logging level name, e.g. ``"INFO"``.
    logfile : pathlib.Path
        Destination of the plain-text run log; parent directories are created.

    Notes
    -----
    The console handler is a :class:`~rich.logging.RichHandler` bound to
    :data:`CONSOLE`; the file handler is a plain formatter so that ``run.log``
    stays greppable and free of ANSI escapes.
    """
    logfile.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(level)
    root.handlers.clear()

    console_handler = RichHandler(**_RICH_HANDLER_KWARGS)
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    root.addHandler(console_handler)

    file_handler = logging.FileHandler(logfile, mode="w", encoding="utf-8")
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)-8s %(name)s: %(message)s")
    )
    root.addHandler(file_handler)


# ---------------------------------------------------------------------------
# input discovery
# ---------------------------------------------------------------------------


def discover_inputs(args: argparse.Namespace, cfg: Config) -> dict[str, Path | None]:
    """Locate each tool's output, honouring explicit overrides first.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    cfg : Config
        Resolved configuration providing the ``input_discovery`` glob patterns.

    Returns
    -------
    dict
        Tool name -> path, or ``None`` when the tool output was not found.
    """
    found: dict[str, Path | None] = {}
    patterns = cfg.raw.get("input_discovery", {})
    for tool in TOOLS:
        override = getattr(args, tool, None)
        if override is not None:
            if not Path(override).exists():
                raise SystemExit(f"--{tool}: path does not exist: {override}")
            found[tool] = Path(override)
            continue
        found[tool] = _glob_first(args.input_dir, patterns.get(tool, []))
    return found


def _glob_first(root: Path | None, patterns: list[str]) -> Path | None:
    """Return the first existing path matching any pattern, deterministically."""
    if root is None:
        return None
    for pattern in patterns:
        matches = sorted(root.glob(pattern))
        if matches:
            return matches[0]
    return None


def file_fingerprint(path: Path) -> dict[str, Any]:
    """Describe an input file: size, line count and SHA-256 checksum."""
    digest = hashlib.sha256()
    n_lines = 0
    if path.is_dir():
        members = sorted(p for p in path.rglob("*") if p.is_file())
        for member in members:
            digest.update(member.name.encode("utf-8"))
        return {
            "size_bytes": sum(p.stat().st_size for p in members),
            "n_lines": len(members),
            "sha256": digest.hexdigest(),
        }
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
            n_lines += block.count(b"\n")
    return {
        "size_bytes": path.stat().st_size,
        "n_lines": n_lines,
        "sha256": digest.hexdigest(),
    }


# ---------------------------------------------------------------------------
# identifier reconciliation
# ---------------------------------------------------------------------------


def build_deepcoil_map(
    deepcoil_ids: set[str], canonical_ids: set[str], cfg: Config
) -> tuple[dict[str, str], list[str]]:
    """Map DeepCoil2 file stems back onto canonical protein IDs.

    DeepCoil2 names its output files after a sanitised protein ID (for the R570
    proteome, the dots are stripped). The mapping is rebuilt by applying the
    same sanitisation to the canonical IDs. The transformation is only usable if
    it is injective over the proteome, which is asserted here rather than
    assumed.

    Returns
    -------
    tuple
        ``(mapping, unmatched)`` where ``mapping`` is DeepCoil stem -> canonical
        ID and ``unmatched`` lists the DeepCoil stems that could not be mapped.

    Raises
    ------
    SystemExit
        If the sanitisation collapses two distinct proteins onto one name.
    """
    form = cfg.raw["ids"].get("deepcoil_canonical_form")
    if form is None:
        direct = {pid: pid for pid in deepcoil_ids & canonical_ids}
        return direct, sorted(deepcoil_ids - canonical_ids)

    lookup: dict[str, str] = {}
    collisions: list[str] = []
    for protein_id in canonical_ids:
        key = normalize_id(protein_id, [form])
        if key in lookup:
            collisions.append(key)
        lookup[key] = protein_id
    if collisions:
        raise SystemExit(
            f"ID normalisation {form!r} is not injective over this proteome "
            f"({len(collisions)} collisions, e.g. {collisions[:3]}); "
            "adjust ids.deepcoil_canonical_form in the configuration"
        )
    mapping = {stem: lookup[stem] for stem in deepcoil_ids if stem in lookup}
    return mapping, sorted(deepcoil_ids - set(mapping))


def reconcile_ids(
    id_sets: dict[str, set[str]], canonical: set[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare each tool's ID set with the canonical proteome.

    Returns
    -------
    tuple
        ``(summary, detail)``. ``summary`` has one row per tool with the counts
        used in the report; ``detail`` lists every unmatched ID with a reason and
        becomes ``unmatched_ids_report.tsv``.
    """
    summary_rows = []
    detail_rows = []
    for tool in sorted(id_sets):
        ids = id_sets[tool]
        missing = sorted(canonical - ids)
        extra = sorted(ids - canonical)
        summary_rows.append(
            {
                "tool": tool,
                "n_ids": len(ids),
                "n_shared_with_proteome": len(ids & canonical),
                "n_absent_from_tool": len(missing),
                "n_not_in_proteome": len(extra),
            }
        )
        detail_rows.extend(
            {
                "tool": tool,
                "protein_id": protein_id,
                "reason": "present in the proteome but absent from this tool's output",
            }
            for protein_id in missing
        )
        detail_rows.extend(
            {
                "tool": tool,
                "protein_id": protein_id,
                "reason": "reported by this tool but absent from the proteome",
            }
            for protein_id in extra
        )
    detail = pd.DataFrame(detail_rows, columns=["tool", "protein_id", "reason"])
    return pd.DataFrame(summary_rows), detail


# ---------------------------------------------------------------------------
# invariants
# ---------------------------------------------------------------------------


def assert_invariants(
    predictions: pd.DataFrame, canonical: set[str], summary: report_mod.Summary
) -> None:
    """Fail loudly if any structural guarantee of the pipeline is violated.

    Checks that every input protein appears exactly once, that every protein
    carries exactly one rule, that the subclass counts add up to the proteome
    size, and that the class definitions were respected (no NLR without NB-ARC,
    no RLK without a kinase feature).
    """
    n_expected = len(canonical)
    if len(predictions) != n_expected:
        raise AssertionError(
            f"output has {len(predictions)} rows for {n_expected} input proteins"
        )
    if predictions["protein_id"].duplicated().any():
        duplicated = predictions.loc[
            predictions["protein_id"].duplicated(), "protein_id"
        ]
        raise AssertionError(
            f"duplicated protein IDs in output: {duplicated.head().tolist()}"
        )
    if set(predictions["protein_id"]) != canonical:
        raise AssertionError("output protein IDs do not match the input proteome")
    if predictions["rule_id"].isna().any():
        raise AssertionError("some proteins carry no rule_id")

    total = int(summary.subclass_counts["n_proteins"].sum())
    if total != n_expected:
        raise AssertionError(f"subclass counts sum to {total}, expected {n_expected}")
    if int(summary.family_counts["n_proteins"].sum()) != n_expected:
        raise AssertionError("family counts do not sum to the proteome size")

    nlr = predictions[predictions["rga_family"] == "NLR"]
    if not nlr.empty and not nlr["features_found"].str.contains("NB-ARC").all():
        raise AssertionError("an NLR was classified without an NB-ARC feature")
    rlk = predictions[predictions["rga_family"] == "RLK"]
    if not rlk.empty and not rlk["features_found"].str.contains("STTK").all():
        raise AssertionError("an RLK was classified without a kinase feature")
    LOGGER.info("All pipeline invariants hold (%d proteins).", n_expected)


def assert_report_consistency(
    summary: report_mod.Summary, summary_tsv: pd.DataFrame, metadata: dict
) -> None:
    """Check that the TSV, the JSON metadata and the report agree on every count."""
    tsv_total = int(
        summary_tsv.loc[summary_tsv["level"] == "subclass", "n_proteins"].sum()
    )
    if tsv_total != summary.n_proteins:
        raise AssertionError(
            f"rga_summary_counts.tsv totals {tsv_total}, summary says {summary.n_proteins}"
        )
    if metadata["counts"]["n_proteins"] != summary.n_proteins:
        raise AssertionError(
            "run_metadata.json disagrees with the summary on n_proteins"
        )
    if metadata["counts"]["n_rga"] != summary.n_rga:
        raise AssertionError("run_metadata.json disagrees with the summary on n_rga")


# ---------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------


def write_tsv(frame: pd.DataFrame, path: Path, cfg: Config) -> None:
    """Write a UTF-8, tab-separated table using the configured missing value."""
    path.parent.mkdir(parents=True, exist_ok=True)
    cleaned = frame.replace({"": None, np.nan: None})
    cleaned.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep=cfg.missing_value,
        encoding="utf-8",
        lineterminator="\n",
    )
    LOGGER.info("wrote %s (%d rows)", path, len(frame))


def locus_summary(predictions: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Collapse isoforms onto their locus.

    The representative isoform is the longest protein of the locus; ties are
    broken by protein ID so the choice is deterministic. ``isoforms_disagree``
    flags loci whose isoforms did not all receive the same subclass, which is a
    direct readout of how much fragmented or alternative gene models perturb the
    counts.

    Notes
    -----
    Implemented with a single pass over plain Python tuples rather than
    ``groupby(...).apply(...)``: on a 300k-row table with ~195k groups the pandas
    route allocates enough intermediate objects to exhaust memory on a modest
    machine, while this loop is linear and allocation-free.

    Parameters
    ----------
    predictions : pandas.DataFrame
        The per-protein prediction table.
    cfg : Config
        Resolved configuration, used for the list separator.

    Returns
    -------
    pandas.DataFrame
        One row per locus, sorted by locus identifier.
    """
    separator = cfg.list_separator
    best: dict[str, tuple[int, str, str, str, str]] = {}
    n_isoforms: dict[str, int] = {}
    n_rga: dict[str, int] = {}
    subclasses: dict[str, set[str]] = {}

    lengths = pd.to_numeric(predictions["sequence_length"], errors="coerce").fillna(0)
    for locus, protein_id, length, is_rga, family, subclass, confidence in zip(
        predictions["locus"],
        predictions["protein_id"],
        lengths,
        predictions["is_rga"],
        predictions["rga_family"],
        predictions["rga_subclass"],
        predictions["confidence"],
    ):
        if not isinstance(locus, str) or not locus:
            continue
        n_isoforms[locus] = n_isoforms.get(locus, 0) + 1
        n_rga[locus] = n_rga.get(locus, 0) + int(bool(is_rga))
        subclasses.setdefault(locus, set()).add(str(subclass))
        candidate = (int(length), protein_id, family, subclass, confidence)
        current = best.get(locus)
        # longest protein wins; ties broken by the smaller protein ID
        if current is None or (-candidate[0], candidate[1]) < (-current[0], current[1]):
            best[locus] = candidate

    loci = sorted(best)
    observed = [separator.join(sorted(subclasses[locus])) for locus in loci]
    return pd.DataFrame(
        {
            "locus": loci,
            "n_isoforms": [n_isoforms[locus] for locus in loci],
            "n_isoforms_rga": [n_rga[locus] for locus in loci],
            "representative_protein_id": [best[locus][1] for locus in loci],
            "rga_family": [best[locus][2] for locus in loci],
            "rga_subclass": [best[locus][3] for locus in loci],
            "subclasses_observed": observed,
            "isoforms_disagree": [len(subclasses[locus]) > 1 for locus in loci],
            "confidence": [best[locus][4] for locus in loci],
        }
    )


def top_rga_table(predictions: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Select the RGA candidates shown in the report, best evidence first."""
    rgas = predictions[predictions["is_rga"]].copy()
    order = {"high": 0, "medium": 1, "low": 2}
    rgas["_rank"] = rgas["confidence"].map(order).fillna(3)
    rgas = rgas.sort_values(
        ["_rank", "rule_priority", "n_lrr", "protein_id"],
        ascending=[True, True, False, True],
        kind="stable",
    )
    columns = [
        "protein_id",
        "rga_family",
        "rga_subclass",
        "domain_architecture",
        "n_lrr",
        "predicted_localization",
        "confidence",
    ]
    return rgas[columns].head(top_n).reset_index(drop=True)


def accession_audit(cfg: Config, counts) -> pd.DataFrame:
    """Report which configured accessions were actually observed in the data."""
    rows = []
    for feature, accessions in cfg.raw["interproscan_features"].items():
        for accession in accessions:
            hits = int(counts.get(str(accession), 0))
            rows.append(
                {
                    "feature": feature,
                    "accession": accession,
                    "n_hits": hits,
                    "status": "used" if hits else "seeded but unused",
                }
            )
    # The domain-level CC channel is configured outside `interproscan_features`
    # (it is not one of the nine features), so it needs its own pass here --
    # otherwise the accessions that now drive the CNL count would be the only
    # evidence accessions missing from the audit.
    for accession in cfg.cc_domain_accessions():
        hits = int(counts.get(accession, 0))
        rows.append(
            {
                "feature": "CC (domain channel)",
                "accession": accession,
                "n_hits": hits,
                "status": "used" if hits else "seeded but unused",
            }
        )
    for accession, note in cfg.raw.get("watch_accessions", {}).items():
        rows.append(
            {
                "feature": "(watch only)",
                "accession": accession,
                "n_hits": int(counts.get(str(accession), 0)),
                "status": note,
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["feature", "accession"], kind="stable")
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------


def resolve_options(args: argparse.Namespace, cfg: Config) -> dict[str, Any]:
    """Merge configuration defaults with command-line overrides."""
    coiled = cfg.coiled_coil
    policies = dict(cfg.policies)
    for channel in ("tm", "sp", "cc"):
        override = getattr(args, f"{channel}_policy")
        if override:
            policies[channel] = override
    return {
        "policies": policies,
        "cc_threshold": _pick(args.cc_threshold, coiled["threshold"]),
        "cc_min_length": _pick(args.cc_min_length, coiled["min_length"]),
        "cc_max_gap": _pick(args.cc_max_gap, coiled["max_gap"]),
        "cc_tm_overlap": coiled["tm_overlap_fraction"],
        "sp_overlap": cfg.raw["transmembrane"]["sp_overlap_fraction"],
        "min_tm_helices": cfg.raw["transmembrane"]["min_helices"],
        "min_lrr_copies": _pick(
            args.min_lrr_copies, cfg.raw["intervals"]["min_lrr_copies"]
        ),
        "merge_min_overlap": cfg.raw["intervals"]["merge_min_overlap"],
        "lrr_repeat_analyses": cfg.raw["intervals"]["lrr_repeat_analyses"],
        "locus_regex": cfg.raw["ids"].get("locus_regex"),
    }


def _pick(override, default):
    """Return ``override`` when it was given on the command line."""
    return default if override is None else override


def load_deepcoil(
    path: Path | None,
    cfg: Config,
    cache: Path,
    refresh: bool,
    workers: int,
    on_progress: ProgressCallback = null_progress,
) -> pd.DataFrame | None:
    """Load the DeepCoil2 raw segment table, using a cache when possible.

    The cached table holds *unfiltered* segments, so changing the threshold, the
    minimum length or the gap parameter never requires re-reading the archives.
    A cache hit finishes the progress bar immediately: there is genuinely no
    work to show.
    """
    if path is None:
        on_progress(0, total=1)
        on_progress(1)
        return None
    if cache.exists() and not refresh:
        LOGGER.info("DeepCoil2: reusing cached segment table %s", cache)
        on_progress(0, total=1)
        frame = pd.read_csv(cache, sep="\t", dtype={"deepcoil_id": str})
        on_progress(1)
        return frame
    segments = parsers.parse_deepcoil(
        path, cfg, workers=workers, on_progress=on_progress
    )
    cache.parent.mkdir(parents=True, exist_ok=True)
    segments.to_csv(cache, sep="\t", index=False)
    LOGGER.info("DeepCoil2: cached segment table written to %s", cache)
    return segments


def collect_id_sets(
    parsed: dict[str, pd.DataFrame | None], ips, deepcoil_path: Path | None, cfg: Config
) -> dict[str, set[str]]:
    """Gather the protein-ID set reported by every available tool."""
    id_sets: dict[str, set[str]] = {"interproscan": set(ips.protein_ids)}
    for tool in ("phobius", "deeptmhmm", "signalp", "deeploc"):
        frame = parsed.get(tool)
        if frame is not None:
            id_sets[tool] = set(frame["protein_id"])
    if deepcoil_path is not None:
        id_sets["deepcoil"] = parsers.deepcoil_protein_ids(deepcoil_path, cfg)
    return id_sets


def environment_record() -> dict[str, str]:
    """Capture the interpreter and package versions used for this run."""
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
    }


def run(args: argparse.Namespace) -> int:
    """Execute the whole pipeline and write every output.

    Returns
    -------
    int
        Process exit status: ``0`` on success.
    """
    cfg = load_config(args.config)
    rules.assert_mutually_exclusive(cfg)
    outdir = args.outdir or Path("results/rgas") / _slug(args.organism_name)
    outdir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.log_level, outdir / "logs" / "run.log")
    LOGGER.info(
        "rgas_prediction.py v%s -- organism: %s", __version__, args.organism_name
    )

    inputs = discover_inputs(args, cfg)
    _report_inputs(inputs)

    options = resolve_options(args, cfg)
    LOGGER.info("policies: %s", options["policies"])

    ips, parsed, segments = parse_all_inputs(
        inputs, cfg, outdir / "cache" / "deepcoil_raw_segments.tsv", args
    )

    id_sets = collect_id_sets(parsed, ips, inputs["deepcoil"], cfg)
    canonical = set().union(
        *(ids for tool, ids in id_sets.items() if tool != "deepcoil")
    )
    deepcoil_map, unmapped = build_deepcoil_map(
        id_sets.get("deepcoil", set()), canonical, cfg
    )
    LOGGER.info(
        "IDs: %d canonical proteins; DeepCoil2 mapped %d, unmatched %d",
        len(canonical),
        len(deepcoil_map),
        len(unmapped),
    )
    id_summary, id_detail = reconcile_ids(
        {tool: ids for tool, ids in id_sets.items() if tool != "deepcoil"}, canonical
    )
    id_summary = _append_deepcoil_row(
        id_summary, id_sets, deepcoil_map, unmapped, canonical
    )
    id_detail = _append_deepcoil_detail(
        id_detail, deepcoil_map, unmapped, canonical, id_sets
    )

    protein_ids = sorted(canonical)
    with tracked("building evidence", total=len(protein_ids)) as report:
        ev = evidence_mod.build_evidence(
            cfg,
            protein_ids,
            ips,
            parsed["phobius"],
            parsed["deeptmhmm"],
            parsed["signalp"],
            parsed["deeploc"],
            segments,
            deepcoil_map,
            options,
            on_progress=report,
        )
    with tracked("classifying proteins", total=len(protein_ids)) as report:
        predictions = rules.classify_proteome(cfg, ev, options, on_progress=report)
    summary = report_mod.summarize(predictions, top_n=20)
    assert_invariants(predictions, canonical, summary)

    context = _build_context(
        args, cfg, options, inputs, ev, predictions, id_summary, outdir
    )
    metadata = _build_metadata(args, cfg, options, context, summary, ev)
    summary_tsv = report_mod.summary_table(summary)
    assert_report_consistency(summary, summary_tsv, metadata)

    _write_outputs(
        cfg,
        outdir,
        args,
        predictions,
        ev,
        summary,
        summary_tsv,
        id_detail,
        context,
        metadata,
        ips,
    )
    LOGGER.info(
        "Done. %d proteins, %d RGA candidates.", summary.n_proteins, summary.n_rga
    )
    return 0


def _report_inputs(inputs: dict[str, Path | None]) -> None:
    """Log which tool outputs were found, and fail if the required one is absent."""
    if inputs["interproscan"] is None:
        raise SystemExit(
            "InterProScan output is required; pass --interproscan or --input-dir"
        )
    for tool, path in sorted(inputs.items()):
        if path is None:
            LOGGER.warning(
                "%s output not found: this evidence channel is unavailable and every call "
                "depending on it will be reported at lower confidence",
                tool,
            )
        else:
            LOGGER.info("%s: %s", tool, path)


def parse_all_inputs(
    inputs: dict[str, Path | None], cfg: Config, cache: Path, args: argparse.Namespace
) -> tuple[parsers.InterProScanResult, dict[str, pd.DataFrame | None], pd.DataFrame]:
    """Parse every supplied tool output, reporting progress per stage.

    Parameters
    ----------
    inputs : dict
        Tool name -> path, as returned by :func:`discover_inputs`.
    cfg : Config
        Resolved configuration.
    cache : pathlib.Path
        Destination of the unfiltered DeepCoil2 segment cache.
    args : argparse.Namespace
        Parsed CLI arguments; ``refresh_deepcoil_cache`` and ``workers`` are read.

    Returns
    -------
    tuple
        ``(interproscan_result, optional_tool_frames, deepcoil_segments)``.

    Notes
    -----
    These three stages dominate the wall time (InterProScan is tens of millions
    of rows, DeepCoil2 is one file per protein), so each gets its own bar:
    InterProScan is measured in bytes of the TSV consumed, DeepCoil2 in source
    directories/archives finished, and the remaining tools one unit each.
    """
    with tracked("InterProScan TSV (bytes)") as report:
        ips = parsers.parse_interproscan(inputs["interproscan"], cfg, report)
    with tracked("TM / SP / localisation", total=len(OPTIONAL_TOOLS)) as report:
        parsed = parse_optional_tools(inputs, cfg, report)
    with tracked("DeepCoil2 segments") as report:
        segments = load_deepcoil(
            inputs["deepcoil"],
            cfg,
            cache,
            args.refresh_deepcoil_cache,
            args.workers,
            report,
        )
    return ips, parsed, segments


#: Optional tools, in parsing order, with the parser each one uses.
OPTIONAL_TOOLS: tuple[tuple[str, str], ...] = (
    ("phobius", "parse_phobius"),
    ("deeptmhmm", "parse_deeptmhmm"),
    ("signalp", "parse_signalp"),
    ("deeploc", "parse_deeploc"),
)


def parse_optional_tools(
    inputs: dict[str, Path | None],
    cfg: Config,
    on_progress: ProgressCallback = null_progress,
) -> dict[str, pd.DataFrame | None]:
    """Parse every optional tool that was supplied, returning ``None`` for the rest.

    Parameters
    ----------
    inputs : dict
        Tool name -> path, as returned by :func:`discover_inputs`.
    cfg : Config
        Resolved configuration.
    on_progress : ProgressCallback, optional
        Advanced by one per tool, whether or not that tool was supplied, so the
        bar reaches 100 % on a partial input set too.

    Returns
    -------
    dict
        Tool name -> parsed frame, or ``None`` when the tool was not supplied.
    """
    parsed: dict[str, pd.DataFrame | None] = {}
    for tool, parser_name in OPTIONAL_TOOLS:
        path = inputs[tool]
        parsed[tool] = getattr(parsers, parser_name)(path, cfg) if path else None
        on_progress(1)
    return parsed


def _slug(text: str) -> str:
    """Turn an organism name into a filesystem-safe directory name."""
    return (
        "".join(c if c.isalnum() or c in "-_" else "_" for c in text).strip("_")
        or "organism"
    )


def _append_deepcoil_row(
    id_summary: pd.DataFrame,
    id_sets: dict[str, set[str]],
    mapping: dict[str, str],
    unmapped: list[str],
    canonical: set[str],
) -> pd.DataFrame:
    """Add the DeepCoil2 line to the ID reconciliation summary."""
    if "deepcoil" not in id_sets:
        return id_summary
    row = {
        "tool": "deepcoil",
        "n_ids": len(id_sets["deepcoil"]),
        "n_shared_with_proteome": len(mapping),
        "n_absent_from_tool": len(canonical - set(mapping.values())),
        "n_not_in_proteome": len(unmapped),
    }
    return (
        pd.concat([id_summary, pd.DataFrame([row])], ignore_index=True)
        .sort_values("tool", kind="stable")
        .reset_index(drop=True)
    )


def _append_deepcoil_detail(
    detail: pd.DataFrame,
    mapping: dict[str, str],
    unmapped: list[str],
    canonical: set[str],
    id_sets: dict[str, set[str]],
) -> pd.DataFrame:
    """Add the DeepCoil2 unmatched IDs to ``unmatched_ids_report.tsv``."""
    if "deepcoil" not in id_sets:
        return detail
    rows = [
        {
            "tool": "deepcoil",
            "protein_id": stem,
            "reason": "DeepCoil2 file could not be mapped back to a canonical protein ID",
        }
        for stem in unmapped
    ]
    rows.extend(
        {
            "tool": "deepcoil",
            "protein_id": protein_id,
            "reason": "present in the proteome but absent from this tool's output",
        }
        for protein_id in sorted(canonical - set(mapping.values()))
    )
    if not rows:
        return detail
    return pd.concat([detail, pd.DataFrame(rows)], ignore_index=True)


def _build_context(
    args: argparse.Namespace,
    cfg: Config,
    options: dict,
    inputs: dict[str, Path | None],
    ev,
    predictions: pd.DataFrame,
    id_summary: pd.DataFrame,
    outdir: Path,
) -> dict[str, Any]:
    """Assemble everything the two report renderers need."""
    input_rows = []
    for tool in TOOLS:
        path = inputs[tool]
        record = {
            "tool": tool,
            "path": str(path) if path else None,
            "available": path is not None,
        }
        record.update(
            file_fingerprint(path)
            if path
            else {"size_bytes": None, "n_lines": None, "sha256": None}
        )
        input_rows.append(record)

    contingency = ev.cc_contingency
    # The 2x2 block compares the two propensity predictors; the two trailing
    # rows report the domain-level channel, which is a different kind of
    # evidence and is therefore counted beside them rather than folded in.
    contingency_table = pd.DataFrame(
        {
            "InterProScan Coils": [
                "CC called",
                "no CC",
                "CC called",
                "no CC",
                "--",
                "no CC",
            ],
            "DeepCoil2": [
                "CC called",
                "CC called",
                "no CC",
                "no CC",
                "--",
                "no CC",
            ],
            "Rx domain": ["--", "--", "--", "--", "CC called", "CC called"],
            "n_proteins": [
                contingency["both"],
                contingency["deepcoil_only"],
                contingency["coils_only"],
                contingency["neither"],
                contingency["rx_domain"],
                contingency["rx_domain_only"],
            ],
        }
    )
    coiled = cfg.coiled_coil
    sensitivity = evidence_mod.cc_sensitivity(
        ev.raw_cc,
        coiled["sensitivity_thresholds"],
        coiled["sensitivity_min_lengths"],
        options["cc_max_gap"],
    )
    policy_counts = rules.cc_policy_sensitivity(cfg, ev).reset_index()

    return {
        "organism": args.organism_name,
        "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
        "version": __version__,
        "config_version": cfg.raw["config_version"],
        # shlex.join, not " ".join: an organism name with spaces has to come back
        # out quoted or the recorded command is not the command that ran.
        "command": shlex.join(
            [_SCRIPT_PATH, *getattr(args, "invoked_with", sys.argv[1:])]
        ),
        "outdir": str(outdir),
        "options": options,
        "inputs": pd.DataFrame(input_rows),
        "availability": pd.DataFrame(
            {
                "channel": list(ev.available),
                "available": [ev.available[c] for c in ev.available],
            }
        ),
        "rules": _rule_table(cfg),
        "cc_contingency": contingency,
        "cc_contingency_table": contingency_table,
        "cc_policy_counts": policy_counts,
        "cc_sensitivity": sensitivity,
        "id_report": id_summary,
        "top_rgas": top_rga_table(
            predictions, cfg.raw["output"]["top_n_rgas_in_report"]
        ),
    }


def _rule_table(cfg: Config) -> pd.DataFrame:
    """Render the ordered rule list exactly as it was applied."""
    separator = cfg.list_separator
    return pd.DataFrame(
        [
            {
                "priority": rule.priority,
                "rule_id": rule.id,
                "family": rule.family,
                "subclass": rule.subclass,
                "requires": separator.join(rule.all_of) or "-",
                "requires_one_of": (
                    " AND ".join("(" + " OR ".join(g) + ")" for g in rule.any_of) or "-"
                ),
                "forbids": separator.join(rule.none_of) or "-",
                "description": rule.description,
            }
            for rule in cfg.rules
        ]
    )


def _build_metadata(
    args: argparse.Namespace,
    cfg: Config,
    options: dict,
    context: dict,
    summary: report_mod.Summary,
    ev,
) -> dict[str, Any]:
    """Assemble ``run_metadata.json``."""
    return {
        "script_version": __version__,
        "config_version": cfg.raw["config_version"],
        "timestamp_utc": context["timestamp"],
        "organism": args.organism_name,
        "command": context["command"],
        "cli_args": {
            k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()
        },
        "resolved_options": options,
        "resolved_config": cfg.raw,
        "inputs": context["inputs"].to_dict(orient="records"),
        "evidence_channels_available": ev.available,
        "environment": environment_record(),
        "counts": summary.as_dict(),
        "cc_contingency_deepcoil_vs_coils": ev.cc_contingency,
        "cc_policy_subclass_counts": context["cc_policy_counts"].to_dict(
            orient="records"
        ),
        "cc_segment_sensitivity": context["cc_sensitivity"].to_dict(orient="records"),
        "id_reconciliation": context["id_report"].to_dict(orient="records"),
    }


def _write_outputs(
    cfg: Config,
    outdir: Path,
    args: argparse.Namespace,
    predictions: pd.DataFrame,
    ev,
    summary: report_mod.Summary,
    summary_tsv: pd.DataFrame,
    id_detail: pd.DataFrame,
    context: dict,
    metadata: dict,
    ips,
) -> None:
    """Write every output file of the run, one progress unit per file.

    The tables are written biggest-first in the sense that matters: the
    172 MB prediction table and the 48 MB evidence table are the two units the
    bar spends most of its time on.
    """
    main_table = (
        predictions if args.keep_non_rga else predictions[predictions["is_rga"]]
    )
    tables: list[tuple[pd.DataFrame, str]] = [
        (main_table, "rga_predictions.tsv"),
        (predictions[predictions["is_rga"]], "rga_predictions_rga_only.tsv"),
        (ev.long, "rga_domain_evidence_long.tsv"),
        (summary_tsv, "rga_summary_counts.tsv"),
        (id_detail, "unmatched_ids_report.tsv"),
        (context["cc_sensitivity"], "cc_segment_sensitivity.tsv"),
        (context["cc_policy_counts"], "cc_policy_sensitivity.tsv"),
        (accession_audit(cfg, ips.accession_counts), "accession_audit.tsv"),
    ]
    documents: list[tuple[str, str]] = [
        ("report.md", report_mod.render_markdown(context, summary)),
        ("report.html", report_mod.render_html(context, summary)),
        ("run_metadata.json", json.dumps(metadata, indent=2, sort_keys=True, default=str)),
    ]
    write_locus = bool(cfg.raw["output"].get("write_locus_summary"))
    total = len(tables) + len(documents) + int(write_locus)

    with tracked("writing outputs", total=total) as report:
        for frame, name in tables:
            write_tsv(frame, outdir / name, cfg)
            report(1)
        if write_locus:
            write_tsv(
                locus_summary(predictions, cfg),
                outdir / "rga_predictions_by_locus.tsv",
                cfg,
            )
            report(1)
        for name, text in documents:
            (outdir / name).write_text(text, encoding="utf-8")
            report(1)
    LOGGER.info("wrote report.md, report.html and run_metadata.json")


#: How this script names itself in the recorded command. Taken from the file
#: rather than from ``sys.argv[0]``, which is the *host* program's path whenever
#: ``main`` is called programmatically -- pytest, a notebook, a Nextflow wrapper.
_SCRIPT_PATH = "code/rgas/rgas_prediction.py"


def main(argv: list[str] | None = None) -> int:
    """Parse the command line and run the pipeline.

    ``argv`` is recorded on the parsed arguments so that the "Reproduce this
    run" block in the report quotes the arguments this call actually received.
    Reading ``sys.argv`` instead would be right in production and wrong whenever
    ``main`` is called programmatically -- from a test, a notebook, or a
    Nextflow wrapper -- which is exactly where a misleading command hurts.
    """
    args = build_parser().parse_args(argv)
    args.invoked_with = list(argv) if argv is not None else sys.argv[1:]
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
