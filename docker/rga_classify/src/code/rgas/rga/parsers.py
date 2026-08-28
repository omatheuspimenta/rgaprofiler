"""One parser per annotation tool, each returning a tidy :class:`pandas.DataFrame`.

Coordinate convention
---------------------
Every interval produced by this module is **1-based and inclusive**, matching
the convention of all six upstream tools:

===================  ==========================================================
InterProScan TSV     columns 7/8, 1-based inclusive
DeepTMHMM GFF3       columns 3/4, 1-based inclusive
Phobius (short)      ``<start>-<end>`` inside the topology string, 1-based
SignalP 6.0          ``CS pos: 30-31`` means the signal peptide is residues 1-30
DeepCoil2 ``.out``   row *i* after the header is residue *i*
===================  ==========================================================
"""

from __future__ import annotations

import csv
import logging
import re
import tarfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

import pandas as pd

from .config import Config, normalize_id
from .progress import ProgressCallback, null_progress

LOGGER = logging.getLogger(__name__)

#: Column names of the 15-column InterProScan TSV (the file carries no header).
INTERPROSCAN_COLUMNS: tuple[str, ...] = (
    "protein_id",
    "sequence_md5",
    "sequence_length",
    "analysis",
    "accession",
    "signature_description",
    "start",
    "end",
    "score",
    "status",
    "date",
    "interpro_accession",
    "interpro_description",
    "go_annotations",
    "pathway_annotations",
)

#: Rows read per chunk when streaming the (potentially multi-GB) InterProScan TSV.
_IPS_CHUNKSIZE = 500_000

_CS_PATTERN = re.compile(r"CS pos:\s*(\d+)-(\d+)\.\s*Pr:\s*([\d.]+)")
_PHOBIUS_TM_PATTERN = re.compile(r"(\d+)-(\d+)")
_PHOBIUS_SP_PATTERN = re.compile(r"^n(\d+)-(\d+)c(\d+)/(\d+)")
_DEEPTMHMM_HEADER = re.compile(
    r"^#\s*(\S+)\s+(Length|Number of predicted TMRs):\s*(\d+)"
)


@dataclass
class InterProScanResult:
    """Filtered InterProScan hits plus the auditing information they carry.

    Attributes
    ----------
    hits : pandas.DataFrame
        Columns ``protein_id, analysis, accession, signature_description,
        start, end, score, feature``. One row per (hit, feature) pair.
    domain_hits : pandas.DataFrame
        All hits from ``integrated_domain_analyses``, used to detect integrated
        domains. Same columns minus ``feature``.
    protein_ids : set of str
        Every protein ID seen in the file, whether or not it produced a hit of
        interest.
    accession_counts : collections.Counter
        Hit counts for every accession referenced by the configuration, plus
        the ``watch_accessions``.
    sequence_lengths : dict
        Protein ID -> sequence length as reported by InterProScan.
    n_rows : int
        Total number of rows read.
    """

    hits: pd.DataFrame
    domain_hits: pd.DataFrame
    protein_ids: set[str]
    accession_counts: Counter[str] = field(default_factory=Counter)
    sequence_lengths: dict[str, int] = field(default_factory=dict)
    n_rows: int = 0
    analyses_seen: Counter[str] = field(default_factory=Counter)
    run_dates: set[str] = field(default_factory=set)


class _CountingReader:
    """Binary file wrapper that reports how many bytes have been consumed.

    ``pandas.read_csv`` accepts any file-like object; wrapping the handle is the
    only way to get a genuine progress bar over a chunked read, because the
    chunk iterator itself knows nothing about the size of the file.

    Every read entry point has to be counted, not just ``read``: pandas wraps a
    binary handle in a :class:`io.TextIOWrapper`, and that calls ``read1``.
    Counting only ``read`` leaves the bar pinned at 0 % for the whole parse,
    which is why ``test_progress.py`` asserts that the reported bytes add up to
    the file size.

    Parameters
    ----------
    handle : io.BufferedReader
        The open binary file.
    on_progress : ProgressCallback
        Called with the number of bytes returned by each ``read``.
    """

    def __init__(self, handle, on_progress: ProgressCallback) -> None:
        """Wrap ``handle``, reporting consumed bytes to ``on_progress``."""
        self._handle = handle
        self._on_progress = on_progress

    def read(self, size: int = -1) -> bytes:
        """Read from the wrapped handle and report the bytes consumed."""
        block = self._handle.read(size)
        if block:
            self._on_progress(len(block))
        return block

    def read1(self, size: int = -1) -> bytes:
        """Read a single block and report the bytes consumed.

        This is the method :class:`io.TextIOWrapper` actually calls, so it is
        the one that drives the bar in practice.
        """
        block = self._handle.read1(size)
        if block:
            self._on_progress(len(block))
        return block

    def readline(self, size: int = -1) -> bytes:
        """Read one line and report the bytes consumed."""
        line = self._handle.readline(size)
        if line:
            self._on_progress(len(line))
        return line

    def __iter__(self):
        """Iterate lines, reporting the bytes consumed by each."""
        for line in self._handle:
            self._on_progress(len(line))
            yield line

    def __getattr__(self, name: str):
        """Delegate everything else (``seekable``, ``close``, ...) to the handle."""
        return getattr(self._handle, name)


def parse_interproscan(
    path: str | Path,
    cfg: Config,
    on_progress: ProgressCallback = null_progress,
) -> InterProScanResult:
    """Stream the InterProScan TSV and keep only the rows the pipeline needs.

    Matching is done on accessions only -- both the signature accession
    (column 5) and the integrated InterPro accession (column 12) -- never on
    description strings.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the InterProScan TSV output (no header, 15 columns).
    cfg : Config
        Resolved pipeline configuration.
    on_progress : ProgressCallback, optional
        Reports the number of bytes of the TSV consumed so far. The total is
        announced as the file size before the first chunk is read.

    Returns
    -------
    InterProScanResult
        Feature hits, domain hits for the integrated-domain scan, and audit
        counters.
    """
    acc_to_features = cfg.accession_to_features()
    excluded = set(cfg.raw["excluded_analyses"])
    domain_analyses = set(cfg.raw["integrated_domain_analyses"])
    watched = set(cfg.raw.get("watch_accessions", {}))
    ops = cfg.raw["ids"]["per_tool"]["interproscan"]
    suffixes = cfg.raw["ids"].get("id_suffixes", [])

    state = _IPSState()
    path = Path(path)
    on_progress(0, total=path.stat().st_size)
    with path.open("rb") as handle:
        reader = pd.read_csv(
            _CountingReader(handle, on_progress),
            sep="\t",
            header=None,
            names=INTERPROSCAN_COLUMNS,
            dtype=str,
            na_filter=False,
            quoting=csv.QUOTE_NONE,
            chunksize=_IPS_CHUNKSIZE,
        )
        for chunk in reader:
            chunk["protein_id"] = [
                normalize_id(value, ops, suffixes) for value in chunk["protein_id"]
            ]
            _accumulate_chunk(
                state, chunk, acc_to_features, excluded, domain_analyses, watched
            )

    protein_ids, counts = state.protein_ids, state.counts
    analyses, dates, n_rows = state.analyses, state.dates, state.n_rows
    lengths = state.lengths
    feature_frames, domain_frames = state.feature_frames, state.domain_frames
    hits = _concat(feature_frames, _feature_columns())
    domain_hits = _concat(domain_frames, _hit_columns())
    lengths = {pid: int(length) for pid, length in lengths.items() if pd.notna(length)}
    LOGGER.info(
        "InterProScan: %d rows, %d proteins, %d feature hits, %d domain hits",
        n_rows,
        len(protein_ids),
        len(hits),
        len(domain_hits),
    )
    return InterProScanResult(
        hits=hits,
        domain_hits=domain_hits,
        protein_ids=protein_ids,
        accession_counts=counts,
        sequence_lengths=lengths,
        n_rows=n_rows,
        analyses_seen=analyses,
        run_dates=dates,
    )


@dataclass
class _IPSState:
    """Mutable accumulator threaded through the InterProScan chunk loop."""

    feature_frames: list[pd.DataFrame] = field(default_factory=list)
    domain_frames: list[pd.DataFrame] = field(default_factory=list)
    protein_ids: set[str] = field(default_factory=set)
    counts: Counter[str] = field(default_factory=Counter)
    analyses: Counter[str] = field(default_factory=Counter)
    lengths: dict[str, float] = field(default_factory=dict)
    dates: set[str] = field(default_factory=set)
    n_rows: int = 0


def _accumulate_chunk(
    state: _IPSState,
    chunk: pd.DataFrame,
    acc_to_features: dict[str, tuple[str, ...]],
    excluded: set[str],
    domain_analyses: set[str],
    watched: set[str],
) -> None:
    """Fold one chunk of the InterProScan TSV into the accumulator."""
    state.n_rows += len(chunk)
    state.protein_ids.update(chunk["protein_id"].unique().tolist())
    state.analyses.update(chunk["analysis"].value_counts().to_dict())
    state.dates.update(chunk["date"].unique().tolist())
    state.lengths.update(
        dict(
            zip(
                chunk["protein_id"],
                pd.to_numeric(chunk["sequence_length"], errors="coerce"),
            )
        )
    )

    chunk = chunk[~chunk["analysis"].isin(excluded)]
    if chunk.empty:
        return

    sig_hit = chunk["accession"].isin(acc_to_features)
    ipr_hit = chunk["interpro_accession"].isin(acc_to_features)
    state.counts.update(chunk.loc[sig_hit, "accession"].value_counts().to_dict())
    state.counts.update(
        chunk.loc[ipr_hit, "interpro_accession"].value_counts().to_dict()
    )
    state.counts.update(
        chunk.loc[chunk["accession"].isin(watched), "accession"]
        .value_counts()
        .to_dict()
    )

    selected = chunk[sig_hit | ipr_hit]
    if not selected.empty:
        state.feature_frames.append(_explode_features(selected, acc_to_features))
    domains = chunk[chunk["analysis"].isin(domain_analyses)]
    if not domains.empty:
        state.domain_frames.append(_tidy_hits(domains))


def _hit_columns() -> list[str]:
    """Column order of a tidy InterProScan hit table."""
    return [
        "protein_id",
        "analysis",
        "accession",
        "signature_description",
        "start",
        "end",
        "score",
    ]


def _feature_columns() -> list[str]:
    """Column order of a tidy InterProScan hit table carrying a feature label."""
    return [*_hit_columns(), "feature"]


def _concat(frames: list[pd.DataFrame], columns: list[str]) -> pd.DataFrame:
    """Concatenate parser chunks, returning a correctly typed empty frame if none."""
    if not frames:
        return pd.DataFrame({name: pd.Series(dtype="object") for name in columns})
    return pd.concat(frames, ignore_index=True)


def _tidy_hits(chunk: pd.DataFrame) -> pd.DataFrame:
    """Reduce an InterProScan chunk to the columns the pipeline stores."""
    out = chunk[
        [
            "protein_id",
            "analysis",
            "accession",
            "signature_description",
            "start",
            "end",
            "score",
        ]
    ].copy()
    out["start"] = pd.to_numeric(out["start"], errors="coerce").astype("Int64")
    out["end"] = pd.to_numeric(out["end"], errors="coerce").astype("Int64")
    return out.dropna(subset=["start", "end"])


def _explode_features(
    chunk: pd.DataFrame, acc_to_features: dict[str, tuple[str, ...]]
) -> pd.DataFrame:
    """Attach the controlled-vocabulary feature(s) implied by each hit.

    A hit may match through its signature accession, its InterPro accession, or
    both; the union of the implied features is used, so a single hit never
    produces duplicated feature rows.
    """
    tidy = _tidy_hits(chunk)
    if tidy.empty:
        return pd.DataFrame(
            {name: pd.Series(dtype="object") for name in _feature_columns()}
        )
    signature = chunk.loc[tidy.index, "accession"]
    integrated = chunk.loc[tidy.index, "interpro_accession"]

    records: list[tuple[int, str]] = []
    for idx, sig_acc, ipr_acc in zip(tidy.index, signature, integrated):
        features = set(acc_to_features.get(sig_acc, ())) | set(
            acc_to_features.get(ipr_acc, ())
        )
        records.extend((idx, feature) for feature in sorted(features))
    if not records:
        return pd.DataFrame(
            {name: pd.Series(dtype="object") for name in _feature_columns()}
        )

    index, features = zip(*records)
    exploded = tidy.loc[list(index)].copy()
    exploded["feature"] = list(features)
    return exploded.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Phobius
# ---------------------------------------------------------------------------


def parse_phobius(path: str | Path, cfg: Config) -> pd.DataFrame:
    """Parse the Phobius *short* output format.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the Phobius short-format file.
    cfg : Config
        Resolved pipeline configuration (used for ID normalisation).

    Returns
    -------
    pandas.DataFrame
        Columns ``protein_id, n_tm_phobius, sp_phobius, tm_intervals_phobius,
        sp_end_phobius``. ``tm_intervals_phobius`` holds a list of
        ``(start, end)`` tuples.

    Notes
    -----
    A Phobius topology string looks like ``n4-15c20/21o396-419i``: the leading
    ``n…c…/…`` block describes the signal peptide and every remaining
    ``<start>-<end>`` pair is a transmembrane helix. The signal-peptide block is
    removed before the helices are read so its coordinates are never mistaken
    for a helix.
    """
    ops = cfg.raw["ids"]["per_tool"]["phobius"]
    suffixes = cfg.raw["ids"].get("id_suffixes", [])
    rows: list[dict[str, object]] = []

    with Path(path).open(encoding="utf-8") as handle:
        next(handle, None)  # header line ("SEQENCE ID  TM SP PREDICTION")
        for line in handle:
            fields = line.split()
            if len(fields) < 4:
                continue
            protein_id = normalize_id(fields[0], ops, suffixes)
            topology = fields[3]
            sp_match = _PHOBIUS_SP_PATTERN.match(topology)
            sp_end = int(sp_match.group(3)) if sp_match else None
            remainder = topology[sp_match.end() :] if sp_match else topology
            intervals = [
                (int(start), int(end))
                for start, end in _PHOBIUS_TM_PATTERN.findall(remainder)
            ]
            rows.append(
                {
                    "protein_id": protein_id,
                    "n_tm_phobius": int(fields[1]),
                    "sp_phobius": fields[2] == "Y",
                    "tm_intervals_phobius": intervals,
                    "sp_end_phobius": sp_end,
                }
            )

    frame = pd.DataFrame(rows)
    LOGGER.info("Phobius: %d proteins", len(frame))
    return frame


# ---------------------------------------------------------------------------
# DeepTMHMM
# ---------------------------------------------------------------------------


def parse_deeptmhmm(path: str | Path, cfg: Config) -> pd.DataFrame:
    """Parse a DeepTMHMM ``TMRs.gff3`` file.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to ``TMRs.gff3``.
    cfg : Config
        Resolved pipeline configuration (used for ID normalisation).

    Returns
    -------
    pandas.DataFrame
        Columns ``protein_id, length_deeptmhmm, n_tm_deeptmhmm,
        tm_intervals_deeptmhmm, sp_end_deeptmhmm``.

    Notes
    -----
    The parser is driven by the data lines, not by the ``#`` comments: the
    legacy implementation keyed on the ``Length:`` comment and silently dropped
    any block that lacked one, as well as the final block after the last ``//``.
    """
    ops = cfg.raw["ids"]["per_tool"]["deeptmhmm"]
    suffixes = cfg.raw["ids"].get("id_suffixes", [])

    tm_intervals: dict[str, list[tuple[int, int]]] = {}
    signal_end: dict[str, int] = {}
    lengths: dict[str, int] = {}
    order: list[str] = []

    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                match = _DEEPTMHMM_HEADER.match(line)
                if match and match.group(2) == "Length":
                    protein_id = normalize_id(match.group(1), ops, suffixes)
                    lengths[protein_id] = int(match.group(3))
                continue
            if line.startswith("//") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            protein_id = normalize_id(fields[0], ops, suffixes)
            if protein_id not in tm_intervals:
                tm_intervals[protein_id] = []
                order.append(protein_id)
            region = fields[1].strip()
            start, end = int(fields[2]), int(fields[3])
            if region == "TMhelix":
                tm_intervals[protein_id].append((start, end))
            elif region == "signal":
                signal_end[protein_id] = max(signal_end.get(protein_id, 0), end)

    frame = pd.DataFrame(
        {
            "protein_id": order,
            "length_deeptmhmm": [lengths.get(pid) for pid in order],
            "n_tm_deeptmhmm": [len(tm_intervals[pid]) for pid in order],
            "tm_intervals_deeptmhmm": [tm_intervals[pid] for pid in order],
            "sp_end_deeptmhmm": [signal_end.get(pid) for pid in order],
        }
    )
    LOGGER.info("DeepTMHMM: %d proteins", len(frame))
    return frame


# ---------------------------------------------------------------------------
# SignalP 6.0
# ---------------------------------------------------------------------------


def parse_signalp(path: str | Path, cfg: Config) -> pd.DataFrame:
    """Parse ``prediction_results.txt`` from SignalP 6.0.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the SignalP 6.0 summary table.
    cfg : Config
        Resolved pipeline configuration.

    Returns
    -------
    pandas.DataFrame
        Columns ``protein_id, signalp_prediction, sp_prob, sp_signalp,
        cleavage_site, sp_end_signalp``. ``sp_prob`` is the probability of the
        predicted class, and ``sp_end_signalp`` is the last residue of the
        signal peptide (the residue before the cleavage site).

    Notes
    -----
    Column 1 of the file is the complete FASTA header, not the protein ID; the
    configured ``strip_after_whitespace`` normalisation recovers the ID.
    """
    ops = cfg.raw["ids"]["per_tool"]["signalp"]
    suffixes = cfg.raw["ids"].get("id_suffixes", [])
    positive = set(cfg.raw["signal_peptide"]["positive_predictions"])
    floor = cfg.raw["signal_peptide"].get("min_probability")

    header_classes: list[str] = []
    rows: list[dict[str, object]] = []

    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                if "\tPrediction\t" in line:
                    header_classes = [c.split("(")[0] for c in line.split("\t")[2:-1]]
                continue
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            protein_id = normalize_id(fields[0], ops, suffixes)
            prediction = fields[1]
            probabilities = [
                _as_float(value) for value in fields[2 : 2 + len(header_classes)]
            ]
            prob = _class_probability(prediction, header_classes, probabilities)
            cs_field = fields[-1] if len(fields) > 2 + len(header_classes) else ""
            match = _CS_PATTERN.search(cs_field)
            sp_end = int(match.group(1)) if match else None
            called = prediction in positive and (
                floor is None or (prob or 0.0) >= floor
            )
            rows.append(
                {
                    "protein_id": protein_id,
                    "signalp_prediction": prediction,
                    "sp_prob": prob,
                    "sp_signalp": called,
                    "cleavage_site": (
                        f"{match.group(1)}-{match.group(2)}" if match else None
                    ),
                    "cleavage_prob": float(match.group(3)) if match else None,
                    "sp_end_signalp": sp_end,
                }
            )

    frame = pd.DataFrame(rows)
    LOGGER.info("SignalP 6.0: %d proteins", len(frame))
    return frame


def _as_float(value: str) -> float | None:
    """Convert a SignalP probability field to ``float``, tolerating blanks."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _class_probability(
    prediction: str, classes: list[str], probabilities: list[float | None]
) -> float | None:
    """Return the probability SignalP assigned to the class it predicted."""
    for name, prob in zip(classes, probabilities):
        if name == prediction:
            return prob
    return max((p for p in probabilities if p is not None), default=None)


# ---------------------------------------------------------------------------
# DeepLoc 2.0
# ---------------------------------------------------------------------------


def parse_deeploc(path: str | Path, cfg: Config) -> pd.DataFrame:
    """Parse a DeepLoc 2.0 ``results*.csv`` file.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the DeepLoc results CSV.
    cfg : Config
        Resolved pipeline configuration.

    Returns
    -------
    pandas.DataFrame
        Columns ``protein_id, predicted_localization, localization_prob,
        all_localizations, deeploc_signals``.

    Notes
    -----
    ``Localizations`` is multi-label and pipe-separated (for example
    ``Cytoplasm|Nucleus``). The first label is reported as the primary
    localisation and its probability is read from the matching per-class column.
    """
    ops = cfg.raw["ids"]["per_tool"]["deeploc"]
    suffixes = cfg.raw["ids"].get("id_suffixes", [])
    frame = pd.read_csv(path, dtype={"Protein_ID": str}, na_filter=False)
    frame = frame.rename(columns={"Protein_ID": "protein_id"})
    frame["protein_id"] = [
        normalize_id(value, ops, suffixes) for value in frame["protein_id"]
    ]

    labels = frame["Localizations"].astype(str)
    primary = labels.str.split("|").str[0]
    probability = [
        _lookup_probability(frame, index, label)
        for index, label in zip(frame.index, primary)
    ]
    out = pd.DataFrame(
        {
            "protein_id": frame["protein_id"],
            "predicted_localization": primary,
            "localization_prob": probability,
            "all_localizations": labels,
            "deeploc_signals": frame.get(
                "Signals", pd.Series([""] * len(frame))
            ).astype(str),
        }
    )
    LOGGER.info("DeepLoc 2.0: %d proteins", len(out))
    return out


def _lookup_probability(frame: pd.DataFrame, index: int, label: str) -> float | None:
    """Read the per-class probability column matching a localisation label."""
    if label in frame.columns:
        try:
            return float(frame.at[index, label])
        except (TypeError, ValueError):
            return None
    return None


# ---------------------------------------------------------------------------
# DeepCoil2
# ---------------------------------------------------------------------------


def _segments_from_lines(lines: Iterable[str]) -> list[tuple[int, int, float]]:
    """Extract raw DeepCoil2 segments from the lines of one ``.out`` file.

    ``cc`` is a per-segment plateau, not a per-residue probability: DeepCoil2
    has already performed peak detection, so a segment is a maximal run of
    residues carrying the *same* non-zero ``cc`` value. Splitting on a change of
    value keeps two adjacent segments with different scores apart, which a
    simple ``cc > threshold`` scan would merge.

    Parameters
    ----------
    lines : iterable of str
        Lines of a DeepCoil2 ``.out`` file, including the header.

    Returns
    -------
    list of tuple
        ``(start, end, cc)`` triples, 1-based inclusive, unfiltered.
    """
    segments: list[tuple[int, int, float]] = []
    position = 0
    current_value = 0.0
    current_start = 0

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("aa\t") or line.startswith("aa "):
            continue
        fields = line.split()
        if len(fields) < 2:
            continue
        position += 1
        try:
            value = float(fields[1])
        except ValueError:
            continue
        if value != current_value:
            if current_value > 0.0:
                segments.append((current_start, position - 1, current_value))
            current_value = value
            current_start = position
    if current_value > 0.0:
        segments.append((current_start, position, current_value))
    return segments


def _segments_from_archive(archive: str) -> list[tuple[str, int, int, float]]:
    """Extract raw segments from every ``.out`` member of one ``.tar.xz``."""
    rows: list[tuple[str, int, int, float]] = []
    with tarfile.open(archive, "r:xz") as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith(".out"):
                continue
            handle = tar.extractfile(member)
            if handle is None:  # pragma: no cover - defensive
                continue
            name = Path(member.name).name[: -len(".out")]
            text = handle.read().decode("utf-8", errors="replace").splitlines()
            rows.extend(
                (name, start, end, value)
                for start, end, value in _segments_from_lines(text)
            )
    return rows


def _segments_from_directory(directory: Path) -> list[tuple[str, int, int, float]]:
    """Extract raw segments from every ``.out`` file in a directory."""
    rows: list[tuple[str, int, int, float]] = []
    for path in sorted(directory.glob("*.out")):
        with path.open(encoding="utf-8", errors="replace") as handle:
            rows.extend(
                (path.stem, start, end, value)
                for start, end, value in _segments_from_lines(handle)
            )
    return rows


def iter_deepcoil_sources(root: str | Path) -> Iterator[Path]:
    """Yield the DeepCoil2 inputs found under ``root``, deterministically.

    Directories of ``.out`` files take precedence over an archive of the same
    name, so an already-extracted part is never parsed twice.
    """
    root = Path(root)
    if root.is_file() and root.name.endswith(".tar.xz"):
        yield root
        return
    if root.is_dir() and any(root.glob("*.out")):
        yield root
        return
    directories = sorted(
        p for p in root.iterdir() if p.is_dir() and any(p.glob("*.out"))
    )
    yield from directories
    extracted = {p.name for p in directories}
    for archive in sorted(root.glob("*.tar.xz")):
        if archive.name[: -len(".tar.xz")] not in extracted:
            yield archive


def parse_deepcoil(
    root: str | Path,
    cfg: Config,
    workers: int = 1,
    on_progress: ProgressCallback = null_progress,
) -> pd.DataFrame:
    """Build the raw DeepCoil2 segment table for a whole proteome.

    The table is intentionally *unfiltered*: it records every segment DeepCoil2
    proposed, with its plateau score. Threshold, minimum-length and gap-merging
    filters are applied downstream, which makes the sensitivity analysis over
    those parameters essentially free.

    Parameters
    ----------
    root : str or pathlib.Path
        Directory holding per-part directories and/or ``.tar.xz`` archives of
        ``.out`` files, a single such directory, or a single archive.
    cfg : Config
        Resolved pipeline configuration.
    workers : int, default 1
        Number of worker processes used to read archives in parallel.
    on_progress : ProgressCallback, optional
        Advanced by one as each directory or archive is finished; the total is
        the number of sources discovered under ``root``.

    Returns
    -------
    pandas.DataFrame
        Columns ``deepcoil_id, start, end, cc``, sorted deterministically.
    """
    sources = list(iter_deepcoil_sources(root))
    LOGGER.info("DeepCoil2: reading %d source(s) from %s", len(sources), root)
    on_progress(0, total=len(sources))
    rows: list[tuple[str, int, int, float]] = []

    archives = [str(p) for p in sources if p.is_file()]
    directories = [p for p in sources if p.is_dir()]
    for directory in directories:
        rows.extend(_segments_from_directory(directory))
        on_progress(1)
    if archives:
        if workers > 1:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                for result in pool.map(_segments_from_archive, archives):
                    rows.extend(result)
                    on_progress(1)
        else:
            for archive in archives:
                rows.extend(_segments_from_archive(archive))
                on_progress(1)

    frame = pd.DataFrame(rows, columns=["deepcoil_id", "start", "end", "cc"])
    frame = frame.sort_values(
        ["deepcoil_id", "start", "end"], kind="stable"
    ).reset_index(drop=True)
    LOGGER.info(
        "DeepCoil2: %d raw segments over %d proteins",
        len(frame),
        frame["deepcoil_id"].nunique(),
    )
    return frame


def deepcoil_protein_ids(root: str | Path, cfg: Config) -> set[str]:
    """List the DeepCoil2 protein IDs without parsing any residue scores."""
    ids: set[str] = set()
    for source in iter_deepcoil_sources(root):
        if source.is_dir():
            ids.update(path.stem for path in source.glob("*.out"))
        else:
            with tarfile.open(source, "r:xz") as tar:
                ids.update(
                    Path(m.name).name[: -len(".out")]
                    for m in tar
                    if m.isfile() and m.name.endswith(".out")
                )
    return ids
