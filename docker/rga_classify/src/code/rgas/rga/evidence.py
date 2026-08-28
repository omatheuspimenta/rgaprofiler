"""Harmonisation of six tool outputs into one controlled feature vocabulary.

The module produces two things:

``long`` evidence
    A DataFrame with one row per protein x feature x supporting hit, written to
    ``rga_domain_evidence_long.tsv`` so that every call can be traced back to
    the exact signature that produced it.

per-protein ``records``
    One dictionary per protein carrying the collapsed feature set plus the
    numeric columns (helix counts, CC statistics, localisation, ...) consumed by
    the rule engine and by the report.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Sequence

import pandas as pd

from .config import CC_DOMAIN_FEATURE, Config
from .progress import (
    PROTEIN_REPORT_INTERVAL,
    ProgressCallback,
    null_progress,
)

LOGGER = logging.getLogger(__name__)

Interval = tuple[int, int]

#: Evidence channels the pipeline knows about, in report order.
CHANNELS: tuple[str, ...] = (
    "interproscan",
    "phobius",
    "deeptmhmm",
    "signalp",
    "deeploc",
    "deepcoil",
)

#: Column order of the long evidence table.
LONG_COLUMNS: tuple[str, ...] = (
    "protein_id",
    "feature",
    "tool",
    "analysis",
    "accession",
    "signature_description",
    "start",
    "end",
    "score_or_evalue",
)


@dataclass
class CCSegments:
    """Coiled-coil segments retained for one protein after filtering.

    Attributes
    ----------
    segments : list of tuple
        ``(start, end, cc)`` triples, 1-based inclusive.
    """

    segments: list[tuple[int, int, float]] = field(default_factory=list)

    @property
    def n(self) -> int:
        """Number of retained segments."""
        return len(self.segments)

    @property
    def max_prob(self) -> float | None:
        """Highest plateau score among the retained segments."""
        return max((s[2] for s in self.segments), default=None)

    @property
    def mean_prob(self) -> float | None:
        """Length-weighted mean plateau score over the retained segments."""
        if not self.segments:
            return None
        total = sum((end - start + 1) for start, end, _ in self.segments)
        weighted = sum((end - start + 1) * value for start, end, value in self.segments)
        return weighted / total

    @property
    def total_length(self) -> int:
        """Total number of residues covered by the retained segments."""
        return sum(end - start + 1 for start, end, _ in self.segments)

    @property
    def intervals(self) -> list[Interval]:
        """Retained segments as plain ``(start, end)`` intervals."""
        return [(start, end) for start, end, _ in self.segments]


# ---------------------------------------------------------------------------
# interval helpers
# ---------------------------------------------------------------------------


def merge_intervals(
    intervals: Iterable[Interval], min_overlap: int = 1
) -> list[Interval]:
    """Merge overlapping 1-based inclusive intervals.

    Parameters
    ----------
    intervals : iterable of tuple
        ``(start, end)`` pairs; order is irrelevant.
    min_overlap : int, default 1
        Number of shared residues required before two intervals are merged.
        ``1`` merges intervals that share at least one residue; a larger value
        keeps intervals that barely touch apart.

    Returns
    -------
    list of tuple
        Disjoint intervals sorted by start coordinate.

    Examples
    --------
    >>> merge_intervals([(10, 20), (18, 30), (50, 60)])
    [(10, 30), (50, 60)]
    """
    ordered = sorted((int(s), int(e)) for s, e in intervals)
    merged: list[Interval] = []
    for start, end in ordered:
        if merged and _overlap(merged[-1], (start, end)) >= min_overlap:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _overlap(first: Interval, second: Interval) -> int:
    """Number of residues shared by two 1-based inclusive intervals."""
    return min(first[1], second[1]) - max(first[0], second[0]) + 1


def overlap_fraction(interval: Interval, others: Sequence[Interval]) -> float:
    """Fraction of ``interval`` covered by any of ``others``.

    Parameters
    ----------
    interval : tuple
        The ``(start, end)`` interval to measure.
    others : sequence of tuple
        Intervals that may cover it.

    Returns
    -------
    float
        Covered residues divided by the length of ``interval``; ``0.0`` when
        there is no overlap.
    """
    length = interval[1] - interval[0] + 1
    if length <= 0:
        return 0.0
    covered = sum(
        max(0, _overlap(interval, other)) for other in merge_intervals(others)
    )
    return min(covered, length) / length


# ---------------------------------------------------------------------------
# coiled-coil segment calling
# ---------------------------------------------------------------------------


def call_cc_segments(
    raw: Sequence[tuple[int, int, float]],
    threshold: float,
    min_length: int,
    max_gap: int,
) -> CCSegments:
    """Turn raw DeepCoil2 segments into called coiled-coil segments.

    The procedure is: (1) drop segments whose plateau score is below
    ``threshold``; (2) merge surviving segments separated by at most
    ``max_gap`` residues; (3) drop merged segments shorter than ``min_length``.
    Gap merging happens *before* the length filter so that a genuine coiled coil
    interrupted by one or two sub-threshold residues is not lost.

    Parameters
    ----------
    raw : sequence of tuple
        ``(start, end, cc)`` triples as emitted by DeepCoil2.
    threshold : float
        Minimum plateau score.
    min_length : int
        Minimum length in residues of a retained segment.
    max_gap : int
        Maximum number of residues separating two segments that are merged.

    Returns
    -------
    CCSegments
        The retained segments and their summary statistics.
    """
    kept = sorted((s for s in raw if s[2] >= threshold), key=lambda s: (s[0], s[1]))
    merged: list[list[float]] = []
    for start, end, value in kept:
        if merged and start - merged[-1][1] - 1 <= max_gap:
            merged[-1][1] = max(merged[-1][1], end)
            merged[-1][2] = max(merged[-1][2], value)
        else:
            merged.append([start, end, value])
    return CCSegments(
        [
            (int(start), int(end), float(value))
            for start, end, value in merged
            if end - start + 1 >= min_length
        ]
    )


def cc_sensitivity(
    raw_by_protein: dict[str, list[tuple[int, int, float]]],
    thresholds: Sequence[float],
    min_lengths: Sequence[int],
    max_gap: int,
) -> pd.DataFrame:
    """Count CC-positive proteins across a grid of segment-calling parameters.

    Parameters
    ----------
    raw_by_protein : dict
        Protein ID -> raw ``(start, end, cc)`` segments.
    thresholds : sequence of float
        Plateau-score cut-offs to evaluate.
    min_lengths : sequence of int
        Minimum segment lengths to evaluate.
    max_gap : int
        Gap-merging parameter, held constant across the grid.

    Returns
    -------
    pandas.DataFrame
        Columns ``threshold, min_length, n_proteins_with_cc, n_segments``.
    """
    rows = []
    for threshold in thresholds:
        for min_length in min_lengths:
            n_proteins = 0
            n_segments = 0
            for raw in raw_by_protein.values():
                called = call_cc_segments(raw, threshold, min_length, max_gap)
                if called.n:
                    n_proteins += 1
                    n_segments += called.n
            rows.append(
                {
                    "threshold": threshold,
                    "min_length": min_length,
                    "n_proteins_with_cc": n_proteins,
                    "n_segments": n_segments,
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# consensus helpers
# ---------------------------------------------------------------------------


def apply_cc_policy(policy: str, channels: dict[str, bool | None]) -> bool:
    """Combine the three coiled-coil channels according to ``--cc-policy``.

    The CC channels are not interchangeable and the policy names say which is
    which: ``rx_domain`` is the domain-level profile HMM channel,
    ``deepcoil`` and ``coils`` are the two propensity predictors. ``union`` and
    ``intersection`` range over every channel that is available.

    As in :func:`apply_policy`, a channel whose tool was not supplied is
    ``None``: it counts as "no evidence" for ``union`` and is skipped for
    ``intersection``, so a missing tool can never silently veto a call. Naming a
    single unavailable channel falls back to the union of the rest rather than
    returning a blanket ``False``.

    Parameters
    ----------
    policy : {'rx_domain', 'deepcoil', 'coils', 'union', 'intersection'}
        The configured or overridden CC policy.
    channels : dict
        ``{'rx_domain': ..., 'deepcoil': ..., 'coils': ...}``, each ``bool`` or
        ``None``.

    Returns
    -------
    bool
        The consensus CC call.
    """
    if policy in channels:
        selected = channels[policy]
        if selected is not None:
            return bool(selected)
        policy = "union"
    available = [value for value in channels.values() if value is not None]
    if not available:
        return False
    if policy == "intersection":
        return all(available)
    return any(available)


def apply_policy(policy: str, first: bool | None, second: bool | None) -> bool:
    """Combine two boolean channels according to a consensus policy.

    Parameters
    ----------
    policy : {'union', 'intersection', 'first', 'second'}
        ``union`` is a logical OR, ``intersection`` a logical AND; ``first`` and
        ``second`` select a single channel. Missing channels (``None``) are
        treated as "no evidence" for ``union`` and are skipped for
        ``intersection`` so that a missing tool never silently vetoes a call.
    first, second : bool or None
        Channel values, ``None`` when the tool was unavailable.

    Returns
    -------
    bool
        The consensus call.
    """
    if policy == "first":
        return bool(first) if first is not None else bool(second)
    if policy == "second":
        return bool(second) if second is not None else bool(first)
    available = [value for value in (first, second) if value is not None]
    if not available:
        return False
    if policy == "intersection":
        return all(available)
    return any(available)


def filter_helices_in_signal(
    helices: Sequence[Interval], sp_end: int | None, fraction: float
) -> tuple[list[Interval], int]:
    """Drop TM helices that lie inside the predicted signal-peptide region.

    Signal peptides are routinely mis-called as transmembrane helices because
    both are hydrophobic. A helix covered by the signal-peptide region
    (residues ``1..sp_end``) by at least ``fraction`` of its own length is
    discarded.

    Parameters
    ----------
    helices : sequence of tuple
        Predicted ``(start, end)`` helices.
    sp_end : int or None
        Last residue of the signal peptide, or ``None`` when no signal peptide
        was predicted by any tool.
    fraction : float
        Overlap fraction above which a helix is discarded.

    Returns
    -------
    tuple
        The retained helices and the number of discarded ones.
    """
    if not sp_end:
        return list(helices), 0
    signal = [(1, int(sp_end))]
    kept = [h for h in helices if overlap_fraction(h, signal) < fraction]
    return kept, len(helices) - len(kept)


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------


@dataclass
class Evidence:
    """Harmonised evidence for a whole proteome.

    Attributes
    ----------
    records : list of dict
        One dictionary per protein. Deliberately *not* a DataFrame: feature keys
        such as ``feat_NB-ARC`` are not valid Python identifiers and would be
        silently renamed by :meth:`pandas.DataFrame.itertuples`, and a
        299k-row object-dtype frame would double the peak memory of the run for
        no benefit.
    long : pandas.DataFrame
        One row per protein x feature x supporting hit.
    available : dict
        Channel name -> whether the corresponding tool output was supplied.
    domain_descriptions : dict
        Accession -> signature description, used to render the integrated-domain
        columns in words as well as accessions.
    cc_contingency : dict
        2x2 counts of DeepCoil2 versus InterProScan Coils CC calls.
    raw_cc : dict
        Protein ID -> raw DeepCoil2 segments, kept for the sensitivity analysis.
    """

    long: pd.DataFrame
    available: dict[str, bool]
    records: list[dict] = field(default_factory=list)
    cc_contingency: dict[str, int] = field(default_factory=dict)
    domain_descriptions: dict[str, str] = field(default_factory=dict)
    raw_cc: dict[str, list[tuple[int, int, float]]] = field(default_factory=dict)


def channel_availability(
    phobius: pd.DataFrame | None,
    deeptmhmm: pd.DataFrame | None,
    signalp: pd.DataFrame | None,
    deeploc: pd.DataFrame | None,
    deepcoil_segments: pd.DataFrame | None,
) -> dict[str, bool]:
    """Report which evidence channels were supplied for this run."""
    return {
        "interproscan": True,
        "phobius": phobius is not None,
        "deeptmhmm": deeptmhmm is not None,
        "signalp": signalp is not None,
        "deeploc": deeploc is not None,
        "deepcoil": deepcoil_segments is not None,
    }


def _group_intervals(
    hits: pd.DataFrame, min_overlap: int, analyses: set[str] | None = None
) -> dict[str, dict[str, list[Interval]]]:
    """Collect merged intervals per protein and per feature from IPS hits.

    Parameters
    ----------
    hits : pandas.DataFrame
        Tidy InterProScan feature hits.
    min_overlap : int
        Residues two intervals must share before they are merged.
    analyses : set of str, optional
        When given, only hits from these signature databases are used. This is
        how the repeat-level LRR copy number is computed separately from the
        all-sources count.
    """
    collected: dict[str, dict[str, list[Interval]]] = {}
    for protein_id, feature, start, end, analysis in zip(
        hits["protein_id"],
        hits["feature"],
        hits["start"],
        hits["end"],
        hits["analysis"],
    ):
        if analyses is not None and analysis not in analyses:
            continue
        collected.setdefault(protein_id, {}).setdefault(feature, []).append(
            (int(start), int(end))
        )
    return {
        protein_id: {
            feature: merge_intervals(intervals, min_overlap)
            for feature, intervals in features.items()
        }
        for protein_id, features in collected.items()
    }


def _group_databases(hits: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Count the distinct signature databases supporting each protein x feature.

    Used by the confidence formula: a domain seen by several independent
    databases is stronger evidence than one seen by a single profile.
    """
    collected: dict[str, dict[str, set[str]]] = {}
    for protein_id, feature, analysis in zip(
        hits["protein_id"], hits["feature"], hits["analysis"]
    ):
        collected.setdefault(protein_id, {}).setdefault(feature, set()).add(analysis)
    return {
        protein_id: {feature: len(analyses) for feature, analyses in features.items()}
        for protein_id, features in collected.items()
    }


def _group_hit_labels(
    hits: pd.DataFrame, max_per_feature: int = 4
) -> dict[str, dict[str, str]]:
    """Build the provenance snippets quoted by the ``reason`` field.

    Each protein x feature is summarised as at most ``max_per_feature``
    ``ACCESSION @ start-end`` labels plus a count of the remaining hits, so the
    reason string stays readable while remaining traceable.
    """
    collected: dict[str, dict[str, list[str]]] = {}
    for protein_id, feature, accession, start, end in zip(
        hits["protein_id"],
        hits["feature"],
        hits["accession"],
        hits["start"],
        hits["end"],
    ):
        collected.setdefault(protein_id, {}).setdefault(feature, []).append(
            f"{accession} @ {int(start)}-{int(end)}"
        )
    labels: dict[str, dict[str, str]] = {}
    for protein_id, features in collected.items():
        labels[protein_id] = {}
        for feature, items in features.items():
            shown = sorted(set(items))
            text = ", ".join(shown[:max_per_feature])
            if len(shown) > max_per_feature:
                text += f", +{len(shown) - max_per_feature} more"
            labels[protein_id][feature] = text
    return labels


def _group_accessions(hits: pd.DataFrame) -> dict[str, dict[str, list[str]]]:
    """Map each protein to the sorted accessions supporting each of its features.

    This is the machine-readable sibling of :func:`_group_hit_labels`: no
    coordinates, no truncation, so a reader can filter ``feature_accessions``
    in a spreadsheet without parsing prose out of ``reason``.
    """
    collected: dict[str, dict[str, set[str]]] = {}
    for protein_id, feature, accession in zip(
        hits["protein_id"], hits["feature"], hits["accession"]
    ):
        collected.setdefault(protein_id, {}).setdefault(feature, set()).add(accession)
    return {
        protein_id: {feature: sorted(values) for feature, values in features.items()}
        for protein_id, features in collected.items()
    }


def _domain_descriptions(domain_hits: pd.DataFrame) -> dict[str, str]:
    """Accession -> signature description, for the integrated-domain columns.

    Built before ``domain_hits`` is released. It is a few thousand entries for a
    whole proteome, so it costs nothing to keep and turns ``PF03106`` in the
    output into ``WRKY DNA-binding domain`` for the reader.
    """
    return {
        str(accession): str(description)
        for accession, description in zip(
            domain_hits["accession"], domain_hits["signature_description"]
        )
        if description and description != "-"
    }


def _noncanonical_pfam(
    domain_hits: pd.DataFrame, canonical: set[str]
) -> dict[str, list[str]]:
    """Map each protein to the sorted non-canonical domain accessions it carries."""
    collected: dict[str, set[str]] = {}
    for protein_id, accession in zip(
        domain_hits["protein_id"], domain_hits["accession"]
    ):
        if accession not in canonical:
            collected.setdefault(protein_id, set()).add(accession)
    return {protein_id: sorted(values) for protein_id, values in collected.items()}


def _lookup_frame(
    frame: pd.DataFrame | None, key: str = "protein_id"
) -> dict[str, dict]:
    """Index a parser DataFrame by protein ID as plain dictionaries.

    A tool that reports the same protein twice (duplicated FASTA entries are
    common) must not produce two output rows, so only the first record of each
    protein is kept and the duplication is logged.
    """
    if frame is None or frame.empty:
        return {}
    duplicated = int(frame[key].duplicated().sum())
    if duplicated:
        LOGGER.warning(
            "%d duplicated protein ID(s) in a tool output; keeping the first",
            duplicated,
        )
        frame = frame.drop_duplicates(subset=[key], keep="first")
    return frame.set_index(key).to_dict(orient="index")


def _raw_cc_by_protein(
    segments: pd.DataFrame | None, id_map: dict[str, str]
) -> dict[str, list[tuple[int, int, float]]]:
    """Group raw DeepCoil2 segments by canonical protein ID."""
    if segments is None or segments.empty:
        return {}
    grouped: dict[str, list[tuple[int, int, float]]] = {}
    for deepcoil_id, start, end, value in zip(
        segments["deepcoil_id"], segments["start"], segments["end"], segments["cc"]
    ):
        protein_id = id_map.get(deepcoil_id)
        if protein_id is None:
            continue
        grouped.setdefault(protein_id, []).append((int(start), int(end), float(value)))
    return grouped


def build_evidence(
    cfg: Config,
    protein_ids: Sequence[str],
    ips,
    phobius: pd.DataFrame | None,
    deeptmhmm: pd.DataFrame | None,
    signalp: pd.DataFrame | None,
    deeploc: pd.DataFrame | None,
    deepcoil_segments: pd.DataFrame | None,
    deepcoil_id_map: dict[str, str],
    options: dict,
    on_progress: ProgressCallback = null_progress,
) -> Evidence:
    """Harmonise every tool output into the per-protein evidence table.

    Parameters
    ----------
    cfg : Config
        Resolved pipeline configuration.
    protein_ids : sequence of str
        Canonical, sorted list of every protein in the proteome.
    ips : InterProScanResult
        Output of :func:`rga.parsers.parse_interproscan`.
    phobius, deeptmhmm, signalp, deeploc : pandas.DataFrame or None
        Parser outputs; ``None`` when the tool was not supplied.
    deepcoil_segments : pandas.DataFrame or None
        Raw DeepCoil2 segment table.
    deepcoil_id_map : dict
        DeepCoil file stem -> canonical protein ID.
    options : dict
        Resolved run options: ``policies``, ``cc_threshold``, ``cc_min_length``,
        ``cc_max_gap``, ``cc_tm_overlap``, ``sp_overlap``, ``min_tm_helices``,
        ``min_lrr_copies``, ``merge_min_overlap``.
    on_progress : ProgressCallback, optional
        Advanced every :data:`~rga.progress.PROTEIN_REPORT_INTERVAL` proteins;
        the total is the size of ``protein_ids``.

    Returns
    -------
    Evidence
        Wide and long evidence tables plus channel availability.
    """
    available = channel_availability(
        phobius, deeptmhmm, signalp, deeploc, deepcoil_segments
    )
    intervals = _group_intervals(ips.hits, options["merge_min_overlap"])
    repeat_intervals = _group_intervals(
        ips.hits, options["merge_min_overlap"], set(options["lrr_repeat_analyses"])
    )
    databases = _group_databases(ips.hits)
    hit_labels = _group_hit_labels(ips.hits)
    hit_accessions = _group_accessions(ips.hits)
    canonical_features = set(cfg.raw["integrated_domain_canonical_features"])
    canonical = {
        accession
        for accession, features in cfg.accession_to_features().items()
        if set(features) & canonical_features
    }
    canonical |= {str(a) for a in cfg.raw["integrated_domain_exclusions"]}
    integrated = _noncanonical_pfam(ips.domain_hits, canonical)
    domain_descriptions = _domain_descriptions(ips.domain_hits)
    # The Pfam hit table is only needed for the integrated-domain scan; release it
    # before building 300k evidence records.
    ips.domain_hits = ips.domain_hits.iloc[0:0]
    raw_cc = _raw_cc_by_protein(deepcoil_segments, deepcoil_id_map)

    lookups = {
        "phobius": _lookup_frame(phobius),
        "deeptmhmm": _lookup_frame(deeptmhmm),
        "signalp": _lookup_frame(signalp),
        "deeploc": _lookup_frame(deeploc),
    }
    on_progress(0, total=len(protein_ids))
    rows: list[dict] = []
    for index, protein_id in enumerate(protein_ids, start=1):
        rows.append(
            _protein_row(
                protein_id,
                cfg,
                options,
                intervals,
                repeat_intervals,
                databases,
                hit_labels,
                hit_accessions,
                integrated,
                raw_cc,
                lookups,
                ips,
                available,
            )
        )
        if index % PROTEIN_REPORT_INTERVAL == 0:
            on_progress(PROTEIN_REPORT_INTERVAL)
    on_progress(len(rows) % PROTEIN_REPORT_INTERVAL)
    contingency = _contingency(rows)
    long = _build_long(ips, rows, available)
    LOGGER.info(
        "Evidence: %d proteins; DeepCoil2/Coils CC agreement %s", len(rows), contingency
    )
    return Evidence(
        long=long,
        available=available,
        records=rows,
        cc_contingency=contingency,
        raw_cc=raw_cc,
        domain_descriptions=domain_descriptions,
    )


def _protein_row(
    protein_id: str,
    cfg: Config,
    options: dict,
    intervals: dict[str, dict[str, list[Interval]]],
    repeat_intervals: dict[str, dict[str, list[Interval]]],
    databases: dict[str, dict[str, int]],
    hit_labels: dict[str, dict[str, str]],
    hit_accessions: dict[str, dict[str, list[str]]],
    integrated: dict[str, list[str]],
    raw_cc: dict[str, list[tuple[int, int, float]]],
    lookups: dict[str, dict],
    ips,
    available: dict[str, bool],
) -> dict[str, object]:
    """Build the evidence record of a single protein."""
    domains = intervals.get(protein_id, {})
    row: dict[str, object] = {
        "protein_id": protein_id,
        "sequence_length": ips.sequence_lengths.get(protein_id),
    }
    row.update(
        _domain_fields(cfg, options, domains, repeat_intervals.get(protein_id, {}))
    )
    row.update(_sp_fields(cfg, options, lookups, protein_id, available))
    row.update(_tm_fields(cfg, options, lookups, protein_id, available, row))
    row.update(_cc_fields(cfg, options, raw_cc, protein_id, domains, available, row))
    row.update(_deeploc_fields(lookups, protein_id, available))
    row["noncanonical_domains"] = integrated.get(protein_id, [])
    row["feature_databases"] = databases.get(protein_id, {})
    row["feature_hit_labels"] = hit_labels.get(protein_id, {})
    # Fold the domain-level CC channel into `CC`: `CC_domain` is an internal
    # pseudo-feature and must not surface in the output as if it were a tenth
    # feature. DeepCoil2 contributes no accession -- it is a predictor, not a
    # signature -- so a CC called by DeepCoil2 alone has coordinates but no
    # accessions, which `cc_source` disambiguates.
    accessions = dict(hit_accessions.get(protein_id, {}))
    cc_accessions = sorted(
        set(accessions.pop(CC_DOMAIN_FEATURE, [])) | set(accessions.get("CC", []))
    )
    if cc_accessions:
        accessions["CC"] = cc_accessions
    row["feature_accessions"] = accessions
    row["feature_intervals"] = {
        **{feature: values for feature, values in domains.items() if feature != "CC"},
        "CC": row["cc_intervals"],
        "TM": row["tm_intervals"],
        "SP": [(1, int(row["sp_end"]))] if row["sp_end"] else [],
    }
    row["features"] = sorted(
        feature for feature in cfg.features if row.get(f"feat_{feature}", False)
    )
    return row


def _domain_fields(
    cfg: Config,
    options: dict,
    domains: dict[str, list[Interval]],
    repeats: dict[str, list[Interval]],
) -> dict:
    """Boolean domain features and the LRR copy number."""
    fields: dict[str, object] = {}
    n_lrr = len(domains.get("LRR", []))
    for feature in cfg.raw["interproscan_features"]:
        if feature == "CC":
            continue
        present = bool(domains.get(feature))
        if feature == "LRR":
            present = n_lrr >= options["min_lrr_copies"]
        fields[f"feat_{feature}"] = present
    fields["n_lrr"] = n_lrr
    fields["n_lrr_repeats"] = len(repeats.get("LRR", []))
    fields["cc_coils"] = bool(domains.get("CC"))
    fields["coils_intervals"] = domains.get("CC", [])
    fields["cc_rx_domain"] = bool(domains.get(CC_DOMAIN_FEATURE))
    fields["rx_domain_intervals"] = domains.get(CC_DOMAIN_FEATURE, [])
    return fields


def _sp_fields(
    cfg: Config,
    options: dict,
    lookups: dict,
    protein_id: str,
    available: dict[str, bool],
) -> dict:
    """Signal-peptide evidence and the configured consensus."""
    signalp = lookups["signalp"].get(protein_id, {}) if available["signalp"] else {}
    phobius = lookups["phobius"].get(protein_id, {}) if available["phobius"] else {}
    sp_signalp = signalp.get("sp_signalp") if available["signalp"] else None
    sp_phobius = phobius.get("sp_phobius") if available["phobius"] else None
    policy = {"signalp": "first", "phobius": "second"}.get(
        options["policies"]["sp"], options["policies"]["sp"]
    )
    consensus = apply_policy(policy, sp_signalp, sp_phobius)
    ends = [
        value
        for value in (signalp.get("sp_end_signalp"), phobius.get("sp_end_phobius"))
        if value is not None and not pd.isna(value)
    ]
    return {
        "sp_signalp": sp_signalp,
        "sp_phobius": sp_phobius,
        "sp_consensus": consensus,
        "feat_SP": consensus,
        "sp_prob": signalp.get("sp_prob"),
        "signalp_prediction": signalp.get("signalp_prediction"),
        "cleavage_site": signalp.get("cleavage_site"),
        "sp_end": max(ends) if ends else None,
    }


def _tm_fields(
    cfg: Config,
    options: dict,
    lookups: dict,
    protein_id: str,
    available: dict[str, bool],
    row: dict,
) -> dict:
    """Transmembrane evidence after removing signal-peptide artefacts."""
    phobius = lookups["phobius"].get(protein_id, {}) if available["phobius"] else {}
    deeptmhmm = (
        lookups["deeptmhmm"].get(protein_id, {}) if available["deeptmhmm"] else {}
    )
    sp_end = row.get("sp_end")
    ends = [sp_end, deeptmhmm.get("sp_end_deeptmhmm"), phobius.get("sp_end_phobius")]
    sp_end = max((e for e in ends if e is not None and not pd.isna(e)), default=None)

    kept_ph, dropped_ph = filter_helices_in_signal(
        phobius.get("tm_intervals_phobius") or [], sp_end, options["sp_overlap"]
    )
    kept_dt, dropped_dt = filter_helices_in_signal(
        deeptmhmm.get("tm_intervals_deeptmhmm") or [], sp_end, options["sp_overlap"]
    )
    minimum = options["min_tm_helices"]
    tm_ph = (len(kept_ph) >= minimum) if available["phobius"] else None
    tm_dt = (len(kept_dt) >= minimum) if available["deeptmhmm"] else None
    policy = {"phobius": "first", "deeptmhmm": "second"}.get(
        options["policies"]["tm"], options["policies"]["tm"]
    )
    consensus = apply_policy(policy, tm_ph, tm_dt)
    contributing = _contributing_helices(options["policies"]["tm"], kept_ph, kept_dt)
    return {
        "n_tm_phobius": len(kept_ph) if available["phobius"] else None,
        "n_tm_deeptmhmm": len(kept_dt) if available["deeptmhmm"] else None,
        "n_tm_phobius_raw": (
            int(phobius.get("n_tm_phobius", 0)) if available["phobius"] else None
        ),
        "n_tm_deeptmhmm_raw": (
            int(deeptmhmm.get("n_tm_deeptmhmm", 0)) if available["deeptmhmm"] else None
        ),
        "n_tm_dropped_in_sp": dropped_ph + dropped_dt,
        "tm_consensus": consensus,
        "feat_TM": consensus,
        "tm_intervals": merge_intervals(contributing) if contributing else [],
    }


def _contributing_helices(
    policy: str, phobius: list[Interval], deeptmhmm: list[Interval]
) -> list[Interval]:
    """Helices contributing coordinates under the requested TM policy."""
    if policy == "phobius":
        return list(phobius)
    if policy == "deeptmhmm":
        return list(deeptmhmm)
    if policy == "intersection":
        return [h for h in phobius if overlap_fraction(h, deeptmhmm) > 0]
    return [*phobius, *deeptmhmm]


def _cc_fields(
    cfg: Config,
    options: dict,
    raw_cc: dict[str, list[tuple[int, int, float]]],
    protein_id: str,
    domains: dict[str, list[Interval]],
    available: dict[str, bool],
    row: dict,
) -> dict:
    """Coiled-coil evidence from all three CC channels.

    The channels are, in descending order of evidential weight: ``rx_domain``
    (a curated profile HMM for a named domain), ``deepcoil`` (a learned
    per-residue propensity predictor) and ``coils`` (the 1991 Lupas algorithm).
    ``cc_source`` names every channel that fired, and the reported coordinates
    come from the most precise channel that did: DeepCoil2 resolves a segment
    per residue, whereas the domain and Coils channels report a signature span.
    """
    called = call_cc_segments(
        raw_cc.get(protein_id, []),
        options["cc_threshold"],
        options["cc_min_length"],
        options["cc_max_gap"],
    )
    cc_deepcoil = called.n > 0 if available["deepcoil"] else None
    cc_coils = bool(row.get("cc_coils"))
    cc_rx_domain = bool(row.get("cc_rx_domain"))
    channels: dict[str, bool | None] = {
        "rx_domain": cc_rx_domain,
        "deepcoil": cc_deepcoil,
        "coils": cc_coils,
    }
    consensus = apply_cc_policy(options["policies"]["cc"], channels)

    # A single contributor keeps the "<channel>_only" spelling the confidence
    # rules and the reports have always used; two or more are joined with "+"
    # in the fixed channel order above, strongest evidence first.
    contributing = [name for name, value in channels.items() if value]
    if not contributing:
        source = None
    elif len(contributing) == 1:
        source = f"{contributing[0]}_only"
    else:
        source = "+".join(contributing)
    if available["deepcoil"] and called.n:
        segments = called.intervals
    elif cc_rx_domain:
        segments = domains.get(CC_DOMAIN_FEATURE, [])
    elif cc_coils:
        segments = domains.get("CC", [])
    else:
        segments = []

    nbarc = domains.get("NB-ARC", [])
    nbarc_start = min((s for s, _ in nbarc), default=None)
    tm_intervals = row.get("tm_intervals") or []
    ambiguous = any(
        overlap_fraction(segment, tm_intervals) > options["cc_tm_overlap"]
        for segment in segments
    )
    return {
        "cc_deepcoil": cc_deepcoil,
        "cc_rx_domain": cc_rx_domain,
        "cc_consensus": consensus,
        "feat_CC": consensus,
        "cc_source": source,
        "n_cc_segments": called.n if available["deepcoil"] else None,
        "cc_max_prob": called.max_prob,
        "cc_mean_prob_in_segments": called.mean_prob,
        "cc_total_length": called.total_length,
        "cc_intervals": segments,
        # The long evidence table must attribute every interval to the channel
        # that produced it. ``cc_intervals`` is whichever channel called this
        # protein, so DeepCoil2's own segments are kept separately; otherwise a
        # CC called by the domain model alone would be exported as a DeepCoil2
        # hit, which is the provenance bug review finding 5 fixed for reasons.
        "deepcoil_intervals": called.intervals if available["deepcoil"] else [],
        "cc_is_n_terminal": (
            all(end < nbarc_start for _, end in segments)
            if segments and nbarc_start is not None
            else None
        ),
        "cc_any_n_terminal": (
            any(end < nbarc_start for _, end in segments)
            if segments and nbarc_start is not None
            else None
        ),
        "cc_tm_ambiguous": ambiguous,
    }


def _deeploc_fields(lookups: dict, protein_id: str, available: dict[str, bool]) -> dict:
    """Localisation columns (supporting evidence only)."""
    record = lookups["deeploc"].get(protein_id, {}) if available["deeploc"] else {}
    return {
        "predicted_localization": record.get("predicted_localization"),
        "localization_prob": record.get("localization_prob"),
        "all_localizations": record.get("all_localizations"),
    }


def _contingency(records: list[dict]) -> dict[str, int]:
    """Agreement between the CC channels.

    The four ``both``/``*_only``/``neither`` cells are the 2x2 table of the two
    *predictors*, kept comparable across configuration versions. The domain
    channel is counted alongside them rather than folded in, because it answers
    a different question -- how often a curated CC domain model fires where the
    predictors do or do not.
    """
    counts = {
        "both": 0,
        "deepcoil_only": 0,
        "coils_only": 0,
        "neither": 0,
        "rx_domain": 0,
        "rx_domain_only": 0,
    }
    for row in records:
        deepcoil = bool(row.get("cc_deepcoil"))
        coils = bool(row.get("cc_coils"))
        if row.get("cc_rx_domain"):
            counts["rx_domain"] += 1
            if not deepcoil and not coils:
                counts["rx_domain_only"] += 1
        key = (
            "both"
            if deepcoil and coils
            else "deepcoil_only"
            if deepcoil
            else "coils_only"
            if coils
            else "neither"
        )
        counts[key] += 1
    return counts


def _build_long(ips, records: list[dict], available: dict[str, bool]) -> pd.DataFrame:
    """Assemble the long evidence table from all contributing channels."""
    frames: list[pd.DataFrame] = []
    if not ips.hits.empty:
        hits = ips.hits.rename(columns={"score": "score_or_evalue"}).copy()
        hits["tool"] = "InterProScan"
        frames.append(hits[list(LONG_COLUMNS)])
    frames.extend(
        _interval_rows(records, column, tool, analysis, feature)
        for column, tool, analysis, feature in _INTERVAL_SOURCES
        if available[_INTERVAL_CHANNEL[column]]
    )
    long = pd.concat([f for f in frames if not f.empty], ignore_index=True)
    return long.sort_values(
        ["protein_id", "feature", "tool", "start"], kind="stable"
    ).reset_index(drop=True)


#: Wide-table interval columns exported to the long evidence table.
_INTERVAL_SOURCES: tuple[tuple[str, str, str, str], ...] = (
    ("deepcoil_intervals", "DeepCoil2", "deepcoil2", "CC"),
    ("tm_intervals", "Phobius/DeepTMHMM", "tm_consensus", "TM"),
)

_INTERVAL_CHANNEL = {"deepcoil_intervals": "deepcoil", "tm_intervals": "deeptmhmm"}


def _interval_rows(
    rows: list[dict], column: str, tool: str, analysis: str, feature: str
) -> pd.DataFrame:
    """Expand one interval field of the evidence records into long evidence rows."""
    records = []
    for row in rows:
        protein_id = row["protein_id"]
        for start, end in row.get(column) or []:
            records.append(
                {
                    "protein_id": protein_id,
                    "feature": feature,
                    "tool": tool,
                    "analysis": analysis,
                    "accession": "-",
                    "signature_description": f"{feature} segment",
                    "start": start,
                    "end": end,
                    "score_or_evalue": "-",
                }
            )
    if not records:
        return pd.DataFrame({name: pd.Series(dtype="object") for name in LONG_COLUMNS})
    return pd.DataFrame(records)[list(LONG_COLUMNS)]
