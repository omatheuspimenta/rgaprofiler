"""Ordered, mutually exclusive classification of proteins into RGA classes.

The rule list lives entirely in the YAML configuration. This module only knows
how to *evaluate* a rule, how to explain the resulting call in plain language,
and how to grade its confidence.

Mutual exclusivity is not a comment here: :func:`find_overlapping_rules`
enumerates every combination of the controlled feature vocabulary and is
asserted by the pipeline and by the test suite.
"""

from __future__ import annotations

import itertools
import logging
import re
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd

from .config import Config, Rule
from .progress import (
    PROTEIN_REPORT_INTERVAL,
    ProgressCallback,
    null_progress,
)

LOGGER = logging.getLogger(__name__)

#: Confidence levels, weakest first.
LEVELS: tuple[str, ...] = ("low", "medium", "high")


@dataclass(frozen=True)
class Call:
    """The outcome of evaluating the rule list for one protein."""

    rule_id: str
    rule_priority: int
    family: str
    subclass: str


# ---------------------------------------------------------------------------
# rule evaluation
# ---------------------------------------------------------------------------


def rule_matches(rule: Rule, features: frozenset[str], core: frozenset[str]) -> bool:
    """Test whether one rule fires for a given feature set.

    Parameters
    ----------
    rule : Rule
        The rule to evaluate.
    features : frozenset of str
        Features present on the protein.
    core : frozenset of str
        Features considered core immune evidence.

    Returns
    -------
    bool
        ``True`` when every clause of the rule is satisfied.
    """
    if rule.any_core and not (features & core):
        return False
    if rule.no_core and (features & core):
        return False
    if not set(rule.all_of) <= features:
        return False
    if set(rule.none_of) & features:
        return False
    return all(set(group) & features for group in rule.any_of)


def classify_features(
    rules: Sequence[Rule], features: frozenset[str], core: frozenset[str]
) -> Call:
    """Return the first rule that fires, evaluating in priority order.

    Raises
    ------
    RuntimeError
        If no rule fires. The configuration must end with a catch-all.
    """
    for rule in rules:
        if rule_matches(rule, features, core):
            return Call(rule.id, rule.priority, rule.family, rule.subclass)
    raise RuntimeError(f"no rule matched feature set {sorted(features)}")


def find_overlapping_rules(cfg: Config) -> list[tuple[str, str, tuple[str, ...]]]:
    """Enumerate feature combinations matched by more than one non-fallback rule.

    Every one of the ``2 ** len(features)`` subsets of the controlled vocabulary
    is tested. Fallback rules (the ordered catch-alls ``Other`` and ``Non-RGA``)
    are excluded because they overlap the specific rules by construction.

    Parameters
    ----------
    cfg : Config
        Resolved pipeline configuration.

    Returns
    -------
    list of tuple
        ``(rule_a, rule_b, features)`` triples; empty when the rule set is
        mutually exclusive.
    """
    specific = [rule for rule in cfg.rules if not rule.fallback]
    core = frozenset(cfg.core_immune_features)
    clashes: list[tuple[str, str, tuple[str, ...]]] = []
    vocabulary = cfg.features
    for size in range(len(vocabulary) + 1):
        for combination in itertools.combinations(vocabulary, size):
            features = frozenset(combination)
            fired = [rule.id for rule in specific if rule_matches(rule, features, core)]
            clashes.extend(
                (first, second, combination)
                for first, second in itertools.combinations(fired, 2)
            )
    return clashes


def assert_mutually_exclusive(cfg: Config) -> None:
    """Fail loudly if two non-fallback rules can fire on the same feature set."""
    clashes = find_overlapping_rules(cfg)
    if clashes:
        preview = "; ".join(f"{a} vs {b} on {list(f)}" for a, b, f in clashes[:5])
        raise AssertionError(
            f"classification rules are not mutually exclusive ({len(clashes)} clashes): {preview}"
        )


# ---------------------------------------------------------------------------
# domain architecture
# ---------------------------------------------------------------------------


def domain_architecture(
    feature_intervals: dict[str, list[tuple[int, int]]], features: Iterable[str]
) -> str:
    """Render the N-terminal to C-terminal domain architecture of a protein.

    Parameters
    ----------
    feature_intervals : dict
        Feature -> merged, sorted intervals.
    features : iterable of str
        Features called on the protein; features without coordinates are
        appended so that they are never silently lost.

    Returns
    -------
    str
        For example ``SP-LRR-TM-STTK`` or ``CC-NB-ARC-LRR``. Consecutive
        repeats of the same feature are collapsed.
    """
    placed: list[tuple[int, str]] = []
    present = set(features)
    for feature, intervals in feature_intervals.items():
        if feature not in present:
            continue
        placed.extend((int(start), feature) for start, _ in intervals or [])
    ordered = [
        feature for _, feature in sorted(placed, key=lambda item: (item[0], item[1]))
    ]
    collapsed = [f for i, f in enumerate(ordered) if i == 0 or f != ordered[i - 1]]
    positioned = {feature for feature in collapsed}
    collapsed.extend(sorted(present - positioned))
    return "-".join(collapsed)


# ---------------------------------------------------------------------------
# confidence
# ---------------------------------------------------------------------------


def _demote(level: str, steps: int) -> str:
    """Lower a confidence level by ``steps``, floored at the weakest level."""
    index = max(0, LEVELS.index(level) - steps)
    return LEVELS[index]


def grade_confidence(
    cfg: Config, call: Call, row: Mapping[str, Any], available: dict[str, bool]
) -> tuple[str, list[str]]:
    """Grade the confidence of one call and collect its triggered demotions.

    Every protein starts at the configured default level and is demoted once per
    triggered rule, floored at ``low``. The formula is fully described in
    ``docs/rga/README.md``.

    Returns
    -------
    tuple
        ``(level, triggered_rule_ids)``.
    """
    settings = cfg.raw["confidence"]
    level = settings["default"]
    uses_cc = call.subclass in settings["classes_using_cc"]
    uses_tm_sp = call.subclass in settings["classes_using_tm_sp"]
    triggered: list[str] = []

    for demotion in settings["demotions"]:
        when = demotion["when"]
        if when.get("classes_using_cc") and not uses_cc:
            continue
        if not _demotion_applies(when, row, cfg, call, available, uses_cc, uses_tm_sp):
            continue
        level = _demote(level, int(demotion["steps"]))
        triggered.append(demotion["id"])
    return level, triggered


#: ``when`` keys evaluated as a plain boolean comparison against the evidence
#: record. Listing them explicitly (rather than treating any unknown key as a
#: column name) means a typo in the configuration raises instead of silently
#: matching everything.
_BOOLEAN_WHEN_KEYS: tuple[str, ...] = (
    "cc_rx_domain",
    "cc_deepcoil",
    "cc_coils",
    "cc_tm_ambiguous",
)


def _demotion_applies(
    when: dict,
    row: Mapping[str, Any],
    cfg: Config,
    call: Call,
    available: dict[str, bool],
    uses_cc: bool,
    uses_tm_sp: bool,
) -> bool:
    """Evaluate the ``when`` clause of a single confidence demotion."""
    if "cc_source" in when and row.get("cc_source") != when["cc_source"]:
        return False
    for key in _BOOLEAN_WHEN_KEYS:
        if key in when and bool(row.get(key)) != when[key]:
            return False
    if "cc_is_n_terminal" in when:
        value = row.get("cc_is_n_terminal")
        if value is None or pd.isna(value) or bool(value) != when["cc_is_n_terminal"]:
            return False
    if "defining_domain_databases_lt" in when:
        databases = row.get("feature_databases") or {}
        defining = _defining_feature(cfg, call)
        if defining is None:
            return False
        if databases.get(defining, 0) >= when["defining_domain_databases_lt"]:
            return False
    if "deeploc_inconsistent" in when:
        if _deeploc_inconsistent(cfg, call, row) != when["deeploc_inconsistent"]:
            return False
    if "call_used_missing_channel" in when:
        missing = (uses_cc and not available["deepcoil"]) or (
            uses_tm_sp
            and not (
                available["phobius"] or available["deeptmhmm"] or available["signalp"]
            )
        )
        if missing != when["call_used_missing_channel"]:
            return False
    return True


def _defining_feature(cfg: Config, call: Call) -> str | None:
    """Return the domain feature that defines a class, if any."""
    for rule in cfg.rules:
        if rule.id != call.rule_id:
            continue
        for feature in rule.all_of:
            if feature in cfg.raw["interproscan_features"] and feature != "CC":
                return feature
    return None


def _deeploc_inconsistent(cfg: Config, call: Call, row: Mapping[str, Any]) -> bool:
    """Test whether the DeepLoc localisation contradicts the assigned class."""
    localization = row.get("predicted_localization")
    if not localization or (isinstance(localization, float) and pd.isna(localization)):
        return False
    settings = cfg.raw["deeploc"]
    labels = set(str(localization).split("|"))
    if call.family in {"NLR", "NLR-associated"}:
        return bool(labels & set(settings["nlr_inconsistent_localizations"]))
    if call.family in {"RLK", "RLP", "TM-CC"}:
        return not (labels & set(settings["membrane_localizations"]))
    return False


# ---------------------------------------------------------------------------
# reason strings
# ---------------------------------------------------------------------------


def build_reason(
    cfg: Config,
    call: Call,
    row: Mapping[str, Any],
    confidence: str,
    warnings: Sequence[str],
    available: dict[str, bool],
) -> str:
    """Compose the human-readable justification of one call.

    The string is generated entirely from the evidence table -- there is no
    per-class template -- so it always reflects what the pipeline actually saw.

    Returns
    -------
    str
        A sentence a biologist can read without opening the code, citing the
        exact signatures, coordinates and tool calls behind the decision.
    """
    rule = next(r for r in cfg.rules if r.id == call.rule_id)
    parts = [
        f"Rule {call.rule_id} (priority {call.rule_priority}): {_positive_clause(row, rule)}"
    ]
    absent = _absent_clause(rule)
    if absent:
        parts.append(absent)
    parts.append(_tm_clause(row, available))
    parts.append(_sp_clause(row, available))
    parts.append(_cc_clause(row, available))
    parts.append(_localization_clause(cfg, call, row, available))
    if warnings:
        parts.append("Warnings: " + "; ".join(warnings) + ".")
    parts.append(f"Confidence: {confidence}.")
    return " ".join(part for part in parts if part).replace("  ", " ")


def _positive_clause(row: Mapping[str, Any], rule: Rule) -> str:
    """List the features that made the rule fire, with their provenance."""
    features = list(rule.all_of) + [f for group in rule.any_of for f in group]
    if rule.any_core:
        features = list(row.get("features") or [])
    shown = [
        _feature_label(row, feature)
        for feature in dict.fromkeys(features)
        if row.get(f"feat_{feature}", False)
    ]
    if shown:
        return ", ".join(shown) + "."
    present = row.get("features") or []
    if present:
        return (
            "no core immune feature; the only features detected were "
            + ", ".join(present)
            + "."
        )
    return "no protein feature of interest was detected."


def _feature_label(row: Mapping[str, Any], feature: str) -> str:
    """Render one feature with the coordinates of the channel that called it."""
    if feature == "CC":
        coords = _format_intervals(row.get("cc_intervals"))
        source = row.get("cc_source") or "unknown source"
        return f"CC [{source} @ {coords}]"
    if feature == "TM":
        return f"TM [{_format_intervals(row.get('tm_intervals'))}]"
    if feature == "SP":
        site = row.get("cleavage_site")
        return f"SP [cleavage site {site}]" if site else "SP"
    labels = row.get("feature_hit_labels") or {}
    return f"{feature} [{labels[feature]}]" if feature in labels else feature


def _absent_clause(rule: Rule) -> str:
    """State the features the rule required to be absent."""
    if not rule.none_of:
        return ""
    return "Excluded: " + ", ".join(f"no {feature}" for feature in rule.none_of) + "."


def _tm_clause(row: Mapping[str, Any], available: dict[str, bool]) -> str:
    """Summarise the transmembrane channel."""
    phobius = _count(row, "n_tm_phobius", available["phobius"])
    deeptmhmm = _count(row, "n_tm_deeptmhmm", available["deeptmhmm"])
    state = "present" if row.get("tm_consensus") else "none"
    dropped = int(row.get("n_tm_dropped_in_sp") or 0)
    extra = (
        f", {dropped} helix/helices discarded inside the signal peptide"
        if dropped
        else ""
    )
    return f"TM: {state} (Phobius {phobius} / DeepTMHMM {deeptmhmm}{extra})."


def _sp_clause(row: Mapping[str, Any], available: dict[str, bool]) -> str:
    """Summarise the signal-peptide channel."""
    if not available["signalp"]:
        prediction = "SignalP6 unavailable"
    else:
        prediction = f"SignalP6 {row.get('signalp_prediction') or 'NO_SP'}"
        probability = row.get("sp_prob")
        if probability is not None and not pd.isna(probability):
            prediction += f" {float(probability):.3f}"
        site = row.get("cleavage_site")
        if site and not (isinstance(site, float) and pd.isna(site)):
            prediction += f", CS {site}"
    state = "present" if row.get("sp_consensus") else "none"
    return f"SP: {state} ({prediction})."


def _cc_clause(row: Mapping[str, Any], available: dict[str, bool]) -> str:
    """Summarise the coiled-coil channel, including the positional check."""
    if not available["deepcoil"]:
        source = "DeepCoil2 unavailable"
    else:
        n_segments = row.get("n_cc_segments") or 0
        best = row.get("cc_max_prob")
        source = f"DeepCoil2 {int(n_segments)} segment(s)"
        if best is not None and not pd.isna(best):
            source += f", max score {float(best):.3f}"
    coils = "Coils yes" if row.get("cc_coils") else "Coils no"
    domain = "Rx domain yes" if row.get("cc_rx_domain") else "Rx domain no"
    coils = f"{domain}; {coils}"
    state = "present" if row.get("cc_consensus") else "none"
    positional = row.get("cc_is_n_terminal")
    if positional is not None and not (
        isinstance(positional, float) and pd.isna(positional)
    ):
        positional = (
            " CC is N-terminal to NB-ARC."
            if positional
            else " CC lies C-terminal to NB-ARC (atypical)."
        )
    else:
        positional = ""
    coords = _format_intervals(row.get("cc_intervals"))
    coords = f" at {coords}" if coords != "NA" else ""
    return f"CC: {state} ({source}; {coils}){coords}.{positional}"


def _localization_clause(
    cfg: Config, call: Call, row: Mapping[str, Any], available: dict[str, bool]
) -> str:
    """Summarise the DeepLoc localisation and whether it fits the class."""
    if not available["deeploc"]:
        return "DeepLoc: unavailable."
    localization = row.get("predicted_localization")
    if not localization or (isinstance(localization, float) and pd.isna(localization)):
        return "DeepLoc: no prediction."
    probability = row.get("localization_prob")
    probability_text = (
        f" ({float(probability):.2f})"
        if probability is not None and not pd.isna(probability)
        else ""
    )
    verdict = "inconsistent" if _deeploc_inconsistent(cfg, call, row) else "consistent"
    return f"DeepLoc: {localization}{probability_text} -- {verdict}."


def _count(row: Mapping[str, Any], column: str, available: bool) -> str:
    """Render a helix count, or ``n/a`` when the tool was not supplied."""
    if not available:
        return "n/a"
    value = row.get(column)
    return "0" if value is None or pd.isna(value) else str(int(value))


def _format_intervals(intervals) -> str:
    """Render a list of intervals as ``start-end,start-end``."""
    if not intervals:
        return "NA"
    return ",".join(f"{int(start)}-{int(end)}" for start, end in intervals)


# ---------------------------------------------------------------------------
# warnings
# ---------------------------------------------------------------------------


def collect_warnings(
    cfg: Config, call: Call, row: Mapping[str, Any], available: dict[str, bool]
) -> list[str]:
    """Collect every caveat attached to one call, in a stable order."""
    warnings: list[str] = []
    uses_cc = call.subclass in cfg.raw["confidence"]["classes_using_cc"]
    if uses_cc and row.get("cc_tm_ambiguous"):
        warnings.append("CC segment overlaps a predicted TM helix (possible artefact)")
    if (
        uses_cc
        and row.get("cc_coils")
        and not row.get("cc_rx_domain")
        and not row.get("cc_deepcoil")
    ):
        warnings.append("CC supported only by InterProScan Coils")
    if call.subclass == "CNL" and row.get("cc_is_n_terminal") is False:
        warnings.append("CC lies C-terminal to the NB-ARC domain")
    if _deeploc_inconsistent(cfg, call, row):
        warnings.append(
            f"DeepLoc localisation ({row.get('predicted_localization')}) "
            f"is inconsistent with class {call.subclass}"
        )
    missing = [name for name, ok in available.items() if not ok]
    if missing:
        warnings.append(
            "evidence channel(s) unavailable: " + ", ".join(sorted(missing))
        )
    if int(row.get("n_tm_dropped_in_sp") or 0) > 0:
        warnings.append(
            f"{int(row['n_tm_dropped_in_sp'])} TM helix/helices discarded as signal peptide"
        )
    return warnings


# ---------------------------------------------------------------------------
# whole-proteome classification
# ---------------------------------------------------------------------------

#: Column order of ``rga_predictions.tsv``.
PREDICTION_COLUMNS: tuple[str, ...] = (
    "protein_id",
    "locus",
    "sequence_length",
    "is_rga",
    "rga_family",
    "rga_subclass",
    "domain_architecture",
    "features_found",
    "feature_coords",
    "feature_accessions",
    "n_lrr",
    "n_lrr_repeats",
    "defining_domain_databases",
    "n_tm_phobius",
    "n_tm_deeptmhmm",
    "n_tm_phobius_raw",
    "n_tm_deeptmhmm_raw",
    "n_tm_dropped_in_sp",
    "n_tm_consensus",
    "tm_consensus",
    "sp_signalp",
    "sp_phobius",
    "sp_consensus",
    "signalp_prediction",
    "sp_prob",
    "cleavage_site",
    "cc_deepcoil",
    "cc_coils",
    "cc_rx_domain",
    "cc_consensus",
    "cc_source",
    "n_cc_segments",
    "cc_max_prob",
    "cc_mean_prob_in_segments",
    "cc_total_length",
    "cc_coords",
    "cc_is_n_terminal",
    "cc_tm_ambiguous",
    "predicted_localization",
    "localization_prob",
    "all_localizations",
    "has_integrated_domain",
    "integrated_domains",
    "integrated_domain_descriptions",
    "rule_id",
    "rule_priority",
    "reason",
    "confidence",
    "confidence_demotions",
    "warnings",
    "evidence_tools_available",
)


def classify_proteome(
    cfg: Config,
    evidence,
    options: dict,
    on_progress: ProgressCallback = null_progress,
) -> pd.DataFrame:
    """Classify every protein and assemble the prediction table.

    Parameters
    ----------
    cfg : Config
        Resolved pipeline configuration.
    evidence : rga.evidence.Evidence
        Harmonised evidence for the whole proteome.
    options : dict
        Resolved run options; ``locus_regex`` and ``list_separator`` are used
        here.
    on_progress : ProgressCallback, optional
        Advanced every :data:`~rga.progress.PROTEIN_REPORT_INTERVAL` proteins;
        the total is the number of evidence records.

    Returns
    -------
    pandas.DataFrame
        One row per protein, columns as in :data:`PREDICTION_COLUMNS`, sorted by
        ``protein_id``.
    """
    core = frozenset(cfg.core_immune_features)
    rules = cfg.rules
    separator = cfg.list_separator
    available = evidence.available
    tools = separator.join(sorted(name for name, ok in available.items() if ok))
    locus_pattern = options.get("locus_regex")
    descriptions = getattr(evidence, "domain_descriptions", {}) or {}

    on_progress(0, total=len(evidence.records))
    records: list[dict[str, object]] = []
    for index, series in enumerate(evidence.records, start=1):
        features = frozenset(series["features"])
        call = classify_features(rules, features, core)
        confidence, demotions = grade_confidence(cfg, call, series, available)
        warnings = collect_warnings(cfg, call, series, available)
        reason = build_reason(cfg, call, series, confidence, warnings, available)
        records.append(
            _prediction_record(
                cfg,
                call,
                series,
                confidence,
                demotions,
                warnings,
                reason,
                tools,
                separator,
                locus_pattern,
                descriptions,
            )
        )
        if index % PROTEIN_REPORT_INTERVAL == 0:
            on_progress(PROTEIN_REPORT_INTERVAL)
    on_progress(len(records) % PROTEIN_REPORT_INTERVAL)
    frame = pd.DataFrame(records, columns=list(PREDICTION_COLUMNS))
    frame = _tidy_dtypes(frame)
    return frame.sort_values("protein_id", kind="stable").reset_index(drop=True)


#: Columns written as integers (nullable) rather than floats.
_INTEGER_COLUMNS: tuple[str, ...] = (
    "sequence_length",
    "n_lrr",
    "n_lrr_repeats",
    "defining_domain_databases",
    "n_tm_phobius",
    "n_tm_deeptmhmm",
    "n_tm_phobius_raw",
    "n_tm_deeptmhmm_raw",
    "n_tm_dropped_in_sp",
    "n_tm_consensus",
    "n_cc_segments",
    "cc_total_length",
    "rule_priority",
)

#: Columns rounded before being written, with the number of decimals.
_ROUNDED_COLUMNS: dict[str, int] = {
    "sp_prob": 6,
    "cc_max_prob": 3,
    "cc_mean_prob_in_segments": 3,
    "localization_prob": 4,
}


def _tidy_dtypes(frame: pd.DataFrame) -> pd.DataFrame:
    """Write counts as integers and probabilities at a sensible precision.

    Without this, a column holding ``None`` for a few proteins is promoted to
    float and a sequence length of 447 is written as ``447.0``.
    """
    for column in _INTEGER_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").astype("Int64")
    for column, decimals in _ROUNDED_COLUMNS.items():
        frame[column] = pd.to_numeric(frame[column], errors="coerce").round(decimals)
    return frame


def _format_feature_map(
    values: Mapping[str, Any], features: Sequence[str], render
) -> str | None:
    """Render a feature -> value mapping as ``FEATURE:value;FEATURE:value``.

    Only features the protein actually carries are emitted, in the same order as
    ``features_found``, so the column can be split on ``;`` and then on the first
    ``:`` without ambiguity. Returns ``None`` (written as ``NA``) when empty.

    Note for ``feature_accessions``: the values are **signature** accessions
    (InterProScan column 5). A hit whose InterPro accession also maps to the
    feature is still recorded once, under the signature that produced it -- the
    same convention ``reason`` has always used.
    """
    parts = []
    for feature in features:
        rendered = render(values.get(feature))
        if rendered and rendered != "NA":
            parts.append(f"{feature}:{rendered}")
    return ";".join(parts) if parts else None


def _prediction_record(
    cfg: Config,
    call: Call,
    row: Mapping[str, Any],
    confidence: str,
    demotions: Sequence[str],
    warnings: Sequence[str],
    reason: str,
    tools: str,
    separator: str,
    locus_pattern: str | None,
    descriptions: Mapping[str, str],
) -> dict[str, object]:
    """Assemble one output row from a call and its evidence."""
    integrated = row.get("noncanonical_domains") or []
    is_nlr = call.family in {"NLR", "NLR-associated"}
    return {
        "protein_id": row["protein_id"],
        "locus": _locus(row["protein_id"], locus_pattern),
        "sequence_length": row.get("sequence_length"),
        "is_rga": call.family in cfg.rga_families,
        "rga_family": call.family,
        "rga_subclass": call.subclass,
        "domain_architecture": domain_architecture(
            row.get("feature_intervals") or {}, row.get("features") or []
        ),
        "features_found": separator.join(row.get("features") or []),
        "feature_coords": _format_feature_map(
            row.get("feature_intervals") or {},
            row.get("features") or [],
            _format_intervals,
        ),
        "feature_accessions": _format_feature_map(
            row.get("feature_accessions") or {},
            row.get("features") or [],
            lambda values: ",".join(values) if values else "",
        ),
        "defining_domain_databases": (row.get("feature_databases") or {}).get(
            _defining_feature(cfg, call)
        ),
        "n_lrr": row.get("n_lrr"),
        "n_lrr_repeats": row.get("n_lrr_repeats"),
        "n_tm_phobius": row.get("n_tm_phobius"),
        "n_tm_deeptmhmm": row.get("n_tm_deeptmhmm"),
        "n_tm_phobius_raw": row.get("n_tm_phobius_raw"),
        "n_tm_deeptmhmm_raw": row.get("n_tm_deeptmhmm_raw"),
        "n_tm_dropped_in_sp": row.get("n_tm_dropped_in_sp"),
        "n_tm_consensus": len(row.get("tm_intervals") or []),
        "tm_consensus": row.get("tm_consensus"),
        "sp_signalp": row.get("sp_signalp"),
        "sp_phobius": row.get("sp_phobius"),
        "sp_consensus": row.get("sp_consensus"),
        "signalp_prediction": row.get("signalp_prediction"),
        "sp_prob": row.get("sp_prob"),
        "cleavage_site": row.get("cleavage_site"),
        "cc_deepcoil": row.get("cc_deepcoil"),
        "cc_coils": row.get("cc_coils"),
        "cc_rx_domain": row.get("cc_rx_domain"),
        "cc_consensus": row.get("cc_consensus"),
        "cc_source": row.get("cc_source"),
        "n_cc_segments": row.get("n_cc_segments"),
        "cc_max_prob": row.get("cc_max_prob"),
        "cc_mean_prob_in_segments": row.get("cc_mean_prob_in_segments"),
        "cc_total_length": row.get("cc_total_length"),
        "cc_coords": _format_intervals(row.get("cc_intervals")),
        "cc_is_n_terminal": row.get("cc_is_n_terminal"),
        "cc_tm_ambiguous": row.get("cc_tm_ambiguous"),
        "predicted_localization": row.get("predicted_localization"),
        "localization_prob": row.get("localization_prob"),
        "all_localizations": row.get("all_localizations"),
        "has_integrated_domain": bool(integrated) and is_nlr,
        "integrated_domains": separator.join(integrated)
        if (integrated and is_nlr)
        else None,
        "integrated_domain_descriptions": separator.join(
            descriptions.get(accession, accession) for accession in integrated
        )
        if (integrated and is_nlr)
        else None,
        "rule_id": call.rule_id,
        "rule_priority": call.rule_priority,
        "reason": reason,
        "confidence": confidence,
        "confidence_demotions": separator.join(demotions) if demotions else None,
        "warnings": separator.join(warnings) if warnings else None,
        "evidence_tools_available": tools,
    }


def _locus(protein_id: str, pattern: str | None) -> str | None:
    """Extract the locus identifier from a protein ID using the configured regex."""
    if not pattern:
        return None
    match = re.match(pattern, protein_id)
    return match.group(1) if match else None


def cc_policy_sensitivity(cfg: Config, evidence) -> pd.DataFrame:
    """Recount every subclass under each ``--cc-policy`` setting.

    Only the CC channel is recomputed; all other evidence is held fixed. This
    answers the question "how sensitive are the CNL/CN/RNL/TM-CC counts to the
    coiled-coil consensus policy?" without re-running the whole pipeline.

    Since config v1.1.0 there are three CC channels and therefore five policies.
    ``union`` and ``intersection`` now range over all three, so their columns are
    **not** comparable with the same-named columns of a v1.0.0 run -- which is
    exactly why the policy is recorded in ``run_metadata.json`` alongside the
    configuration version.

    Parameters
    ----------
    cfg : Config
        Resolved pipeline configuration.
    evidence : rga.evidence.Evidence
        Harmonised evidence for the whole proteome.

    Returns
    -------
    pandas.DataFrame
        Subclasses as rows, policies as columns, counts as values.
    """
    from .evidence import apply_cc_policy

    core = frozenset(cfg.core_immune_features)
    base_features = [set(row["features"]) for row in evidence.records]
    channels = [
        {
            "rx_domain": bool(row.get("cc_rx_domain")),
            "deepcoil": row.get("cc_deepcoil"),
            "coils": bool(row.get("cc_coils")),
        }
        for row in evidence.records
    ]

    counts: dict[str, dict[str, int]] = {}
    for policy in ("rx_domain", "deepcoil", "coils", "union", "intersection"):
        tally: dict[str, int] = {}
        for features, row_channels in zip(base_features, channels):
            has_cc = apply_cc_policy(policy, row_channels)
            adjusted = (features | {"CC"}) if has_cc else (features - {"CC"})
            call = classify_features(cfg.rules, frozenset(adjusted), core)
            tally[call.subclass] = tally.get(call.subclass, 0) + 1
        counts[policy] = tally

    subclasses = sorted({s for tally in counts.values() for s in tally})
    return pd.DataFrame(
        {policy: [counts[policy].get(s, 0) for s in subclasses] for policy in counts},
        index=pd.Index(subclasses, name="rga_subclass"),
    )
