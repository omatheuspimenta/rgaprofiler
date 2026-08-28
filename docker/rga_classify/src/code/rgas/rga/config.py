"""Loading and validation of the RGA pipeline configuration.

The configuration is the single source of truth for every accession, threshold
and rule used by the pipeline. Nothing in the code may hard-code a biological
constant; this module is what makes that promise enforceable.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

#: Token that may appear inside a rule's ``any_of`` group and is expanded to the
#: configured ``ectodomain_features`` list at load time.
ECTODOMAIN_TOKEN = "ECTODOMAIN_FEATURES"

#: Normalisation operations understood by :func:`normalize_id`.
KNOWN_ID_OPS = frozenset(
    {"strip_after_whitespace", "strip_dots", "strip_after_pipe", "rstrip_suffixes"}
)

#: Reserved pseudo-feature carrying the domain-level coiled-coil channel.
#: It is deliberately *not* a member of the nine-feature controlled vocabulary:
#: it never reaches a rule, it only feeds the ``cc_rx_domain`` evidence channel.
#: Keeping it out of ``features`` is what stops it from enlarging the 2**9
#: exclusivity proof or creating a ``feat_`` column of its own.
CC_DOMAIN_FEATURE = "CC_domain"

#: Consensus policies accepted per channel.
VALID_POLICIES: dict[str, frozenset[str]] = {
    "tm": frozenset({"union", "intersection", "deeptmhmm", "phobius"}),
    "sp": frozenset({"union", "intersection", "signalp", "phobius"}),
    "cc": frozenset(
        {"rx_domain", "deepcoil", "coils", "union", "intersection"}
    ),
}


class ConfigError(ValueError):
    """Raised when the configuration file is missing a key or is inconsistent."""


@dataclass(frozen=True)
class Rule:
    """One classification rule.

    Parameters
    ----------
    id : str
        Stable rule identifier, written to ``rule_id`` in the output.
    priority : int
        Evaluation order; lower is evaluated first.
    family : str
        RGA family the rule assigns.
    subclass : str
        RGA subclass the rule assigns.
    description : str
        Human-readable statement of the rule, reused in the report.
    all_of : tuple of str
        Features that must all be present.
    none_of : tuple of str
        Features that must all be absent.
    any_of : tuple of tuple of str
        Groups of features; each group must contribute at least one feature.
    fallback : bool
        Ordered catch-all rule, excluded from the pairwise-disjointness proof.
    any_core : bool
        Fires when at least one core immune feature is present.
    no_core : bool
        Fires when no core immune feature is present.
    """

    id: str
    priority: int
    family: str
    subclass: str
    description: str
    all_of: tuple[str, ...] = ()
    none_of: tuple[str, ...] = ()
    any_of: tuple[tuple[str, ...], ...] = ()
    fallback: bool = False
    any_core: bool = False
    no_core: bool = False


@dataclass
class Config:
    """Fully resolved pipeline configuration.

    Attributes
    ----------
    raw : dict
        The resolved configuration mapping, written verbatim to
        ``run_metadata.json`` for reproducibility.
    rules : list of Rule
        Classification rules, sorted by ``priority``.
    """

    raw: dict[str, Any]
    rules: list[Rule] = field(default_factory=list)

    # -- convenience accessors -------------------------------------------------

    @property
    def features(self) -> list[str]:
        """Controlled feature vocabulary, in canonical display order."""
        return list(self.raw["features"])

    @property
    def core_immune_features(self) -> list[str]:
        """Features that make a protein an RGA candidate."""
        return list(self.raw["core_immune_features"])

    @property
    def ectodomain_features(self) -> list[str]:
        """Features used to sub-classify RLK/RLP receptors."""
        return list(self.raw["ectodomain_features"])

    @property
    def policies(self) -> dict[str, str]:
        """Consensus policy per channel (``tm``, ``sp``, ``cc``)."""
        return dict(self.raw["policies"])

    @property
    def coiled_coil(self) -> dict[str, Any]:
        """DeepCoil2 segment-calling parameters."""
        return dict(self.raw["coiled_coil"])

    @property
    def rga_families(self) -> set[str]:
        """Families reported with ``is_rga = True``."""
        return set(self.raw["rga_families"])

    @property
    def missing_value(self) -> str:
        """Placeholder written for missing values in every output file."""
        return str(self.raw["output"]["missing_value"])

    @property
    def list_separator(self) -> str:
        """Separator used inside multi-valued output cells."""
        return str(self.raw["output"]["list_separator"])

    def cc_domain_accessions(self) -> tuple[str, ...]:
        """Accessions of the domain-level coiled-coil channel."""
        return tuple(str(a) for a in self.raw.get("cc_domain_accessions", []))

    def accession_to_features(self) -> dict[str, tuple[str, ...]]:
        """Invert the feature -> accessions mapping.

        Returns
        -------
        dict
            Mapping from accession to the tuple of features it supports. An
            accession may legitimately support more than one feature.
        """
        inverted: dict[str, list[str]] = {}
        for feature, accessions in self.raw["interproscan_features"].items():
            for accession in accessions:
                inverted.setdefault(str(accession), []).append(feature)
        for accession in self.raw.get("cc_domain_accessions", []):
            inverted.setdefault(str(accession), []).append(CC_DOMAIN_FEATURE)
        return {acc: tuple(sorted(feats)) for acc, feats in sorted(inverted.items())}


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------

_REQUIRED_TOP_LEVEL = (
    "config_version",
    "ids",
    "features",
    "core_immune_features",
    "ectodomain_features",
    "interproscan_features",
    "excluded_analyses",
    "integrated_domain_analyses",
    "integrated_domain_canonical_features",
    "integrated_domain_exclusions",
    "cc_domain_accessions",
    "intervals",
    "coiled_coil",
    "policies",
    "signal_peptide",
    "transmembrane",
    "deeploc",
    "rules",
    "rga_families",
    "confidence",
    "output",
)


def load_config(path: str | Path) -> Config:
    """Read, validate and resolve the YAML configuration.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to ``rga_config.yaml``.

    Returns
    -------
    Config
        The resolved configuration.

    Raises
    ------
    ConfigError
        If a required key is missing, a policy is unknown, a rule references an
        unknown feature, or rule priorities/ids are not unique.
    """
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"configuration file not found: {path}")
    with path.open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ConfigError(f"configuration must be a mapping, got {type(raw).__name__}")

    missing = [key for key in _REQUIRED_TOP_LEVEL if key not in raw]
    if missing:
        raise ConfigError(f"missing configuration keys: {', '.join(missing)}")

    _validate_ids(raw["ids"])
    _validate_policies(raw["policies"])
    _validate_features(raw)
    _validate_cc_domain_accessions(raw)

    rules = _build_rules(raw)
    return Config(raw=copy.deepcopy(raw), rules=rules)


def _validate_ids(ids_cfg: dict[str, Any]) -> None:
    """Check that every declared ID normalisation operation is implemented."""
    for tool, ops in ids_cfg.get("per_tool", {}).items():
        unknown = sorted(set(ops) - KNOWN_ID_OPS)
        if unknown:
            raise ConfigError(f"unknown ID normalisation op(s) for {tool}: {unknown}")
    form = ids_cfg.get("deepcoil_canonical_form")
    if form is not None and form not in KNOWN_ID_OPS:
        raise ConfigError(f"unknown deepcoil_canonical_form: {form!r}")


def _validate_policies(policies: dict[str, str]) -> None:
    """Check that each consensus policy is one of the accepted values."""
    for channel, allowed in VALID_POLICIES.items():
        if channel not in policies:
            raise ConfigError(f"missing policy for channel {channel!r}")
        if policies[channel] not in allowed:
            raise ConfigError(
                f"invalid {channel} policy {policies[channel]!r}; "
                f"expected one of {sorted(allowed)}"
            )


def _validate_cc_domain_accessions(raw: dict[str, Any]) -> None:
    """Check the domain-level CC channel is a flat list of accession strings.

    The list may be empty -- that simply disables the channel and restores the
    two-predictor behaviour of config v1.0.0 -- but it may not overlap the
    accessions already assigned to a feature, because an accession that fed two
    CC channels at once would make ``--cc-policy intersection`` meaningless.
    """
    accessions = raw["cc_domain_accessions"]
    if not isinstance(accessions, list):
        raise ConfigError(
            "cc_domain_accessions must be a list, got "
            f"{type(accessions).__name__}"
        )
    assigned = {
        str(accession)
        for accessions_of_feature in raw["interproscan_features"].values()
        for accession in accessions_of_feature
    }
    clash = sorted({str(a) for a in accessions} & assigned)
    if clash:
        raise ConfigError(
            "cc_domain_accessions may not repeat an accession already used as "
            f"feature evidence: {clash}"
        )


def _validate_features(raw: dict[str, Any]) -> None:
    """Check that feature lists only reference the controlled vocabulary."""
    vocabulary = set(raw["features"])
    for key in (
        "core_immune_features",
        "ectodomain_features",
        "integrated_domain_canonical_features",
    ):
        unknown = sorted(set(raw[key]) - vocabulary)
        if unknown:
            raise ConfigError(f"{key} references unknown feature(s): {unknown}")
    unknown = sorted(set(raw["interproscan_features"]) - vocabulary)
    if unknown:
        raise ConfigError(
            f"interproscan_features references unknown feature(s): {unknown}"
        )


def _build_rules(raw: dict[str, Any]) -> list[Rule]:
    """Turn the YAML rule list into validated, sorted :class:`Rule` objects."""
    vocabulary = set(raw["features"])
    ectodomains = tuple(raw["ectodomain_features"])
    rules: list[Rule] = []

    for entry in raw["rules"]:
        for key in ("id", "priority", "family", "subclass"):
            if key not in entry:
                raise ConfigError(f"rule {entry.get('id', '?')!r} is missing {key!r}")

        any_of = tuple(
            tuple(ectodomains) if group == ECTODOMAIN_TOKEN else tuple(group)
            for group in entry.get("any_of", [])
        )
        rule = Rule(
            id=str(entry["id"]),
            priority=int(entry["priority"]),
            family=str(entry["family"]),
            subclass=str(entry["subclass"]),
            description=str(entry.get("description", "")),
            all_of=tuple(entry.get("all_of", [])),
            none_of=tuple(entry.get("none_of", [])),
            any_of=any_of,
            fallback=bool(entry.get("fallback", False)),
            any_core=bool(entry.get("any_core", False)),
            no_core=bool(entry.get("no_core", False)),
        )
        referenced = (
            set(rule.all_of) | set(rule.none_of) | {f for g in rule.any_of for f in g}
        )
        unknown = sorted(referenced - vocabulary)
        if unknown:
            raise ConfigError(
                f"rule {rule.id!r} references unknown feature(s): {unknown}"
            )
        rules.append(rule)

    _check_unique(rules)
    return sorted(rules, key=lambda r: r.priority)


def _check_unique(rules: list[Rule]) -> None:
    """Fail if rule ids or priorities are duplicated."""
    ids = [rule.id for rule in rules]
    priorities = [rule.priority for rule in rules]
    if len(set(ids)) != len(ids):
        duplicates = sorted({i for i in ids if ids.count(i) > 1})
        raise ConfigError(f"duplicate rule id(s): {duplicates}")
    if len(set(priorities)) != len(priorities):
        duplicates = sorted({p for p in priorities if priorities.count(p) > 1})
        raise ConfigError(f"duplicate rule priorities: {duplicates}")


# ---------------------------------------------------------------------------
# identifier normalisation
# ---------------------------------------------------------------------------


def normalize_id(raw_id: str, ops: list[str], suffixes: list[str] | None = None) -> str:
    """Apply the configured normalisation operations to a raw protein ID.

    Parameters
    ----------
    raw_id : str
        Identifier exactly as emitted by a tool.
    ops : list of str
        Operation names, applied in order. See :data:`KNOWN_ID_OPS`.
    suffixes : list of str, optional
        Suffixes removed by the ``rstrip_suffixes`` operation.

    Returns
    -------
    str
        The normalised identifier.

    Examples
    --------
    >>> normalize_id("PROT.1.p pacid=1 locus=x", ["strip_after_whitespace"])
    'PROT.1.p'
    >>> normalize_id("PROT.1.p", ["strip_dots"])
    'PROT1p'
    """
    value = raw_id.strip()
    for op in ops:
        if op == "strip_after_whitespace":
            value = value.split()[0] if value.split() else value
        elif op == "strip_after_pipe":
            value = value.split("|", 1)[0]
        elif op == "strip_dots":
            value = value.replace(".", "")
        elif op == "rstrip_suffixes":
            for suffix in suffixes or []:
                if value.endswith(suffix):
                    value = value[: -len(suffix)]
                    break
        else:  # pragma: no cover - guarded by _validate_ids
            raise ConfigError(f"unknown ID normalisation op: {op!r}")
    return value
