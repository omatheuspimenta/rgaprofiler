"""Genome-wide prediction and classification of Resistance Gene Analogs (RGAs).

The package turns pre-computed protein-annotation outputs (InterProScan,
Phobius, DeepTMHMM, SignalP 6.0, DeepLoc 2.0, DeepCoil2) into a per-protein
RGA call with an explicit, human-readable justification.

Modules
-------
config
    Loading and validation of the YAML configuration.
parsers
    One parser per tool, each returning a tidy :class:`pandas.DataFrame`.
evidence
    Raw annotations -> controlled feature vocabulary.
rules
    Ordered, mutually exclusive classification rules.
report
    HTML and Markdown report generation.
"""

__version__ = "0.0.1"

__all__ = ["__version__"]
