"""The progress-reporting contract between pipeline stages and the CLI.

The stages in this package must not import ``rich``: ``parsers``, ``evidence``
and ``rules`` are meant to be usable from a notebook, a test or another program
with no console attached. They therefore report progress through an opaque
callback, and the CLI in ``rgas_prediction.py`` is the only place that binds
that callback to a :class:`rich.progress.Progress` task.

A stage calls its callback in one of two ways::

    on_progress(0, total=len(sources))   # I now know how much work there is
    on_progress(n)                       # I have just completed n units

Both are optional and both are cheap: when no callback is supplied the stage
uses :func:`null_progress`, which does nothing.
"""

from __future__ import annotations

from typing import Protocol


class ProgressCallback(Protocol):
    """Callable used by a stage to report how far along it is.

    Parameters
    ----------
    advance : int
        Units of work completed since the last call. ``0`` is legitimate and is
        how a stage announces a total without claiming progress.
    total : int, optional
        Total units of work, when the stage discovers it at run time (the
        number of DeepCoil2 archives, say). ``None`` leaves the total unchanged.
    """

    def __call__(self, advance: int = 0, total: int | None = None) -> None:
        """Report progress."""


def null_progress(advance: int = 0, total: int | None = None) -> None:
    """Discard a progress report.

    The default for every ``on_progress`` parameter, so that stages need no
    ``if callback is not None`` guard in their hot loops.
    """


#: How many proteins a per-protein loop processes between progress reports.
#: Small enough that the bar moves smoothly on a 300k-protein proteome, large
#: enough that the callback is not measurable against the work it reports on.
PROTEIN_REPORT_INTERVAL = 2_000
