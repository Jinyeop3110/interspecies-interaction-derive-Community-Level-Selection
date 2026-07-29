"""Shared analysis library for the microbial community coalescence study.

Every figure script in ``figures/`` imports from here.  Nothing in this
package writes files or draws figures; keep plotting in the figure scripts.
"""

from . import (
    decomposition, io, lv, outcomes, processing, selection_correlation,
)

__version__ = "1.0.0"

__all__ = [
    "decomposition",
    "io",
    "lv",
    "outcomes",
    "processing",
    "selection_correlation",
]
