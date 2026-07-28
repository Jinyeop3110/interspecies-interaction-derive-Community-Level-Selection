"""Shared setup for the figure scripts.

Import this first in every ``make_panels.py``; it puts ``src/`` on the path so
the scripts run from their own directory without installation.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

#: Panel colours, matched to the published figures.
OUTCOME_COLORS = {"Dominance": "#bc5090", "Mixture": "#003f5c", "Restructuring": "#ffa600"}
MEDIUM_COLORS = {"Nutr-": "#bc5090", "Base": "#7a4a2b", "Nutr+": "#ffa600"}

MM = 1 / 25.4  # matplotlib works in inches; the figures are specified in mm


def use_paper_style():
    """Apply the rcParams used for the published panels."""
    mpl.rcParams.update(
        {
            "svg.fonttype": "none",  # keep text editable for figure assembly
            "font.size": 7,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.frameon": False,
        }
    )


def panel_dir():
    """Output directory for the calling figure, created if needed."""
    caller = Path(sys.argv[0]).resolve().parent
    out = caller / "panels"
    out.mkdir(exist_ok=True)
    return out


def save(fig, name):
    """Write one panel as SVG and report where it went."""
    path = panel_dir() / f"{name}.svg"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path.relative_to(REPO_ROOT)}")
    return path
