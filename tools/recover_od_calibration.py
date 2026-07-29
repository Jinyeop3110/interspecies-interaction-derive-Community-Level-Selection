#!/usr/bin/env python3
"""Recover the OD calibration curve without MATLAB.

Plate-reader absorbance saturates at high cell density, so raw readings are
passed through a calibration curve before use.  In the original pipeline that
curve is a Curve Fitting Toolbox ``cfit`` object stored in
``OD_correction_function_2202.mat``.  That file holds a serialized MATLAB class
instance (MCOS): its bytes survive a round trip through scipy, but the object
cannot be reconstructed outside MATLAB, so the curve is not directly portable.

It does not need to be.  The curve is a smooth, deterministic, monotone
transform, and the archived data contains both sides of it — raw plate readings
on one side, corrected OD in ``Metadata.xlsx`` on the other.  Fitting a
polynomial through those pairs recovers the transform to well below plate-reader
precision, with no MATLAB dependency.

    python tools/recover_od_calibration.py \
        --raw-root /path/to/working/tree \
        --degree 9

Prints coefficients suitable for
``coalescence.processing.measurements.RECOVERED_OD_CORRECTION`` and reports the
accuracy achieved on held-out points.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from coalescence.processing.measurements import (  # noqa: E402
    read_plate_block, reshape_parental_plate,
)

#: Parental-community OD workbooks, by medium code.
OD_FILES = {"L": "LN_SC_OD.xlsx", "M": "MN_SC_OD.xlsx", "H": "HN_SC_OD.xlsx"}

#: Where the OD workbooks sit under the working tree.
OD_SUBDIR = Path("ExperimentalResult/Data/2208_Coalescence_processed/OD")


def collect_pairs(raw_root):
    """Gather (raw noise-subtracted, corrected) reading pairs.

    Raw readings come from the plate workbooks with the same noise subtraction
    the MATLAB applies; corrected readings come from the archived metadata.
    """
    metadata = pd.read_excel(raw_root / "Postprocessed" / "Metadata.xlsx")
    raw_values, corrected_values = [], []

    for medium, filename in OD_FILES.items():
        path = raw_root / OD_SUBDIR / filename
        for sheet in range(7):
            block = read_plate_block(path, sheet)
            noise = np.nanmean(block[block < 1.05 * np.nanmin(block)])
            reshaped = reshape_parental_plate(block - noise)

            for replicate in (1, 2):
                raw = np.concatenate(
                    [reshaped[k][:, replicate - 1] for k in ("s6", "s12", "s24")]
                )
                selection = metadata[
                    (metadata.Timepoint == "F")
                    & (metadata.CommunityOrigin == "S")
                    & (metadata.Medium == medium)
                    & (metadata.CoalescenceType == "S")
                    & (metadata.Replicate == replicate)
                ]
                corrected = selection[f"fieldOD{sheet + 1}"].to_numpy(dtype=float)
                count = min(len(raw), len(corrected))
                raw_values.append(raw[:count])
                corrected_values.append(corrected[:count])

    raw = np.concatenate(raw_values)
    corrected = np.concatenate(corrected_values)
    usable = np.isfinite(raw) & np.isfinite(corrected)
    return raw[usable], corrected[usable]


def fit(raw, corrected, degree):
    """Fit the curve and report held-out accuracy.

    Alternate points are held out.  Because the underlying transform is
    deterministic rather than noisy, held-out error falls monotonically with
    degree; it is reported so the choice of degree is justified rather than
    assumed.
    """
    order = np.argsort(raw)
    raw, corrected = raw[order], corrected[order]

    train = np.zeros(len(raw), dtype=bool)
    train[::2] = True

    holdout_fit = np.poly1d(np.polyfit(raw[train], corrected[train], degree))
    residual = corrected[~train] - holdout_fit(raw[~train])

    coefficients = np.polyfit(raw, corrected, degree)
    return coefficients, np.abs(residual).max(), float(np.sqrt((residual ** 2).mean()))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw-root", required=True, type=Path,
                        help="working tree containing ExperimentalResult/ and Postprocessed/")
    parser.add_argument("--degree", type=int, default=9)
    args = parser.parse_args()

    raw, corrected = collect_pairs(args.raw_root)
    print(f"{len(raw)} paired readings; "
          f"raw {raw.min():.4f}..{raw.max():.4f}, "
          f"corrected {corrected.min():.4f}..{corrected.max():.4f}")

    coefficients, worst, rms = fit(raw, corrected, args.degree)
    print(f"\ndegree {args.degree}: held-out max error {worst:.3e}, RMS {rms:.3e}")
    print("\nRECOVERED_OD_CORRECTION = [")
    for value in coefficients:
        print(f"    {value!r},")
    print("]")


if __name__ == "__main__":
    main()
