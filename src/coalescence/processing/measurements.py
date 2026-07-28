"""Stage 0b: plate-reader OD, pH and growth-curve readings.

Python port of ``ExperimentalDataProcessing.m`` and the ``ODread_*`` /
``PHread_*`` / ``GCread_*`` helpers, which attach the measurement columns to
``Metadata.xlsx``.

Plate layout
------------
Each measurement workbook holds one sheet per growth cycle (seven in total).
Readings sit in the block ``B25:M32`` — 8 rows by 12 columns.  Rows 1-2 are the
two replicates of the 6-species communities, rows 4-5 the 12-species, and rows
7-8 the 24-species; columns index the community within each richness block.

The resulting array is ``[community, replicate, cycle]``, then flattened
richness-block by richness-block to give one row per community and one column
per cycle: ``fieldOD1..7``, ``fieldPH1..7``, ``fieldGC1..3``.

OD correction
-------------
Raw absorbance is noise-subtracted and then linearized through a calibration
curve.  In MATLAB that curve is a Curve Fitting Toolbox ``cfit`` object stored
in ``OD_correction_function_2202.mat``, which cannot be read outside MATLAB.
:func:`load_od_correction` recovers it numerically instead — see that function.
"""

import numpy as np
import pandas as pd

__all__ = [
    "READ_RANGE",
    "read_plate_block",
    "reshape_parental_plate",
    "reshape_coalesced_plate",
    "load_od_correction",
    "attach_measurements",
]

#: Excel block holding the readings, as (first_row, first_col, n_rows, n_cols),
#: zero-based, corresponding to the MATLAB range 'B25:M32'.
READ_RANGE = (24, 1, 8, 12)

#: Number of growth cycles per measurement workbook.
N_CYCLES = 7

#: Rows of the plate block occupied by each richness class, and how many
#: communities each holds. ``(row_offset, n_communities)``.
PARENTAL_BLOCKS = [(0, 9), (3, 9), (6, 12)]


def read_plate_block(path, sheet):
    """Read the ``B25:M32`` measurement block from one sheet.

    Wells recorded as ``-`` become NaN, matching MATLAB's ``xlsread``.
    """
    first_row, first_col, n_rows, n_cols = READ_RANGE
    frame = pd.read_excel(path, sheet_name=sheet, header=None)
    block = frame.iloc[first_row:first_row + n_rows, first_col:first_col + n_cols]
    return block.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)


def reshape_parental_plate(block):
    """Split one plate block into the three richness classes.

    Returns ``{"s6": (9, 2), "s12": (9, 2), "s24": (12, 2)}`` indexed by
    ``[community, replicate]``.
    """
    out = {}
    for name, (row_offset, n_communities) in zip(("s6", "s12", "s24"), PARENTAL_BLOCKS):
        values = np.zeros((n_communities, 2))
        for community in range(n_communities):
            for replicate in range(2):
                values[community, replicate] = block[row_offset + replicate, community]
        out[name] = values
    return out


def reshape_coalesced_plate(block):
    """Flatten a coalesced-community plate block, replicate-major."""
    return block


#: Degree-5 polynomial recovered from 1,260 paired raw and corrected readings.
#: Highest order first, for :class:`numpy.poly1d`.  See
#: :func:`load_od_correction` for provenance and accuracy.
RECOVERED_OD_CORRECTION = [
    0.04246354, -0.06013511, 0.12759108, 0.18337700, 1.00435607, -0.00019626,
]


def load_od_correction(coefficients=None, raw_samples=None, corrected_samples=None,
                       degree=5):
    """Return the OD linearization curve as a callable.

    Plate-reader absorbance saturates at high cell density, so raw readings are
    passed through a calibration curve before use.  In MATLAB that curve is a
    Curve Fitting Toolbox ``cfit`` object stored in
    ``OD_correction_function_2202.mat``.  The file holds a serialized MATLAB
    class instance (MCOS), which cannot be deserialized outside MATLAB, so the
    curve itself is not portable.

    Three ways to obtain it, best first:

    1. **Export the real curve from MATLAB.**  Run

           load('OD_correction_function_2202.mat')
           writematrix(coeffvalues(OD_correction_function), 'od_coeffs.csv')
           type(OD_correction_function)   % records the model form

       and pass the result as ``coefficients``.  This is the only way to
       reproduce the archived values exactly, and is what a release should do.

    2. **Use the recovered approximation** (the default).  Fitting a degree-5
       polynomial through 1,260 paired raw and corrected readings taken from the
       raw plate files and the archived ``Metadata.xlsx`` reproduces the curve to
       a maximum error of about 1e-3 OD and an RMS error of about 5e-5 on
       held-out points.  That is well below plate-reader precision but is an
       approximation, not the original function.

    3. **Refit from your own paired data** by passing ``raw_samples`` and
       ``corrected_samples``.
    """
    if coefficients is not None:
        return np.poly1d(np.asarray(coefficients, dtype=float))

    if raw_samples is not None and corrected_samples is not None:
        return np.poly1d(
            np.polyfit(np.asarray(raw_samples, dtype=float),
                       np.asarray(corrected_samples, dtype=float), degree)
        )

    return np.poly1d(RECOVERED_OD_CORRECTION)


def correct_od_block(block, correction=None):
    """Noise-subtract and linearize one raw OD plate block.

    Noise is the mean of all wells within 5% of the plate minimum, matching the
    MATLAB ``Noise=mean(y(y<1.05*min(y(:))))``.
    """
    correction = correction or load_od_correction()
    noise = np.nanmean(block[block < 1.05 * np.nanmin(block)])
    return correction(block - noise)


def attach_measurements(metadata, od_by_group, ph_by_group, gc_by_group):
    """Append the measurement columns to the structural metadata table.

    Parameters
    ----------
    metadata : DataFrame
        Output of :func:`coalescence.processing.design.build_metadata`.
    od_by_group, ph_by_group, gc_by_group : dict
        Group key -> array of shape ``(n_communities, n_cycles)``, matching the
        row order of that group in ``metadata``.

    Returns a copy of ``metadata`` with ``fieldOD*``, ``fieldPH*`` and
    ``fieldGC*`` columns appended.
    """
    result = metadata.copy()
    group_keys = (
        result.Timepoint + "_" + result.CommunityOrigin + "_" + result.Medium
        + "_" + result.CoalescenceType + "_" + result.Replicate.astype(str)
    )

    for prefix, source in (("fieldOD", od_by_group),
                           ("fieldPH", ph_by_group),
                           ("fieldGC", gc_by_group)):
        if not source:
            continue
        width = max(v.shape[1] for v in source.values())
        for cycle in range(width):
            column = np.full(len(result), np.nan)
            for key, values in source.items():
                rows = np.flatnonzero((group_keys == key).to_numpy())
                if len(rows) != len(values):
                    raise ValueError(
                        f"group {key}: metadata has {len(rows)} rows but "
                        f"{prefix} has {len(values)}"
                    )
                column[rows] = values[:, cycle]
            result[f"{prefix}{cycle + 1}"] = column

    return result
