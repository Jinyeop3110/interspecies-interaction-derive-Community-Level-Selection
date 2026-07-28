"""Stage 1: ASV tables to per-sample composition tables.

Python port of ``PostProcessingSequences_synthetic.m`` and
``PostProcessingSequences_natural.m``.

Takes the DADA2 output (OTU table, taxonomy table, unique sequences) and
produces ``processed_Sequences_*.xlsx``: one row per sample, relative abundance
across the retained ASVs.

Two things happen here that are specific to this study:

*Cluster merging.* The 54-isolate library contains isolates whose 16S V4 region
DADA2 splits into several ASVs.  Those splits are merged back using
correspondences established from sequencing the isolates in monoculture.  Most
are simple sums; two (CL3, CL4) are ratio-based splits where a single ASV
carries reads from two different isolates in a known proportion.

*Filtering.* ASVs beyond the known library, or with too few total reads, or
present in too few samples, are dropped, and trace abundances are zeroed.

The MATLAB indices were 1-based; the merge tables below keep the original
1-based numbering and convert on use, so they can be checked against the
MATLAB source line by line.
"""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

__all__ = ["SYNTHETIC", "NATURAL", "SequenceProcessingConfig", "process_sequences"]


@dataclass
class SequenceProcessingConfig:
    """Per-dataset constants from the MATLAB scripts."""

    #: Keep only this many ASV columns before any merging.
    initial_cutoff_asv_num: int = 200
    #: ASV indices above this are outside the known isolate library.
    max_asv_num: int = 56
    #: Minimum total reads for an ASV to be retained.
    min_asv_reads: int = 500
    #: Relative abundances below this are set to zero.
    trace_abundance: float = 0.001
    #: Shannon-diversity floor across samples, in nats.
    min_shannon: float = np.log(2)
    #: ``target -> [sources]``: sum the source ASVs into the target, 1-based.
    merge_sums: dict = field(default_factory=dict)
    #: Ratio splits, 1-based ``(a, b, ratio)``; see ``_apply_ratio_split``.
    merge_ratio_splits: list = field(default_factory=list)
    #: ASVs with no match in the isolate library, 1-based.
    drop_not_found: list = field(default_factory=list)
    #: ASVs made redundant by the merges above, 1-based.
    drop_clustered: list = field(default_factory=list)


#: ``PostProcessingSequences_synthetic.m``
SYNTHETIC = SequenceProcessingConfig(
    merge_sums={
        18: [20, 28, 30],  # CL1, 3:2:1:1
        12: [23],          # CL2, 3:1
        38: [51],          # CL5, 6:1
        11: [36],          # CL6, 8:1
    },
    merge_ratio_splits=[
        (16, 33, 6),  # CL3 / CL3': ASV33 carries 1/6 of the reads assigned to ASV16
        (34, 50, 9),  # CL4 / CL4': ratio 8:1, see note in _apply_ratio_split
    ],
    drop_not_found=[43, 44, 45, 48, 49, 53, 55, 57],
    drop_clustered=[20, 28, 30, 23, 51, 36],
)

#: ``PostProcessingSequences_natural.m``
#:
#: Natural sample-derived communities are not built from the isolate library,
#: so no cluster merging applies and every ASV surviving the read-count and
#: abundance filters is kept.  The MATLAB declares ``MAX_ASV_num=50`` but never
#: uses it; the effective ceiling is the 200-column initial cutoff, which is
#: what is reproduced here.
NATURAL = SequenceProcessingConfig(
    max_asv_num=200,
    merge_sums={},
    merge_ratio_splits=[],
    drop_not_found=[],
    drop_clustered=[],
)


def _apply_ratio_split(counts, first, second, factor):
    """Split reads that DADA2 pooled into one ASV across two isolates.

    ``first`` and ``second`` are 1-based ASV indices.  The second ASV is a
    marker present at a known fixed ratio, so its reads are scaled up by
    ``factor`` to recover the true abundance of one isolate, and the remainder
    of the first ASV is attributed to the other.

    Note the asymmetry carried over from the MATLAB: the scale-up uses
    ``factor`` while the subtraction uses ``factor - 1`` for CL4 (9 and 8
    respectively).  This is reproduced exactly rather than normalized, since it
    is what generated the archived tables.
    """
    i, j = first - 1, second - 1
    subtract = factor if factor == 6 else factor - 1
    scaled = factor * counts[:, j]
    remainder = np.maximum(counts[:, i] - subtract * counts[:, j], 0)
    counts[:, i] = scaled
    counts[:, j] = remainder
    return counts


def _filter_asvs(counts, asv_index, config):
    """Port of the MATLAB ``Filtering`` helper.

    Returns the filtered count matrix and the surviving ASV indices.

    Carries over a quirk of the original: the Shannon-diversity criterion drops
    entries from the index vector without dropping the corresponding columns
    from the count matrix.  On the published data the criterion selects nothing,
    so the two stay aligned; the behaviour is preserved so that any future input
    reproduces the archived output rather than a silently corrected version.
    A warning fires if it ever does select a column.
    """
    input_reads = counts.sum()

    keep = ~((counts.sum(axis=0) < config.min_asv_reads) | (asv_index > config.max_asv_num))
    counts = counts[:, keep]
    asv_index = asv_index[keep]

    column_fraction = counts / counts.sum(axis=0)
    shannon = -np.sum(np.log(column_fraction + 1e-20) * column_fraction, axis=0)
    low_diversity = shannon < config.min_shannon
    if low_diversity.any():
        import warnings

        warnings.warn(
            f"{low_diversity.sum()} ASV(s) fall below the Shannon floor. The original "
            "MATLAB drops these from the index but not from the abundance matrix, "
            "which desynchronizes the two. Inspect before trusting the output.",
            RuntimeWarning,
        )
    asv_index = asv_index[~low_diversity]

    counts = np.where(counts / (counts.sum(axis=1, keepdims=True) + 1e-6)
                      < config.trace_abundance, 0, counts)

    print(f"  retained {counts.sum() / input_reads:.6f} of reads")
    return counts, asv_index


def process_sequences(otu_table, taxonomy_table, unique_sequences, config=SYNTHETIC):
    """Run stage 1.

    Parameters
    ----------
    otu_table : DataFrame
        DADA2 OTU table; first column is the sample name, remaining columns are
        per-ASV counts ordered by abundance.
    taxonomy_table, unique_sequences : DataFrame
        DADA2 taxonomy assignments and representative sequences, row-aligned to
        the ASV columns of ``otu_table``.
    config : SequenceProcessingConfig

    Returns
    -------
    compositions : DataFrame
        ``SampleIDX`` plus ``NormalizedAbundance1..n``; sheet 1 of the output.
    annotations : DataFrame
        Taxonomy and representative sequence per retained ASV; sheet 2.
    """
    table = otu_table.iloc[:, : config.initial_cutoff_asv_num + 1]
    sample_names = table.iloc[:, 0].astype(str)
    counts = table.iloc[:, 1:].to_numpy(dtype=float)
    asv_index = np.arange(1, counts.shape[1] + 1)

    for target, sources in config.merge_sums.items():
        for source in sources:
            counts[:, target - 1] += counts[:, source - 1]

    for first, second, factor in config.merge_ratio_splits:
        counts = _apply_ratio_split(counts, first, second, factor)

    to_drop = sorted(set(config.drop_not_found) | set(config.drop_clustered))
    if to_drop:
        mask = np.ones(counts.shape[1], dtype=bool)
        mask[[i - 1 for i in to_drop]] = False
        counts = counts[:, mask]
        asv_index = asv_index[mask]

    counts, asv_index = _filter_asvs(counts, asv_index, config)

    abundance = counts / (counts.sum(axis=1, keepdims=True) + 1e-6)

    compositions = pd.DataFrame(
        abundance,
        columns=[f"NormalizedAbundance{i + 1}" for i in range(abundance.shape[1])],
    )
    compositions.insert(
        0, "SampleIDX", sample_names.str.replace("_F_filt.fastq.gz", "", regex=False).values
    )

    rows = [i - 1 for i in asv_index]
    annotations = taxonomy_table.iloc[rows].reset_index(drop=True).copy()
    annotations.columns = ["ASV"] + list(annotations.columns[1:])
    annotations["UniqueSeuquences"] = unique_sequences.iloc[rows, 1].values

    return compositions, annotations
