#!/usr/bin/env python3
"""Run the processing pipeline: ASV tables to processed data tables.

Python replacement for ``Main_synthetic.m`` and ``Main_natural.m``.

    python pipeline/02_processing/run_pipeline.py \
        --asv-dir SEQanalysis/excludeNatural \
        --design-dir Postprocessed \
        --out-dir Analyzed \
        --source synthetic

Reproduces the archived tables to floating-point tolerance; run
``tests/test_pipeline_parity.py`` to confirm against a known-good set.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from coalescence.processing.coalescence_events import analyze_coalescence_events  # noqa: E402
from coalescence.processing.communities import analyze_communities  # noqa: E402
from coalescence.processing.sequences import NATURAL, SYNTHETIC, process_sequences  # noqa: E402

# Per-dataset settings that differ between the two MATLAB entry points.
SOURCES = {
    "synthetic": {
        "config": SYNTHETIC,
        "community_origin": "S",
        "recipe_sheet": 0,
        "recipe_columns": (
            "CommunityIDX_Coal", "CommunityIDX_Sub1", "CommunityIDX_Sub2",
        ),
    },
    "natural": {
        "config": NATURAL,
        "community_origin": "N",
        "recipe_sheet": 1,
        "recipe_columns": (
            "CommunityIDX_Coal_natural",
            "CommunityIDX_Sub1_natural",
            "CommunityIDX_Sub2_natural",
        ),
    },
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--asv-dir", required=True, type=Path,
                        help="directory holding the DADA2 output CSVs")
    parser.add_argument("--design-dir", required=True, type=Path,
                        help="directory holding Metadata.xlsx and CoalescenceRecipe.xlsx")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="where to write the processed tables")
    parser.add_argument("--source", choices=sorted(SOURCES), default="synthetic")
    args = parser.parse_args()

    settings = SOURCES[args.source]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/3] sequences ({args.source})")
    compositions, annotations = process_sequences(
        pd.read_csv(args.asv_dir / "M_OTUtableGreenGenes.csv"),
        pd.read_csv(args.asv_dir / "M_TAXAtableGreenGenes.csv"),
        pd.read_csv(args.asv_dir / "M_UniqueSequences.csv"),
        settings["config"],
    )
    sequences_path = args.out_dir / f"processed_Sequences_{args.source}.xlsx"
    with pd.ExcelWriter(sequences_path) as writer:
        compositions.to_excel(writer, sheet_name="Sheet1", index=False)
        annotations.to_excel(writer, sheet_name="Sheet2", index=False)
    print(f"      -> {sequences_path}  {compositions.shape}")

    metadata = pd.read_excel(args.design_dir / "Metadata.xlsx")
    recipe = pd.read_excel(args.design_dir / "CoalescenceRecipe.xlsx",
                           sheet_name=settings["recipe_sheet"])
    sequences = compositions.set_index("SampleIDX")
    sequences.index = sequences.index.astype(str)

    print(f"[2/3] coalescence events ({args.source})")
    events = analyze_coalescence_events(
        metadata, sequences, recipe,
        community_origin=settings["community_origin"],
        recipe_columns=settings["recipe_columns"],
    )
    events_path = args.out_dir / f"processed_CoalescenceEvent_{args.source}.xlsx"
    events.to_excel(events_path, index=False)
    print(f"      -> {events_path}  {events.shape}")

    print(f"[3/3] communities ({args.source})")
    community_table = analyze_communities(
        metadata, sequences, community_origin=settings["community_origin"],
    )
    communities_path = args.out_dir / f"processed_Communities_{args.source}.xlsx"
    community_table.to_excel(communities_path, index=False)
    print(f"      -> {communities_path}  {community_table.shape}")


if __name__ == "__main__":
    main()
