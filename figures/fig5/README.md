# Figure 5 — Predictability of Dominance direction

## Panels

| Panel | Content | Status |
|---|---|---|
| a | Schematic: can dominant-species competition predict the winner? | Illustration |
| b | Relative abundance of dominant species across media | Not yet ported |
| c | Dominant-species competitive success vs event PDI | Not yet ported |

## Inputs

- `data/processed_CoalescenceEvent_synthetic.xlsx`
- `data/processed_Sequences_synthetic.xlsx`
- Pairwise invasion assay outcomes — **not currently in the Dryad archive**

## Blocker

Panel c regresses event PDI against the outcome of the pairwise competition
between the two parental communities' dominant species. The assay results are
not among the processed tables, so the panel cannot be regenerated from the
archive as it stands. PDI itself is available from
`coalescence.outcomes.outcome_table`.

Panel c applies two filters that are not in the shared library and must be
ported with it:

- a mixing filter, keeping only events with `u² + v² > 0.5` (that is, excluding
  Restructuring);
- a bistability filter on the monoculture reference measurements.

## Source

`Main_Fig_5/Generate_Fig5_1..5_5.ipynb` (regression panels),
`Main_Fig_6/Generate_Fig6_2_MostAbundant.ipynb` (dominant-species abundance).
Note both working-tree directories 5 and 6 feed manuscript Figure 5.
