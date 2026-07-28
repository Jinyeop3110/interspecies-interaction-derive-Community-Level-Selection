# Figure 3 — Interaction strength controls the transition between outcome types

## Panels

| Panel | Content | Status |
|---|---|---|
| a | Simulated outcome maps at μ = 0.3, 0.6, 0.8 | Not yet ported |
| b | Outcome fractions across μ = 0 to 1.2 | Not yet ported |

## Inputs

- `data/processed_CoalescenceEvent_simulation.xlsx` (`I` indexes interaction
  strength, `S` the parental-community richness)
- Simulation composition vectors — **not currently in the Dryad archive**

## Blocker

Same as Fig. 2: the classification needs composition vectors that the archived
simulation summary table does not carry. See
[`../fig2/README.md`](../fig2/README.md).

Each reported fraction is computed from 1,200 simulated coalescence events per
value of μ (200 independently sampled species pools × 6 pairwise combinations
within each pool).

## Source

`Main_Fig_3/Generate_Fig3_1..3_3.ipynb`, `Main_Fig_3/Impt_Generate_Fig3_4.ipynb`.
