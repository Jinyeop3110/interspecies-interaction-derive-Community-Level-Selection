# Figure 2 — Generalized Lotka-Volterra model of community coalescence

## Panels

| Panel | Content | Status |
|---|---|---|
| a | Simulation workflow schematic | Illustration |
| b | Simulated outcomes in similarity space at μ = 0.6 | Not yet ported |
| c | Interaction coefficients before and after assembly | Not yet ported |
| d | Pairwise selection correlation, within vs across parental communities | Not yet ported |

## Inputs

- `data/processed_CoalescenceEvent_simulation.xlsx` — summary table
- Simulation composition vectors — **not currently in the Dryad archive**

## Blocker

Panels b and d need the per-species composition vectors of the simulated
communities, not just the per-event summaries. The published simulation table
carries precomputed similarity columns but not the composition vectors that
`coalescence.decomposition` operates on, so these panels cannot be regenerated
from the archive as it stands.

Either the simulation outputs must be added to the archive, or
`pipeline/03_simulation/` must be able to regenerate them deterministically from
a seed. The latter is preferable — it makes the model results reproducible
rather than merely re-plottable.

## Source

`Main_Fig_2/Generate_Fig2_1..2_4.ipynb`, `Main_Fig_2/Impt_Generate_Fig2_5.ipynb`.
