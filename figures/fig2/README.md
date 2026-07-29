# Figure 2 — Generalized Lotka-Volterra model of community coalescence

The simulation is in
[`../../pipeline/03_simulation`](../../pipeline/03_simulation). It is seeded per
replicate, so the figures regenerate from code and no simulation output needs to
be archived.

## Panels

| Panel | Content | Status |
|---|---|---|
| a | Simulation workflow schematic | Illustration |
| b | Simulated outcomes in similarity space at μ = 0.6 | Runs; see configuration below |
| c | Interaction coefficients before and after assembly | Not yet ported |
| d | Pairwise selection correlation, within vs across parental communities | Not yet ported |

## Configuration

The sweep is parameterized, and the outcome fractions depend on it. Measured at
μ = 0.6:

| Communities | Pool | Replicates | n | Dominance | Mixture | Restructuring |
|---|---|---|---|---|---|---|
| 2 | 24 | 1,200 | 1,200 | 59.7% | 15.0% | 25.3% |
| 4 | 48 | 200 | 1,200 | 61.3% | 16.2% | 22.5% |
| 4 | 54 | 200 | 1,200 | 57.4% | 17.0% | 25.6% |
| — | — | — | — | *61%* | *13%* | *26%* |

The last row is the published panel. The two-community configuration is closest
across all three classes, and is the configuration recorded in every cached
simulation session in the working tree.

Reproduce it with:

```bash
python pipeline/03_simulation/coalescence_simulation.py \
    --mu-min 0.6 --mu-max 0.7 --mu-step 0.1 \
    --replicates 1200 --communities 2 --species 12
```

Residual differences of one to two percentage points are not yet accounted for.
Candidates are the extinction threshold, the ODE solver tolerance, and the exact
value of μ — several cached sessions used μ = 0.55, some with a spread term.

## Source

`Main_Fig_2/Generate_Fig2_1..2_4.ipynb`, `Main_Fig_2/Impt_Generate_Fig2_5.ipynb`.
