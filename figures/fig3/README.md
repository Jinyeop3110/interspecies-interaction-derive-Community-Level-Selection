# Figure 3 — Interaction strength controls the transition between outcome types

The simulation is in
[`../../pipeline/03_simulation`](../../pipeline/03_simulation), seeded per
replicate.

## Panels

| Panel | Content | Status |
|---|---|---|
| a | Simulated outcome maps at μ = 0.3, 0.6, 0.8 | Runs from the sweep |
| b | Outcome fractions across μ = 0 to 1.2 | Runs from the sweep |

## Running the sweep

```bash
python pipeline/03_simulation/coalescence_simulation.py \
    --mu-min 0 --mu-max 1.2 --mu-step 0.1 \
    --replicates 1200 --communities 2 --species 12 \
    --out fig3_sweep.csv
```

`outcome_fractions()` in that module gives the panel b curves directly.

The sweep configuration determines the reported event count; see
[`../fig2/README.md`](../fig2/README.md) for how the outcome fractions vary with
it and which configuration the cached sessions used.

## Verification

The simulation reproduces the cached notebook output exactly — five decimal
places on the similarity-map coordinates across the range of μ. See
[`../../pipeline/03_simulation/README.md`](../../pipeline/03_simulation/README.md).

## Source

`Main_Fig_3/Generate_Fig3_1..3_3.ipynb`, `Main_Fig_3/Impt_Generate_Fig3_4.ipynb`.
