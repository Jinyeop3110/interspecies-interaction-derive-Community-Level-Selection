# Pipeline: raw reads to processed tables

You do not need to run any of this to reproduce the figures. Start from the
processed tables on Dryad; see [`../data/README.md`](../data/README.md).

This directory documents how those tables were produced.

## 01_sequencing — amplicon processing (R)

`Merging-SeqWorkflow_ForCoalescence.R` runs DADA2 over demultiplexed 16S V4
fastq files: filter and trim (truncation at 145 bases both directions, maximum
2 expected errors), learn error profiles, denoise, merge pairs, remove chimeras
by consensus, and assign taxonomy against SILVA 138.

Outputs the OTU table, taxonomy table and unique-sequence table consumed by
step 02.

Set the input and reference paths at the top of the script before running.

## 02_processing — experiment tables

**Python is the reference implementation.** `run_pipeline.py` replaces
`Main_synthetic.m` and `Main_natural.m`:

```bash
python pipeline/02_processing/run_pipeline.py \
    --asv-dir SEQanalysis/excludeNatural \
    --design-dir Postprocessed \
    --out-dir Analyzed \
    --source synthetic
```

It reproduces every archived table to floating-point tolerance. `tests/
test_pipeline_parity.py` checks this stage by stage; set `COALESCENCE_RAW` to a
tree containing `SEQanalysis/`, `Postprocessed/` and `Analyzed/` to run it.

| Table | Rows x cols | Largest difference from archived |
|---|---|---|
| `processed_Sequences_synthetic` | 590 x 44 | 0 |
| `processed_Sequences_natural` | 144 x 131 | 0 |
| `processed_CoalescenceEvent_synthetic` | 282 x 143 | 7e-16 |
| `processed_CoalescenceEvent_natural` | 90 x 143 | 7e-16 |
| `processed_Communities_synthetic` | 462 x 46 | 3e-15 |
| `processed_Communities_natural` | 126 x 46 | 2e-15 |

The MATLAB originals are kept alongside for reference. They were driven by
`Main_synthetic.m` and `Main_natural.m`, which call the other scripts in order:

| Script | Produces |
|---|---|
| `PostProcessingSequences_*.m` | `processed_Sequences_*.xlsx` |
| `MetadataGenerator.m` | `Metadata.xlsx` |
| `RecipeGenerator.m`, `CoalescenceRecipe.m` | `CoalescenceRecipe.xlsx` |
| `ExperimentalDataProcessing.m` | OD and pH readings |
| `AnalyzeCoalescence_*.m` | `processed_CoalescenceEvent_*.xlsx` |
| `AnalyzeCommunities_*.m` | `processed_Communities_*.xlsx` |
| `BringSimulationResults.m` | collects simulation runs into `processed_simulations.xlsx` |

### Design and measurement tables

`Metadata.xlsx` and `CoalescenceRecipe.xlsx` encode the plate layouts, the
parental pairings and the plate-reader readings. Both are also ported:

| Module | Replaces | Parity |
|---|---|---|
| `processing/design.py` | `MetadataGenerator.m`, `RecipeGenerator.m`, `CoalescenceRecipe.m` | exact |
| `processing/measurements.py` | `ExperimentalDataProcessing.m` and the `ODread_*`/`PHread_*`/`GCread_*` helpers | pH exact; OD to 6e-7 |

### The OD calibration curve

The MATLAB applies a Curve Fitting Toolbox `cfit` object stored in
`OD_correction_function_2202.mat`. That is a serialized MATLAB class instance
(MCOS) and cannot be deserialized outside MATLAB — but it does not need to be.

The curve is a smooth, monotone, deterministic transform, and the archived data
contains both sides of it: raw plate readings on one side, corrected OD in
`Metadata.xlsx` on the other. Fitting a polynomial through those pairs recovers
it directly:

```bash
python tools/recover_od_calibration.py --raw-root /path/to/working/tree --degree 9
```

Held-out accuracy is 9e-6 OD maximum (RMS 4e-7), and end-to-end the pipeline
reproduces the archived `fieldOD` columns to **6e-7 across 1,260 readings** —
four orders of magnitude below plate-reader precision. The coefficients are
baked into `measurements.RECOVERED_OD_CORRECTION`; the tool regenerates them.

No MATLAB is required at any stage.

## 03_simulation — generalized Lotka-Volterra model (Python)

Species grow logistically and compete pairwise:

    dn_i/dt = r_i n_i (1 - sum_j alpha_ij n_j)

with `r_i = 1`, self-interaction `alpha_ii = 1`, and off-diagonal competition
coefficients drawn from `U(0, 2μ)`, so μ sets both the mean and the spread of
competitive effects.

Each replicate draws a pool of 54 species, assembles four non-overlapping
12-species parental communities, equilibrates them (species below 1e-3 relative
abundance set to zero), and coalesces all six pairs at equal proportions —
mirroring the experimental protocol of parental assembly followed by pairwise
coalescence. 200 replicate pools × 6 pairs = 1,200 coalescence events per value
of μ.

`coalescence_simulation.py` implements the model and the sampling design above.
It is seeded per replicate, so sweeps are fully determined by their parameters
and no simulation output needs to be archived. See
[`03_simulation/README.md`](03_simulation/README.md) for verification against
the cached notebook outputs, and
[`../figures/fig2/README.md`](../figures/fig2/README.md) for how the outcome
fractions depend on the sweep configuration.
