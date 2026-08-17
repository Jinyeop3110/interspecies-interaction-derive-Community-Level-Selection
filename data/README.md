# Data

Processed data are archived on Dryad: **https://doi.org/10.5061/dryad.2z34tmq0z**

This directory is empty in git. Unpack the Dryad archive here, or set
`COALESCENCE_DATA` to wherever you put it.

## Files

Every figure script reads from this set. `source` is `synthetic` (isolate-derived
communities), `natural` (environmental sample-derived communities) or
`simulation` (generalized Lotka-Volterra).

| File | Rows | Contents |
|---|---|---|
| `processed_Sequences_synthetic.xlsx` | one per sample | Relative abundance across ASVs, indexed by `SampleIDX`. The composition vectors the similarity decomposition operates on. |
| `processed_Sequences_natural.xlsx` | one per sample | As above, natural communities. |
| `processed_CoalescenceEvent_synthetic.xlsx` | one per event | Event metadata (`Medium`, `Replicate`, parent sample IDs `SampleIDX_Sub1/2`), OD and pH readings, and precomputed similarity/diversity summaries at five abundance thresholds. |
| `processed_CoalescenceEvent_natural.xlsx` | one per event | As above, natural communities. |
| `processed_CoalescenceEvent_simulation.xlsx` | one per event | Simulated coalescence events, indexed by species-pool size `S` and interaction strength `I`. |
| `processed_Communities_synthetic.xlsx` | one per community | Parental and coalesced communities with diversity metrics. |
| `processed_Communities_natural.xlsx` | one per community | As above, natural communities. |
| `Metadata.xlsx` | one per community | Experimental design factors (timepoint, origin, medium, coalescence type, replicate, community index) plus endpoint OD, pH and growth-curve readings. Loaded by `io.load_metadata()`. |
| `CoalescenceRecipe.xlsx` | one per coalescence | Which two parental community indices were mixed to produce each coalesced community. Sheet 1 synthetic, sheet 2 natural. Loaded by `io.load_coalescence_recipe()`. |

| `processed_CoalescenceEvent_simulation_uniform.xlsx` | one per event | Simulated events under the uniform interaction ensemble. Loaded by the Extended Data / Supplementary harness in `pipeline/04_figure_generation/extended_supplementary/`. |
| `M_OTUtableGreenGenes.csv` | one per sample | Raw DADA2 ASV counts before processing. Needed by analyses that work from read counts rather than relative abundance, such as the additive null model behind Extended Data Fig. 3. |

The last four entries are not optional. `Metadata.xlsx` and
`CoalescenceRecipe.xlsx` are required by `pipeline/04_figure_generation/fig5`,
`pipeline/04_figure_generation/supplementary` and the processing pipeline; the simulation-uniform
table and the raw count table are required by `pipeline/04_figure_generation/extended_supplementary/`. All of them must be
in the archive alongside the six processed tables.

### Key columns

- `SampleIDX` — primary key, joins the event tables to the sequence tables.
- `SampleIDX_Sub1`, `SampleIDX_Sub2` — the two parental communities of an event.
- `Medium` — `L` (Nutr−, no added glucose/urea), `M` (Base), `H` (Nutr+, high
  supplementation).
- `NormalizedAbundance*` — per-ASV relative abundance, summing to 1 per sample.
- `Threshold_level_1..5` — the abundance thresholds (0.1, 0.033, 0.01, 0.0033,
  0.001) at which the precomputed summaries were evaluated.
- `SimilarityTo1_*`, `SimilarityTo2_*`, `Assymetricity_*` — precomputed
  per-event summaries under several similarity measures (`BC` Bray-Curtis, `J`
  Jaccard, `JS` Jensen-Shannon, `DP` inner product). The figures do **not** read
  these; they recompute coordinates from the composition vectors via
  `coalescence.decomposition`. They are retained for cross-checking.

## How the processed data are generated

You do not need to run any of this to reproduce the figures — start from the
Dryad tables. This is the provenance chain, and the code is in
[`../pipeline`](../pipeline).

**1. Raw reads → ASV tables.** `pipeline/01_sequencing/` runs the DADA2 workflow
on demultiplexed 16S V4 fastq files: filter and trim, learn error profiles,
denoise, merge pairs, remove chimeras, assign taxonomy against SILVA 138.
Produces the OTU table, taxonomy table and unique-sequence table.

**2. ASV tables → per-sample and per-event tables.** `pipeline/02_processing/`
(MATLAB) is driven by `Main_synthetic.m` and `Main_natural.m`:

```
PostProcessingSequences_*.m   ASV table  -> processed_Sequences_*.xlsx
MetadataGenerator.m           plate maps -> Metadata.xlsx
RecipeGenerator.m             plate maps -> CoalescenceRecipe.xlsx
ExperimentalDataProcessing.m  OD and pH readings
AnalyzeCoalescence_*.m        -> processed_CoalescenceEvent_*.xlsx
AnalyzeCommunities_*.m        -> processed_Communities_*.xlsx
```

**3. Simulations.** `pipeline/03_simulation/` runs the generalized
Lotka-Volterra model: sample a pool of 54 species, assemble four non-overlapping
12-species parental communities, equilibrate, then coalesce all six pairs.
Repeating over independently sampled pools gives the simulated event set, which
`AnalyzeCoalescence_synthetic_simulations.m` summarizes into
`processed_CoalescenceEvent_simulation.xlsx`.

## Raw sequencing reads

The raw demultiplexed 16S V4 fastq files are in the same Dryad deposit as the
processed tables: **https://doi.org/10.5061/dryad.2z34tmq0z**

You only need them to re-run step 1 of the provenance chain above. The figure
scripts start from the processed tables and never touch the reads. To re-run
DADA2, unpack the reads and point `COALESCENCE_FASTQ` at the directory:

```bash
export COALESCENCE_FASTQ=/path/to/demultiplexed_fastq
Rscript pipeline/01_sequencing/Merging-SeqWorkflow_ForCoalescence.R
```

Read `pipeline/01_sequencing/README.md` first: it records three unresolved
issues in that script, one of which changes the resulting ASV tables.
