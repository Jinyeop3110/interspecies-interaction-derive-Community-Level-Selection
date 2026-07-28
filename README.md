# Interspecies interactions drive community-level selection in microbial coalescence

Analysis and figure-generation code for the manuscript.

Coalescence experiments mix two separately assembled bacterial communities and
ask which one shapes the result. This repository contains the shared analysis
library that classifies coalescence outcomes, the scripts that generate the
main-text figure panels from processed data, and the upstream pipeline that
turns raw 16S reads into that processed data.

## Data

The processed data are archived on Dryad:

**https://doi.org/10.5061/dryad.2z34tmq0z**

Nothing in this repository ships data. Download the archive, unpack it into
`data/`, and every figure script will find it:

```bash
# from the repository root
unzip ~/Downloads/dryad_download.zip -d data/
```

Or point the code somewhere else:

```bash
export COALESCENCE_DATA=/path/to/processed/tables
```

See [`data/README.md`](data/README.md) for the file list, the schema of each
table, and how the tables are generated from raw reads.

## Quick start

```bash
conda env create -f environment.yml
conda activate coalescence

python figures/fig4/make_panels.py     # writes figures/fig4/panels/*.svg
pytest tests/                          # checks published values still reproduce
```

Each figure script prints the quantities it plots, so you can compare against
the manuscript without opening the SVGs.

## Layout

```
data/                 processed tables (downloaded from Dryad; not in git)
pipeline/             raw reads -> processed tables
  01_sequencing/      DADA2 amplicon workflow (R)
  02_processing/      community and coalescence-event tables (Python)
  03_simulation/      generalized Lotka-Volterra simulations (Python)
src/coalescence/      shared analysis library, imported by every figure script
figures/fig1..fig6/   one directory per main-text figure
tests/                reproduction checks against published values
tools/                release manifest and drift check
```

## How the analysis fits together

The central quantity is where a coalescence event lands in the **two-parent
similarity map**. Given the coalesced community and its two parents as
relative-abundance vectors, `coalescence.decomposition` expresses the coalesced
community in the basis spanned by the two parents. That yields two coordinates
plus a residual, normalized so they lie on the unit sphere:

- **retention magnitude** `r` — how much of the coalesced composition the two
  parents jointly explain. Low `r` means the community reorganized into
  something neither parent predicts.
- **angular asymmetry** — where the event sits between the diagonal (both
  parents contribute equally) and the axes (one parent explains everything).
  The parental dominance index reported in the paper is a monotone
  reparameterization of the same angle.

Outcomes are then classified as **Restructuring** (`r² ≤ 0.5`), **Mixture**
(high retention, near the diagonal) or **Dominance** (high retention, near an
axis).

The parental basis is oblique: whenever the two parental communities share
species, their composition vectors are not orthogonal, so the coordinates come
from solving the 2×2 Gram system rather than from two independent cosine
similarities. For disjoint parents the two agree.

## Figures

One directory per main-text figure, each with its own README recording which
panels are generated here, which inputs they read, and where any not-yet-ported
panels come from.

| Figure | Content | Directory |
|---|---|---|
| 1 | Coalescence of synthetic communities frequently yields Dominance | [`figures/fig1`](figures/fig1) |
| 2 | Generalized Lotka-Volterra model of coalescence | [`figures/fig2`](figures/fig2) |
| 3 | Interaction strength controls the outcome transition | [`figures/fig3`](figures/fig3) |
| 4 | Nutrient concentration modulates invasion resistance and outcomes | [`figures/fig4`](figures/fig4) |
| 5 | Predictability of Dominance direction | [`figures/fig5`](figures/fig5) |
| 6 | Dominance in natural sample-derived communities | [`figures/fig6`](figures/fig6) |

Extended Data and Supplementary figures are not covered here; they are
generated from the same library and the same processed tables.

Scripts emit **individual panels as SVG**. The composed figures in the
manuscript were assembled from these panels in a vector editor, so a panel will
match its published counterpart in content and scaling but not in layout,
lettering or annotation.

## Reproduction status

`tests/` checks the library against numbers printed in the manuscript.

| Check | Published | Reproduced |
|---|---|---|
| Fig. 6c Dominance, natural (Nutr−/Base/Nutr+) | 37% / 70% / 77% | 36.7% / 70.0% / 76.7% |
| Fig. 4d Dominance, synthetic (Nutr−/Base/Nutr+) | 39% / 65% / 76% | 38.0% / 63.6% / 73.4% |
| Fig. 4d Mixture, synthetic | 53% / 4% / 6% | 53.3% / 3.4% / 5.3% |
| Fig. 1e Base (Dominance / Restructuring / Mixture) | 65% / 31% / 4% | 63.6% / 33.0% / 3.4% |
| Fig. 5b dominant-species abundance | 44 / 51 / 67% | 42.4 / 59.6 / 66.7% |

The natural-community numbers reproduce exactly. The synthetic outcome
fractions sit within about 3 percentage points because a per-figure
quality-control filter drops a small number of events (published *n* = 90/83/90
against 92/88/94 from the unfiltered tables); that filter is documented in
[`figures/fig4/README.md`](figures/fig4/README.md) and is not yet part of the
shared library.

Figure 5b reproduces Nutr− and Nutr+ to within 1.6 percentage points but gives
59.6% for Base against a published 51%. That gap is confined to one medium and
is not explained by richness block; see
[`figures/fig5/README.md`](figures/fig5/README.md). Treat Base in panel 5b as
unverified.

## Citing

See [`CITATION.cff`](CITATION.cff).

## License

MIT — see [`LICENSE`](LICENSE).
