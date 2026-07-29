# Supplementary figures

Only one supplementary figure currently has a script here. The remaining 45 are
not covered; their generating code exists but has not been prepared for release.

| Figure | Content | Script |
|---|---|---|
| 14 | Parental OD difference vs PDI, and parental endpoint OD vs endpoint pH | `make_supp_fig14_od_ph.py` |

```bash
python figures/supplementary/make_supp_fig14_od_ph.py
```

Panels are written to `panels/`, like every other figure directory.

The script previously wrote straight into the manuscript tree at a path that
only exists on the author's machine, which meant a fresh clone silently created
a stray `latex/` directory beside the repository. To regenerate the manuscript
asset, point it there explicitly:

```bash
COALESCENCE_FIGURE_OUT=../../../latex/supplementary_figs \
    python figures/supplementary/make_supp_fig14_od_ph.py
```

The output basename is `supp_fig14_od_vs_ph`; the manuscript includes
`Fig_R1_1B_OD_vs_PDI.pdf`, so rename on copy or adjust `OUT_BASE` if you want
the two to match.

## Inputs

- `data/Metadata.xlsx` — endpoint OD and pH per community
- `data/processed_Sequences_synthetic.xlsx`
- `data/processed_CoalescenceEvent_synthetic.xlsx`

The script reports Spearman correlations for each medium, both for the signed
parental OD difference against PDI and for parental endpoint OD against
endpoint pH.
