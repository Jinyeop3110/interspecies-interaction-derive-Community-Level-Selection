# Stage 4: figures

The last stage of the pipeline. Everything here reads the processed tables from
[`../../data`](../../data) — or runs the simulation in
[`../03_simulation`](../03_simulation) — and emits figure panels.

```
fig1/ .. fig6/            one directory per main-text figure
supplementary/            Supplementary figures
extended_supplementary/   Extended Data and Supplementary analysis scripts,
                          published as they were run
```

Each figure directory has its own README recording which panels it generates,
which inputs they read, and how the output compares with the published figure.
Panels are written to `panels/` inside each directory and are not tracked.

Scripts emit **individual panels as SVG**. The composed figures in the
manuscript were assembled from these panels in a vector editor, so a panel
matches its published counterpart in content and scaling but not in layout,
lettering or annotation.

## Coverage

| | Panels with code |
|---|---|
| Main text | 9 of 15 data panels |
| Extended Data | 1 of 8 figures |
| Supplementary | 1 of 46 figures |

`fig1/README.md` through `fig6/README.md` state per panel what is and is not
covered, and which published numbers reproduce.
