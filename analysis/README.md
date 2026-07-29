# Extended Data and Supplementary analyses

These are the working analysis scripts, published as they were run. They are
not written to the standard of `src/coalescence` — they carry their own
plotting setup, some duplicate each other, and their variable naming is
inconsistent. They are here because they are the code behind the figures, and
availability matters more than polish.

For new work, prefer the library in [`../src/coalescence`](../src/coalescence).

## Running them

```bash
python analysis/analyze_additive_null.py
```

Each script writes its figures and tables to `output/`, which is not tracked.
Data is read through `coalescence.io`, so `COALESCENCE_DATA` works here too.

## The shared harness

`common_setup.py` is what every script imports. It provides the plotting
defaults, the loaded data tables, the community-lookup helpers, and the
similarity-decomposition functions.

It is a port of the working-tree module of the same name, changed so the
scripts run off this machine:

- data is loaded through `coalescence.io` instead of hardcoded `../../Analyzed/`
  paths;
- nothing is created or written at import time;
- `exception_list` is `coalescence.io.EXCLUDED_SAMPLES`, so the quality-control
  exclusions have one definition in the repository rather than an inline copy
  per script.

Helper signatures are unchanged, so scripts port with only their data paths
adjusted.

`COLORMAP.py` is the palette module the scripts import directly.

## Coverage

| Figure | Script | Status |
|---|---|---|
| Extended Data Fig. 3 | `analyze_additive_null.py` | Runs; reproduces the published Dominance counts (35 / 54 / 68 for Nutr−/Base/Nutr+) |

The remaining Extended Data and Supplementary figures are not yet ported. Their
generating code exists in the working tree; porting each one is a matter of
repointing its data paths at `coalescence.io`, as
`analyze_additive_null.py` shows.

## Porting recipe

Each working-tree script computes its own `ROOT` by walking up from `__file__`
and joins `Analyzed/`, `Postprocessed/` or `SEQanalysis/` onto it. Replace that
block with:

```python
from coalescence import io as _io

DATA = _io.data_dir()
coalescence_path = DATA / "processed_CoalescenceEvent_synthetic.xlsx"
```

and replace any inline `exception_list` with `set(_io.EXCLUDED_SAMPLES)`.
Everything else generally runs unchanged.
