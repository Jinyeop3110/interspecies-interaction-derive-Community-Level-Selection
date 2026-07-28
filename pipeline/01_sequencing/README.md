# Amplicon sequence processing

`Merging-SeqWorkflow_ForCoalescence.R` runs the DADA2 workflow that turns
demultiplexed 16S V4 fastq files into the ASV tables consumed by step 02.

## Running

Set these at the top of the script before running:

- `miseq_path` — directory of demultiplexed fastq files
- the SILVA training set path passed to `assignTaxonomy`

## Known issues to resolve before release

These are recorded here rather than silently patched, because changing them
changes the ASV tables behind the published results and so requires a re-run
and a comparison against the archived tables.

1. **Reverse error model is learned from forward reads.** `errR` is fitted with
   `learnErrors(filtFs[exists], ...)` where it should use `filtRs`. The result
   feeds `dada(derepRs, err=errR, ...)`. Denoising of the reverse reads
   therefore used a forward-read error profile. Re-run with the corrected call
   and confirm the resulting ASV table matches the archived one before
   publishing either version.

2. **Misleading GreenGenes naming.** Variables and output files are named
   `GreenGenesTaxM` / `M_GreenGenesTaxa.csv` and a comment refers to GreenGenes,
   but taxonomy is assigned against SILVA 138, which is what the manuscript
   Methods states. Renaming is cosmetic and safe, but touches the filenames that
   step 02 reads.

3. **Hardcoded absolute paths.** Input and reference paths point into a local
   machine and must be parameterized.
