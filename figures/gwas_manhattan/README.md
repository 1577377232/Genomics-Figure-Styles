# SNP + TE Joint Manhattan Style — v2.0

## Goals
Use one reproducible Nature Portfolio-inspired (not official journal) system for 89 chicken traits. This repository **does not contain the 89 GWAS input data or their completed plots**.

## Files
- [Style rules](STYLE_GUIDE.md)
- [Figure tokens](theme.json)
- [Three named palettes](palettes.json)
- [Color semantics](COLOR_SEMANTICS.md)
- [Codex implementation brief](CODEX_TASK.md)
- [Literature rationale](LITERATURE_REFERENCES.md)
- [Trait metadata template](trait_manifest.template.csv)

Default theme: **Pangenome Spectrum**. Alternative palettes: **Botanical Genome** and **Chromatic NC**.

## Expected Codex outputs (in research worktree, not tracked here)
89 individual joint SNP+TE Manhattan plots (PDF and PNG 600 dpi), grouped contact sheets, provenance and QC metadata, optional fixed-y comparison versions. Real result paths, genome FAI, sample models and 89 trait IDs must be supplied/discovered before the plot batch can run.

SNPs: per-chromosome six-color cycle, small translucent circles, selectively rasterized in PDF. TE insertions: coral diamonds, deletions: violet triangles, other/unknown: amber with an explicit class. Highlight only TE that surpass a *documented* threshold.

The palettes are independently composed, inspired by publications; not copied or endorsed.
