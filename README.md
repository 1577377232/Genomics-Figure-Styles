# Genomics Figure Styles

A reusable, publication-oriented visual style library for chicken pangenomics, transposable elements (TEs), structural variants (SVs), GWAS and molecular-QTL figures.

> These independently designed, literature-inspired palettes are **not** official Nature Genetics, Nature Communications or Advanced Science templates.

## Start here

- [SNP + TE Manhattan style guide](figures/gwas_manhattan/STYLE_GUIDE.md)
- [Machine-readable theme](figures/gwas_manhattan/theme.json)
- [Three coordinated palettes](figures/gwas_manhattan/palettes.json)
- [Codex implementation task](figures/gwas_manhattan/CODEX_TASK.md)
- [Palette preview comparison](figures/gwas_manhattan/STYLE_COMPARISON_simulated.png)
- [Research inspiration / paper links](figures/gwas_manhattan/LITERATURE_REFERENCES.md)

## Design principles

One genome reference; reproducible chromosome coordinates; consistent color semantics; clear foreground TE markers and subdued SNP backgrounds; transparent model provenance; vector-preserving publication exports; honest interpretation of association versus fine-mapping/colocalization.

## Scope

The v2.0 style kit provides visual configuration, documentation, preview generator and **simulated** reference previews. It does **not** contain actual GWAS results or a completed 89-trait plotting application. Codex should implement that pipeline from the real local input data using the instructions in this repository.

## Agent instructions

Read [AGENTS.md](AGENTS.md) first. Further module-specific instructions are provided under `figures/gwas_manhattan/AGENTS.md`.

## Privacy

Do not commit individual-level genotypes, phenotypes, expression profiles or raw genome-wide association matrices. Use a private GitHub repository if you wish to keep your evolving research visualization rules private.
