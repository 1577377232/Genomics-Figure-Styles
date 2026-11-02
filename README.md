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

The v2.0 style kit provides visual configuration and documentation. The [executable R workflow](figures/gwas_manhattan/R/README.md) adds tested joint SNP/TE Manhattan plotting for 89 quantitative traits and seven multivariate modules, with the user-approved two-colour preset and joint-sex Z provenance checks. Actual GWAS inputs and completed biological figures remain in the research workspace.

## Agent instructions

Read [AGENTS.md](AGENTS.md) first. Further module-specific instructions are provided under `figures/gwas_manhattan/AGENTS.md`.

## Privacy

Do not commit individual-level genotypes, phenotypes, expression profiles or raw genome-wide association matrices. Use a private GitHub repository if you wish to keep your evolving research visualization rules private.
