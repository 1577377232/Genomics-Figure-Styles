# Genotype–phenotype boxplot | 基因型–表型箱线图

Publication-oriented R template for raw phenotype distributions stratified by diploid genotype and dominant insertion-carrier status.

## Run

Install ggplot2 and optionally svglite. Cairo graphics and Liberation Sans are recommended.

```bash
Rscript 02_boxplots.R /path/to/input_directory
```

Without an argument, the script uses its own directory (or the working directory when sourced interactively). Outputs are written alongside the three CSV inputs and overwrite matching filenames.

## Input contracts

- `plot_data.csv`: numeric `GT` (0, 1, 2), numeric `abdominal_fat_percent` (already in percent).
- `pairwise_LMM.csv`: `comparison` (`GT0_vs_GT1`, `GT1_vs_GT2`, `GT0_vs_GT2`) and `P_Holm`; one row per contrast.
- `dominant_GWAS.csv`: one row with `P` from the original dominant-model association analysis.

Change `TRAIT_COL`, `Y_LABEL`, `EXPECTED_N`, and output stems for other verified traits. The current default expects 822 observations and all three genotype groups. Do not invent missing group comparisons or replace source LMM P values with unadjusted tests.

## Style

85 × 85 mm, white background, no grid, blue `#386C8C`, lavender `#9689B9`, coral `#E75A46` (coordinated with the repository variant-class palette), ink `#253246`, opaque 78% color strength. Separate whisker and box geometry prevents internal vertical lines. Bottom padding separates whiskers from the horizontal axis. REF = 0/0; INS = 0/1 + 1/1. These group labels apply to insertion variants only; do not confuse the genotype groups with SNP/SV/TEV variant classes.

P values use literal e notation with three significant digits (e.g. `4.92e-8`), with italic P and subscript adj for Holm-adjusted pairwise tests. Boxes show type-7 quartiles, medians and observed whisker endpoints within 1.5 IQR. Outliers are shown by default. PDF, 600 dpi PNG, optional SVG, descriptive statistics, exact box statistics and R session information are exported.

## Scientific integrity and QA

Raw phenotype distributions and covariate/kinship-adjusted LMM P values are distinct quantities. State this in the figure legend. If comparing panels directly, set `SAME_Y_LIMIT <- TRUE`. Do not commit individual-level data or generated biological plots. This template has been reviewed statically; execution and visual QA with real CSV inputs remain required.

[Complete R script](02_boxplots.R)
