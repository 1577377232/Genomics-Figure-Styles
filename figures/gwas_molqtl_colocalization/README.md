# GWAS–molQTL colocalization

Editable **ggplot2 + patchwork** regional association figures with separate annotated and clean outputs. This module defines plotting style and input contracts; biological datasets and generated figures remain in the user's analysis workspace.

## Run

Install `ggplot2`, `patchwork`, `jsonlite`, `scales`, and `svglite`; use a Cairo-capable R build and Liberation Sans.

```bash
Rscript plot_coloc_template.R /path/to/case_directory
```

Each case directory contains `data/` and optionally `case_config.json`. Figures and audits are written to `figures/` within that case. Edit the style constants at the top of the script; selected labels can also be overridden in the case configuration. Inputs are never changed. No GWAS, molQTL, SuSiE, or coloc fit is run.

## Style

- 105 × 100 mm, white background, Liberation Sans, ink `#253246`, left and bottom axes.
- Two vertically aligned regional association tracks, sharing physical coordinates and panel width. Default: focal position ±1 Mb, total 2 Mb; x labels only on the bottom track.
- Continuous fill represents **population-specific dosage LD r²**, anchored at the same focal biological SV. Stops: 0 `#386C8C`, .25 `#6681A2`, .5 `#9689B9`, .75 `#C6A8A0`, 1 `#E5C65F`. Missing LD stays NA, displayed light gray.
- Shapes: SNP 21, SV 22, TEV 23. Sizes: 1.15, 1.40, 1.65; focal coral diamond 3.0. Thin dark outlines. Main scatter fill never encodes variant class.
- Explicit layer order: SNP → SV → TEV → focal. Within each type, missing LD first and finite LD ascending. The focal is excluded from its ordinary layer and drawn once, last.
- Right-side shared legends: Variant Type above vertical LD r² colorbar.
- Annotated: truthful trait/population/model labels, original P values, regional H4, focal joint shared posterior, and available original PIP with its fitted window. A faint focal line is enabled by default.
- Clean: axes identifying GWAS and molQTL, legends, and focal highlight; sample and statistic text removed, focal line off. No overall title.
- H0–H4 bar chart is a separate supplement. No significance line unless its actual threshold and source are supplied.

## Inputs

UTF-8 TSVs may be gzip compressed. Default filenames are shown below. Keep original analysis records, including variants outside the shared colocalization set.

| File | Required fields / role |
|---|---|
| `data/GWAS.tsv.gz` | `variant_id`, `variant_key`, `position`, `P`, `type` |
| `data/molQTL.tsv.gz` | `variant_id`, `variant_key`, `position_mol`, `pval_g1`, `type_mol` |
| `data/GWAS_dosage_LD.tsv.gz` | `variant_id`, `variant_key`, `LD_r2`, `LD_N_complete`, `LD_status` |
| `data/C367_dosage_LD.tsv.gz` | molecular-population LD; same columns; historical default name, actual population comes from case metadata |
| `data/SV_TE_coverage.tsv.gz` | `marker`, `TE_coverage`; exact complete allele identity must have been verified upstream; TEV coverage ≥.5 |
| `data/配对信息.tsv` | one row: `coordinate` (`chrom:pos`), `GWAS_model`, `molecular_model`, `scope`, `kind`, `N_GWAS`, `N_mol`, `focal_GWAS_P`, `focal_mol_P`, `PPH4`, `focal_conditional_shared_PP`, `focal_joint_shared_PP` |
| `data/case.json` | existing result, including `focal_variant` and `PPH0`…`PPH4`; posterior values must sum to 1 |
| `data/focal_PIP.json` | optional existing result: `variant_id`, `PIP`, `model`, `half_window_bp`, `original_GWAS_P`; omission means unavailable, never zero |

GWAS and molecular LD must use their **own actual analysis individuals** and their own actual genotype version. Compute Pearson r² from 0/1/2 dosage using pairwise complete individuals; separately document any haploid 0/1 or 0/2 coding. Do not call REC/DOM predictor correlations dosage LD. Retain sample lists, sample-order hashes, genotype provenance, and per-variant effective N outside the public repository. Exact IDs or chromosome/position/complete alleles are required; never merge on position alone.

A focal absent from a particular original analysis is not added as a fake point; LD unavailable in its native genotype panel remains missing. A display window change does not refit colocalization or generate association records beyond the supplied input range.

## Statistical labels

`PPH4` is the regional shared-signal posterior. `focal_conditional_shared_PP` is the focal allocation conditional on H4. `focal_joint_shared_PP = PPH4 × focal_conditional_shared_PP`. SuSiE PIP is a separate result with its original model and input window. Never substitute these quantities for one another. Preserve any original model mismatch and sparse variant coverage in the case notes.

## Outputs and checks

`colocalization_annotated` and `colocalization_clean`: PDF with vector text/points, SVG, 600 dpi PNG, and small previews. A separate `*_H0_H4_supplement` is also exported. Original input P values and posterior values are not modified; invalid P values are flagged rather than capped.

Audits contain per-variant LD and effective N, exact plotting order, input/plotted/type/focal/missing counts for both versions, original focal statistics, panel alignment, and R session information. All original association rows in the requested window are retained, independent of LD availability or membership in the coloc overlap set.

## Batch organization

Recommended local layout: `locus/gene__gene_id/population_kind_GWAS-model_mol-model_feature_id/`. Preserve every original comparison in the gene table. Give each distinct model/population/feature combination a separate plot pair; deduplicate only identical combinations. Best regional, best same-model, and all passing posterior pairs are **selection categories**, not independent replications or causal proof. Report selection rules and coverage gaps.

## References

The implementation is independent; it does not import these packages or query human genomes or 1000 Genomes:

- [gassocplot2](https://github.com/jrs95/gassocplot2): stacked regional association layout.
- [locuscomparer](https://github.com/boxiangliu/locuscomparer): paired regional association presentation and focal marker.
- [locuszoomr](https://github.com/myles-lewis/locuszoomr): LD point ordering and locus display.
- Repository [genotype–phenotype boxplot](../genotype_phenotype_boxplot): typography, P-value e notation, and export conventions.

The two-version template was executed on a real chicken CR1L–COG3 case; individual-level inputs and generated biological figures are deliberately not committed here.
