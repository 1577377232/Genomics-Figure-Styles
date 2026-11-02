# Codex implementation brief — 89 joint SNP + TE Manhattan plots

Use this style repository `1577377232/Genomics-Figure-Styles` as the source of truth. Read `AGENTS.md`, `STYLE_GUIDE.md`, `theme.json`, `palettes.json` and `trait_manifest.template.csv` before coding. Implement and execute the workflow in the **user's data/workspace repository**, not here; do not upload raw GWAS matrices.

1. Inventory the real 89 trait SNP/TE GWAS results, all file paths, columns, TE subtype annotations, exact reference genome .fai, tested variant counts, model/software, sample size and covariate provenance. Never invent these.
2. Build a robust Python CLI `plot_joint_manhattan.py` with `--manifest`, `--fai`, `--theme`, `--palette`, `--outdir`, `--trait-id`, `--preflight`, `--jobs`, `--strict`, `--overwrite`. Support common TSV, CSV, compressed and parquet formats where supported.
3. Validate 89 *real* unique trait IDs and accurate trait labels/groups. If absent, generate only an incomplete draft manifest and report missing inputs. Never fabricate 89 trait results.
4. Build common cumulative genome coordinates from verified FAI (chr1–39, Z, W only when present). Never change real physical chromosome width. Avoid microchromosome tick collisions without dropping data.
5. SNP layer: six-color cycle from selected palette, small translucent dots, rasterized in PDF. TE layer: larger INS coral diamond, DEL violet triangle and explicitly labeled unknown/other amber. Highlight large TE only if a validated per-family or joint threshold is supplied. Make TE/SNP legends accurate.
6. Preserve source P/negative-log P (no cap; handle zero underflow as a QC failure unless original precise source score is available). Log rejected invalid inputs. Do not downsample final GWAS data.
7. Never assume 5e-8 significance. Document SNP, TE and optional joint Bonferroni/multiple-testing family as actually selected, with tested variant counts. If unknown, omit horizontal line.
8. Render 89 PDF + 89 PNG at 600 dpi, grouped overview and optional fixed-y comparison versions. All axes and text publication quality; ensure only dense SNP scatter is rasterized in PDF.
9. Write `outputs/audit/trait_qc.tsv`, `plot_manifest.json`, `failures.tsv`, and reproducible execution log; output files must stay ignored by Git.
10. Add tests (synthetic test-only data) for chromosome order, coordinate mapping, invalid P, threshold families, variant type styling, correct output counts and file provenance.
11. Preflight 1–2 *real* traits, visually verify PDF and PNG, then run 89 only when inputs and QC pass.

Final deliverables: reproducible CLI and tests, 89 real plots if data are present, scientific audit reports and clear explanation of any missing data. Do not misrepresent synthetic previews as real findings.
