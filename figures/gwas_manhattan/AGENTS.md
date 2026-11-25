

## Current compact SNP/SV/TEV profile (2026-11-25)
For joint SNP/SV/TEV plots, use the `compact_snp_sv_tev` profile in `theme.json` and [CLEAN_PROFILE.md](CLEAN_PROFILE.md). It takes precedence over legacy chromosome-cycle colors, TE subtype shapes, fixed batch size, and TE-vector rendering for this profile. It uses 183 × 68 mm, blue/lavender/coral variant classes, outer chromosome padding, English trait titles and no peak annotations or explanatory footers. TEV means TE-only coverage **≥50%**, joined to the actual GWAS marker inventory by exact complete REF/ALT identity. Unknown coverage is separately audited. The earlier specifications below remain available for other plot styles.
# Codex instructions — chicken pangenomics figures

Read STYLE_GUIDE.md, theme.json, palettes.json, COLOR_SEMANTICS.md, LITERATURE_REFERENCES.md and CODEX_TASK.md before writing any figures. theme.json is the source of layout and the named palettes in palettes.json are the source of colors. Default to pangenome_spectrum. Preserve 89 real phenotype names; do not fabricate inputs, reference genome, association results, thresholds or provenance. No arbitrary peak capping or GWAS marker downsampling in final results. Keep data in the user's working data repository; this repository holds only style and code specifications. Do not change unrelated files.
