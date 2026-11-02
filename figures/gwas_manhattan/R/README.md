# Joint SNP + TE Manhattan plots in R

Executable R implementation for 89 quantitative traits and seven frozen multivariate modules. Common packages: `ggplot2`, `data.table`, `patchwork`, `ragg`, `png`, `jsonlite`, `digest`. Tested with R 4.1.2 and ggplot2 4.0.3 on Linux; use ggplot2 >=3.5 for inside-panel legends. Scientific inputs and completed biological figures stay in the research workspace.

## Explicit user preset

The supplied preset overrides the parent style defaults at the user's request:

- SNP `#386C8C`, TE `#E75A46`, from `pangenome_spectrum`; all TE types share one colour and circle shape.
- Equal 1.10 mm circles, chromosome gaps of 12 Mb, clear axes, white background.
- Biological chromosomes 1–37 and Z; chromosome 38 and W are excluded from these figures. Exclusion counts remain in the audit.
- One dashed reference line at P=1e-5, without a numerical label. This is a plotting reference, not a universal significance or Bonferroni threshold.
- 260 x 75 mm per panel, PNG 600 dpi, editable PDF text and TE points. Every valid requested SNP point reaches `geom_point`; only that dense layer is rasterized. No marker sampling, bin maxima, min-P merge or arbitrary peak cap.

## Model and Z provenance

Single-trait panels require native GEMMA additive Wald `p_wald`; module panels require native multivariate `p_score`. A module is a joint multivariate phenotype matrix, not an average or principal-component score. The plotting program never fits association models or substitutes FarmCPU/nonadditive P values.

Set `Z_analysis` to `joint` only after fitting an actual joint-sex model. Each trait/marker pair must then have exactly one `A` and one `Z_joint` source, with identical complete-case N. Female/male stratified P values cannot be concatenated, minimized or combined and described as a joint fit.

The research Z rerun retained corrected DNA sex, PC1–10, numeric batch and the full autosomal GRM; the conditional carcass module also retained log live weight. Its explicit primary coding was female hemizygous 0/2 standardized dosage and male 0/1/2. This estimates a common standardized dosage effect; it does not establish biological complete dosage compensation. Female heterozygous calls were masked before QC, never randomly resolved. No validated project PAR interval was supplied, so that masking policy is conservative throughout Z. SNP/TE pooled-panel MAF thresholds were .05/.01 and maximum missingness .1; phenotype-specific native filters can still give different tested Z counts. Such differences must be reported rather than filled with synthetic points.

## Run with local inputs

Install R packages into your normal library, or set `GWAS_R_LIB` to a custom library directory. For the optional workbook/gallery verification, install Python `pandas` and `openpyxl` and the Poppler commands `pdfinfo`, `pdffonts`, `pdfimages`, `pdftotext`.

Create an output directory with `config/analysis.json`, `config/traits.tsv`, `config/sources.tsv` and `config/user_style.json`. Start with the header-only templates and the supplied preset. Supply actual files and hashes; templates contain no fabricated research results.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript plot_joint_manhattan.R --outdir /path/to/local/output --trait-id REAL_ID --jobs 1 --no-overview
Rscript plot_joint_manhattan.R --outdir /path/to/local/output --trait-id all --jobs 30
python verify_plot_delivery.py --outdir /path/to/local/output
python tests/test_joint_manhattan.py
```

Each worker is single-threaded. Thirty workers were used for the research batch. Budget RAM for a complete SNP table per worker; reduce concurrency for smaller machines. CLI also accepts `--config`, `--manifest`, `--sources`, `--style`, `--fai`, `--threshold-p`, `--assemble-only`, and `--no-overview`. If changing the reference line with the CLI, update the style file too before running delivery verification.

`traits.tsv` contains unique real IDs, true bilingual labels, categories, `kind=single|module`, module dimensions and QC flags. `sources.tsv` contains one source row per trait, marker class and panel; it records N, P column, reported/valid counts, minimum P, source SHA256, method and model. Required association columns are `chr`, `ps`, and the declared P column; TE also needs `rs`, and univariate TE needs finite `beta` and positive `se`. Source files must be tab-delimited text readable by `data.table::fread`, optionally gzip compressed.

Supply the exact reference FAI and SHA256, biological chromosome order, exclusions and any explicitly verified aliases. Numeric chromosome codes are not universally mapped to Z: biological chromosome 39 and PLINK code 39 must be distinguished for each project. Record model, covariates, GRM and reference provenance in `analysis.json`.

## Outputs and checks

Outputs include `traits/` and `modules/` PDFs/PNGs, separate 89-page and seven-page PDF booklets, grouped overview pages, source/chromosome counts, source hashes, true minimum P, per-source and combined-marker Bonferroni metadata, R session information and failure receipts. Per-trait Bonferroni metadata uses all valid source tests, including explicitly hidden chromosome 38 tests; only the user reference line is drawn. There is no correction across all traits. Resume fingerprints include code, style, configuration, threshold override, source rows and file metadata; changed plotted sources are rehashed on rendering. The audit does not claim that GWAS association proves causality or colocalization.

`verify_plot_delivery.py` checks source counts, one line, equal circle sizes, no chromosome 38, pooled-Z structure, no P underflow or clipping, DPI, output hashes, PDF text/fonts and SNP-only rasterization. It writes a local Excel catalogue, searchable HTML gallery and final PASS receipt. Use a full batch when requesting booklet/page verification.

The end-to-end tests generate tiny, conspicuously labelled `TEST_ONLY` fixtures in temporary directories. They check univariate/multivariate column handling, biological chromosome order, cumulative coordinates, explicit exclusions, counts, minimum P, styling, resume and P=0 rejection. They are not biological previews and are never used in the real figures.

Keep association matrices, individual-level data, actual source manifests, genotype files and biological output figures outside this public repository.
