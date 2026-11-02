# Scientific figure style and statistical integrity — v2.0

## Visual language
White background, no grids, no chromosome stripe shading, 3D, gradients or drop shadows. Default `pangenome_spectrum` in `palettes.json`; can switch to `botanical_genome` or `chromatic_nc`. Six distinct but coordinated SNP colors cycle by biological chromosome order. Color maps **chromosome identity**, never P value. SNP scatter: small alpha ~0.52 circular markers. TE insertion: coral diamond; deletion: violet triangle; other/unknown: amber. Draw foreground TEs above rasterized SNP cloud. Render nonsignificant TE small and semitransparent, enlarge *only* above validated thresholds. Use shape as redundant encoding for accessibility.

## Dimensions and typography
Double-column 183 × 76 mm. Arial/Helvetica with Liberation Sans/DejaVu Sans fallback; final figure text ~6.5–7 pt, minimum 5 pt. Export vector PDF with editable text and TE, selectively rasterize dense SNP collections; PNG 600 dpi. Use x-axis Chromosome and y-axis −log₁₀(P), concise sentence-case true trait title. Small nonoverlapping legends, left/bottom spines only, no giant 'A' on single panels.

## Genome positions and labels
Require validated .fai for the actual reference build; default biological preferred order chr1–39, Z, W **only when present**. Use FAI chromosome lengths for physical cumulative positions and stable offsets across all 89 phenotypes. Do not give equal widths to chromosomes. Microchromosome tick labels may be omitted by collision detection, but points must all remain. Reject unknown chromosomes, invalid positions or mixed reference builds until resolved.

## Data/statistics safeguards
Visually combined association results do not prove shared causality or colocalization. Record sample size, phenotype name/units, traits, SNP and TE GWAS software/model, covariates and GRM. Handle original P or -log10P. Reject invalid P; for P=0 caused by machine underflow use source negative-log P if available, else flag (do not cap). Every threshold must be explicitly supplied or derived under a declared alpha/correction family and tested-variant counts. Distinguish SNP and TE family Bonferroni alpha/M from optional deliberately defined joint alpha/(M_snp+M_te) family. Never default to P=5e-8 and never copy another paper's dashed line.

## Batching and audit
Exactly 89 **real** phenotype IDs in completed manifest, not fabricated codes. Maintain true labels (muscle weight versus muscle percentage must not be interchanged). Preserve one joint PDF and one PNG per trait, plus grouped overview pages by real category, optional separately labeled fixed-y comparison pages. Under ignored outputs/audit/ include trait_qc.tsv, plot_manifest.json, failures.tsv, software/model provenance and input hashes. Keep raw genotype, association datasets and phenotype files outside this style repository.

## Mandatory QA
Validate exact plotted marker counts and min P against input; no missing Z/W, wrong chromosome order, clipped peaks, legend obstruction or microchromosome label overlap; verify TE shapes/colors in grayscale; inspect PDF text vectors and actual PNG DPI. Preflight 1–2 real phenotypes before all 89.

This is a self-developed publication-inspired design system, not an official NG/NC/Advanced Science template.
