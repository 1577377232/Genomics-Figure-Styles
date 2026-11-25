# Compact SNP / SV / TEV profile (2026-11-25)

This profile is the current default for joint SNP/SV/TEV Manhattan plots. It follows the existing qualitative-trait figures, with a shorter 183 × 68 mm canvas. Earlier chromosome-cycle and TE-subtype styles remain available as legacy profiles.

## Rendering
- White background; left and bottom spines; English trait title, axes and SNP/SV/TEV legend only.
- No footer, counts, software/model caption, candidate circle, peak position label or leader line. Keep provenance in the adjacent JSON/TSV and encode software/model in filenames.
- SNP blue #386C8C; SV lavender #9689B9; TEV coral #E75A46. All circles, without selectively enlarging lead points.
- Use actual FAI lengths and a common biological chromosome order. Preserve every valid point; tick omission prevents label collisions and never removes markers.
- Preserve the existing 12 Mb interchromosome spacing; add 12 Mb beyond the left edge of chromosome 1 and the right edge of the terminal chromosome (Z in this dataset).
- PDF text remains editable. All dense variant collections are rasterized at 600 dpi, matching the qualitative-plot treatment. This profile overrides the legacy TE-vector requirement.
- Automatic y limits include every peak, with space above the maximum. No peak capping.
- Thresholds must be explicitly supplied and documented. The project's P=1e-5 line is a candidate display reference, not an automatically calibrated genome-wide significance threshold. The generic renderer accepts threshold_p=None to omit it.

## TEV classification
A variant is TEV when the union of TE-derived aligned bases covers **at least 50%** of its trimmed changed allele sequence. Merge overlapping intervals before dividing, remove shared REF/ALT prefix and suffix, and count genuine TE families separately from simple repeats, satellites, RNA and unknown repeats.

Join annotations using chromosome, coordinate and **complete REF and tested ALT SHA256**, not coordinate alone. Match the annotation inventory to the actual unfiltered GWAS marker set: reusing a smaller GQ/DP-filtered annotation table silently loses TEVs. Supplement with independently verified candidate annotations where available. Keep unknown/unassessed coverage as missing; it remains in the SV display layer but must be counted separately and is not evidence that a variant is non-TE.

Counts are per association result set and model, not a universal property of the entire catalogue. SV and TEV plotted counts are mutually exclusive and sum to the valid structural-variant rows.

## Trait labels
Translate the recorded measurement and denominator, preserving weight versus yield. Examples: Eviscerated weight (g), Semi-eviscerated weight (g), Eviscerated yield, Breast muscle yield, Leg muscle yield, Abdominal fat percentage, Shank and foot weight (g), and Breast angle (°). Preserve the actual denominator in the manifest even when omitted from the short title. An organ sum including proventriculus must explicitly include it rather than using an ambiguous shortened label. Confirm ages from local phenotype definitions before using “at 7 days” or “at 84 days”.

Terminology references: [chicken carcass definitions](https://pmc.ncbi.nlm.nih.gov/articles/PMC9495107/), [carcass traits and chest angle](https://pmc.ncbi.nlm.nih.gov/articles/PMC12138444/). Definitions and denominators in the user's source data take precedence over conventions in these papers.

## Implementation and reuse
Use python/manhattan_clean.py::draw with the actual SNP chromosome codes (Z=40), physical positions, P values, SV table (chrom, position, P, TE_coverage), verified FAI-derived layout, true trait title and output PDF path. Pass the documented threshold explicitly. Optional regional=(chromosome, center_position) produces an unannotated ±2.5 Mb regional panel. Generate PNG from PDF at 600 dpi. Reuse identical trait/model panels across candidate folders and regenerate combined atlases from the clean source panels; never just cover old annotations with white rectangles.

Audit exact plotted counts, min P, unknown annotations, dimensions, padding, threshold provenance and input hashes. Require the actual manifest size, whether 22 candidate records or 89 phenotypes; never fabricate missing records. No individual-level or raw association data belong in this repository.
