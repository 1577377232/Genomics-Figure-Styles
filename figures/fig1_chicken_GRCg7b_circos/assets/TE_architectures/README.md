# Chicken Figure 1 — editable TE structure callouts

These are **original vector schematics** intended for the right side of the existing **V15 GRCg7b Circos** panel. They borrow the visual language (modular pastel segments and simple arrows) of the published human ATEC figure, but **not** its artwork or code. Every segment, arrow and text item remains an editable SVG element. Structures are schematic, not sequence-to-scale.

## Assets

- `LINE_CR1.svg` — representative full-length chicken CR1-like LINE: 5′ UTR, ORF1, ORF2, 3′ UTR, terminal repeats. Chicken CR1 does not have the universal human L1-style Poly(A) tail.
- `LTR_ERV.svg` — an ERV-like LTR element; Gag/Pol, with `Env*` shown as an **optional** component. The full LTR annotation class also contains Gypsy and other subtypes.
- `SINE_generic.svg` — a **schematic tRNA-type example** with an internal Pol III A/B-box promoter. Chicken SINE annotations also include MIR/5S/other families, which do not necessarily share this exact arrangement. This is **not** the two-monomer human Alu cartoon.
- `DNA_TIR.svg` — a representative TIR-type DNA transposon. The source project's RepeatMasker normalization keeps Helitron/RC in a separate class, not this DNA track.
- `TE_architectures_4up.svg` — a convenient, still-editable four-row composition with four individual connector-arrow groups; place, ungroup and adjust in Adobe Illustrator.

All four correspond to the confirmed **V15** track palette:
`LINE #386C8C`, `LTR #7BB6AC`, `SINE #9689B9`, `DNA #B4AD74`; text/outline `#253246`. White/transparent artwork; editable Liberation Sans/Arial/Helvetica SVG text; no embedded font or third-party image.

## Recommended Illustrator workflow

1. Open the **V15 SVG or vector PDF** from `../../07_figures/R_v15/` in Illustrator as the immutable left-hand Circos panel.
2. Set a wider artboard (roughly **1.55–1.75:1**). Keep the Circos circle on the left and the existing 48° opening facing right. Do not stretch/distort the circle.
3. Place the four independent SVGs at the right, align all structural bars, distribute vertically; alternatively open the `TE_architectures_4up.svg` stack and ungroup.
4. On a separate **Connectors** layer, draw four simple dark arrows from the circle's opening region to the four structures. Adjust the existing stack connectors to align with the final figure; they are **schematic**, not genomic position mappings.
5. Keep at least 2–3 mm white clearance around the Circos right-edge tick labels and density axes. Do not hide underlying data.
6. Save the editable `.ai` master; export vector PDF/SVG and 600 dpi PNG. Use a consistent font substitution if Liberation Sans is not installed. Do not rasterize SVG structures or redistribute fonts.

## Automated assembly alternative

R can create these diagrams using `grid.rect`, `grid.roundrect`, `grid.text` and `grid.lines` within a right-hand viewport after drawing the V15 circlize panel. The V15 PDF/SVG can also be combined with these SVGs through an SVG compositor. **Do not change V15 genomics inputs or track calculations merely to add decorative callouts.** For one final publication panel, Illustrator is more practical; for repeated automated export, use code-based composition.

## Scientific limitations and attribution

The shapes show **representative TE architectures**, not intact structures of all copies annotated in the reference. Cite chicken CR1 structural work and the Fourth Report on Chicken Genes and Chromosomes for biological context. The original human paper is only visual inspiration. No existing V15 file, figure or statistics are replaced.
