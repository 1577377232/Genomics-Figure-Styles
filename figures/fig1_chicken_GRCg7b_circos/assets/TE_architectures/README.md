# Chicken Figure 1 — editable TE architecture schematics

These four **editable vector SVGs** are designed as the right-hand structural callouts for the V15 GRCg7b Circos figure. The visual style follows the supplied clean, modular schematic preview: bold class title, rounded pastel modules, consistent dark text, and no footnotes beneath the structures.

## Files

- `LINE_CR1.svg` — representative chicken CR1-like LINE architecture: 5′ UTR, ORF1, ORF2, 3′ UTR and a terminal 3′ repeat region. This is not a human L1 schematic.
- `LTR_ERV.svg` — representative ERV-like LTR architecture: LTR–Gag–Pol–Env*–LTR. The dashed `Env*` module indicates that Env is optional / lineage-dependent.
- `SINE_generic.svg` — simplified generic chicken SINE architecture: 5′ head, internal Pol III promoter, 3′ tail and a variable terminal region. It is not intended to represent every chicken SINE family or the human Alu structure.
- `DNA_TIR.svg` — representative TIR-type DNA transposon architecture: TIR–Transposase–TIR. It does not represent Helitron rolling-circle elements.

- `TE_architectures_4up.svg` — all four editable structures arranged vertically for preview or reuse.

## Palette

The colors follow the repository's `pangenome_spectrum` mapping:
- LINE: `#386C8C`
- LTR: `#7BB6AC`
- SINE: `#9689B9`
- DNA: `#B4AD74`
- Text: `#253246`
- Neutral terminal modules: `#A1A8BE`

Module fills are light tints of the corresponding palette color; outlines and titles use the full class color. The diagrams are schematic and are **not drawn to sequence scale**.

## Editing and assembly

All titles, labels, module rectangles, and outlines are native SVG elements grouped by function. They can be edited in R, Inkscape, or Adobe Illustrator. Use the four individual SVGs with the independent V15 R composition script; do not modify the V15 Circos data or plotting script just to alter these schematic callouts.

No external image, rasterized text, or embedded font is used.
