# TE-mediated structural variation schematic

This directory contains an editable SVG schematic illustrating the project's operational TE-overlap rule for classifying structural variants.

## Classification rule

For each called SV interval:

`TE coverage (%) = 100 × (number of unique bases in the SV interval overlapped by one or more TE annotations) / SV length`

- **TE-mediated SV:** TE coverage is strictly **greater than 80%**.
- **Does not meet the criterion:** TE coverage is **80% or less**.

The numerator is the **union** of TE-annotated bases within the SV interval. Overlapping TE annotations must not be double-counted. The denominator is the SV length as defined by the project's SV callset and must be applied consistently across variant classes.

## Figure contents

- Panel A: 85 kb TE-covered out of a 100 kb SV interval (85%); meets the rule.
- Panel B: 65 kb TE-covered out of a 100 kb SV interval (65%); does not meet the rule.
- TE colors follow the repository's pangenome spectrum palette: LINE `#386C8C`, LTR `#7BB6AC`, SINE `#9689B9`, DNA `#B4AD74`.
- Segment order is illustrative. Segment widths represent proportions and do not imply genomic order, orientation, or a particular SV mechanism.

## Interpretation caveat

This is an operational **sequence-overlap classification**. A high TE-overlap fraction supports describing an SV as TE-mediated under the project's rule, but does not by itself establish the molecular mechanism that generated the SV. Mechanistic claims require additional evidence (for example, breakpoint structure, target-site duplications, TE insertion signatures, or homology between breakpoints).

The SVG is composed of native editable vector shapes and text. No raster artwork or embedded fonts are used.
