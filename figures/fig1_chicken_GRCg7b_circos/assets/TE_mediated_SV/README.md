# TE-mediated SV (TEV) schematic

A minimal one-row schematic showing the project's operational classification rule:

**For a called SV interval, if more than 80% of the SV sequence is covered by TE sequence, classify that SV as a TE-mediated variant (TEV).**

- The outlined bar represents the **SV interval / SV sequence**.
- The coral segment represents **TE-covered sequence (>80%)**.
- The small neutral segment represents **other sequence (<20%)**.
- The arrow indicates classification of the qualifying SV as **TEV**.

For quantitative implementation, calculate TE coverage from the union of TE-annotated bases within the SV interval so overlapping TE annotations are not double-counted. The SVG is a conceptual illustration, not a depiction of exact genomic coordinates or a specific observed SV.

The figure uses the repository's established coral accent (`#E85B4A`), dark slate typography (`#253246`), neutral grey, rounded modules, and a 1400 × 250 canvas matching the TE architecture schematics.
