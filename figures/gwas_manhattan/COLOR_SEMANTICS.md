# Color semantics across a multi-panel genetics paper

- **SNP Manhattan:** six categorical hues distinguish chromosomes; use real FAI-defined physical widths.
- **TE:** INS = coral diamond; DEL = purple triangle; other/unknown = amber marker. Maintain consistent identity in regional GWAS, molecular QTL, LD, fine-mapping and colocalization panels. Distinguish marker types by shape as well as color.
- **Associations:** effect significance is encoded only by vertical -log10(P) position, not by hue. No falsely enhanced lead points.
- **Statistical evidence:** GWAS peaks, colocalization H4, credible set PIP, and functional validation are separate evidence layers; do not apply identical glyphs to imply equal confidence.
- **SNP background:** small, translucent and below TE in z-order; foreground TE must remain legible without masking the genome-wide point density.
- **Thresholds:** gray dashed, lightly styled, and only shown with documented provenance. SNP-specific, TE-specific and joint tests are different families and must be labeled.
- **Special loci:** optional dark outline only for explicitly declared candidates; log criteria and do not recolor by desired biological outcome.
- **Journal:** neutral white background, clean typography, no rainbow continuous scales, no non-data-driven decorative panels.
