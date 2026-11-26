# GWAS–molQTL colocalization module

Follow README.md and the editable style constants in plot_coloc_template.R. This module's approved continuous LD color semantics, 105 × 100 mm two-track layout, and two output versions apply here; chromosome colors and Manhattan dimensions from other modules do not apply.

Never rerun statistical models as a plotting side effect. Preserve original P, model, population, PIP, and posterior values. Separate regional H4, conditional focal probability, joint focal probability, and SuSiE PIP. Do not invent missing LD/PIP or apply a default significance threshold.

Use real population-specific dosage LD and record actual samples and effective N. Exact variant identity is mandatory. Keep all genuine association rows in the window. SNP < SV < TEV < focal, LD ascending within type; focal drawn once. Preserve missing/unavailable focal records explicitly.

Keep all biological inputs, sample IDs, genotypes, generated biological plots, and local absolute paths outside this public repository. Commit only reusable code, input contracts, and style documentation. Do not edit unrelated modules.
