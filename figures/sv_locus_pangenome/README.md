# SV locus pangenome visualization | ODGI + gggenes

## Scope
Original method specification for Codex; no upstream library source is vendored. Produce a reference-anchored regional overview and a breakpoint-scale structural haplotype view for a selected chicken SV/TEV. This module is independent of Manhattan plotting. Do not modify other figure modules or publish actual research data in this public repository.

## Required inputs
- One graph: assembly-path GFA (GFA1 P-lines preferred), ODGI .og, VG .gbz, or .xg. Prefer original assembly graph, not a VCF-only graph, when comparing assemblies.
- Reference build and FASTA .fai; exact reference path name and mapping from VCF/GFF chromosome to graph path.
- Target: VCF record ID, CHROM, POS, REF, ALT, SVTYPE, END/SVLEN as applicable. Resolve multiallelic ALT explicitly. Default flank = 500000 bp per side; insertion uses POS anchor, spanning event uses both breakpoints. Clip to chromosome bounds. VCF POS is 1-based; BED and ODGI ranges are 0-based half-open.
- Reference gene GFF3/GTF or BED6 (unique feature IDs, gene names, strand); optional assembly-specific gene annotation, TE BED and ALT sequence/repeat annotations.
- Path manifest: exact path name, assembly, haplotype, population, source (assembly / phased observed genotype / synthetic / unknown), and coordinate offset. Original FASTA/assembly alignments required if actual assembly paths are absent. No genotype inference from graph topology alone.

## 1. Inspect and prepare graph
Run on a compute node under the local scheduler, not on login/storage node. Record software versions, help text, input checksums, reference build and all commands. Check free scratch and memory; reuse cached per-chromosome OG where available. The following are templates, not hardcoded research paths.

```bash
# If GFA already has compatible P-lines, skip conversion.
vg convert -fW input.gbz > graph.gfa
# Same conversion can be attempted for XG; absent paths cannot be recovered by conversion.
odgi build -g graph.gfa -o graph.og -t 16
odgi paths -i graph.og -L > all_paths.txt
```

Check vg/ODGI version and help before using options. Verify that graph path names, lengths, reference sequence and provenance match the original assemblies. GBZ may contain synthetic or panel-derived haplotypes; XG may omit haplotype threads and GBZ may lack unsupported ALT branches. Do not claim the graph contains 51 assembly paths without verification. If only VCF graph exists, label figures as variant-graph/reference/alternate paths, not assembly haplotypes. Use original graph plus reliable assembly-path mappings for genuine assembly plots.

## 2. Extract two windows
Set REF to exact graph path, START0 and END0 to validated coordinates from the target plus flanks. For a single-base anchor at POS, start0=max(0,POS-1-500000); end0=min(chr_len,POS+500000). For long DEL/DUP/INV include both breakpoints plus flanks. For BND handle each linked locus separately. Record requested and actual extracted ranges; node edges may extend beyond requested window.

```bash
odgi extract -i graph.og -r "${REF}:${START0}-${END0}" -E -o region.og -t 16
odgi paths -i region.og -L > region.paths.txt
odgi view -i region.og -g > region.gfa
```

-E includes nodes between selected nodes in graph sort order; only use with suitable ordering, preferably reference-guided. It does not guarantee that every alternate allele is retained. Check actual ALT sequence/traversal, both breakpoint anchors, graph connectivity, path continuity, and reference coverage. If missing, use validated path lacing (-R), graph context (-c or -L) or explicit target nodes with a documented strategy; do not silently drop branches or merge disconnected fragments. Preserve an unmodified region graph. For the breakpoint panel extract a second region with 2–10 kb flanks (adjust to event length and repeats). Do not draw all SNP bubbles across 1 Mb.

## 3. Untangle against actual local reference
Select the exact extracted reference path from region.paths.txt; path may have a coordinate suffix. Prepare query_paths.txt with exact verified local path names (no annotation paths, no synthetic paths in assembly mode). Keep original offsets; local starts are not genomic coordinates.
```bash
odgi untangle -i region.og -r "$LOCAL_REF" -Q query_paths.txt -t 16 > region.untangle.bedpe
```
Inspect column schema: query name/start/end, reference name/start/end, score, relative orientation, self coverage, rank (verify installed version). Preserve orientation and ambiguous mappings. Reference-relative blocks are segment homology, not base-resolution alignments. Unmapped gaps are not automatically insertions/deletions. Do not treat self coverage as exact copy number without independent evidence. Use query path coordinates for true linear gene arrows; for the shared reference axis project validated intervals/anchors explicitly and distinguish off-reference sequence.

## 4. Project genes to paths (optional but recommended)
Create region.genes.bed with full reference graph path names, 0-based intervals and collision-free feature IDs in column 4. Preserve biological strand in a metadata table. Select features overlapping the region; flag partial genes at boundaries. Convert GFF gene start to start-1; do not use transcript/exon rows as separate genes unless desired.
```bash
odgi procbed -i region.og -b region.genes.bed > region.genes.local.bed
odgi inject -i region.og -b region.genes.local.bed -o region.genes.og
cut -f4 region.genes.local.bed | sort -u > gene_paths.txt
odgi untangle -i region.genes.og -Q query_paths.txt -R gene_paths.txt -j 0.5 -t 16 -g > region.gggenes.tsv
```
Check installed ODGI supports -g (gggenes output), procbed and inject. The official injection tutorial documents -g; if unavailable, use BEDPE and implement a validated converter, never assume BEDPE columns are gggenes columns. Inspect actual output header. Confirm annotation target IDs exist and query names still match after injection. Strand in graph projection is relative orientation: combine with reference gene strand if annotation path was injected in reference-forward orientation; test on known minus-strand genes. If paths are reversed for visual consistency, retain original coordinates and recompute feature and target positions; never flip a true local inversion away. The 0.5 Jaccard cutoff is a starting heuristic, not a biological confidence threshold; inspect ambiguous paralogs and sensitivity. Projected reference annotations are not independently verified assembly annotations.

## 5. Identify target allele
Match the exact VCF record to reference and alternate sequence/traversal and breakpoint anchors. Verify sample/assembly membership from real path traversal or assembly alignments. For insertion, identify the actual ALT nodes/sequence and its query-local coordinates; reference-only BED cannot describe novel sequence. For deletion, show breakpoint junction and verified missing reference interval, not a fabricated block. Retain unknown/partial/ambiguous paths as unknown; no absence call from missing annotation or no hit.
Write target_alleles.tsv: path_id, allele, query_start0, query_end0, ref_anchor0, evidence, status; for a deletion junction use equal query start/end. Validate spanning flanks. If TEV, include family, orientation, length, and annotation evidence only if supported. Distinguish number of assembly paths vs individuals vs population genotype frequencies.

## 6. Plot design
Panel A — default overview ±500 kb: shared reference genomic coordinate in Mb; reference gene track; verified homologous blocks or gene tracks for selected assembly paths; soft accent around target interval and explicit breakpoint marker. Show target insertion as an anchored marker, not a reference-length interval. Optional faint synteny ribbons (gggenomes).
Panel B — breakpoint zoom: independently validated query-local structural paths, with shared flanking anchors and accurate sequence lengths; REF grey-blue, target ALT terracotta; inversion orientation visible; TE substructure only from actual annotations. Label local coordinate axes clearly; if schematic spacing used, label not to scale. No unverified haplotype labels.
Panel C optional — compact graph topology from Bandage/ODGI, no fake topology. A gene-arrow map alone is not a graph topology visualization.

Main text: choose representative *verified* structural patterns (typically 3–8), state grouping rule, show n and denominator only when computed. Complete assembly rows in supplement, with adaptive height. If gene content is identical, do not force artificial gene arrows to distinguish SV alleles: use actual breakpoints/variant blocks.

### R rendering contract
Use ggplot2 + gggenes, optional gggenomes/patchwork. Normalize validated arrow TSV to columns: molecule, start, end, feature_id, feature_label, forward (TRUE/FALSE biological direction in displayed coordinates), feature_class. Normalize path backbone TSV: molecule, start, end. Normalize target markers separately, joined by exact path IDs; never paint a target on every path by default. In overview use shared reference coordinate projection; in zoom use oriented query coordinates with explicit anchor alignment/offset. No mixing coordinate systems.
Minimal original plotting idiom (after validation and coordinate handling):
```r
library(ggplot2)
library(gggenes)
library(grid)
genes <- read.delim("genes.normalized.tsv", check.names=FALSE)
backbone <- read.delim("backbone.normalized.tsv")
stopifnot(all(genes$forward %in% c(TRUE,FALSE)), all(genes$end > genes$start))
p <- ggplot() +
  geom_segment(data=backbone,
    aes(x=start,xend=end,y=molecule,yend=molecule),
    linewidth=0.35,colour="#9DA9B2") +
  geom_gene_arrow(data=genes,
    aes(xmin=start,xmax=end,y=molecule,fill=feature_class,forward=forward),
    arrowhead_height=unit(2,"mm"),arrowhead_width=unit(1.1,"mm"),
    arrow_body_height=unit(1.6,"mm"),colour="#273746",linewidth=0.2) +
  scale_fill_manual(values=c(Gene="#578A9E",Other="#83A895",Target="#C36B4B")) +
  theme_genes() +
  theme(text=element_text(family="sans",size=8),legend.position="top",
        panel.grid=element_blank()) +
  labs(x="Validated display coordinate (bp)",y=NULL,fill=NULL)
# Add verified target interval/junction layer separately, annotate reference anchor.
# Choose axis labels and limits based on actual coordinate system, not this generic example.
ggsave("locus.pdf",p,width=183,height=max(65,16+5*nrow(backbone)),units="mm",bg="white")
ggsave("locus.png",p,width=183,height=max(65,16+5*nrow(backbone)),units="mm",dpi=600,bg="white")
```
Above is a rendering fragment, not an end-to-end executable pipeline. Codex must implement data adapters, explicit factor order, target markers, axis offsets, and QC based on actual inputs. Use vector PDF and SVG (svglite if installed), 600 dpi PNG. Suggested colors: background #FFFFFF, ink #273746, backbone #9DA9B2, gene #578A9E, contextual feature #83A895, focal SV #C36B4B, unknown #BFC6CB. Use 7–9 pt final labels; avoid overlabeling; no large title, no major grid, no excessive decoration.

## 7. Acceptance checks
- Exact record ID/ALT, build, path and window coordinates recorded; reference extracted sequence matches FASTA; coordinate round-trip including minus strand and offset checked.
- Confirm both breakpoint anchors, target sequence/traversal and path coverage; distinguish REF / ALT / UNKNOWN. Compare representative sequences with original FASTA or independent alignments when available.
- Preserve ambiguous/partial/repeated features; no unsupported copy number or assembly labels. Report missing paths/genes and extraction truncation.
- Target appears in overview and zoom, supported rows only; labels and axes consistent; no clipped features, unreadable labels or false linear scaling.
- Save versions, commands, path manifest, locus metadata, original BEDPE, normalized plotting tables, QC TSV, and figures locally. Record tests and failures. Never upload genotypes/assemblies/results to this public style repository.
- Pilot one simple insertion and one deletion/inversion before batch use; test a synthetic graph with known REF/ALT and a negative-strand gene. This specification has not been executed on the user's graph.

## References
- ODGI extraction: https://odgi.readthedocs.io/en/latest/rst/commands/odgi_extract.html
- ODGI untangle: https://odgi.readthedocs.io/en/latest/rst/commands/odgi_untangle.html
- Official gene projection tutorial: https://odgi.readthedocs.io/en/latest/rst/tutorials/injecting_gene_arrows.html
- Official C4 untangle tutorial: https://odgi.readthedocs.io/en/latest/rst/tutorials/untangling_the_pangenome.html
- gggenes: https://github.com/wilkox/gggenes
- VG graph formats: https://github.com/vgteam/vg
- Chicken graph pangenome, K locus Fig. 3: https://doi.org/10.1186/s12915-023-01758-0
- HPRC complex loci: https://doi.org/10.1038/s41586-023-05896-x

No upstream repository is copied; dependencies are installed by the research environment.
