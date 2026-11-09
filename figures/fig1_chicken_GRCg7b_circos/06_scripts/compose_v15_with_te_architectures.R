#!/usr/bin/env Rscript
# ==============================================================================
# Figure 1: immutable V15 SVG + four repository TE architecture SVGs.
# R-only compositor: original vector layers are embedded intact; grid/svglite
# generates a separate connector layer; librsvg exports vector PDF and PNG.
# No V15 plot code is sourced, executed or modified. No genomics is recomputed.
#
# Why not redraw the V15 SVG with grImport2::grid.picture()?
# Its default ext="none" drops complex clip paths; ext="clipbbox" only clips
# their bounding boxes. V15 relies on 164 real annular clipping components.
# Native SVG embedding preserves those clips and editable architecture text.
# ==============================================================================
options(stringsAsFactors = FALSE, scipen = 20)

script_project_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) return(dirname(dirname(normalizePath(sub("^--file=", "", arg[1])))))
  if (dir.exists("04_window_statistics")) return(normalizePath("."))
  normalizePath("chicken_7b_TE_circos", mustWork = TRUE)
}

composition_config <- function(project_dir = script_project_dir()) {
  stem <- paste0("chicken_GRCg7b_R_v15_1000kb_adaptive_discovery_",
                 "TIP-extended_SV-all_palette-pangenome_spectrum")
  list(
    project_dir = project_dir,
    source_dir = file.path(project_dir, "07_figures", "R_v15"),
    source_stem = stem,
    asset_dir = file.path(project_dir, "assets", "TE_architectures"),
    output_dir = file.path(project_dir, "07_figures", "R_v15_with_TE_architectures"),
    output_stem = "chicken_GRCg7b_R_v15_with_TE_architectures",
    page_width_in = 14, page_height_in = 10,
    dpi = 600, preview_dpi = 150,
    preview_only = FALSE,
    # These are physical inches, top-left origin. Scaling is uniform.
    main = list(left_in = 0.15, top_in = 0.15, height_in = 9.70),
    architectures = data.frame(
      track = c("LINE", "LTR", "SINE", "DNA"),
      file = c("LINE_CR1.svg", "LTR_ERV.svg", "SINE_generic.svg", "DNA_TIR.svg"),
      left_in = 10.0, top_in = c(1.45, 3.40, 5.35, 7.30), width_in = 3.95,
      bar_center_source_y = 125,
      color = c("#386C8C", "#7BB6AC", "#9689B9", "#B4AD74")),
    # A short common stem originates immediately inside the LOWER opening edge,
    # near the TE tracks. Four schematic branches fan out without crossing data.
    # They identify representative architectures, not genomic loci/annotations.
    connectors = list(
      anchor_radius = 0.72, anchor_inside_opening_deg = 1.5,
      stem_bend_in = c(7.65, 5.75), junction_in = c(8.00, 5.75),
      branch_bends_in = list(LINE = c(8.70, 3.05), LTR = c(8.70, 4.30),
                             SINE = c(8.75, 5.84), DNA = c(8.70, 6.95)),
      target_bend_x_in = 9.35, target_clearance_in = 0.18,
      color = "#111111", line_width = 1.25, arrow_length_in = 0.10,
      arrow_angle = 25),
    formats = c("svg", "pdf", "png")
  )
}

stopf <- function(...) stop(sprintf(...), call. = FALSE)
assert_file <- function(path) {
  if (!file.exists(path)) stopf("Missing input: %s", path)
  normalizePath(path, mustWork = TRUE)
}
fmt <- function(x) format(x, digits = 14, scientific = FALSE, trim = TRUE)
read_tsv <- function(path) utils::read.delim(assert_file(path), check.names = FALSE)
write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

# Minimal, lossless SVG envelope handling; no SVG paths or styles are rebuilt.
svg_attribute <- function(root, key) {
  rx <- paste0("\\b", key, "\\s*=\\s*(['\"])(.*?)\\1")
  hit <- regmatches(root, regexec(rx, root, perl = TRUE))[[1]]
  if (length(hit) < 3) return(NULL)
  hit[3]
}
read_svg <- function(path) {
  path <- assert_file(path)
  text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  hit <- regexpr("<svg\\b[^>]*>", text, perl = TRUE)
  if (hit[1] < 1) stopf("Invalid SVG root: %s", path)
  root <- regmatches(text, hit)
  box <- as.numeric(strsplit(svg_attribute(root, "viewBox"), "[ ,]+")[[1]])
  if (length(box) != 4 || any(!is.finite(box)) || any(box[3:4] <= 0))
    stopf("A finite positive viewBox is required: %s", path)
  body <- substr(text, hit[1] + attr(hit, "match.length"), nchar(text))
  body <- sub("</svg>\\s*$", "", body, perl = TRUE)
  if (grepl("<(image|script)\\b", body, perl = TRUE, ignore.case = TRUE))
    stopf("Inputs must be self-contained vectors without raster/script: %s", path)
  list(path = path, body = body, viewbox = box)
}
namespace_body <- function(body, prefix) {
  # Empty IDs in the supplied titles are removed only in the composed layer.
  # Original files remain byte-identical; all nonempty IDs/references get a prefix.
  body <- gsub("\\s+id\\s*=\\s*(['\"])\\1", "", body, perl = TRUE)
  body <- gsub("(\\bid\\s*=\\s*['\"])([^'\"]+)(['\"])",
               paste0("\\1", prefix, "\\2\\3"), body, perl = TRUE)
  body <- gsub("(\\b(?:xlink:)?href\\s*=\\s*['\"]#)",
               paste0("\\1", prefix), body, perl = TRUE)
  gsub("(url\\(\\s*['\"]?#)", paste0("\\1", prefix), body, perl = TRUE)
}
embed_svg <- function(source, id, left_pt, top_pt, width_pt, height_pt) {
  paste0('<svg id="', id, '" x="', fmt(left_pt), '" y="', fmt(top_pt),
    '" width="', fmt(width_pt), '" height="', fmt(height_pt),
    '" viewBox="', paste(fmt(source$viewbox), collapse = " "),
    '" preserveAspectRatio="xMinYMin meet" overflow="hidden">\n',
    namespace_body(source$body, paste0(id, "__")), "\n</svg>")
}

prepare_layout <- function(cfg, sources) {
  main <- sources$circos
  scale <- cfg$main$height_in * 72 / main$viewbox[4]
  main_place <- data.frame(id = "V15-circos", left_pt = cfg$main$left_in * 72,
    top_pt = cfg$main$top_in * 72, width_pt = main$viewbox[3] * scale,
    height_pt = main$viewbox[4] * scale, scale_x = scale, scale_y = scale)
  arch <- cfg$architectures
  arch$left_pt <- arch$left_in * 72
  arch$top_pt <- arch$top_in * 72
  arch$width_pt <- arch$width_in * 72
  arch$height_pt <- vapply(arch$track, function(t) {
    src <- sources[[t]]; arch$width_pt[arch$track == t] * src$viewbox[4]/src$viewbox[3]
  }, numeric(1))
  arch$scale_x <- arch$width_pt/vapply(arch$track, function(t) sources[[t]]$viewbox[3], numeric(1))
  arch$scale_y <- arch$height_pt/vapply(arch$track, function(t) sources[[t]]$viewbox[4], numeric(1))
  arch$target_y_in <- arch$top_in + arch$bar_center_source_y * arch$scale_y/72
  # Target lies before the first visible bar, not on top of an SVG label/shape.
  arch$target_x_in <- arch$left_in + 30 * arch$scale_x/72 - cfg$connectors$target_clearance_in
  if (any(abs(arch$scale_x - arch$scale_y) > 1e-10)) stopf("Architecture aspect ratio changed.")
  center <- read_tsv(file.path(cfg$source_dir, paste0(cfg$source_stem, ".center_legend.tsv")))[1, ]
  geometry <- read_tsv(file.path(cfg$source_dir, paste0(cfg$source_stem, ".figure_geometry.tsv")))[1, ]
  mapping <- read_tsv(file.path(cfg$source_dir, paste0(cfg$source_stem, ".chromosome_display_mapping.tsv")))
  if (isTRUE(geometry$clockwise)) stopf("This layout expects unchanged counterclockwise V15.")
  cx <- main_place$left_pt + center$circle_center_x_in/geometry$width_in * main$viewbox[3] * scale
  cy <- main_place$top_pt + (1-center$circle_center_y_in/geometry$height_in) * main$viewbox[4] * scale
  radius <- geometry$circle_diameter_in/geometry$height_in * main$viewbox[4]/2 * scale
  lower_edge <- mapping$end_degrees[mapping$chromosome == "W"]
  theta <- (lower_edge + cfg$connectors$anchor_inside_opening_deg) * pi/180
  anchor <- c(cx + radius*cfg$connectors$anchor_radius*cos(theta),
              cy - radius*cfg$connectors$anchor_radius*sin(theta))/72
  stem <- rbind(anchor, cfg$connectors$stem_bend_in, cfg$connectors$junction_in)
  branches <- lapply(arch$track, function(t) {
    a <- arch[arch$track == t, ]
    rbind(cfg$connectors$junction_in, cfg$connectors$branch_bends_in[[t]],
      c(cfg$connectors$target_bend_x_in, a$target_y_in), c(a$target_x_in, a$target_y_in))
  })
  names(branches) <- arch$track
  # Check every connector segment against the actual annular data region.
  samples <- function(m) do.call(rbind, lapply(seq_len(nrow(m)-1), function(i)
    cbind(seq(m[i,1],m[i+1,1],length.out=501),seq(m[i,2],m[i+1,2],length.out=501))))
  for (m in c(list(stem), branches)) {
    points <- samples(m)*72
    dx <- points[,1]-cx;dy <- cy-points[,2]
    rr <- sqrt(dx^2+dy^2)/radius;aa <- (atan2(dy,dx)*180/pi) %% 360
    inside_data <- rr >= 0.3185 & rr <= 0.8895 & aa >= 24 & aa <= 336
    if (any(inside_data)) stopf("A connector crosses a V15 data panel; adjust composition_config().")
    if (any(points[,1] < 0 | points[,1] > cfg$page_width_in*72 |
            points[,2] < 0 | points[,2] > cfg$page_height_in*72)) stopf("Connector outside page.")
  }
  list(main = main_place, architectures = arch, stem = stem, branches = branches,
    circle_center_pt = c(cx,cy), circle_radius_pt = radius, lower_edge_deg = lower_edge)
}

draw_connector_layer <- function(cfg, layout, filename) {
  svglite::svglite(filename, width = cfg$page_width_in, height = cfg$page_height_in,
                   bg = "transparent", pointsize = 10)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  polyline <- function(points, arrow = FALSE) grid::grid.polyline(
    x = grid::unit(points[,1]/cfg$page_width_in, "npc"),
    y = grid::unit(1-points[,2]/cfg$page_height_in, "npc"),
    gp = grid::gpar(col = cfg$connectors$color, fill = cfg$connectors$color,
      lwd = cfg$connectors$line_width, lineend = "butt", linejoin = "mitre"),
    arrow = if (arrow) grid::arrow(type = "closed", ends = "last",
      angle = cfg$connectors$arrow_angle,
      length = grid::unit(cfg$connectors$arrow_length_in, "inches")) else NULL)
  polyline(layout$stem)
  for (points in layout$branches) polyline(points, arrow = TRUE)
}

write_png <- function(svg, path, width_in, height_in, dpi) {
  image <- rsvg::rsvg_nativeraster(svg, width = as.integer(round(width_in*dpi)),
                                  height = as.integer(round(height_in*dpi)))
  png::writePNG(image, path, dpi = c(dpi,dpi))
  rm(image); invisible(gc())
}

compose_v15 <- function(cfg = composition_config()) {
  lib <- file.path(cfg$project_dir, "06_scripts", "composition_R_library")
  if (dir.exists(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
  packages <- c("svglite", "rsvg", "png")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stopf("Missing R packages: %s. Run 06_scripts/install_composition_R_dependencies.sh",
                            paste(missing,collapse=", "))
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  cfg$output_dir <- normalizePath(cfg$output_dir, mustWork = TRUE)
  protected_dir <- normalizePath(cfg$source_dir, mustWork = TRUE)
  if (cfg$output_dir == protected_dir || startsWith(cfg$output_dir, paste0(protected_dir,"/")))
    stopf("Composition outputs must not be written into the immutable R_v15 directory.")
  sources <- list(circos = read_svg(file.path(cfg$source_dir, paste0(cfg$source_stem, ".svg"))))
  for (i in seq_len(nrow(cfg$architectures)))
    sources[[cfg$architectures$track[i]]] <- read_svg(file.path(cfg$asset_dir,cfg$architectures$file[i]))
  # Protect all finalized V15 files and the main R/shell; archive files are not read.
  protected <- list.files(cfg$source_dir, full.names = TRUE, recursive = FALSE)
  protected <- protected[file.info(protected)$isdir == FALSE & !grepl("\\.zip$",protected)]
  protected <- c(protected, file.path(cfg$project_dir,"06_scripts",c("plot_chicken_circos_v15.R","run_circos_R_v15.sh")),
    vapply(sources,`[[`,character(1),"path"))
  protected <- c(protected, file.path(cfg$project_dir,"06_scripts",c(
    "compose_v15_with_te_architectures.R", "run_compose_v15_with_te_architectures.sh",
    "install_composition_R_dependencies.sh")))
  provenance_files <- file.path(cfg$asset_dir,c("README.md","source_manifest.json"))
  protected <- c(protected,provenance_files[file.exists(provenance_files)])
  protected <- unique(protected); before <- tools::md5sum(protected)
  layout <- prepare_layout(cfg, sources)
  prefix <- file.path(cfg$output_dir,cfg$output_stem)
  overlay_path <- file.path(cfg$output_dir,"connectors_only.svg")
  draw_connector_layer(cfg,layout,overlay_path)
  overlay <- read_svg(overlay_path)
  parts <- embed_svg(sources$circos,"V15-circos",layout$main$left_pt,layout$main$top_pt,
                     layout$main$width_pt,layout$main$height_pt)
  for (i in seq_len(nrow(layout$architectures))) {
    a <- layout$architectures[i,]
    parts <- c(parts,embed_svg(sources[[a$track]],paste0("architecture-",a$track),
      a$left_pt,a$top_pt,a$width_pt,a$height_pt))
  }
  overlay_body <- namespace_body(overlay$body,"connector__")
  svg <- paste0('<?xml version="1.0" encoding="UTF-8"?>\n',
    '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"',
    ' width="',fmt(cfg$page_width_in),'in" height="',fmt(cfg$page_height_in),
    'in" viewBox="0 0 ',fmt(cfg$page_width_in*72),' ',fmt(cfg$page_height_in*72),'">\n',
    '<title>Chicken GRCg7b V15 Circos with representative TE architectures</title>\n',
    '<desc>Unchanged V15 genomic panel; four original repository vector schematics; ',
    'four schematic arrows fan from the lower-right opening. Not sequence to scale.</desc>\n',
    '<rect width="100%" height="100%" fill="#FFFFFF"/>\n',paste(parts,collapse="\n"),
    '\n<g id="connector-layer" class="svglite">\n',overlay_body,'\n</g>\n</svg>')
  svg_path <- paste0(prefix,".svg")
  writeLines(svg,svg_path,useBytes=TRUE)
  write_png(svg_path,paste0(prefix,".preview.png"),cfg$page_width_in,cfg$page_height_in,cfg$preview_dpi)
  if (!cfg$preview_only) {
    if ("pdf" %in% cfg$formats) rsvg::rsvg_pdf(svg_path,paste0(prefix,".pdf"))
    if ("png" %in% cfg$formats) write_png(svg_path,paste0(prefix,".png"),
      cfg$page_width_in,cfg$page_height_in,cfg$dpi)
  }
  after <- tools::md5sum(protected)
  if (!identical(before,after)) stopf("A protected V15/source file changed during composition.")
  manifests <- data.frame(path=names(before),bytes=file.info(names(before))$size,
                          md5=unname(before),unchanged=unname(before==after))
  write_tsv(manifests,file.path(cfg$output_dir,"input_preservation.tsv"))
  write_tsv(layout$main,file.path(cfg$output_dir,"circos_placement.tsv"))
  write_tsv(layout$architectures,file.path(cfg$output_dir,"architecture_placements.tsv"))
  paths <- c(list(stem=layout$stem),layout$branches)
  connector_rows <- do.call(rbind,lapply(names(paths),function(n) data.frame(
    connector=n,vertex=seq_len(nrow(paths[[n]])),x_in=paths[[n]][,1],y_from_top_in=paths[[n]][,2],
    arrowhead=n!="stem" & seq_len(nrow(paths[[n]]))==nrow(paths[[n]]))))
  write_tsv(connector_rows,file.path(cfg$output_dir,"connector_vertices.tsv"))
  dput(cfg,file=file.path(cfg$output_dir,"composition_settings.R"))
  dput(layout[c("circle_center_pt","circle_radius_pt","lower_edge_deg")],
       file=file.path(cfg$output_dir,"source_geometry.R"))
  writeLines(capture.output(sessionInfo()),file.path(cfg$output_dir,"sessionInfo.txt"))
  outputs <- list.files(cfg$output_dir,full.names=TRUE)
  outputs <- outputs[file.info(outputs)$isdir==FALSE & !grepl("(^output_manifest[.]tsv$|[.]log$)",basename(outputs))]
  write_tsv(data.frame(path=basename(outputs),bytes=file.info(outputs)$size,
    md5=unname(tools::md5sum(outputs))),file.path(cfg$output_dir,"output_manifest.tsv"))
  message("Done: ",cfg$output_dir,"; unchanged protected inputs: ",length(before))
  invisible(list(configuration=cfg,layout=layout,prefix=prefix))
}

parse_args <- function(args,cfg) {
  i <- 1L
  while(i<=length(args)) {
    a <- args[i]
    if(a=="--help") {
      cat("Usage: Rscript 06_scripts/compose_v15_with_te_architectures.R [--project-dir PATH] [--out-dir PATH] [--preview-only] [--dpi 600]\n")
      quit(status=0)
    }
    if(a=="--preview-only") {cfg$preview_only<-TRUE;i<-i+1L;next}
    if(i==length(args)) stopf("Missing argument for %s",a)
    value<-args[i+1L]
    if(a=="--project-dir") cfg<-composition_config(normalizePath(value,mustWork=TRUE))
    else if(a=="--out-dir") cfg$output_dir<-value
    else if(a=="--dpi") {cfg$dpi<-as.integer(value);if(is.na(cfg$dpi)||cfg$dpi<72)stopf("Invalid dpi")}
    else stopf("Unknown argument: %s",a)
    i<-i+2L
  }
  cfg
}
if(sys.nframe()==0L) compose_v15(parse_args(commandArgs(TRUE),composition_config()))
