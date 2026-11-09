#!/usr/bin/env Rscript
# ==============================================================================
# Figure 1 v3_r2: immutable V15 + five linked modules; TEV illustration >80%, immutable V15 data >50%.
# R compositor: original vector layers remain intact; R writes rounded native
# SVG connector paths, and librsvg exports the same vector scene as PDF/PNG.
# No V15 plot code is sourced, executed or modified. No genomics is recomputed.
#
# Why not redraw the V15 SVG with grImport2::grid.picture()?
# Its default ext="none" drops complex clip paths; ext="clipbbox" only clips
# their bounding boxes. V15 relies on complex annular clipping components.
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
    source_dir = file.path(project_dir, "07_figures", "R_v15"), source_stem = stem,
    asset_dir = file.path(project_dir, "assets"),
    output_dir = file.path(project_dir,"07_figures","R_v15_with_TE_architectures_v3_r2"),
    output_stem = "chicken_GRCg7b_R_v15_with_TE_architectures_v3_r2",
    page_width_in = 17, page_height_in = 10, dpi = 600, preview_dpi = 150,
    preview_only = FALSE,
    main = list(left_in = 0.15, top_in = 0.15, height_in = 9.70),
    architectures = data.frame(
      track = c("LINE","LTR","SINE","DNA","TEV"),
      file = c("TE_architectures/LINE_CR1.svg","TE_architectures/LTR_ERV.svg",
        "TE_architectures/SINE_generic.svg","TE_architectures/DNA_TIR.svg",
        "TE_mediated_SV/TE_mediated_SV_classification_v3_r2.svg"),
      left_in = 11.75, top_in = c(1.50,2.95,4.40,5.85,7.30), width_in = 4.9375,
      bar_center_source_y = rep(167,5),
      bar_left_source_x = c(63.5,138.5,167,148,110),
      connector = rep(TRUE,5),
      color = c("#386C8C","#7BB6AC","#9689B9","#B4AD74","#E75A46")),
    # Exact upper Chr1 sector boundary for all five true track origins.
    # Native cubic launches avoid the immutable axis labels, then feed five
    # independent smooth fan-out curves in the clear corridor to the right.
    connectors = list(
      tracks = c("LINE","LTR","SINE","DNA","TEV"),
      anchor_track_fraction = c(LINE=0.65,LTR=0.65,SINE=0.65,DNA=0.65,TEV=0.25),
      launch_end_degrees = c(LINE=21.5,LTR=21,SINE=20,DNA=20,TEV=0),
      launch_tangent_in = c(LINE=0.08,LTR=0.08,SINE=0.08,DNA=0.08,TEV=0.26),
      launch_extension_in = c(LINE=0.15,LTR=0.15,SINE=0.15,DNA=0.15,TEV=0.35),
      fan_start_x_in = c(LINE=9.50,LTR=9.40,SINE=9.30,DNA=9.20,TEV=8.90),
      fan_end_x_in = c(LINE=11.30,LTR=11.30,SINE=11.30,DNA=11.30,TEV=11.30),
      fan_control_fraction = 0.35,
      target_clearance_in = 0.15, color = "#253246", line_width = 0.80,
      arrow_length_in = 0.075, arrow_angle = 25,
      annotation_clearance_pt = 1.2, minimum_pair_clearance_pt = 3.5),
    debug = list(zoom_x_in=5.35,zoom_y_in=3.05,zoom_width_in=3.55,zoom_height_in=2.30),
    formats = c("svg","pdf","png")
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

# Geometric tests operate on the final page coordinates, not genome data.
segment_crosses <- function(a,b,c,d) {
  cross <- function(x,y,z) (y[1]-x[1])*(z[2]-x[2])-(y[2]-x[2])*(z[1]-x[1])
  boxes <- max(min(a[1],b[1]),min(c[1],d[1])) <= min(max(a[1],b[1]),max(c[1],d[1])) &&
           max(min(a[2],b[2]),min(c[2],d[2])) <= min(max(a[2],b[2]),max(c[2],d[2]))
  boxes && cross(a,b,c)*cross(a,b,d)<=0 && cross(c,d,a)*cross(c,d,b)<=0
}
segment_distance <- function(a,b,c,d) {
  if(segment_crosses(a,b,c,d)) return(0)
  point_distance <- function(p,a,b) {
    v<-b-a; t<-max(0,min(1,sum((p-a)*v)/sum(v*v)))
    sqrt(sum((p-(a+t*v))^2))
  }
  min(point_distance(a,c,d),point_distance(b,c,d),point_distance(c,a,b),point_distance(d,a,b))
}
annotation_boxes <- function(cfg,main_place,source) {
  files<-paste0(cfg$source_stem,".",c("chromosome_labels.tsv","coordinate_labels.tsv",
    "TE_axis_labels.tsv","variant_axis_labels.tsv"))
  xx<-list(); ux<-source$viewbox[3]/9.8*main_place$scale_x
  uy<-source$viewbox[4]/10*main_place$scale_y
  for(file in files) {
    d<-read_tsv(file.path(cfg$source_dir,file))
    if(!"center_x_in" %in% names(d))next
    if("label_status" %in% names(d))d<-d[!d$label_status %in% c("hidden","omitted"),]
    if(!"tangent_x" %in% names(d)){d$tangent_x<-1;d$tangent_y<-0;d$radial_x<-0;d$radial_y<-1}
    xx[[file]]<-data.frame(label=d$label,
      x=main_place$left_pt+d$center_x_in*ux,y=main_place$top_pt+(10-d$center_y_in)*uy,
      hw=d$width_in*ux/2+cfg$connectors$annotation_clearance_pt,
      hh=d$height_in*uy/2+cfg$connectors$annotation_clearance_pt,
      tx=d$tangent_x,ty=-d$tangent_y,rx=d$radial_x,ry=-d$radial_y)
  }
  do.call(rbind,xx)
}
prepare_layout <- function(cfg,sources) {
  main<-sources$circos;scale<-cfg$main$height_in*72/main$viewbox[4]
  mp<-data.frame(id="V15-circos",left_pt=cfg$main$left_in*72,top_pt=cfg$main$top_in*72,
    width_pt=main$viewbox[3]*scale,height_pt=main$viewbox[4]*scale,scale_x=scale,scale_y=scale)
  arch<-cfg$architectures
  arch$left_pt<-arch$left_in*72;arch$top_pt<-arch$top_in*72;arch$width_pt<-arch$width_in*72
  arch$scale_x<-arch$width_pt/vapply(arch$track,function(t)sources[[t]]$viewbox[3],numeric(1))
  arch$height_pt<-arch$scale_x*vapply(arch$track,function(t)sources[[t]]$viewbox[4],numeric(1))
  arch$scale_y<-arch$scale_x
  arch$target_y_in<-arch$top_in+arch$bar_center_source_y*arch$scale_y/72
  arch$target_x_in<-arch$left_in+arch$bar_left_source_x*arch$scale_x/72-cfg$connectors$target_clearance_in
  center<-read_tsv(file.path(cfg$source_dir,paste0(cfg$source_stem,".center_legend.tsv")))[1,]
  geometry<-read_tsv(file.path(cfg$source_dir,paste0(cfg$source_stem,".figure_geometry.tsv")))[1,]
  mapping<-read_tsv(file.path(cfg$source_dir,paste0(cfg$source_stem,".chromosome_display_mapping.tsv")))
  tracks<-read_tsv(file.path(cfg$source_dir,paste0(cfg$source_stem,".track_geometry.tsv")))
  if(isTRUE(geometry$clockwise))stopf("Expected unchanged counterclockwise V15.")
  cx<-mp$left_pt+center$circle_center_x_in/geometry$width_in*main$viewbox[3]*scale
  cy<-mp$top_pt+(1-center$circle_center_y_in/geometry$height_in)*main$viewbox[4]*scale
  radius<-geometry$circle_diameter_in/geometry$height_in*main$viewbox[4]/2*scale
  upper<-mapping$start_degrees[mapping$chromosome=="1"]
  lower<-mapping$end_degrees[mapping$chromosome=="W"]
  branches<-list();anchors<-list();curves<-list()
  for(t in cfg$connectors$tracks) {
    source_track<-if(t=="TEV")"TIP" else t
    tr<-tracks[tracks$track==source_track,]
    if(nrow(tr)!=1L)stopf("Expected exactly one geometry row for %s",source_track)
    f<-cfg$connectors$anchor_track_fraction[[t]]
    r<-tr$inner_radius+f*(tr$outer_radius-tr$inner_radius)
    # This is the actual x=0 edge of Chr1, not a manually shifted gap point.
    theta<-upper*pi/180
    start<-c(cx+radius*r*cos(theta),cy-radius*r*sin(theta))/72
    endtheta<-cfg$connectors$launch_end_degrees[[t]]*pi/180
    launch_polar<-c(cx+radius*r*cos(endtheta),cy-radius*r*sin(endtheta))/72
    launch_end<-launch_polar+c(cfg$connectors$launch_extension_in[[t]],0)
    control1<-start+cfg$connectors$launch_tangent_in[[t]]*c(sin(theta),cos(theta))
    control2<-launch_polar
    a<-arch[arch$track==t,]
    fan_start<-c(cfg$connectors$fan_start_x_in[[t]],launch_end[2])
    fan_end<-c(cfg$connectors$fan_end_x_in[[t]],a$target_y_in)
    target<-c(a$target_x_in,a$target_y_in)
    span<-fan_end[1]-fan_start[1]
    fc1<-fan_start+c(span*cfg$connectors$fan_control_fraction,0)
    fc2<-fan_end-c(span*cfg$connectors$fan_control_fraction,0)
    branches[[t]]<-rbind(start,launch_end,fan_start,fan_end,target)
    curves[[t]]<-smooth_connector(start,control1,control2,launch_end,fan_start,fc1,fc2,fan_end,target)
    anchors[[t]]<-data.frame(track=t,source_track=source_track,source_track_index=tr$track_index,
      source_chromosome="1",source_position_bp=0,
      source_inner_radius=tr$inner_radius,source_outer_radius=tr$outer_radius,
      anchor_radius=r,anchor_track_fraction=f,anchor_degrees=upper,anchor_side="upper",
      reference_edge_degrees=upper,angular_offset_from_edge_deg=0,offset_from_track_edge_mm=0,
      x_in=start[1],y_from_top_in=start[2],origin_above_circle_center=start[2]*72<cy,
      exact_sector_endpoint=TRUE)
  }
  labels<-annotation_boxes(cfg,mp,main)
  for(t in names(curves)) {
    p<-curves[[t]]$samples
    rr<-sqrt((p[,1]-cx)^2+(cy-p[,2])^2)/radius
    aa<-(atan2(cy-p[,2],p[,1]-cx)*180/pi)%%360
    if(any(rr>=0.3185 & rr<=0.9975 & aa>upper+1e-6 & aa<lower-1e-6))stopf("%s rounded connector crosses main panel.",t)
    for(j in seq_len(nrow(labels))) {
      dx<-p[,1]-labels$x[j];dy<-p[,2]-labels$y[j]
      hit<-abs(dx*labels$tx[j]+dy*labels$ty[j])<labels$hw[j] &
           abs(dx*labels$rx[j]+dy*labels$ry[j])<labels$hh[j]
      if(any(hit))stopf("%s rounded connector covers label %s",t,labels$label[j])
    }
  }
  pairs<-list()
  for(i in seq_len(length(curves)-1L))for(j in (i+1L):length(curves)) {
    a<-curves[[i]]$samples;b<-curves[[j]]$samples
    dist<-min(sqrt(outer(a[,1],b[,1],"-")^2+outer(a[,2],b[,2],"-")^2))
    if(dist<cfg$connectors$minimum_pair_clearance_pt)stopf("%s/%s curve clearance %.3f pt is too small.",names(curves)[i],names(curves)[j],dist)
    pairs[[length(pairs)+1L]]<-data.frame(track1=names(curves)[i],track2=names(curves)[j],sampled_curve_clearance_pt=dist,crossing=FALSE,shared_vertex=FALSE)
  }
  list(main=mp,architectures=arch,branches=branches,curves=curves,anchors=do.call(rbind,anchors),
    connector_pair_QC=do.call(rbind,pairs),circle_center_pt=c(cx,cy),circle_radius_pt=radius,
    upper_edge_deg=upper,lower_edge_deg=lower)
}
# Continuous tangent launch + horizontal feeder + smooth fan-out + arrow lead.
# All positions are in physical inches; path geometry is emitted in SVG points.
smooth_connector <- function(start,c1,c2,launch,fan_start,f1,f2,fan_end,target) {
  coords<-list(start,c1,c2,launch,fan_start,f1,f2,fan_end,target)
  coords<-lapply(coords,function(z)z*72)
  pt<-function(q)paste(fmt(q),collapse=" ")
  d<-paste("M",pt(coords[[1]]),"C",pt(coords[[2]]),pt(coords[[3]]),pt(coords[[4]]),
    "L",pt(coords[[5]]),"C",pt(coords[[6]]),pt(coords[[7]]),pt(coords[[8]]),"L",pt(coords[[9]]))
  cubic<-function(a,b,c,d) {
    t<-seq(0,1,length.out=401)
    outer((1-t)^3,a)+outer(3*(1-t)^2*t,b)+outer(3*(1-t)*t^2,c)+outer(t^3,d)
  }
  line<-function(a,b) {t<-seq(0,1,length.out=201);outer(1-t,a)+outer(t,b)}
  samples<-rbind(cubic(coords[[1]],coords[[2]],coords[[3]],coords[[4]]),
    line(coords[[4]],coords[[5]])[-1,],cubic(coords[[5]],coords[[6]],coords[[7]],coords[[8]])[-1,],
    line(coords[[8]],coords[[9]])[-1,])
  list(d=d,samples=samples)
}

draw_connector_layer <- function(cfg,layout,filename) {
  parts<-character()
  for(t in names(layout$branches)) {
    points<-layout$branches[[t]]; route<-layout$curves[[t]]
    # Preserve a horizontal final segment; use a small closed triangular arrow.
    end<-points[nrow(points),]*72; len<-cfg$connectors$arrow_length_in*72
    half<-len*tan(cfg$connectors$arrow_angle*pi/180)
    triangle<-rbind(end,end+c(-len,-half),end+c(-len,half))
    parts<-c(parts,sprintf('<path id="connector-%s" d="%s" fill="none" stroke="%s" stroke-width="%s" stroke-linecap="round" stroke-linejoin="round"/>',t,route$d,cfg$connectors$color,fmt(cfg$connectors$line_width)),
      sprintf('<polygon id="arrow-%s" points="%s" fill="%s"/>',t,paste(apply(triangle,1,function(x)paste(fmt(x),collapse=",")),collapse=" "),cfg$connectors$color))
  }
  writeLines(paste0('<svg xmlns="http://www.w3.org/2000/svg" width="',cfg$page_width_in,'in" height="',cfg$page_height_in,'in" viewBox="0 0 ',cfg$page_width_in*72,' ',cfg$page_height_in*72,'">\n',paste(parts,collapse="\n"),'\n</svg>'),filename)
}

make_debug_svg <- function(svg,cfg,layout) {
  parts<-character()
  for(i in seq_len(nrow(layout$anchors))) {
    a<-layout$anchors[i,];xy<-c(a$x_in,a$y_from_top_in)*72
    color<-cfg$architectures$color[cfg$architectures$track==a$track]
    parts<-c(parts,sprintf('<circle cx="%s" cy="%s" r="4.2" fill="%s" stroke="#FFFFFF" stroke-width="1"/>',fmt(xy[1]),fmt(xy[2]),color),
      sprintf('<text x="%s" y="%s" text-anchor="middle" dominant-baseline="central" fill="#FFFFFF" font-family="Arial, sans-serif" font-size="5.5" font-weight="700">%s</text>',fmt(xy[1]),fmt(xy[2]),i))
  }
  # A key in the blank page header gives explicit track IDs without obscuring anchors.
  parts<-c(parts,'<g font-family="Arial, sans-serif" font-size="8" fill="#273746">',
    '<text x="680" y="23" font-weight="700">DEBUG: true upper Chr1 track endpoints (x=0 bp)</text>')
  for(i in seq_len(nrow(layout$anchors))) {
    a<-layout$anchors[i,]
    parts<-c(parts,sprintf('<text x="680" y="%s">%s  %s — source track %s (%s), angle 24°, r=%s</text>',
      24+i*12,i,a$track,a$source_track_index,a$source_track,fmt(a$anchor_radius)))
  }
  parts<-c(parts,'</g>')
  sub('</svg>\\s*$',paste0('<g id="DEBUG-anchor-markers">',paste(parts,collapse="\n"),'</g>\n</svg>'),svg,perl=TRUE)
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
  packages <- c("rsvg", "png")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stopf("Missing R packages: %s. Run 06_scripts/install_composition_R_dependencies.sh",
                            paste(missing,collapse=", "))
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  cfg$output_dir <- normalizePath(cfg$output_dir, mustWork = TRUE)
  protected_dir <- normalizePath(cfg$source_dir, mustWork = TRUE)
  if (cfg$output_dir == protected_dir || startsWith(cfg$output_dir,paste0(protected_dir,"/")) ||
      cfg$output_dir %in% file.path(dirname(protected_dir),c("R_v15_with_TE_architectures","R_v15_with_TE_architectures_v2","R_v15_with_TE_architectures_v3")))
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
    "compose_v15_with_te_architectures_v3_r2.R", "run_compose_v15_with_te_architectures_v3_r2.sh",
    "install_composition_R_dependencies.sh")))
  provenance_files <- file.path(cfg$asset_dir,c("TE_architectures/README.md","TE_architectures/source_manifest.json",
    "TE_architectures/TE_architectures_4up.svg","TE_mediated_SV/README.md",
    "TE_mediated_SV/TE_mediated_SV_classification_v3_r2.svg","TE_mediated_SV/TE_mediated_SV_classification_V15_gt50.provenance.json"))
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
    '<title>Chicken GRCg7b V15 Circos with four TE architectures and TEV classification</title>\n',
    '<desc>Unchanged V15 genomic panel; five repository-derived vector schematics; ',
    'five independent smooth arrows originate at exact upper Chr1 endpoints of LINE, LTR, SINE, DNA and innermost TEV tracks. TEV schematic shows the user-requested >80% rule; V15 data retain >50%. Not sequence to scale.</desc>\n',
    '<rect width="100%" height="100%" fill="#FFFFFF"/>\n',paste(parts,collapse="\n"),
    '\n<g id="connector-layer" class="svglite">\n',overlay_body,'\n</g>\n</svg>')
  svg_path <- paste0(prefix,".svg")
  writeLines(svg,svg_path,useBytes=TRUE)
  write_png(svg_path,paste0(prefix,".preview.png"),cfg$page_width_in,cfg$page_height_in,cfg$preview_dpi)
  debug_svg<-make_debug_svg(svg,cfg,layout)
  writeLines(debug_svg,paste0(prefix,".DEBUG.svg"),useBytes=TRUE)
  write_png(paste0(prefix,".DEBUG.svg"),paste0(prefix,".DEBUG.preview.png"),
    cfg$page_width_in,cfg$page_height_in,cfg$preview_dpi)
  crop<-cfg$debug
  zoom<-sub('viewBox="0 0 [^"]+"',paste0('viewBox="',paste(fmt(c(crop$zoom_x_in,crop$zoom_y_in,crop$zoom_width_in,crop$zoom_height_in)*72),collapse=" "),'"'),debug_svg)
  writeLines(zoom,paste0(prefix,".DEBUG.anchor_zoom.svg"),useBytes=TRUE)
  write_png(paste0(prefix,".DEBUG.anchor_zoom.svg"),paste0(prefix,".DEBUG.anchor_zoom.png"),
    7.1,4.6,300)
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
  write_tsv(layout$anchors,file.path(cfg$output_dir,"connector_anchors.tsv"))
  write_tsv(layout$connector_pair_QC,file.path(cfg$output_dir,"connector_pair_QC.tsv"))
  write_tsv(layout$main,file.path(cfg$output_dir,"circos_placement.tsv"))
  write_tsv(layout$architectures,file.path(cfg$output_dir,"architecture_placements.tsv"))
  paths <- layout$branches
  connector_rows <- do.call(rbind,lapply(names(paths),function(n) data.frame(
    connector=n,vertex=seq_len(nrow(paths[[n]])),x_in=paths[[n]][,1],y_from_top_in=paths[[n]][,2],
    arrowhead=n!="stem" & seq_len(nrow(paths[[n]]))==nrow(paths[[n]]))))
  write_tsv(connector_rows,file.path(cfg$output_dir,"connector_vertices.tsv"))
  write_tsv(data.frame(connector=names(layout$curves),svg_path=vapply(layout$curves,`[[`,character(1),"d")),
    file.path(cfg$output_dir,"connector_paths.tsv"))
  samples<-do.call(rbind,lapply(names(layout$curves),function(t) {
    m<-layout$curves[[t]]$samples
    data.frame(connector=t,sample=seq_len(nrow(m)),x_pt=m[,1],y_from_top_pt=m[,2])
  }))
  write_tsv(samples,file.path(cfg$output_dir,"connector_curve_samples.tsv"))
  dput(cfg,file=file.path(cfg$output_dir,"composition_settings.R"))
  dput(layout[c("circle_center_pt","circle_radius_pt","upper_edge_deg","lower_edge_deg")],
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
      cat("Usage: Rscript 06_scripts/compose_v15_with_te_architectures_v3_r2.R [--project-dir PATH] [--out-dir PATH] [--preview-only] [--dpi 600]\n")
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
