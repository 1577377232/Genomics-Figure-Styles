#!/usr/bin/env Rscript
# ==============================================================================
# Figure 1 v3_r3: immutable V15 + five linked modules; TEV illustration >80%, immutable V15 data >50%.
# R compositor: original vector layers remain intact; R writes straight native
# SVG M/L connector paths, and librsvg exports the same vector scene as PDF/PNG.
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
    output_dir = file.path(project_dir,"07_figures","R_v15_with_TE_architectures_v3_r3"),
    output_stem = "chicken_GRCg7b_R_v15_with_TE_architectures_v3_r3",
    page_width_in = NA_real_, page_height_in = 10, dpi = 600, preview_dpi = 150,
    preview_only = FALSE,
    main = list(left_in = 0.15, top_in = 0.15, height_in = 9.70),
    architectures = data.frame(
      track = c("LINE","LTR","SINE","DNA","TEV"),
      file = c("TE_architectures/LINE_CR1_composition_1_5x.svg","TE_architectures/LTR_ERV_composition_1_5x.svg",
        "TE_architectures/SINE_generic_composition_1_5x.svg","TE_architectures/DNA_TIR_composition_1_5x.svg",
        "TE_mediated_SV/TEV_classification_composition_1_5x.svg"),
      left_in = NA_real_, top_in = NA_real_, width_in = 7.40625,
      bar_center_source_y = rep(167,5),
      bar_left_source_x = c(63.5,138.5,167,148,110),
      bar_outline_source_width = rep(5,5),
      connector = rep(TRUE,5),
      color = c("#386C8C","#7BB6AC","#9689B9","#B4AD74","#E75A46")),
    # All true anchors remain on the upper Chr1 endpoints of real tracks.
    # User clarified: offset into the gap along the END-FACE NORMAL by 0.14in,
    # not radially from the circle center (which would touch adjacent tracks).
    layout = list(module_gap_in=0.70,module_canvas_gap_in=0.08,
      max_column_to_circos_height=0.92,right_margin_in=0.20,stroke_allowance_pt=0.5),
    connectors = list(
      tracks = c("LINE","LTR","SINE","DNA","TEV"),
      anchor_track_fraction = c(LINE=0.65,LTR=0.65,SINE=0.65,DNA=0.65,TEV=0.25),
      start_clearance_in = 0.14,clearance_direction="end_face_normal_into_opening",
      first_elbow_offset_from_bound_in = c(LINE=0.12,LTR=0.17,SINE=0.12,DNA=0.02,TEV=-0.13),
      first_leg_drop_in = c(LINE=0,LTR=0,SINE=0.30,DNA=0.60,TEV=0.80),
      elbow_offset_from_bound_in = c(LINE=0.42,LTR=0.49,SINE=0.56,DNA=0.63,TEV=0.70),
      target_clearance_in = 0.09, color = "#253246", line_width = 0.85,
      arrow_length_in = 0.065, arrow_angle = 25,
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
resolve_canvas <- function(cfg,sources) {
  main<-sources$circos
  scale<-cfg$main$height_in*72/main$viewbox[4]
  mp<-data.frame(left_pt=cfg$main$left_in*72,top_pt=cfg$main$top_in*72,scale_x=scale,scale_y=scale)
  center<-read_tsv(file.path(cfg$source_dir,paste0(cfg$source_stem,".center_legend.tsv")))[1,]
  g<-read_tsv(file.path(cfg$source_dir,paste0(cfg$source_stem,".figure_geometry.tsv")))[1,]
  mapping<-read_tsv(file.path(cfg$source_dir,paste0(cfg$source_stem,".chromosome_display_mapping.tsv")))
  tracks<-read_tsv(file.path(cfg$source_dir,paste0(cfg$source_stem,".track_geometry.tsv")))
  cx<-mp$left_pt+center$circle_center_x_in/g$width_in*main$viewbox[3]*scale
  cy<-mp$top_pt+(1-center$circle_center_y_in/g$height_in)*main$viewbox[4]*scale
  radius<-g$circle_diameter_in/g$height_in*main$viewbox[4]/2*scale
  max_cos<-vapply(seq_len(nrow(mapping)),function(i) {
    a<-mapping$start_degrees[i];b<-mapping$end_degrees[i]
    if(b<a)stopf("Expected counterclockwise sectors without wrap.")
    candidates<-c(a,b,seq(-720,720,by=360))
    candidates<-candidates[candidates>=a & candidates<=b]
    max(cos(candidates*pi/180))
  },numeric(1))
  sector_right<-cx+radius*max(tracks$outer_radius)*max(max_cos)
  labels<-annotation_boxes(cfg,mp,main)
  pad<-cfg$connectors$annotation_clearance_pt
  label_right<-max(labels$x+abs(labels$tx)*(labels$hw-pad)+abs(labels$rx)*(labels$hh-pad))
  bound<-(max(sector_right,label_right)+cfg$layout$stroke_allowance_pt)/72
  cfg$architectures$left_in<-bound+cfg$layout$module_gap_in
  cfg$page_width_in<-max(cfg$architectures$left_in+cfg$architectures$width_in)+cfg$layout$right_margin_in
  extrema_sin<-vapply(seq_len(nrow(mapping)),function(i) {
    a<-mapping$start_degrees[i];b<-mapping$end_degrees[i]
    candidates<-c(a,b,seq(-630,810,by=180))
    candidates<-candidates[candidates>=a & candidates<=b]
    range(sin(candidates*pi/180))
  },numeric(2))
  outer<-radius*max(tracks$outer_radius)
  label_support_y<-abs(labels$ty)*(labels$hw-pad)+abs(labels$ry)*(labels$hh-pad)
  circos_top<-min(cy-outer*max(extrema_sin),min(labels$y-label_support_y))/72
  circos_bottom<-max(cy-outer*min(extrema_sin),max(labels$y+label_support_y))/72
  module_heights<-cfg$architectures$width_in*vapply(cfg$architectures$track,
    function(t)sources[[t]]$viewbox[4]/sources[[t]]$viewbox[3],numeric(1))
  column_height<-sum(module_heights)+(length(module_heights)-1L)*cfg$layout$module_canvas_gap_in
  if(column_height>(circos_bottom-circos_top)*cfg$layout$max_column_to_circos_height)
    stopf("Enlarged module column does not fit the requested Circos height ratio.")
  column_top<-(circos_top+circos_bottom-column_height)/2
  cfg$architectures$top_in<-column_top+c(0,cumsum(head(module_heights,-1)+cfg$layout$module_canvas_gap_in))
  cfg$resolved_geometry<-data.frame(circle_center_x_in=cx/72,circle_center_y_from_top_in=cy/72,
    circle_radius_in=radius/72,virtual_circle_right_in=(cx+radius)/72,
    sector_envelope_right_in=sector_right/72,annotation_right_in=label_right/72,
    stroke_allowance_pt=cfg$layout$stroke_allowance_pt,circos_right_bound_in=bound,
    module_left_in=cfg$architectures$left_in[1],module_gap_in=cfg$layout$module_gap_in,
    module_width_in=cfg$architectures$width_in[1],page_width_in=cfg$page_width_in,
    page_height_in=cfg$page_height_in,right_margin_in=cfg$layout$right_margin_in,
    circos_top_in=circos_top,circos_bottom_in=circos_bottom,circos_height_in=circos_bottom-circos_top,
    module_canvas_height_in=module_heights[1],module_canvas_gap_in=cfg$layout$module_canvas_gap_in,
    column_top_in=column_top,column_bottom_in=column_top+column_height,column_height_in=column_height,
    column_to_circos_height=column_height/(circos_bottom-circos_top))
  cfg
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
  arch$target_x_in<-arch$left_in+(arch$bar_left_source_x-arch$bar_outline_source_width/2)*arch$scale_x/72-cfg$connectors$target_clearance_in
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
  branches<-list();anchors<-list();polylines<-list()
  bound<-cfg$resolved_geometry$circos_right_bound_in
  for(t in cfg$connectors$tracks) {
    source_track<-if(t=="TEV")"TIP" else t
    tr<-tracks[tracks$track==source_track,]
    if(nrow(tr)!=1L)stopf("Expected one geometry row for %s",source_track)
    f<-cfg$connectors$anchor_track_fraction[[t]]
    r<-tr$inner_radius+f*(tr$outer_radius-tr$inner_radius)
    theta<-upper*pi/180
    true_anchor<-c(cx+radius*r*cos(theta),cy-radius*r*sin(theta))/72
    # For a CCW upper sector edge, this unit normal points INTO the white gap.
    normal<-c(sin(theta),cos(theta))
    start<-true_anchor+normal*cfg$connectors$start_clearance_in
    a<-arch[arch$track==t,]
    first<-c(bound+cfg$connectors$first_elbow_offset_from_bound_in[[t]],start[2]+cfg$connectors$first_leg_drop_in[[t]])
    second<-c(bound+cfg$connectors$elbow_offset_from_bound_in[[t]],a$target_y_in)
    target<-c(a$target_x_in,a$target_y_in)
    if(first[1]<=start[1] || second[1]<=first[1] || target[1]-second[1]<cfg$connectors$arrow_length_in)
      stopf("Invalid or too-short polyline segment for %s",t)
    branches[[t]]<-rbind(start,first,second,target)
    polylines[[t]]<-polyline_connector(branches[[t]])
    anchors[[t]]<-data.frame(track=t,source_track=source_track,source_track_index=tr$track_index,
      source_chromosome="1",source_position_bp=0,source_inner_radius=tr$inner_radius,
      source_outer_radius=tr$outer_radius,anchor_radius=r,anchor_track_fraction=f,
      anchor_degrees=upper,anchor_side="upper",reference_edge_degrees=upper,
      true_x_in=true_anchor[1],true_y_from_top_in=true_anchor[2],
      connector_start_x_in=start[1],connector_start_y_from_top_in=start[2],
      clearance_normal_x=normal[1],clearance_normal_y=normal[2],
      clearance_in=cfg$connectors$start_clearance_in,
      clearance_direction=cfg$connectors$clearance_direction,
      origin_above_circle_center=start[2]*72<cy,exact_true_sector_endpoint=TRUE)
  }
  labels<-annotation_boxes(cfg,mp,main)
  for(t in names(polylines)) {
    p<-polylines[[t]]$samples
    rr<-sqrt((p[,1]-cx)^2+(cy-p[,2])^2)/radius
    aa<-(atan2(cy-p[,2],p[,1]-cx)*180/pi)%%360
    if(any(rr>=0.3185 & rr<=0.9975 & aa>upper+1e-6 & aa<lower-1e-6))stopf("%s polyline crosses main panel.",t)
    for(j in seq_len(nrow(labels))) {
      dx<-p[,1]-labels$x[j];dy<-p[,2]-labels$y[j]
      hit<-abs(dx*labels$tx[j]+dy*labels$ty[j])<labels$hw[j] &
           abs(dx*labels$rx[j]+dy*labels$ry[j])<labels$hh[j]
      if(any(hit))stopf("%s polyline covers label %s",t,labels$label[j])
    }
  }
  # No connector may run across another SVG canvas (own white margin is allowed).
  for(t in names(polylines))for(j in which(arch$track!=t)) {
    p<-polylines[[t]]$samples/72
    if(any(p[,1]>=arch$left_in[j] & p[,1]<=arch$left_in[j]+arch$width_in[j] &
           p[,2]>=arch$top_in[j] & p[,2]<=arch$top_in[j]+arch$height_pt[j]/72))
      stopf("%s connector crosses %s module canvas.",t,arch$track[j])
  }
  pairs<-list()
  for(i in seq_len(length(branches)-1L))for(j in (i+1L):length(branches)) {
    a<-branches[[i]]*72;b<-branches[[j]]*72;dist<-Inf
    for(k in 1:3)for(l in 1:3)dist<-min(dist,segment_distance(a[k,],a[k+1,],b[l,],b[l+1,]))
    if(dist<cfg$connectors$minimum_pair_clearance_pt)stopf("%s/%s polyline clearance %.3f pt is too small.",names(branches)[i],names(branches)[j],dist)
    pairs[[length(pairs)+1L]]<-data.frame(track1=names(branches)[i],track2=names(branches)[j],minimum_clearance_pt=dist,crossing=FALSE,shared_segment=FALSE)
  }
  list(main=mp,architectures=arch,branches=branches,polylines=polylines,anchors=do.call(rbind,anchors),
    connector_pair_QC=do.call(rbind,pairs),circle_center_pt=c(cx,cy),circle_radius_pt=radius,
    upper_edge_deg=upper,lower_edge_deg=lower)
}
# Exactly three straight segments; the exported path contains M and L only.
polyline_connector <- function(points) {
  p<-points*72;pt<-function(x)paste(fmt(x),collapse=" ")
  d<-paste("M",pt(p[1,]),paste(vapply(2:4,function(i)paste("L",pt(p[i,])),character(1)),collapse=" "))
  tt<-seq(0,1,length.out=801)
  samples<-do.call(rbind,lapply(1:3,function(i)outer(1-tt,p[i,])+outer(tt,p[i+1,])))
  list(d=d,samples=samples)
}

draw_connector_layer <- function(cfg,layout,filename) {
  parts<-character()
  for(t in names(layout$branches)) {
    points<-layout$branches[[t]]; route<-layout$polylines[[t]]
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
    a<-layout$anchors[i,];xy<-c(a$true_x_in,a$true_y_from_top_in)*72
    start<-c(a$connector_start_x_in,a$connector_start_y_from_top_in)*72
    color<-cfg$architectures$color[cfg$architectures$track==a$track]
    parts<-c(parts,sprintf('<path class="DEBUG-gap" d="M %s L %s" fill="none" stroke="%s" stroke-width="0.5" stroke-dasharray="1.5 1.5"/>',paste(fmt(xy),collapse=" "),paste(fmt(start),collapse=" "),color),
      sprintf('<circle class="true-anchor" cx="%s" cy="%s" r="2.7" fill="none" stroke="%s" stroke-width="1"/>',fmt(xy[1]),fmt(xy[2]),color),
      sprintf('<circle class="connector-start" cx="%s" cy="%s" r="3.6" fill="%s" stroke="#FFFFFF" stroke-width="0.7"/>',fmt(start[1]),fmt(start[2]),color),
      sprintf('<text x="%s" y="%s" text-anchor="middle" dominant-baseline="central" fill="#FFFFFF" font-family="Arial, sans-serif" font-size="5.3" font-weight="700">%s</text>',fmt(start[1]),fmt(start[2]),i))
  }
  parts<-c(parts,'<g font-family="Arial, sans-serif" font-size="8" fill="#273746">',
    '<text x="680" y="20" font-weight="700">DEBUG: hollow = true endpoint; solid = start after 0.14in normal gap</text>')
  for(i in seq_len(nrow(layout$anchors))) {
    a<-layout$anchors[i,]
    parts<-c(parts,sprintf('<text x="680" y="%s">%s  %s: source track %s (%s), upper Chr1/0bp, 24°</text>',23+i*11,i,a$track,a$source_track_index,a$source_track))
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

# This rsvg version quantizes PDF dimensions to whole points. Pad the page up
# by less than 1pt, keeping the original SVG point coordinates at 1:1 scale.
write_pdf <- function(svg_path,pdf_path,width_in,height_in) {
  wp<-ceiling(width_in*72);hp<-ceiling(height_in*72)
  svg<-paste(readLines(svg_path,warn=FALSE),collapse="\n")
  svg<-sub('width="[^"]+"',paste0('width="',wp,'pt"'),svg)
  svg<-sub('height="[^"]+"',paste0('height="',hp,'pt"'),svg)
  svg<-sub('viewBox="[^"]+"',paste0('viewBox="0 0 ',wp,' ',hp,'"'),svg)
  tmp<-tempfile("fig1_r3_pdf_",fileext=".svg");on.exit(unlink(tmp),add=TRUE)
  writeLines(svg,tmp,useBytes=TRUE)
  rsvg::rsvg_pdf(tmp,pdf_path,width=wp,height=hp)
  invisible(data.frame(nominal_width_pt=width_in*72,nominal_height_pt=height_in*72,
    page_width_pt=wp,page_height_pt=hp,right_padding_pt=wp-width_in*72,
    bottom_padding_pt=hp-height_in*72,content_scale=1))
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
      cfg$output_dir %in% file.path(dirname(protected_dir),c("R_v15_with_TE_architectures","R_v15_with_TE_architectures_v2","R_v15_with_TE_architectures_v3","R_v15_with_TE_architectures_v3_r2")))
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
    "compose_v15_with_te_architectures_v3_r3.R", "run_compose_v15_with_te_architectures_v3_r3.sh",
    "install_composition_R_dependencies.sh")))
  provenance_files <- file.path(cfg$asset_dir,c("TE_architectures/README.md","TE_architectures/README_composition_1_5x.md","TE_architectures/source_manifest.json",
    "TE_architectures/TE_architectures_4up.svg","TE_mediated_SV/README.md",
    "TE_mediated_SV/TEV_classification_composition_1_5x.svg","TE_mediated_SV/TE_mediated_SV_classification_V15_gt50.provenance.json"))
  protected <- c(protected,provenance_files[file.exists(provenance_files)])
  protected <- unique(protected); before <- tools::md5sum(protected)
  cfg <- resolve_canvas(cfg,sources)
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
    'five independent three-segment M/L arrows begin 0.14in along the endpoint outward normal into the upper opening from the true LINE, LTR, SINE, DNA and innermost TEV anchors. TEV schematic shows the user-requested >80% rule; V15 data retain >50%. Not sequence to scale.</desc>\n',
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
    if ("pdf" %in% cfg$formats) {
      pdf_geometry<-write_pdf(svg_path,paste0(prefix,".pdf"),cfg$page_width_in,cfg$page_height_in)
      write_tsv(pdf_geometry,file.path(cfg$output_dir,"PDF_geometry.tsv"))
    }
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
  write_tsv(cfg$resolved_geometry,file.path(cfg$output_dir,"layout_geometry.tsv"))
  write_tsv(layout$architectures,file.path(cfg$output_dir,"architecture_placements.tsv"))
  paths <- layout$branches
  connector_rows <- do.call(rbind,lapply(names(paths),function(n) data.frame(
    connector=n,vertex=seq_len(nrow(paths[[n]])),x_in=paths[[n]][,1],y_from_top_in=paths[[n]][,2],
    arrowhead=n!="stem" & seq_len(nrow(paths[[n]]))==nrow(paths[[n]]))))
  write_tsv(connector_rows,file.path(cfg$output_dir,"connector_vertices.tsv"))
  write_tsv(data.frame(connector=names(layout$polylines),svg_path=vapply(layout$polylines,`[[`,character(1),"d")),
    file.path(cfg$output_dir,"connector_paths.tsv"))
  samples<-do.call(rbind,lapply(names(layout$polylines),function(t) {
    m<-layout$polylines[[t]]$samples
    data.frame(connector=t,sample=seq_len(nrow(m)),x_pt=m[,1],y_from_top_pt=m[,2])
  }))
  write_tsv(samples,file.path(cfg$output_dir,"connector_line_samples.tsv"))
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
      cat("Usage: Rscript 06_scripts/compose_v15_with_te_architectures_v3_r3.R [--project-dir PATH] [--out-dir PATH] [--preview-only] [--dpi 600]\n")
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
