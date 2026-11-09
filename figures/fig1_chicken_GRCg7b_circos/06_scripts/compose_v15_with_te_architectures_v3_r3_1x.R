#!/usr/bin/env Rscript
# Figure1 R3 1x: original-size SVG modules in a compact, centered column.
# Original display baseline: 4.9375in = 7.40625in / 1.5.
# Keep page, main Circos, true anchors, connector starts and horizontal placement.
# Only resize/reposition modules and retarget the final two connector vertices.
options(stringsAsFactors=FALSE,scipen=20)
script_project_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) return(dirname(dirname(normalizePath(sub("^--file=", "", arg[1])))))
  if (dir.exists("04_window_statistics")) return(normalizePath("."))
  normalizePath("chicken_7b_TE_circos", mustWork = TRUE)
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


replace_asset_block <- function(text,id,replacement) {
  pattern<-paste0('(?s)<svg\\b[^>]*\\bid="',id,'"[^>]*>.*?</svg>')
  # The five supplied assets contain no nested SVG root of their own.
  hits<-gregexpr(pattern,text,perl=TRUE)[[1]]
  if(length(hits)!=1L || hits[1]<0)stopf("Expected one block for %s",id)
  at<-hits[1];len<-attr(hits,"match.length")
  paste0(substr(text,1,at-1L),replacement,substr(text,at+len,nchar(text)))
}


# Update only the geometry attribute of one native SVG element. All styles stay.
replace_element_attribute <- function(text, tag, id, attribute, value) {
  pattern<-paste0('<',tag,'\\b[^>]*\\bid="',id,'"[^>]*>')
  hit<-regexpr(pattern,text,perl=TRUE)
  if(hit[1]<0)stopf("Missing %s: %s",tag,id)
  element<-regmatches(text,hit)
  updated<-sub(paste0('\\b',attribute,'="[^"]*"'),paste0(attribute,'="',value,'"'),element,perl=TRUE)
  if(identical(updated,element))stopf("Failed to update %s for %s",attribute,id)
  regmatches(text,hit)<-updated
  text
}
retarget_connectors <- function(text, vertices, old_lanes, new_lanes, id_prefix="connector__") {
  for(track in names(new_lanes)) {
    v<-vertices[vertices$connector==track,]
    d<-paste(paste(c("M","L","L","L"),fmt(v$x_in*72),fmt(v$y_from_top_in*72)),collapse=" ")
    text<-replace_element_attribute(text,"path",paste0(id_prefix,"connector-",track),"d",d)
    id<-paste0(id_prefix,"arrow-",track)
    hit<-regexpr(paste0('<polygon\\b[^>]*\\bid="',id,'"[^>]*>'),text,perl=TRUE)
    element<-regmatches(text,hit)
    points<-svg_attribute(element,"points")
    xy<-matrix(as.numeric(strsplit(points,"[ ,]+")[[1]]),ncol=2,byrow=TRUE)
    xy[,2]<-xy[,2]+(new_lanes[track]-old_lanes[track])*72
    value<-paste(paste(fmt(xy[,1]),fmt(xy[,2]),sep=","),collapse=" ")
    text<-replace_element_attribute(text,"polygon",id,"points",value)
  }
  text
}
segment_crosses <- function(a,b,c,d) {
  cross<-function(x,y,z)(y[1]-x[1])*(z[2]-x[2])-(y[2]-x[2])*(z[1]-x[1])
  boxes<-max(min(a[1],b[1]),min(c[1],d[1]))<=min(max(a[1],b[1]),max(c[1],d[1])) &&
    max(min(a[2],b[2]),min(c[2],d[2]))<=min(max(a[2],b[2]),max(c[2],d[2]))
  boxes && cross(a,b,c)*cross(a,b,d)<=0 && cross(c,d,a)*cross(c,d,b)<=0
}
segment_distance <- function(a,b,c,d) {
  if(segment_crosses(a,b,c,d))return(0)
  pd<-function(p,a,b){v<-b-a;t<-max(0,min(1,sum((p-a)*v)/sum(v*v)));sqrt(sum((p-a-t*v)^2))}
  min(pd(a,c,d),pd(b,c,d),pd(c,a,b),pd(d,a,b))
}

compose_original_size <- function(project_dir=script_project_dir(),out_dir=NULL,preview_only=FALSE) {
  project_dir<-normalizePath(project_dir,mustWork=TRUE)
  lib<-file.path(project_dir,"06_scripts","composition_R_library")
  if(dir.exists(lib)).libPaths(c(normalizePath(lib),.libPaths()))
  for(pkg in c("rsvg","png"))if(!requireNamespace(pkg,quietly=TRUE))
    stopf("Missing %s; run 06_scripts/install_composition_R_dependencies.sh",pkg)
  frozen<-file.path(project_dir,"07_figures","R_v15_with_TE_architectures_v3_r3")
  if(is.null(out_dir))out_dir<-file.path(project_dir,"07_figures","R_v15_with_TE_architectures_v3_r3_1x")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  out_dir<-normalizePath(out_dir,mustWork=TRUE)
  if(out_dir==normalizePath(frozen) || startsWith(out_dir,paste0(normalizePath(frozen),"/")))
    stopf("The frozen R3 input directory cannot be an output directory")
  parent<-basename(dirname(out_dir));base<-basename(out_dir)
  if(parent=="07_figures" && base %in% c("R_v15","R_v15_with_TE_architectures","R_v15_with_TE_architectures_v2","R_v15_with_TE_architectures_v3","R_v15_with_TE_architectures_v3_r2"))
    stopf("Do not overwrite any older figure version")
  old_stem<-"chicken_GRCg7b_R_v15_with_TE_architectures_v3_r3"
  stem<-paste0(old_stem,"_1x");prefix<-file.path(out_dir,stem)
  a<-read_tsv(file.path(frozen,"architecture_placements.tsv"))
  a$source_R3_file<-a$file;a$source_R3_width_in<-a$width_in
  a$file<-c("TE_architectures/LINE_CR1.svg","TE_architectures/LTR_ERV.svg",
    "TE_architectures/SINE_generic.svg","TE_architectures/DNA_TIR.svg",
    "TE_mediated_SV/TE_mediated_SV_classification_v3_r2.svg")
  protected<-list.files(frozen,full.names=TRUE,recursive=TRUE)
  protected<-protected[file.info(protected)$isdir==FALSE]
  protected<-unique(c(protected,file.path(project_dir,"assets",a$file),
    file.path(project_dir,"06_scripts",c("compose_v15_with_te_architectures_v3_r3_1x.R","run_compose_v15_with_te_architectures_v3_r3_1x.sh"))))
  before<-tools::md5sum(protected)
  sources<-lapply(file.path(project_dir,"assets",a$file),read_svg)
  # Preserve visible left edges; center a tighter original-size module column.
  g<-read_tsv(file.path(frozen,"layout_geometry.tsv"))[1,]
  old_lanes<-setNames(a$target_y_in,a$track)
  edge_pt<-a$left_pt+(a$bar_left_source_x-a$bar_outline_source_width/2)*a$scale_x
  a$width_in<-4.9375;a$width_pt<-a$width_in*72
  a$scale_x<-a$width_pt/vapply(sources,function(z)z$viewbox[3],numeric(1))
  a$scale_y<-a$scale_x
  a$height_pt<-a$scale_x*vapply(sources,function(z)z$viewbox[4],numeric(1))
  module_gap_in<-0.08
  column_height_in<-sum(a$height_pt/72)+(nrow(a)-1)*module_gap_in
  column_top_in<-mean(c(g$circos_top_in,g$circos_bottom_in))-column_height_in/2
  a$top_in<-column_top_in+c(0,cumsum(head(a$height_pt/72,-1)+module_gap_in))
  center_pt<-a$top_in*72+a$bar_center_source_y*a$scale_y
  a$target_y_in<-center_pt/72
  new_lanes<-setNames(a$target_y_in,a$track)
  vertices<-read_tsv(file.path(frozen,"connector_vertices.tsv"))
  vertices$y_from_top_in[vertices$vertex>=3]<-new_lanes[vertices$connector[vertices$vertex>=3]]
  pairs<-list()
  for(i in 1:4)for(j in (i+1):5) {
    p<-as.matrix(vertices[vertices$connector==a$track[i],c("x_in","y_from_top_in")])*72
    q<-as.matrix(vertices[vertices$connector==a$track[j],c("x_in","y_from_top_in")])*72
    distance<-min(vapply(1:3,function(k)min(vapply(1:3,function(l)segment_distance(p[k,],p[k+1,],q[l,],q[l+1,]),numeric(1))),numeric(1)))
    if(distance<3.5)stopf("Connector clearance failed: %s / %s: %.3f pt",a$track[i],a$track[j],distance)
    pairs[[length(pairs)+1]]<-data.frame(track1=a$track[i],track2=a$track[j],minimum_clearance_pt=distance,crossing=FALSE,shared_segment=FALSE)
  }
  a$left_pt<-edge_pt-(a$bar_left_source_x-a$bar_outline_source_width/2)*a$scale_x
  a$top_pt<-center_pt-a$bar_center_source_y*a$scale_y
  a$left_in<-a$left_pt/72;a$top_in<-a$top_pt/72
  a$connection_edge_x_pt<-edge_pt;a$connection_center_y_pt<-center_pt
  replacements<-vapply(seq_len(nrow(a)),function(i)embed_svg(sources[[i]],paste0("architecture-",a$track[i]),a$left_pt[i],a$top_pt[i],a$width_pt[i],a$height_pt[i]),character(1))
  changed<-list()
  for(suffix in c(".svg",".DEBUG.svg",".DEBUG.anchor_zoom.svg")) {
    input<-file.path(frozen,paste0(old_stem,suffix))
    text<-paste(readLines(assert_file(input),warn=FALSE,encoding="UTF-8"),collapse="\n")
    for(i in seq_len(nrow(a)))text<-replace_asset_block(text,paste0("architecture-",a$track[i]),replacements[i])
    text<-retarget_connectors(text,vertices,old_lanes,new_lanes)
    output<-paste0(prefix,suffix);writeLines(text,output,useBytes=TRUE);changed[[suffix]]<-output
  }
  W<-g$page_width_in;H<-g$page_height_in
  write_png(changed[[".svg"]],paste0(prefix,".preview.png"),W,H,150)
  write_png(changed[[".DEBUG.svg"]],paste0(prefix,".DEBUG.preview.png"),W,H,150)
  write_png(changed[[".DEBUG.anchor_zoom.svg"]],paste0(prefix,".DEBUG.anchor_zoom.png"),7.1,4.6,300)
  if(!preview_only) {
    write_png(changed[[".svg"]],paste0(prefix,".png"),W,H,600)
    pdf_geometry<-write_pdf(changed[[".svg"]],paste0(prefix,".pdf"),W,H)
    write_tsv(pdf_geometry,file.path(out_dir,"PDF_geometry.tsv"))
  }
  write_tsv(a,file.path(out_dir,"architecture_placements.tsv"))
  write_tsv(vertices,file.path(out_dir,"connector_vertices.tsv"))
  write_tsv(do.call(rbind,pairs),file.path(out_dir,"connector_pair_QC.tsv"))
  paths<-do.call(rbind,lapply(a$track,function(track) {
    v<-vertices[vertices$connector==track,]
    data.frame(connector=track,svg_path=paste(paste(c("M","L","L","L"),fmt(v$x_in*72),fmt(v$y_from_top_in*72)),collapse=" "))
  }))
  write_tsv(paths,file.path(out_dir,"connector_paths.tsv"))
  only<-paste(readLines(file.path(frozen,"connectors_only.svg"),warn=FALSE),collapse="\n")
  writeLines(retarget_connectors(only,vertices,old_lanes,new_lanes,id_prefix=""),file.path(out_dir,"connectors_only.svg"),useBytes=TRUE)
  for(file in c("connector_anchors.tsv","circos_placement.tsv")) {
    if(!file.copy(file.path(frozen,file),file.path(out_dir,file),overwrite=TRUE))stopf("Copy failed: %s",file)
  }
  file.copy(file.path(frozen,"layout_geometry.tsv"),file.path(out_dir,"frozen_R3_layout_geometry.tsv"),overwrite=TRUE)
  after<-tools::md5sum(protected)
  if(!identical(before,after))stopf("A frozen R3/input file changed")
  write_tsv(data.frame(path=names(before),md5=unname(before),unchanged=unname(before==after)),file.path(out_dir,"input_preservation.tsv"))
  write_tsv(data.frame(setting=c("original_display_width_in","source_R3_width_in","relative_scale","page_width_in","page_height_in","main_and_connector_starts_frozen","module_canvas_gap_in","column_height_in","column_top_in","column_center_in"),
    value=c("4.9375","7.40625","0.6666666666666667",fmt(W),fmt(H),"TRUE",fmt(module_gap_in),fmt(column_height_in),fmt(column_top_in),fmt(column_top_in+column_height_in/2))),file.path(out_dir,"resize_settings.tsv"))
  writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))
  message("Done: ",out_dir,"; compact original-size modules and retargeted line ends; frozen inputs unchanged: ",length(protected))
}

if(sys.nframe()==0L) {
  args<-commandArgs(TRUE);project_dir<-script_project_dir();out_dir<-NULL;preview<-FALSE;i<-1L
  while(i<=length(args)) {
    if(args[i]=="--preview-only") {preview<-TRUE;i<-i+1L;next}
    if(args[i]=="--help") {cat("Usage: Rscript compose_v15_with_te_architectures_v3_r3_1x.R [--project-dir PATH] [--out-dir PATH] [--preview-only]\n");quit(status=0)}
    if(i==length(args))stopf("Missing value for %s",args[i])
    if(args[i]=="--project-dir")project_dir<-args[i+1L]
    else if(args[i]=="--out-dir")out_dir<-args[i+1L]
    else stopf("Unknown argument: %s",args[i])
    i<-i+2L
  }
  compose_original_size(project_dir,out_dir,preview)
}
