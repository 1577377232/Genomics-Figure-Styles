#!/usr/bin/env Rscript
# GWAS–molQTL: original associations, population-specific dosage LD.
# Run: Rscript plot_coloc_template.R [template_directory]
# Design references (no dependency; no copied code or human reference data):
# https://github.com/jrs95/gassocplot2 (stacked regional association)
# https://github.com/boxiangliu/locuscomparer (paired association visualisation)
# https://github.com/myles-lewis/locuszoomr (LD colour, layering, index symbol)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})

# ====================== EDIT STYLE / INPUT SETTINGS HERE ======================
WIDTH_MM <- 105
HEIGHT_MM <- 100
SUPPLEMENT_WIDTH_MM <- 75
SUPPLEMENT_HEIGHT_MM <- 45
FONT_FAMILY <- 'Liberation Sans'
BASE_SIZE <- 7
LABEL_SIZE_PT <- 6.2
LABEL_WRAP_WIDTH <- 28 # wrap long verified trait/gene labels without changing their identity
INK <- '#253246'
LD_COLOURS <- c('#386C8C','#6681A2','#9689B9','#C6A8A0','#E5C65F')
LD_MISSING_COLOUR <- '#D9DDE2'
FOCAL_COLOUR <- '#E75A46'
POINT_SIZES <- c(SNP=1.15,SV=1.40,TEV=1.65)
POINT_SHAPES <- c(SNP=21,SV=22,TEV=23)
POINT_STROKE <- .07
POINT_ALPHA <- .84
FOCAL_SIZE <- 3.0
FOCAL_STROKE <- .38
HALF_WINDOW_BP <- 1000000
SHOW_FOCAL_LINE_ANNOTATED <- TRUE
SHOW_FOCAL_LINE_CLEAN <- FALSE
SHOW_CLEAN_FOCAL_LABEL <- FALSE
FOCAL_LINE_ALPHA <- .25
FOCAL_LINE_WIDTH <- .25
FOCAL_LINE_COLOUR <- '#D7B661'
TE_COVERAGE_CUTOFF <- .50
TRACK_HEIGHTS <- c(1,1)
HEADROOM_ANNOTATED <- 1.40 # space for true-statistic labels
HEADROOM_CLEAN <- 1.12
TRACK_MARGIN_PT <- .6
LD_ANCHORS <- c(0,.25,.50,.75,1)
TYPE_BASE_COLOURS <- c(SNP='#386C8C',SV='#9689B9',TEV='#E75A46') # supplemental/legend use; not main fills
TEXT_LEFT_X_FRAC <- .025
TEXT_RIGHT_X_FRAC <- .985
TEXT_Y_FRAC <- .98
GWAS_THRESHOLD_P <- NULL # only enter a real, documented analysis threshold
MOL_THRESHOLD_P <- NULL
GWAS_THRESHOLD_SOURCE <- NULL
MOL_THRESHOLD_SOURCE <- NULL
FOCAL_LABEL <- 'Focal SV'
TRAIT_LABEL <- 'GWAS trait'
GWAS_POP_LABEL <- 'GWAS population'
GENE_LABEL <- 'Molecular phenotype'
OUTPUT_PREFIX <- 'colocalization' # retain original output prefix and directory
GWAS_INPUT <- 'data/GWAS.tsv.gz'
MOL_INPUT <- 'data/molQTL.tsv.gz'
GWAS_LD_INPUT <- 'data/GWAS_dosage_LD.tsv.gz'
MOL_LD_INPUT <- 'data/C367_dosage_LD.tsv.gz'
PIP_INPUT <- 'data/focal_PIP.json' # optional; absent PIP is not replaced by zero
# =============================================================================
args <- commandArgs(trailingOnly=TRUE)
self <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
root <- if(length(args)) normalizePath(args[1]) else dirname(normalizePath(self))
# Optional per-case JSON overrides: style and labels only, never statistical results.
config_file <- file.path(root,'case_config.json')
if(file.exists(config_file)) {
  config <- fromJSON(config_file)
  allowed <- c('FOCAL_LABEL','TRAIT_LABEL','GWAS_POP_LABEL','GENE_LABEL','OUTPUT_PREFIX','WIDTH_MM','HEIGHT_MM','HALF_WINDOW_BP','SHOW_FOCAL_LINE_ANNOTATED','SHOW_FOCAL_LINE_CLEAN')
  if(length(setdiff(names(config),allowed)))stop('Unknown configuration field')
  for(nm in names(config))assign(nm,config[[nm]])
}
out <- file.path(root,'figures');dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_tsv <- function(rel) read.delim(gzfile(file.path(root,rel)),check.names=FALSE,stringsAsFactors=FALSE)
row <- read.delim(file.path(root,'data/配对信息.tsv'),check.names=FALSE)
meta <- fromJSON(file.path(root,'data/case.json'))
stopifnot(nrow(row)==1)
focal_id <- meta$focal_variant
coord <- strsplit(row$coordinate,':',fixed=TRUE)[[1]]
chr <- coord[1];focal_pos <- as.numeric(coord[2]);limits <- pmax(1,focal_pos+c(-HALF_WINDOW_BP,HALF_WINDOW_BP))/1e6
coverage <- read_tsv('data/SV_TE_coverage.tsv.gz')
stopifnot(!anyDuplicated(coverage$marker))
fmt_p <- function(p) if(length(p)!=1 || !is.finite(p)) 'unavailable' else sub('e([+-])0+','e\\1',sprintf('%.2e',p))
fmt_pp <- function(p) if(length(p)!=1 || !is.finite(p)) 'unavailable' else sprintf('%.3f',p)
pip <- if(file.exists(file.path(root,PIP_INPUT))) fromJSON(file.path(root,PIP_INPUT)) else NULL
if(!is.null(pip)) stopifnot(pip$variant_id==focal_id,pip$PIP>=0,pip$PIP<=1)

prepare <- function(input,ld_file,pos_col,p_col,type_col,track) {
  d <- read_tsv(input);ld <- read_tsv(ld_file)
  stopifnot(!anyDuplicated(d$variant_id),!anyDuplicated(ld$variant_id))
  idx <- match(d$variant_id,ld$variant_id) # full ID, never POS-only
  matched <- !is.na(idx)
  if(any(d$variant_key[matched]!=ld$variant_key[idx[matched]])) stop('Complete variant key mismatch in ',track)
  d$position_bp <- d[[pos_col]];d$association_P <- d[[p_col]]
  d$LD_r2 <- ld$LD_r2[idx];d$LD_N_complete <- ld$LD_N_complete[idx]
  d$LD_status <- ld$LD_status[idx];d$LD_status[!matched] <- 'unmatched_full_variant_ID'
  stopifnot(all(is.na(d$LD_r2)|(d$LD_r2>=0&d$LD_r2<=1)))
  d$TE_coverage <- coverage$TE_coverage[match(d$variant_id,coverage$marker)]
  d$plot_class <- ifelse(d[[type_col]]=='SNP','SNP',ifelse(!is.na(d$TE_coverage)&d$TE_coverage>=TE_COVERAGE_CUTOFF,'TEV','SV'))
  d$plot_class <- factor(d$plot_class,levels=c('SNP','SV','TEV'))
  d$within_window <- is.finite(d$position_bp)&d$position_bp>=limits[1]*1e6&d$position_bp<=limits[2]*1e6
  d$valid_P <- is.finite(d$association_P)&d$association_P>0&d$association_P<=1
  if(any(d$within_window&!d$valid_P)) warning(track,': invalid P rows exported, not assigned artificial values.')
  write.table(d[,c('variant_id','variant_key','position_bp','association_P','plot_class','TE_coverage','LD_r2','LD_N_complete','LD_status','within_window','valid_P')],file.path(out,paste0(track,'_LD_and_record_audit.tsv')),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
  plot <- d[d$within_window&d$valid_P,,drop=FALSE]
  # All original records retained, independent of coloc overlap or LD availability.
  plot <- plot[order(as.integer(plot$plot_class),!is.na(plot$LD_r2),plot$LD_r2,na.last=FALSE),,drop=FALSE]
  plot$x <- plot$position_bp/1e6;plot$y <- -log10(plot$association_P)
  stopifnot(sum(plot$variant_id==focal_id)<=1)
  # Explicit layer audit: focal is removed from its base class and drawn once, last.
  draw <- plot[plot$variant_id!=focal_id,,drop=FALSE]
  draw$draw_layer <- as.character(draw$plot_class)
  ff <- plot[plot$variant_id==focal_id,,drop=FALSE];ff$draw_layer <- rep('focal',nrow(ff))
  draw <- rbind(draw,ff);draw$draw_order <- seq_len(nrow(draw))
  stopifnot(nrow(draw)==nrow(plot),!anyDuplicated(draw$variant_id))
  if(nrow(ff))stopifnot(tail(draw$variant_id,1)==focal_id)
  for(cl in c('SNP','SV','TEV')) {
    ld_values <- draw$LD_r2[draw$draw_layer==cl];ld_values <- ld_values[is.finite(ld_values)]
    stopifnot(!is.unsorted(ld_values))
  }
  write.table(draw[,c('variant_id','plot_class','draw_layer','draw_order','LD_r2')],file.path(out,paste0(track,'_drawing_order_audit.tsv')),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
  list(data=plot,total=nrow(d),in_window=sum(d$within_window),invalid=sum(d$within_window&!d$valid_P),track=track)
}
g <- prepare(GWAS_INPUT,GWAS_LD_INPUT,'position','P','type','GWAS')
m <- prepare(MOL_INPUT,MOL_LD_INPUT,'position_mol','pval_g1','type_mol','molQTL')
fg <- g$data[g$data$variant_id==focal_id,,drop=FALSE]
fm <- m$data[m$data$variant_id==focal_id,,drop=FALSE]
same_number <- function(x,y) {
  if(length(x)==0)return(is.na(y))
  (is.na(x)&&is.na(y)) || (is.finite(x)&&is.finite(y)&&abs(x-y)<1e-12)
}
stopifnot(same_number(fg$association_P,row$focal_GWAS_P),same_number(fm$association_P,row$focal_mol_P),
          abs(meta$PPH4-row$PPH4)<1e-12,
          same_number(row$PPH4*row$focal_conditional_shared_PP,row$focal_joint_shared_PP))
if(!is.null(pip))stopifnot(same_number(pip$original_GWAS_P,row$focal_GWAS_P))

theme_base <- theme_classic(base_size=BASE_SIZE,base_family=FONT_FAMILY) +
  theme(text=element_text(colour=INK),axis.text=element_text(colour=INK),
        axis.title=element_text(size=BASE_SIZE),axis.line=element_line(colour=INK,linewidth=.30),
        axis.ticks=element_line(colour=INK,linewidth=.25),axis.ticks.length=grid::unit(1,'mm'),
        legend.position='right',legend.box='vertical',legend.box.just='left',
        legend.title=element_text(size=BASE_SIZE),legend.text=element_text(size=BASE_SIZE-1),
        legend.key.height=grid::unit(2.7,'mm'),legend.key.width=grid::unit(2.7,'mm'),
        legend.spacing.y=grid::unit(1.6,'mm'),legend.margin=margin(0,0,0,1),
        plot.margin=margin(TRACK_MARGIN_PT,1,TRACK_MARGIN_PT,1),panel.grid=element_blank())

panel <- function(case,left_label,right_label,bottom=FALSE,threshold=NULL,threshold_source=NULL,version='annotated') {
  d <- case$data;f <- d[d$variant_id==focal_id,,drop=FALSE]
  annotated <- identical(version,'annotated')
  ytop <- max(d$y)*if(annotated) HEADROOM_ANNOTATED else HEADROOM_CLEAN
  p <- ggplot(d,aes(x,y))
  # Separate geom layers enforce SNP < SV < TEV < focal, independent of LD.
  for(cl in c('SNP','SV','TEV')) {
    layer_data <- d[d$plot_class==cl & d$variant_id!=focal_id,,drop=FALSE]
    p <- p + geom_point(data=layer_data,aes(fill=LD_r2,shape=plot_class,size=plot_class),colour=INK,stroke=POINT_STROKE,alpha=POINT_ALPHA)
  }
  p <- p +
    scale_fill_gradientn(colours=LD_COLOURS,values=LD_ANCHORS,limits=c(0,1),breaks=c(0,.5,1),na.value=LD_MISSING_COLOUR,name=expression(LD~(r^2)),
      guide=guide_colourbar(order=2,title.position='top',barheight=grid::unit(22,'mm'),barwidth=grid::unit(2.4,'mm'),ticks=TRUE)) +
    scale_shape_manual(values=POINT_SHAPES,drop=FALSE,name='Variant Type',
      guide=guide_legend(order=1,override.aes=list(fill='#D5DBE2',colour=INK,size=1.7,alpha=1,stroke=.20))) +
    scale_size_manual(values=POINT_SIZES,guide='none') +
    scale_x_continuous(limits=limits,expand=expansion(mult=.015),breaks=scales::breaks_pretty(4)) +
    scale_y_continuous(limits=c(0,ytop),expand=expansion(mult=c(.025,0)),breaks=scales::breaks_pretty(4)) +
    labs(x=paste0('Chromosome ',chr,' position (Mb)'),y=expression(-log[10](P))) + theme_base
  if(annotated) {
    p <- p + annotate('text',x=limits[1]+diff(limits)*TEXT_LEFT_X_FRAC,y=ytop*TEXT_Y_FRAC,label=left_label,
      hjust=0,vjust=1,size=LABEL_SIZE_PT/ggplot2::.pt,family=FONT_FAMILY,colour=INK,lineheight=1.15) +
      annotate('text',x=limits[1]+diff(limits)*TEXT_RIGHT_X_FRAC,y=ytop*TEXT_Y_FRAC,label=right_label,
      hjust=1,vjust=1,size=LABEL_SIZE_PT/ggplot2::.pt,family=FONT_FAMILY,colour=INK,lineheight=1.15)
  } else {
    # Track identity is retained only as an axis label; all statistic/sample text is omitted.
    p <- p + labs(y=bquote(atop(.(if(bottom) 'molQTL' else 'GWAS'),-log[10](P))))
    if(SHOW_CLEAN_FOCAL_LABEL) p <- p + annotate('text',x=f$x,y=ytop*.98,label=FOCAL_LABEL,vjust=1,
      family=FONT_FAMILY,colour=INK,size=LABEL_SIZE_PT/ggplot2::.pt)
  }
  if(if(annotated) SHOW_FOCAL_LINE_ANNOTATED else SHOW_FOCAL_LINE_CLEAN) p <- p + geom_vline(xintercept=focal_pos/1e6,colour=FOCAL_LINE_COLOUR,alpha=FOCAL_LINE_ALPHA,linewidth=FOCAL_LINE_WIDTH)
  if(annotated && !is.null(threshold)) {
    if(is.null(threshold_source)||!nzchar(threshold_source)||threshold<=0||threshold>=1)stop('A real threshold and source are required.')
    p <- p + geom_hline(yintercept=-log10(threshold),linetype='dashed',colour='#A0A8B0',linewidth=.3)
  }
  # Last layer, independent of the colourbar; same focal biological allele in both populations.
  p <- p + geom_point(data=f,shape=23,fill=FOCAL_COLOUR,colour=INK,size=FOCAL_SIZE,stroke=FOCAL_STROKE,show.legend=FALSE)
  if(!bottom)p <- p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  p
}
pip_label <- if(is.null(pip)) 'PIP unavailable' else sprintf('PIP (+/-%g kb) = %s',pip$half_window_bp/1000,fmt_pp(pip$PIP))
g_label <- paste0(FOCAL_LABEL,'\nP = ',fmt_p(fg$association_P),'\n',pip_label)
m_label <- paste0('P = ',fmt_p(fm$association_P),'\nPPH4 = ',fmt_pp(row$PPH4),'\nSV joint PP = ',fmt_pp(row$focal_joint_shared_PP))
wrap_label <- function(x) paste(strwrap(x,width=LABEL_WRAP_WIDTH),collapse='\n')
make_main <- function(version) {
  p1 <- panel(g,paste0(wrap_label(TRAIT_LABEL),'\nGWAS | ',row$GWAS_model,' | ',GWAS_POP_LABEL,' (n=',row$N_GWAS,')'),g_label,FALSE,GWAS_THRESHOLD_P,GWAS_THRESHOLD_SOURCE,version)
  p2 <- panel(m,paste0(wrap_label(GENE_LABEL),'\n',row$scope,' | ',row$molecular_model,' (n=',row$N_mol,')'),m_label,TRUE,MOL_THRESHOLD_P,MOL_THRESHOLD_SOURCE,version)
  ((p1/p2)+plot_layout(heights=TRACK_HEIGHTS,guides='collect')) & theme(legend.position='right')
}
main_annotated <- make_main('annotated')
main_clean <- make_main('clean')

save_formats <- function(plot,stem,w,h) {
  ggsave(file.path(out,paste0(stem,'.pdf')),plot,width=w,height=h,units='mm',device=cairo_pdf,bg='white')
  ggsave(file.path(out,paste0(stem,'.svg')),plot,width=w,height=h,units='mm',device=svglite::svglite,bg='white')
  ggsave(file.path(out,paste0(stem,'.png')),plot,width=w,height=h,units='mm',dpi=600,bg='white')
  ggsave(file.path(out,paste0(stem,'_preview.png')),plot,width=w,height=h,units='mm',dpi=160,bg='white')
}
save_formats(main_annotated,paste0(OUTPUT_PREFIX,'_annotated'),WIDTH_MM,HEIGHT_MM)
save_formats(main_clean,paste0(OUTPUT_PREFIX,'_clean'),WIDTH_MM,HEIGHT_MM)
post <- as.numeric(unlist(meta[paste0('PPH',0:4)]));stopifnot(length(post)==5,abs(sum(post)-1)<1e-6)
pp <- data.frame(H=factor(paste0('H',0:4),levels=paste0('H',0:4)),posterior=post)
supplement <- ggplot(pp,aes(H,posterior,fill=H))+geom_col(width=.62)+
  scale_fill_manual(values=c('#D9DDE2','#6681A2','#9689B9','#E5C65F',FOCAL_COLOUR),guide='none')+
  scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=c(0,.03)))+
  labs(x=NULL,y='Posterior probability')+theme_base+theme(legend.position='none')
save_formats(supplement,paste0(OUTPUT_PREFIX,'_H0_H4_supplement'),SUPPLEMENT_WIDTH_MM,SUPPLEMENT_HEIGHT_MM)

count <- function(case) {
  d<-case$data;valid<-is.finite(d$LD_r2)
  data.frame(track=case$track,input_records=case$total,within_window=case$in_window,invalid_P=case$invalid,
    plotted_records=nrow(d),SNP=sum(d$plot_class=='SNP'),SV=sum(d$plot_class=='SV'),TEV=sum(d$plot_class=='TEV'),
    LD_estimable=sum(valid),LD_missing=sum(!valid),LD_N_min=if(any(valid)) min(d$LD_N_complete[valid]) else NA,
    LD_N_max=if(any(valid)) max(d$LD_N_complete[valid]) else NA,focal_records=sum(d$variant_id==focal_id))
}
counts <- rbind(count(g),count(m))
counts <- rbind(transform(counts,version='annotated'),transform(counts,version='clean'))
counts$window_start1 <- limits[1]*1e6;counts$window_end1 <- limits[2]*1e6
counts$layering <- 'SNP < SV < TEV < focal; within class LD ascending, NA first'
write.table(counts,file.path(out,'plotted_input_counts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(focal_id=focal_id,GWAS_P=row$focal_GWAS_P,molQTL_P=row$focal_mol_P,
  PIP=if(is.null(pip)) NA else pip$PIP,PIP_model=if(is.null(pip)) NA else pip$model,
  PPH4=row$PPH4,focal_conditional_shared_PP=row$focal_conditional_shared_PP,
  focal_joint_shared_PP=row$focal_joint_shared_PP,GWAS_model=row$GWAS_model,molecular_model=row$molecular_model,
  GWAS_N=row$N_GWAS,molecular_N=row$N_mol),file.path(out,'focal_statistics_audit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
# Verify equal physical-coordinate panels for BOTH versions, including collected legends.
alignment <- list()
for(version in c('annotated','clean')) {
  audit_device <- tempfile(fileext='.pdf')
  grDevices::cairo_pdf(audit_device,family=FONT_FAMILY,width=WIDTH_MM/25.4,height=HEIGHT_MM/25.4)
  gt <- patchwork::patchworkGrob(if(version=='annotated') main_annotated else main_clean)
  invisible(grDevices::dev.off());unlink(audit_device)
  panel_cols <- gt$layout[grepl('^panel-[12]$',gt$layout$name),c('name','l','r')]
  stopifnot(nrow(panel_cols)==2,length(unique(panel_cols$l))==1,length(unique(panel_cols$r))==1)
  panel_cols$version <- version;alignment[[version]] <- panel_cols
}
write.table(do.call(rbind,alignment),file.path(out,'panel_alignment_audit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
message('PASS: annotated + clean; ',nrow(g$data),' GWAS records; ',nrow(m$data),' molQTL records. Outputs: ',out)
