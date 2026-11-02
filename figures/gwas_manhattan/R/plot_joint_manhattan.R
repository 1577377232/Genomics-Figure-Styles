#!/usr/bin/env Rscript
# Common R packages; every valid source point is passed to the SNP raster layer.
# SNP-only raster cache avoids drawing 12 million points again for each export.
# TE markers, thresholds, axes and text remain editable vectors in the PDF.
.libPaths(unique(c(Sys.getenv('GWAS_R_LIB', Sys.getenv('R_LIBS_USER')), .libPaths())))
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(jsonlite)
  library(ragg); library(png); library(digest); library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
opt <- function(key, default=NULL) {
  i <- match(key, args); if (is.na(i)) return(default)
  if (i == length(args)) stop('Missing value for ', key); args[i+1L]
}
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1])
out <- normalizePath(opt('--outdir', file.path(dirname(script), '..')), mustWork=TRUE)
for(subdir in c('audit','traits','modules','overview','cache','logs'))dir.create(file.path(out,subdir),showWarnings=FALSE,recursive=TRUE)
cfg <- fromJSON(opt('--config',file.path(out,'config/analysis.json')))
# Explicit user edits take precedence over the downloaded style defaults.
style_file <- opt('--style',file.path(out,cfg$display_style_file))
style <- fromJSON(style_file)
style_sha <- digest(style_file,algo='sha256',file=TRUE)
script_sha <- digest(script,algo='sha256',file=TRUE)
palette_name <- style$version
display_threshold <- as.numeric(opt('--threshold-p', as.character(style$threshold_p)))
stopifnot(is.finite(display_threshold),display_threshold>0,display_threshold<1,
          length(style$colours)==2L,style$chromosome_gap_bp>0)
cfg$gap_bp <- style$chromosome_gap_bp
traits <- fread(opt('--manifest', file.path(out, 'config/traits.tsv')), encoding='UTF-8')
sources <- fread(opt('--sources',file.path(out,'config/sources.tsv')),encoding='UTF-8')
stopifnot(nrow(traits)==cfg$expected_count, uniqueN(traits$trait_id)==cfg$expected_count, setequal(unique(sources$trait_id),traits$trait_id))
stopifnot(all(traits$kind %in% c('single','module')))
expected_p<-ifelse(traits$kind[match(sources$trait_id,traits$trait_id)]=='module','p_score','p_wald')
stopifnot(all(sources$p_column==expected_p))
if (identical(cfg$Z_analysis,'joint')) {
  stopifnot(all(sources$panel %in% c('A','Z_joint')))
  pairs<-sources[,.(n=.N,panels=paste(sort(panel),collapse=','),Ns=uniqueN(N)),by=.(trait_id,marker)]
  stopifnot(all(pairs$n==2L),all(pairs$panels=='A,Z_joint'),all(pairs$Ns==1L))
}
selection <- opt('--trait-id','all')
ids <- if (identical(selection,'all')) traits$trait_id else strsplit(selection,',',fixed=TRUE)[[1]]
stopifnot(all(ids %in% traits$trait_id), !anyDuplicated(ids))
jobs <- as.integer(opt('--jobs','1')); stopifnot(is.finite(jobs),jobs>=1L,jobs<=30L)
setDTthreads(1)
font <- style$font
pt <- 72.27/25.4
fai_file <- opt('--fai', cfg$fai)
stopifnot(digest(fai_file,algo='sha256',file=TRUE)==cfg$fai_sha256)
fai <- fread(fai_file,header=FALSE,select=1:2,col.names=c('chr','length_bp'),colClasses=list(character=1,numeric=2))
layout <- fai[match(cfg$chromosomes,chr)]
stopifnot(!anyNA(layout),!anyDuplicated(layout$chr),all(layout$length_bp>0))
biological_order<-match(layout$chr,c(as.character(1:39),'Z','W'))
stopifnot(!anyNA(biological_order),all(diff(biological_order)>0))
layout[, offset_bp:=c(0,head(cumsum(length_bp+cfg$gap_bp),-1L))]
layout[, midpoint_mb:=(offset_bp+length_bp/2)/1e6]
layout[, snp_colour:=style$colours$SNP]
xmax <- (tail(layout$offset_bp,1)+tail(layout$length_bp,1))/1e6
# Greedy spacing in physical millimetres; keep Z and all points, hide only colliding tick labels.
axis_width_mm <- style$width_mm-22
tick_spacing <- style$tick_min_spacing_mm
ticks <- 1L
for (i in 2:(nrow(layout)-1L)) {
  if ((layout$midpoint_mb[i]-layout$midpoint_mb[tail(ticks,1L)])/xmax*axis_width_mm>=tick_spacing &&
      (tail(layout$midpoint_mb,1L)-layout$midpoint_mb[i])/xmax*axis_width_mm>=tick_spacing) ticks<-c(ticks,i)
}
ticks <- c(ticks,nrow(layout))
for (i in 2:length(ticks)) stopifnot(diff(layout$midpoint_mb[ticks[(i-1):i]])/xmax*axis_width_mm>=tick_spacing)
fwrite(layout,file.path(out,'audit/chromosome_layout.tsv'),sep='\t')
fwrite(layout[ticks],file.path(out,'audit/chromosome_ticks.tsv'),sep='\t')

map_chr <- function(x) {
  x <- sub('^chr','',as.character(x),ignore.case=TRUE)
  # Numeric PLINK codes are mapped only when declared; biological chr39 is distinct.
  for(alias in names(cfg$chromosome_aliases))x[x==alias]<-cfg$chromosome_aliases[[alias]]
  x
}
te_type <- function(id) fifelse(grepl('\\.INS\\.',id,ignore.case=TRUE),'TE insertion',
                       fifelse(grepl('\\.DEL\\.',id,ignore.case=TRUE),'TE deletion','TE other'))

load_marker <- function(tid, marker) {
  wanted_marker <- marker
  ss <- sources[trait_id==tid & marker==wanted_marker]
  stopifnot(nrow(ss)>0L,uniqueN(ss$p_column)==1L,ss$p_column[1L] %in% c('p_wald','p_score'))
  chunks<-list(); checks<-list(); bychr<-list()
  for (i in seq_len(nrow(ss))) {
    s<-ss[i]; message(tid,' ',marker,' ',s$panel,': reading full association file')
    hash<-digest(s$file,algo='sha256',file=TRUE)
    stopifnot(hash==s$expected_sha256)
    pcol<-s$p_column;is_single<-pcol=='p_wald'
    cols<-c('chr','ps',pcol);if(marker=='TE')cols<-c(cols,'rs')
    if(marker=='TE' && is_single)cols<-c(cols,'beta','se')
    d<-fread(s$file,select=cols,colClasses=list(character=c('chr',if(marker=='TE')'rs')),
             na.strings=c('NA','NaN','nan','-nan','NAN','-NaN',''),showProgress=FALSE)
    setnames(d,pcol,'P')
    reported<-nrow(d);stopifnot(reported==s$expected_reported)
    valid<-is.finite(d$P)&d$P>=0&d$P<=1
    if(marker=='TE' && is_single)valid<-valid&is.finite(d$beta)&is.finite(d$se)&d$se>0
    invalid<-sum(!valid);d<-d[valid]
    stopifnot(nrow(d)==s$expected_valid,all(is.finite(d$ps)),all(d$ps>0),all(d$ps==floor(d$ps)))
    if(any(d$P==0))stop('P=0 underflow requires original precise source score; no arbitrary cap: ',s$file)
    stopifnot(isTRUE(all.equal(min(d$P),s$expected_minP,tolerance=2e-6)))
    d[,chr:=map_chr(chr)]
    if(s$panel=='A')stopifnot(all(d$chr %in% as.character(1:38))) else stopifnot(all(d$chr=='Z'))
    all_j<-match(d$chr,fai$chr);stopifnot(!anyNA(all_j),all(d$ps<=fai$length_bp[all_j]))
    valid_count<-nrow(d);native_minP<-min(d$P);native_bonf<-cfg$alpha/valid_count
    native_significant<-sum(d$P<=native_bonf)
    # Explicit user chromosome exclusion is separate from invalid tests or downsampling.
    excluded_count<-sum(d$chr %in% cfg$excluded_plot_chromosomes)
    d<-d[!chr %in% cfg$excluded_plot_chromosomes]
    j<-match(d$chr,layout$chr);stopifnot(!anyNA(j))
    d[,`:=`(x_mb=(ps+layout$offset_bp[j])/1e6,y=-log10(P),source_panel=s$panel)]
    if(marker=='TE')d[,variant_type:=te_type(rs)]
    if(marker=='TE' && is_single)d[,c('beta','se'):=NULL]
    bychr[[i]]<-d[,.(valid_tests=.N,minP=min(P)),by=.(chr)][,`:=`(trait_id=tid,marker=marker,panel=s$panel)]
    checks[[i]]<-data.table(trait_id=tid,marker=marker,panel=s$panel,N=s$N,file=s$file,
      source_sha256=hash,reported=reported,valid=valid_count,plotted=nrow(d),excluded_chr38=excluded_count,
      invalid_or_unestimable=invalid,minP=native_minP,minP_plotted=min(d$P),P_zero=0L,
      native_Bonferroni=native_bonf,native_Bonferroni_tests=native_significant)
    chunks[[i]]<-d
  }
  list(data=rbindlist(chunks,use.names=TRUE),qc=rbindlist(checks),bychr=rbindlist(bychr))
}

publication_theme <- function() {
  theme_classic(base_family=font,base_size=10)+
    theme(text=element_text(colour='#202830'),
      axis.line=element_line(colour='#202830',linewidth=.7/pt),
      axis.ticks=element_line(colour='#202830',linewidth=.6/pt),
      axis.ticks.length=grid::unit(1.3,'mm'),
      axis.text=element_text(size=style$axis_tick_pt,colour='#202830'),
      axis.title=element_text(size=style$axis_title_pt),
      axis.title.x=element_text(margin=margin(t=4)),axis.title.y=element_text(margin=margin(r=5)),
      plot.caption=element_text(size=7.1,hjust=1,margin=margin(t=2),colour='#68737C'),
      plot.margin=margin(3,4,2,3,'mm'),legend.position='inside',
      legend.position.inside=c(.995,.985),legend.justification=c(1,1),legend.direction='horizontal',
      legend.background=element_rect(fill='white',colour=NA),
      legend.key=element_blank(),legend.title=element_blank(),legend.text=element_text(size=style$legend_pt),
      legend.key.height=grid::unit(4,'mm'),legend.key.width=grid::unit(4,'mm'),
      legend.spacing.x=grid::unit(2,'mm'),legend.margin=margin(0,0,0,0,'mm'))
}

panel_dimensions_mm <- function(ymax,caption) {
  # Measure the actual ggplot panel, rather than estimating its margins. This
  # keeps raster SNP circles round and the same physical diameter as vector TEs.
  probe<-ggplot()+geom_blank(aes(x=0,y=0))+
    scale_x_continuous(breaks=layout$midpoint_mb[ticks],labels=layout$chr[ticks],expand=c(0,0),limits=c(0,xmax))+
    scale_y_continuous(breaks=seq(0,ymax,by=if(ymax<=14)2 else 4),limits=c(0,ymax),expand=c(0,0))+
    labs(x='Chromosome',y=expression(-log[10](P)),caption=caption)+publication_theme()
  ragg::agg_capture(width=style$width_mm,height=style$height_mm,units='mm',res=600,background='white')
  on.exit(grDevices::dev.off())
  g<-ggplotGrob(probe);grid::grid.newpage();grid::grid.draw(g);grid::grid.force()
  cell<-g$layout[g$layout$name=='panel',c('t','l','b','r')]
  stopifnot(nrow(cell)==1L)
  vp<-paste0('panel.',paste(unlist(cell),collapse='-'))
  grid::seekViewport(vp)
  c(width=grid::convertWidth(grid::unit(1,'npc'),'mm',valueOnly=TRUE),
    height=grid::convertHeight(grid::unit(1,'npc'),'mm',valueOnly=TRUE))
}

plot_trait <- function(tid) {
  tryCatch({
    begin<-proc.time()[['elapsed']];t<-traits[trait_id==tid]
    task_sources<-sources[trait_id==tid]
    stamp<-file.info(task_sources$file)
    fingerprint<-digest(list(style=style_sha,config=cfg,code=script_sha,threshold=display_threshold,trait=as.list(t),sources=as.data.frame(task_sources),
                             source_size=stamp$size,source_mtime=as.numeric(stamp$mtime)),algo='sha256')
    receipt_file<-file.path(out,'audit',paste0(tid,'__PASS.json'))
    if(file.exists(receipt_file)) {
      old<-fromJSON(receipt_file)
      if(identical(old$run_fingerprint,fingerprint) && file.exists(old$png) && file.exists(old$pdf) &&
         file.exists(file.path(out,'cache',paste0(tid,'__plot.rds'))) &&
         digest(old$png,algo='sha256',file=TRUE)==old$png_sha256 && digest(old$pdf,algo='sha256',file=TRUE)==old$pdf_sha256) {
        message(tid,': REUSE PASS');return(list(trait_id=tid,status='PASS'))
      }
    }
    snp<-load_marker(tid,'SNP');te<-load_marker(tid,'TE')
    M_snp<-sum(snp$qc$valid);M_te<-sum(te$qc$valid)
    plotted_snp<-nrow(snp$data);plotted_te<-nrow(te$data)
    snp_line<-cfg$alpha/M_snp;te_line<-cfg$alpha/M_te
    peak<-max(snp$data$y,te$data$y,-log10(display_threshold))
    ymax<-ceiling((peak+max(2,peak*.14))/2)*2
    stopifnot(is.finite(ymax),ymax>peak)
    N<-snp$qc[panel=='A',N][1];stopifnot(N==te$qc[panel=='A',N][1])
    method_label<-if(t$kind=='module')'GEMMA mvGWAS score' else 'GEMMA'
    caption<-paste0(tid,'  |  ',method_label,'  |  N = ',N)
    panel_mm<-panel_dimensions_mm(ymax,caption)
    raster_file<-file.path(out,'cache',paste0(tid,'__all_SNP_points_600dpi.png'))
    message(tid,': rasterizing ALL ',plotted_snp,' SNP tests on requested chromosomes at 600 dpi; chr38 excluded')
    # No sampling, bin maxima, pixel filtering, or min-P merging is performed.
    # Only this dense layer is rasterized; all original points reach geom_point.
    cloud<-ggplot(snp$data,aes(x=x_mb,y=y))+
      geom_point(shape=style$shared_point_shape,size=style$snp_size_mm,colour=style$colours$SNP,alpha=style$snp_alpha,stroke=0)+
      scale_x_continuous(limits=c(0,xmax),expand=c(0,0))+
      scale_y_continuous(limits=c(0,ymax),expand=c(0,0))+
      theme_void()+theme(plot.margin=margin(0,0,0,0),plot.background=element_rect(fill='transparent',colour=NA),
                         panel.background=element_rect(fill='transparent',colour=NA))
    # Panel-only raster exceeds the final panel's pixel dimensions; avoid repeat rendering.
    ggsave(raster_file,cloud,device=ragg::agg_png,width=panel_mm['width'],height=panel_mm['height'],units='mm',dpi=600,bg='transparent',limitsize=FALSE)
    stopifnot(nrow(cloud$data)==plotted_snp)
    rm(cloud);gc(FALSE)
    image<-png::readPNG(raster_file,native=TRUE)
    te$data[,`:=`(significant=P<=te_line,display_threshold_pass=P<=display_threshold)]
    cols<-unlist(style$colours)
    legend_data<-data.table(x_mb=c(0,0),y=c(0,0),marker=factor(c('SNP','TE'),levels=c('SNP','TE')))
    N<-snp$qc[panel=='A',N][1];stopifnot(N==te$qc[panel=='A',N][1])
    z_snp<-sum(snp$qc[panel!='A',valid]);z_te<-sum(te$qc[panel!='A',valid])
    # One plotting reference line, exactly as requested; it is not labelled Bonferroni.
    caption<-paste0(tid,'  |  ',method_label,'  |  N = ',N)
    display_title<-if(tid=='D056')'Breast muscle proportion' else t$display_name_en
    p<-ggplot()+annotation_raster(image,xmin=0,xmax=xmax,ymin=0,ymax=ymax,interpolate=FALSE)+
      geom_hline(yintercept=-log10(display_threshold),colour='#838B91',linewidth=.35,linetype='dashed',show.legend=FALSE)+
      geom_point(data=te$data,aes(x=x_mb,y=y),colour=style$colours$TE,fill=style$colours$TE,
                 shape=style$shared_point_shape,size=style$te_size_mm,stroke=style$te_stroke_mm,alpha=style$te_alpha)+
      geom_point(data=legend_data,aes(x=x_mb,y=y,colour=marker),alpha=0,size=0,shape=16,show.legend=TRUE)+
      scale_colour_manual(values=cols,breaks=c('SNP','TE'),drop=FALSE)+
      scale_x_continuous(breaks=layout$midpoint_mb[ticks],labels=layout$chr[ticks],expand=c(0,0),limits=c(0,xmax))+
      scale_y_continuous(breaks=seq(0,ymax,by=if(ymax<=14)2 else 4),limits=c(0,ymax),expand=c(0,0))+
      labs(x='Chromosome',y=expression(-log[10](P)),caption=caption)+
      guides(colour=guide_legend(override.aes=list(alpha=1,size=1.8,shape=16),nrow=1))+
      publication_theme()
    p<-p+annotate('text',x=xmax*.006,y=ymax*.965,label=display_title,hjust=0,vjust=1,
                   family=font,size=style$title_pt/pt,colour='#202830')
    # Legend fixed upper right must not hide any input peak; check its intended genomic rectangle.
    legend_xmin<-xmax*.84;legend_ymin<-ymax*.88
    overlap<-sum(snp$data$x_mb>=legend_xmin&snp$data$y>=legend_ymin)+
             sum(te$data$x_mb>=legend_xmin&te$data$y>=legend_ymin)
    if(overlap>0L)p<-p+theme(legend.position='top',legend.justification='right',legend.direction='horizontal')
    # These plots contain only ~10k vector TE points plus one full-SNP raster.
    built<-ggplot_build(p)
    stopifnot(nrow(built$data[[3]])==plotted_te,all(is.finite(built$data[[3]]$y)))
    stopifnot(max(built$data[[3]]$y)<=ymax,min(snp$data$y)>=0,max(snp$data$y)<ymax)
    stopifnot(nrow(built$data[[2]])==1L,built$data[[2]]$yintercept==-log10(display_threshold))
    stopifnot(uniqueN(built$data[[3]]$colour)==1L,uniqueN(built$data[[3]]$shape)==1L)
    stopifnot(style$snp_size_mm==style$te_size_mm,style$te_stroke_mm==0,style$shared_point_shape==16)
    stopifnot(!'38' %in% layout$chr,!any(snp$data$chr=='38'),!any(te$data$chr=='38'))
    prefix<-file.path(out,if(t$kind=='module')'modules' else 'traits',paste0(tid,'__joint_manhattan'))
    ggsave(paste0(prefix,'.png'),p,device=ragg::agg_png,width=style$width_mm,height=style$height_mm,units='mm',dpi=600,bg='white')
    ggsave(paste0(prefix,'.pdf'),p,device=grDevices::cairo_pdf,width=style$width_mm,height=style$height_mm,units='mm',bg='white')
    # Serialize the completed grob, not a ggplot quosure environment containing full raw SNP tables.
    saveRDS(ggplotGrob(p),file.path(out,'cache',paste0(tid,'__plot.rds')),compress=FALSE)
    checks<-rbind(snp$qc,te$qc);checks[,`:=`(display_name_en=t$display_name_en,palette=palette_name)]
    fwrite(checks,file.path(out,'audit',paste0(tid,'__sources.tsv')),sep='\t')
    fwrite(rbind(snp$bychr,te$bychr),file.path(out,'audit',paste0(tid,'__chromosome_counts.tsv')),sep='\t')
    top<-te$data[order(P)][1:min(20,.N),.(rs,chr,ps,P,source_panel,variant_type,significant,display_threshold_pass)]
    fwrite(top,file.path(out,'audit',paste0(tid,'__top_TE.tsv')),sep='\t')
    result<-list(trait_id=tid,display_name_zh=t$display_name_zh,display_name_en=t$display_name_en,
      Z_analysis=if(is.null(cfg$Z_analysis))'unspecified' else cfg$Z_analysis,
      analysis_kind=t$kind,p_column=if(t$kind=='module')'p_score' else 'p_wald',module_dimensions=t$module_dimensions,
      category=t$category,N_autosomes=N,M_SNP=M_snp,M_TE=M_te,
      minP_SNP=min(snp$qc$minP),minP_TE=min(te$qc$minP),
      minP_plotted_SNP=min(snp$data$P),minP_plotted_TE=min(te$data$P),
      Bonferroni_SNP=snp_line,Bonferroni_TE=te_line,
      significant_SNP_tests=sum(snp$data$P<=snp_line),significant_TE_tests=sum(te$data$significant),
      display_threshold_p=display_threshold,threshold_line_count=1L,
      SNP_below_display_threshold=sum(snp$data$P<=display_threshold),
      TE_below_display_threshold=sum(te$data$display_threshold_pass),
      SNP_colour=style$colours$SNP,TE_colour=style$colours$TE,TE_type_distinction=FALSE,Z_sex_distinction=FALSE,
      TE_display_shape=style$shared_point_shape,SNP_display_shape=style$shared_point_shape,
      SNP_point_size_mm=style$snp_size_mm,TE_point_size_mm=style$te_size_mm,threshold_label=FALSE,
      chr38_plotted=FALSE,SNP_excluded_chr38=sum(snp$qc$excluded_chr38),TE_excluded_chr38=sum(te$qc$excluded_chr38),
      chromosome_gap_bp=style$chromosome_gap_bp,width_mm=style$width_mm,height_mm=style$height_mm,
      plot_style=style$version,style_config_sha256=style_sha,run_fingerprint=fingerprint,
      SNP_raster_panel_width_mm=unname(panel_mm['width']),SNP_raster_panel_height_mm=unname(panel_mm['height']),
      SNP_raster_measured_to_match_vector_panel=TRUE,
      SNP_Z_tests=z_snp,TE_Z_tests=z_te,SNP_plotted=plotted_snp,TE_plotted=plotted_te,
      ymax=ymax,clipped_peaks=0L,P_zero=0L,downsampling=FALSE,
      PDF_text_vector=TRUE,PDF_TE_vector=TRUE,SNP_only_raster=TRUE,
      legend_relocated=(overlap>0L),BW42_source_flag=t$BW42_source_flag,
      existing_QC_HOLD_TE=ifelse(is.na(t$existing_QC_HOLD_TE),'',t$existing_QC_HOLD_TE),palette=palette_name,reference=cfg$reference_build,
      fai_sha256=cfg$fai_sha256,style_commit=cfg$style_commit,dpi=600,
      png=paste0(prefix,'.png'),pdf=paste0(prefix,'.pdf'),elapsed_seconds=proc.time()[['elapsed']]-begin,
      png_sha256=digest(paste0(prefix,'.png'),algo='sha256',file=TRUE),pdf_sha256=digest(paste0(prefix,'.pdf'),algo='sha256',file=TRUE),
      status='PASS')
    unlink(file.path(out,'audit',paste0(tid,'__FAIL.json')))
    write_json(result,file.path(out,'audit',paste0(tid,'__PASS.json')),pretty=TRUE,auto_unbox=TRUE,na='null',digits=16)
    message(tid,': PASS; ',round(result$elapsed_seconds,1),' seconds')
    list(trait_id=tid,status='PASS')
  },error=function(e){
    write_json(list(trait_id=tid,status='FAIL',error=conditionMessage(e)),
      file.path(out,'audit',paste0(tid,'__FAIL.json')),pretty=TRUE,auto_unbox=TRUE)
    message(tid,': FAIL: ',conditionMessage(e));list(trait_id=tid,status='FAIL',error=conditionMessage(e))
  })
}

if (!('--assemble-only' %in% args)) {
  results<-if(jobs>1L)parallel::mclapply(ids,plot_trait,mc.cores=min(jobs,length(ids)),mc.preschedule=FALSE) else lapply(ids,plot_trait)
  if(any(vapply(results,function(x)x$status!='PASS',logical(1))))quit(status=1)
}
qc_files<-file.path(out,'audit',paste0(ids,'__PASS.json'))
stopifnot(all(file.exists(qc_files)))
qc<-lapply(qc_files,fromJSON)
fwrite(rbindlist(qc,fill=TRUE),file.path(out,'audit/trait_qc.tsv'),sep='\t')
write_json(list(phase='Selected actual traits and frozen multivariate modules',
               selected_ids=ids,real_trait_count=nrow(traits),style_commit=cfg$style_commit,provenance=cfg,user_style=style,
               packages=lapply(c('ggplot2','data.table','patchwork','ragg','png','jsonlite','digest'),function(x)list(name=x,version=as.character(packageVersion(x)))),
               figures=qc),file.path(out,'audit/plot_manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
fwrite(data.table(trait_id=character(),error=character()),file.path(out,'audit/failures.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'audit/R_sessionInfo.txt'))
if(length(ids)>1L && !('--no-overview' %in% args)) {
  # Each analysis kind gets its own booklet and readable category overview pages.
  selected<-traits[match(ids,trait_id)]
  for(kind in unique(selected$kind)) {
    # Explicit extraction avoids confusing data.table column names with loop values.
    kid<-ids[traits$kind[match(ids,traits$trait_id)]==kind]
    groups<-split(kid,traits$category[match(kid,traits$trait_id)])
    pages<-unlist(lapply(groups,function(v)split(v,ceiling(seq_along(v)/5L))),recursive=FALSE)
    label<-if(kind=='module')'modules' else 'quantitative_traits'
    for(page in seq_along(pages)) {
      these<-pages[[page]]
      plots<-lapply(these,function(tid)wrap_elements(full=readRDS(file.path(out,'cache',paste0(tid,'__plot.rds')))))
      composite<-wrap_plots(plots,ncol=1)+plot_annotation(tag_levels='A')&
        theme(plot.tag=element_text(family=font,size=12,face='bold'),plot.tag.position=c(.005,.995))
      name<-sprintf('SNP_TE_%s_overview_%02d',label,page)
      prefix<-file.path(out,'overview',name)
      ggsave(paste0(prefix,'.pdf'),composite,device=grDevices::cairo_pdf,width=style$width_mm,height=style$height_mm*length(these),units='mm',limitsize=FALSE,bg='white')
      ggsave(paste0(prefix,'.png'),composite,device=ragg::agg_png,width=style$width_mm,height=style$height_mm*length(these),units='mm',dpi=600,limitsize=FALSE,bg='white')
      ggsave(paste0(prefix,'_screen.png'),composite,device=ragg::agg_png,width=style$width_mm,height=style$height_mm*length(these),units='mm',dpi=180,limitsize=FALSE,bg='white')
    }
    book<-file.path(out,'overview',paste0('SNP_TE_',length(kid),'_',label,'_booklet.pdf'))
    grDevices::cairo_pdf(book,width=style$width_mm/25.4,height=style$height_mm/25.4,onefile=TRUE)
    for(tid in kid) {grid::grid.newpage();grid::grid.draw(readRDS(file.path(out,'cache',paste0(tid,'__plot.rds'))))}
    grDevices::dev.off()
  }
}
message('Completed ',length(ids),' joint Manhattan plots. Output: ',out)
