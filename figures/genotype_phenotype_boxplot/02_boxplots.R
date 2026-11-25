#!/usr/bin/env Rscript
# Genotype–phenotype boxplots V4 | publication main figure
# Run: Rscript 02_boxplots.R [input_directory]
# Inputs: plot_data.csv, pairwise_LMM.csv, dominant_GWAS.csv
# Raw phenotype distributions; P values are imported, never recomputed.
suppressPackageStartupMessages(library(ggplot2))

# ---- User settings ----
COLORS_GT <- c("0/0"="#386C8C", "0/1"="#9689B9", "1/1"="#E75A46")
COLORS_DOM <- c("REF"="#386C8C", "INS"="#E75A46")
INK <- "#253246"
FILL_STRENGTH <- 0.78
FONT_FAMILY <- "Liberation Sans"
FIG_WIDTH_MM <- 85
FIG_HEIGHT_MM <- 85
DPI <- 600
BASE_FONT_PT <- 9
P_TEXT_PT <- 9
BOX_WIDTH <- 0.44
CAP_WIDTH <- 0.22
BOX_LW <- 0.40
MEDIAN_LW <- 0.60
AXIS_LW <- 0.45
BRACKET_LW <- 0.35
BOTTOM_PADDING <- 0.045
BRACKET_START <- 0.055
BRACKET_STEP <- 0.105
BRACKET_TIP <- 0.022
BRACKET_TEXT <- 0.018
TOP_PADDING <- 0.085
SAME_Y_LIMIT <- FALSE
SHOW_OUTLIERS <- TRUE
SHOW_JITTER <- FALSE
SHOW_MEAN <- FALSE
SHOW_SAMPLE_SIZE <- TRUE
SAVE_SVG <- TRUE

TRAIT_COL <- "abdominal_fat_percent"
Y_LABEL <- "Abdominal fat (%)"
X_LABEL_GT <- "Genotype"
EXPECTED_N <- 822L  # Set NULL for other datasets.
STEM_GT <- "01_三种基因型_腹脂率箱线图"
STEM_DOM <- "02_显性模型_插入有无_腹脂率箱线图"

# ---- Input and output directory ----
args <- commandArgs(trailingOnly=TRUE)
a <- commandArgs(trailingOnly=FALSE)
f <- sub("^--file=", "", a[grepl("^--file=", a)])
D <- normalizePath(if(length(args)) args[1] else
  if(length(f)) dirname(f[1]) else getwd(), mustWork=TRUE)
message("Input/output directory: ", D)
d <- read.csv(file.path(D, "plot_data.csv"), na.strings="NA")
pairs <- read.csv(file.path(D, "pairwise_LMM.csv"))
dom <- read.csv(file.path(D, "dominant_GWAS.csv"))
stopifnot(all(c("GT", TRAIT_COL) %in% names(d)),
          is.numeric(d$GT), all(d$GT %in% 0:2),
          is.numeric(d[[TRAIT_COL]]), all(is.finite(d[[TRAIT_COL]])),
          all(c("comparison", "P_Holm") %in% names(pairs)),
          "P" %in% names(dom), nrow(dom)==1L)
if(!is.null(EXPECTED_N)) stopifnot(nrow(d)==EXPECTED_N)
d$y <- d[[TRAIT_COL]]
d$genotype <- factor(c("0/0","0/1","1/1")[d$GT+1L],
                     levels=c("0/0","0/1","1/1"))
d$dominant <- factor(ifelse(d$GT==0,"REF","INS"), levels=c("REF","INS"))
stopifnot(all(table(d$genotype)>0), all(table(d$dominant)>0))

# ---- P-value format: literal e notation, three significant digits ----
p_label <- function(p, adjusted=FALSE) {
  stopifnot(length(p)==1L, is.finite(p), p>0, p<=1)
  s <- strsplit(sprintf("%.2e", p), "e", fixed=TRUE)[[1]]
  sci <- paste0(s[1], "e", as.integer(s[2]))
  symbol <- if(adjusted) "italic(P)[adj]" else "italic(P)"
  paste0(symbol, ' == "', sci, '"')
}
comparisons <- c("GT0_vs_GT1","GT1_vs_GT2","GT0_vs_GT2")
pv <- vapply(comparisons, function(z) {
  i <- which(pairs$comparison==z)
  stopifnot(length(i)==1L)
  pairs$P_Holm[i]
}, numeric(1))
p3_labels <- vapply(pv, p_label, character(1), adjusted=TRUE)
p2_label <- p_label(dom$P[1])

# ---- Opaque fills; standard Tukey statistics ----
blend_white <- function(colors) {
  m <- grDevices::col2rgb(colors)/255
  m <- m*FILL_STRENGTH + (1-FILL_STRENGTH)
  setNames(grDevices::rgb(m[1,],m[2,],m[3,]), names(colors))
}
box_stats <- function(data, variable) {
  lev <- levels(data[[variable]])
  b <- lapply(seq_along(lev), function(i) {
    y <- data$y[data[[variable]]==lev[i]]
    q <- quantile(y,c(.25,.5,.75),type=7,names=FALSE)
    fence <- c(q[1]-1.5*(q[3]-q[1]), q[3]+1.5*(q[3]-q[1]))
    inside <- y[y>=fence[1] & y<=fence[2]]
    out <- y[y<fence[1] | y>fence[2]]
    list(box=data.frame(group=lev[i],x=i,N=length(y),
      q1=q[1],median=q[2],q3=q[3],ymin=min(inside),ymax=max(inside),
      mean=mean(y)),
      out=data.frame(x=rep(i,length(out)),y=out))
  })
  list(boxes=do.call(rbind,lapply(b,function(z) z$box)),
       outliers=do.call(rbind,lapply(b,function(z) z$out)),
       levels=lev)
}

# ---- Axis range: retain zero, add visible space below whiskers ----
axis_base <- min(0,min(d$y))
ymax <- max(d$y)
yspan <- ymax-axis_base
stopifnot(is.finite(yspan),yspan>0)
y_lower <- axis_base-yspan*BOTTOM_PADDING
h3 <- ymax+yspan*(BRACKET_START+(0:2)*BRACKET_STEP)
h2 <- ymax+yspan*BRACKET_START
upper3 <- max(h3)+yspan*TOP_PADDING
upper2 <- h2+yspan*TOP_PADDING
if(SAME_Y_LIMIT) upper2 <- upper3
yb <- pretty(c(axis_base,ymax),n=5)
yb <- yb[yb>=axis_base & yb<=ymax]

theme_publication <- function() {
  theme_classic(base_size=BASE_FONT_PT,base_family=FONT_FAMILY)+
    theme(legend.position="none",panel.border=element_blank(),
      axis.line=element_line(colour=INK,linewidth=AXIS_LW),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_line(colour=INK,linewidth=.40),
      axis.ticks.length=grid::unit(1.2,"mm"),
      axis.text.x=element_text(colour=INK,size=9,lineheight=1.08,
                              margin=margin(t=4)),
      axis.text.y=element_text(colour=INK,size=9,margin=margin(r=3)),
      axis.title.x=element_text(colour=INK,size=10,margin=margin(t=7)),
      axis.title.y=element_text(colour=INK,size=10,margin=margin(r=7)),
      plot.margin=margin(5,6,5,5))
}

# ---- Separate whiskers, caps, opaque box and median ----
base_plot <- function(data,variable,palette,labels,xlab,upper) {
  s <- box_stats(data,variable)
  b <- s$boxes
  w <- BOX_WIDTH/2
  cap <- CAP_WIDTH/2
  p <- ggplot()
  if(SHOW_JITTER) {
    j <- data.frame(x=as.integer(data[[variable]]),y=data$y,
                    group=data[[variable]])
    p <- p+geom_point(data=j,aes(x,y,colour=group),
       position=position_jitter(width=.14,height=0,seed=167184720),
       size=.55,alpha=.18)+scale_colour_manual(values=palette)
  }
  p <- p+
    geom_segment(data=b,aes(x=x,xend=x,y=ymin,yend=q1),
                 colour=INK,linewidth=BOX_LW,lineend="butt")+
    geom_segment(data=b,aes(x=x,xend=x,y=q3,yend=ymax),
                 colour=INK,linewidth=BOX_LW,lineend="butt")+
    geom_segment(data=b,aes(x=x-cap,xend=x+cap,y=ymin,yend=ymin),
                 colour=INK,linewidth=BOX_LW)+
    geom_segment(data=b,aes(x=x-cap,xend=x+cap,y=ymax,yend=ymax),
                 colour=INK,linewidth=BOX_LW)+
    geom_rect(data=b,aes(xmin=x-w,xmax=x+w,ymin=q1,ymax=q3,fill=group),
              colour=INK,linewidth=BOX_LW,alpha=1)+
    geom_segment(data=b,aes(x=x-w,xend=x+w,y=median,yend=median),
                 colour=INK,linewidth=MEDIAN_LW,lineend="butt")
  if(SHOW_OUTLIERS && !SHOW_JITTER && nrow(s$outliers)>0) {
    p <- p+geom_point(data=s$outliers,aes(x,y),shape=16,
                      size=.9,colour="#909090",alpha=.65)
  }
  if(SHOW_MEAN) {
    p <- p+geom_point(data=b,aes(x=x,y=mean),shape=23,
                      size=1.8,fill="white",colour=INK,stroke=.4)
  }
  p+scale_fill_manual(values=blend_white(palette))+
    scale_x_continuous(breaks=seq_along(s$levels),
                       labels=unname(labels[s$levels]))+
    scale_y_continuous(breaks=yb)+
    coord_cartesian(xlim=c(.45,length(s$levels)+.55),
                    ylim=c(y_lower,upper),expand=FALSE,clip="on")+
    labs(x=xlab,y=Y_LABEL)+theme_publication()
}
add_bracket <- function(p,x1,x2,h,label) {
  tip <- yspan*BRACKET_TIP
  p+annotate("segment",x=c(x1,x1,x2),xend=c(x1,x2,x2),
             y=c(h-tip,h,h),yend=c(h,h,h-tip),
             colour=INK,linewidth=BRACKET_LW,lineend="square")+
    annotate("text",x=(x1+x2)/2,y=h+yspan*BRACKET_TEXT,
             label=label,parse=TRUE,vjust=0,size=P_TEXT_PT/2.845,
             family=FONT_FAMILY,colour=INK)
}
make_labels <- function(x) {
  n <- table(x)
  setNames(if(SHOW_SAMPLE_SIZE)
    paste0(names(n),"\n(n = ",as.integer(n),")") else names(n), names(n))
}

# ---- Construct both panels before export ----
p3 <- base_plot(d,"genotype",COLORS_GT,make_labels(d$genotype),
                X_LABEL_GT,upper3)
for(i in seq_along(comparisons)) {
  p3 <- add_bracket(p3,c(1,2,1)[i],c(2,3,3)[i],h3[i],p3_labels[i])
}
p2 <- base_plot(d,"dominant",COLORS_DOM,make_labels(d$dominant),NULL,upper2)
p2 <- add_bracket(p2,1,2,h2,p2_label)

# ---- Export: same original directory and filenames ----
save_plot <- function(p,stem) {
  file <- function(ext) file.path(D,paste0(stem,".",ext))
  ggsave(file("pdf"),p,device=grDevices::cairo_pdf,
         width=FIG_WIDTH_MM,height=FIG_HEIGHT_MM,units="mm",bg="white")
  ggsave(file("png"),p,device=grDevices::png,type="cairo",dpi=DPI,
         width=FIG_WIDTH_MM,height=FIG_HEIGHT_MM,units="mm",bg="white")
  if(SAVE_SVG) {
    if(requireNamespace("svglite",quietly=TRUE)) {
      ggsave(file("svg"),p,device=svglite::svglite,
             width=FIG_WIDTH_MM,height=FIG_HEIGHT_MM,units="mm",bg="white")
    } else warning("svglite is unavailable; SVG skipped.")
  }
  message("Saved: ",file.path(D,stem))
}
save_plot(p3,STEM_GT)
save_plot(p2,STEM_DOM)

# ---- Descriptive and exact boxplot statistics ----
summary_one <- function(x) c(N=length(x),mean=mean(x),SD=sd(x),
  median=median(x),Q1=unname(quantile(x,.25)),Q3=unname(quantile(x,.75)),
  min=min(x),max=max(x))
write_summary <- function(variable,group_name,filename,boxfile) {
  z <- do.call(rbind,lapply(split(d$y,d[[variable]]),summary_one))
  out <- data.frame(group=rownames(z),z,row.names=NULL)
  names(out)[1] <- group_name
  write.csv(out,file.path(D,filename),row.names=FALSE)
  write.csv(box_stats(d,variable)$boxes,file.path(D,boxfile),row.names=FALSE)
  print(out)
}
write_summary("genotype","genotype","三种基因型_表型描述统计.csv",
              "三种基因型_箱线图统计量.csv")
write_summary("dominant","group","显性模型_表型描述统计.csv",
              "显性模型_箱线图统计量.csv")
capture.output(sessionInfo(),file=file.path(D,"R_sessionInfo.txt"))
if(interactive()) {print(p3);print(p2)}
message("Finished V4. Output directory: ",D)
