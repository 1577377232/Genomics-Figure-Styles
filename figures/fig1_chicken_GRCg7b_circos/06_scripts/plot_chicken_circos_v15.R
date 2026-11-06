#!/usr/bin/env Rscript
# ==============================================================================
# V15 鸡 GRCg7b：紧凑无引导轴、统一核型、TE 扇区裁剪
# ==============================================================================
# 作者说明：基于项目 V13 精修；参考用户提供的人类图与水稻图图例，独立实现。
# 依赖：circlize + R 自带 grid；项目内 R_library 优先，不修改系统 R 库。
# 输入：已有真实窗口统计表、主要染色体长度和近似着丝粒坐标。
# 输出：SVG、PDF、600 dpi PNG、预览图、尺度表、输入校验记录和运行参数。
#
# 坐标：GRCg7b，0-based、左闭右开；所有轨道共用同一窗口表。
# 外 -> 内：LINE / LTR / SINE / DNA 覆盖比例 + SNP / SV / TIP 位点密度。
# 前四条轨道：类内合并重叠区间后的覆盖比例，不是多态性。
# V11 前四类 TE 使用独立显示上限（LTR q95，其余 q99）；三窗口平滑仅用于图形展示。
# 后三条轨道：独立位点数 / 实际窗口长度；缺失值不替换为 0。
# TIP extended 包含 simple repeat / low complexity / satellite，且包括 INS/DEL。
# adaptive_discovery 为探索性 SNP 集合，不代表所有样本均可分型。
#
# 常见修改入口：default_config() 的 colors / geometry / typography / scales。
# RStudio 中可 source() 本文件，再运行：
#   cfg <- default_config("/path/to/chicken_7b_TE_circos")
#   cfg$snp_panel <- "population95"
#   render_circos(cfg)
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 20)


# ------------------------------------------------------------------------------
# 1. 集中配置：建议优先在这里修改图形，而不改动下面的绘图函数
# ------------------------------------------------------------------------------

detect_project_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
    return(dirname(dirname(script)))
  }
  if (dir.exists(file.path(getwd(), "04_window_statistics"))) return(getwd())
  file.path(getwd(), "chicken_7b_TE_circos")
}

default_config <- function(project_dir = detect_project_dir()) {
  list(
    project_dir = project_dir,
    window_bp = 1000000,
    snp_panel = "adaptive_discovery",  # population95 / discovery50 / adaptive_discovery
    tip_set = "extended",             # strict / extended
    sv_set = "all",                   # all / nonTE_assessed / nonrepeat_assessed
    table = NULL,
    chrom_sizes = NULL,
    centromeres = NULL,
    out_dir = NULL,                    # NULL: 07_figures/R_v15；保留旧版结果
    formats = c("svg", "pdf", "png"),
    dpi = 600,
    preview_dpi = 160,
    validate_only = FALSE,
    all_windows = FALSE,
    te_classes = c("LINE", "LTR", "SINE", "DNA"),
    palette = "pangenome_spectrum",
    colors = NULL,
    geometry = list(
      width_in = 9.8, height_in = 10.0,
      opening_deg = 48, opening_center_deg = 0, clockwise = FALSE,
      gap_deg = 1.1,
      min_micro_deg = 2.1,
      min_W_deg = 3.0,
      micro_chromosomes = c(as.character(10:39), "W"),
      xlim = c(-1.08, 1.20), ylim = c(-1.06, 1.06),
      figure_region = c(0.025, 0.99, 0.085, 0.925),
      coordinate_height = 0.050,
      chromosome_label_height = 0.044,
      ideogram_height = 0.022,
      ideogram_track_margin = c(0.0065, 0.0065),
      te_height = 0.080,
      te_track_margin = c(0.0045, 0.0045),
      ideogram_fill = "#F1F4F6",
      heat_height = 0.058,
      track_margin = c(0.0025, 0.0025),
      cell_padding = c(0, 0, 0, 0),
      ideogram_lwd = 1.10,
      track_lwd = 0.65,
      annotation_lwd = 0.38,
      label_padding_mm = 0.40,
      unit_circle_segments = 1440
    ),
    typography = list(
      family = "Liberation Sans", base_pt = 10,
      title_pt = 11.5, subtitle_pt = 7.5,
      chromosome_pt = 7.2,
      tick_macro_pt = 5.0, tick_micro_pt = 4.6,
      legend_pt = 8.2, coverage_axis_pt = 4.2
    ),
    scales = list(
      te_quantile = 0.99,
      te_quantiles = c(LINE = 0.99, LTR = 0.95, SINE = 0.99, DNA = 0.99),
      heat_quantile = 0.99,
      snp_cap_population95 = 25000,
      snp_cap_discovery = 50000,
      snp_display_divisor = 1,
      snp_transform = "linear",
      sv_transform = "log1p", tip_transform = "log1p"
    ),
    te_smoothing = list(kernel = c(0.25, 0.5, 0.25), minimum_windows = 3L,
      points_per_window = 10L, method = "length-weighted 3-window triangular mean; shape-preserving cubic Hermite"),
    centromere_color = "#E75A46",
    na_fill = "#F1F3F5", na_hatch = "#9AA6AF",
    border_color = "#253246", text_color = "#273746"
  )
}

get_palette_colors <- function(palette) {
  if (!identical(palette, "pangenome_spectrum")) stopf("V15 固定使用 pangenome_spectrum。")
  # Explicit user-defined Circos mapping; no hue adaptation or alternate theme.
  c(LINE = "#386C8C", LTR = "#7BB6AC", SINE = "#9689B9", DNA = "#B4AD74",
    SNP = "#5899B2", SV = "#A1A8BE", TIP = "#E75A46")
}

apply_palette <- function(cfg) {
  if (is.null(cfg$colors)) cfg$colors <- get_palette_colors(cfg$palette)
  required <- c("LINE", "LTR", "SINE", "DNA", "SNP", "SV", "TIP")
  if (!all(required %in% names(cfg$colors))) {
    stopf("自定义颜色必须包含：%s", paste(required, collapse = ", "))
  }
  cfg$colors <- cfg$colors[required]
  if (!identical(toupper(cfg$colors), get_palette_colors(cfg$palette))) {
    stopf("V15 七轨颜色必须使用用户指定的固定 HEX 映射。")
  }
  cfg
}

ideogram_fill_table <- function(chromosomes, cfg) {
  # One neutral, from the existing center background; no chromosome color cycle.
  data.frame(chromosome = chromosomes, fill_color = cfg$geometry$ideogram_fill,
    fill_source = "repository theme border_soft #E3E9ED mixed with 50% white")
}

# ------------------------------------------------------------------------------
# 2. 输入读取与数据检查：禁止混入其他参考坐标，禁止将 NA 当成 0
# ------------------------------------------------------------------------------

stopf <- function(...) stop(sprintf(...), call. = FALSE)

canonical_chromosome <- function(x) {
  sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE)
}

read_tsv <- function(path) {
  if (!file.exists(path)) stopf("输入文件不存在：%s", path)
  read.delim(path, check.names = FALSE, na.strings = c("NA", ""),
             quote = "", comment.char = "", stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA", fileEncoding = "UTF-8")
}

resolve_paths <- function(cfg) {
  cfg$project_dir <- normalizePath(cfg$project_dir, mustWork = TRUE)
  if (is.null(cfg$table)) {
    cfg$table <- file.path(cfg$project_dir, "04_window_statistics", "v11",
                           sprintf("windows_%skb.tsv", cfg$window_bp / 1000))
  }
  if (is.null(cfg$chrom_sizes)) {
    cfg$chrom_sizes <- file.path(cfg$project_dir, "01_reference", "main.chrom.sizes")
  }
  if (is.null(cfg$centromeres)) {
    cfg$centromeres <- file.path(cfg$project_dir, "01_reference", "v14",
                                "GRCg7b.centromeres.approximate.tsv")
  }
  if (is.null(cfg$out_dir)) cfg$out_dir <- file.path(cfg$project_dir, "07_figures", "R_v15")
  cfg
}

track_columns <- function(cfg) {
  snp <- switch(cfg$snp_panel,
    population95 = "SNP_per_Mb",
    discovery50 = "SNP_discovery50_per_Mb",
    adaptive_discovery = "SNP_adaptive_discovery_per_Mb",
    stopf("未知 SNP 集合：%s", cfg$snp_panel)
  )
  sv <- switch(cfg$sv_set,
    all = "SV_per_Mb",
    nonTE_assessed = "SV_nonTE_assessed_per_Mb",
    nonrepeat_assessed = "SV_nonrepeat_assessed_per_Mb",
    stopf("未知 SV 集合：%s", cfg$sv_set)
  )
  if (!cfg$tip_set %in% c("strict", "extended")) stopf("未知 TIP 集合：%s", cfg$tip_set)
  c(SNP = snp, SV = sv, TIP = paste0("TIP_", cfg$tip_set, "_per_Mb"))
}

read_and_validate_inputs <- function(cfg) {
  expected <- c(as.character(1:39), "Z", "W")
  sizes <- read.table(cfg$chrom_sizes, sep = "\t", header = FALSE,
                       col.names = c("chrom", "length_bp"), colClasses = c("character", "numeric"))
  sizes$chrom <- canonical_chromosome(sizes$chrom)
  if (anyDuplicated(sizes$chrom) || !setequal(sizes$chrom, expected)) {
    stopf("染色体长度文件必须恰好包含 Chr1-39、Z、W，共 41 条且没有重复。")
  }
  sizes <- sizes[match(expected, sizes$chrom), , drop = FALSE]
  if (any(!is.finite(sizes$length_bp) | sizes$length_bp <= 0 |
          sizes$length_bp != floor(sizes$length_bp))) stopf("染色体长度必须是正整数。")
  lengths <- setNames(sizes$length_bp, sizes$chrom)

  d <- read_tsv(cfg$table)
  coverage_columns <- paste0(c("LINE", "LTR", "SINE", "DNA"), "_coverage_fraction")
  heat_columns <- track_columns(cfg)
  count_columns <- sub("_per_Mb$", "_count", heat_columns)
  required <- c("Chromosome", "Start", "End", "Midpoint", "effective_window_bp",
                coverage_columns, heat_columns, count_columns)
  absent <- setdiff(required, names(d))
  if (length(absent)) stopf("窗口表缺少列：%s", paste(absent, collapse = ", "))
  d$Chromosome <- canonical_chromosome(d$Chromosome)
  if (!setequal(unique(d$Chromosome), expected)) stopf("窗口表染色体集合与 41 条主要染色体不一致。")
  numeric_columns <- setdiff(required, "Chromosome")
  if (!all(vapply(d[numeric_columns], is.numeric, logical(1)))) {
    stopf("必需的坐标、覆盖度、计数和密度列必须是数值或 NA。")
  }
  d <- d[order(match(d$Chromosome, expected), d$Start, d$End), , drop = FALSE]
  key <- paste(d$Chromosome, d$Start, d$End, sep = ":")
  if (anyDuplicated(key)) stopf("窗口表存在重复坐标。")
  if (anyNA(d[c("Start", "End", "Midpoint", "effective_window_bp")])) {
    stopf("窗口坐标或实际窗口长度存在 NA。")
  }
  if (any(d$Start < 0 | d$End <= d$Start | d$End > lengths[d$Chromosome]) ||
      any(d$Start != floor(d$Start) | d$End != floor(d$End))) stopf("窗口坐标越界或不是整数。")
  if (any(abs(d$effective_window_bp - (d$End - d$Start)) > 1e-7) ||
      any(abs(d$Midpoint - (d$Start + d$End) / 2) > 1e-7)) {
    stopf("窗口长度或中点与 Start/End 不一致。")
  }
  for (chrom in expected) {
    z <- d[d$Chromosome == chrom, , drop = FALSE]
    starts <- seq(0, lengths[chrom] - 1, by = cfg$window_bp)
    ends <- pmin(starts + cfg$window_bp, lengths[chrom])
    if (nrow(z) != length(starts) || any(z$Start != starts) || any(z$End != ends)) {
      stopf("Chr%s 的窗口有缺失、重叠，或不是指定的 %s bp 窗口。", chrom, cfg$window_bp)
    }
  }
  for (column in coverage_columns) {
    x <- d[[column]]
    if (anyNA(x) || any(!is.finite(x) | x < 0 | x > 1)) {
      stopf("%s 必须是 0 到 1 之间的真实覆盖比例，且不能缺失。", column)
    }
  }
  for (i in seq_along(heat_columns)) {
    density <- d[[heat_columns[i]]]
    counts <- d[[count_columns[i]]]
    if (any(is.na(density) != is.na(counts))) stopf("%s 的计数和密度缺失标记不一致。", names(heat_columns)[i])
    ok <- !is.na(density)
    if (any(!is.finite(density[ok]) | density[ok] < 0) ||
        any(!is.finite(counts[ok]) | counts[ok] < 0 | counts[ok] != floor(counts[ok]))) {
      stopf("%s 存在无效计数或密度。", names(heat_columns)[i])
    }
    calculated <- counts[ok] * 1e6 / d$effective_window_bp[ok]
    if (any(abs(density[ok] - calculated) > 1e-6 * pmax(1, calculated))) {
      stopf("%s 密度与 count / 实际窗口长度不一致。", names(heat_columns)[i])
    }
  }

  centers <- NULL
  if (!identical(cfg$centromeres, "none")) {
    centers <- read_tsv(cfg$centromeres)
    needed <- c("chrom", "start", "end", "evidence")
    if (!all(needed %in% names(centers))) stopf("着丝粒表缺少 chrom/start/end/evidence 列。")
    centers$chrom <- canonical_chromosome(centers$chrom)
    if (any(!centers$chrom %in% expected) ||
        anyNA(centers[needed]) ||
        any(centers$start < 0 | centers$end <= centers$start |
            centers$end > lengths[centers$chrom])) stopf("着丝粒坐标缺失或超出 GRCg7b 边界。")
    if (!all(centers$evidence %in% c("high_coverage_transfer", "approximate_transfer", "approximate_partial_transfer"))) {
      message("着丝粒 evidence 保留在 QC；主图全部使用同一种红色矩形。")
    }
  }
  list(data = d, sizes = sizes, lengths = lengths, chromosomes = expected,
       centers = centers, heat_columns = heat_columns)
}


# ------------------------------------------------------------------------------
# 3. 扇区与尺度：小染色体可展开；所有轨道共用同一坐标映射
# ------------------------------------------------------------------------------

make_layout <- function(inputs, cfg) {
  g <- cfg$geometry
  chromosomes <- inputs$chromosomes
  lengths <- inputs$lengths
  available <- 360 - g$opening_deg - g$gap_deg * (length(chromosomes) - 1)
  if (!is.finite(available) || available <= 0 || g$opening_deg < 0 || g$opening_deg >= 360 ||
      g$gap_deg < 0 || g$min_micro_deg < 0) stopf("开口、间隔或小染色体角度配置无效。")
  micro <- intersect(g$micro_chromosomes, chromosomes)
  minimum <- setNames(rep(0, length(chromosomes)), chromosomes)
  minimum[micro] <- g$min_micro_deg
  minimum["W"] <- max(minimum["W"], g$min_W_deg)
  if (sum(minimum) >= available) stopf("小染色体最小角度之和超过可用角度。")
  arcs <- setNames(rep(NA_real_, length(chromosomes)), chromosomes)
  remaining <- chromosomes
  repeat {
    scale <- (available - sum(arcs, na.rm = TRUE)) / sum(lengths[remaining])
    newly_fixed <- remaining[lengths[remaining] * scale < minimum[remaining]]
    if (!length(newly_fixed)) break
    arcs[newly_fixed] <- minimum[newly_fixed]
    remaining <- setdiff(remaining, newly_fixed)
  }
  arcs[remaining] <- lengths[remaining] * scale
  # 右侧开口；用户选择 Chr1 右上起始、逆时针。sector.width 仅控制显示角度。
  direction <- if (g$clockwise) -1 else 1
  first <- g$opening_center_deg + if (g$clockwise) -g$opening_deg / 2 else g$opening_deg / 2
  start <- first + direction * c(0, head(cumsum(arcs + g$gap_deg), -1))
  layout <- data.frame(
    chrom = chromosomes, length_bp = unname(lengths[chromosomes]),
    chromosome = chromosomes, real_length_bp = unname(lengths[chromosomes]),
    displayed_sector_width = unname(arcs / sum(arcs)), displayed_angle = unname(arcs),
    expansion_factor = unname(arcs / (available * lengths[chromosomes] / sum(lengths))),
    start_degrees = unname(start), end_degrees = unname(start + direction * arcs),
    arc_degrees = unname(arcs), degrees_per_bp = unname(arcs / lengths[chromosomes]),
    enlargement_factor = unname(arcs / lengths[chromosomes] / scale)
  )
  rownames(layout) <- chromosomes
  expected_end <- first + direction * (360 - g$opening_deg)
  if (abs(tail(layout$end_degrees, 1) - expected_end) > 1e-7) {
    stopf("扇区角度总和异常。")
  }
  layout
}

positive_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(1)
  upper <- unname(quantile(x, probability, type = 7))
  if (upper <= 0) upper <- max(x)
  if (upper <= 0) upper <- 1
  upper
}

make_scales <- function(inputs, cfg) {
  d <- inputs$data
  classes <- cfg$te_classes
  result <- vector("list", 7)
  names(result) <- c(classes, "SNP", "SV", "TIP")

  # Four independent coverage scales. Only the LTR display cap changes to q95.
  # Raw values and the V10 bounded display smoothing remain unchanged.
  for (class in classes) {
    probability <- cfg$scales$te_quantiles[[class]]
    raw <- d[[paste0(class, "_coverage_fraction")]]
    quantile_value <- positive_quantile(raw, probability)
    quantum <- 10^(floor(log10(quantile_value)) - 1)
    cap <- min(1, ceiling(quantile_value / quantum - 1e-10) * quantum)
    if (!is.finite(cap) || cap <= 0) stopf("%s 独立 TE 上限无效。", class)
    result[[class]] <- list(column = paste0(class, "_coverage_fraction"),
      cap = cap, class_q99 = positive_quantile(raw, 0.99),
      selected_quantile = probability, quantile_value = quantile_value,
      transform = "linear", display_divisor = 0.01, unit = "coverage (%)",
      scale_basis = sprintf("class-specific raw q%g, rounded upward to 2 significant digits", probability * 100))
  }

  for (name in c("SNP", "SV", "TIP")) {
    column <- inputs$heat_columns[[name]]
    if (name == "SNP") {
      cap <- if (cfg$snp_panel == "population95") cfg$scales$snp_cap_population95 else cfg$scales$snp_cap_discovery
      divisor <- cfg$scales$snp_display_divisor
      unit <- if (divisor == 1000) "SNP/kb" else "SNP/MB"
      transform <- cfg$scales$snp_transform
    } else {
      cap <- positive_quantile(d[[column]], cfg$scales$heat_quantile)
      divisor <- 1
      unit <- paste0(name, "/MB")
      transform <- cfg$scales[[paste0(tolower(name), "_transform")]]
    }
    if (!is.finite(cap) || cap <= 0 || !transform %in% c("linear", "log1p")) {
      stopf("%s 的显示尺度无效。", name)
    }
    result[[name]] <- list(column = column, cap = cap, transform = transform,
                           display_divisor = divisor, unit = unit, has_missing = anyNA(d[[column]]))
  }
  result
}

scale_fraction <- function(x, scale) {
  fraction <- if (scale$transform == "log1p") log1p(x) / log1p(scale$cap) else x / scale$cap
  pmax(0, pmin(1, fraction))
}


# ------------------------------------------------------------------------------
# 4. circlize 绘图引擎：sector 为染色体，track 为数据层；grid 仅排版图例
# ------------------------------------------------------------------------------

load_circlize <- function(cfg) {
  local_library <- file.path(cfg$project_dir, "06_scripts", "R_library")
  if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stopf("缺少 circlize；请先运行 06_scripts/install_circos_R_dependencies.sh。")
  }
}

density_colors <- function(fraction, color) {
  # Only the colour transfer changes; density transforms/caps remain unchanged.
  base <- as.numeric(col2rgb(color)) / 255
  light <- grDevices::rgb(0.68 + 0.32 * base[1], 0.68 + 0.32 * base[2], 0.68 + 0.32 * base[3])
  dark <- grDevices::rgb(0.70 * base[1], 0.70 * base[2], 0.70 * base[3])
  ramp <- grDevices::colorRamp(c("#FFFFFF", light, color, dark), space = "Lab")
  mapped <- approx(c(0, 0.25, 0.65, 1), c(0, 1/3, 2/3, 1),
    xout = pmax(0, pmin(1, fraction)), rule = 2)$y
  rgb <- pmin(pmax(ramp(mapped), 0), 255)
  grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)
}

format_value <- function(value) {
  if (abs(value) < 1e-10) return("0")
  if (value >= 1000) return(format(round(value), big.mark = ",", scientific = FALSE, trim = TRUE))
  format(signif(value, 3), scientific = FALSE, trim = TRUE)
}

chromosome_font <- function(chrom, cfg) {
  cfg$typography$chromosome_pt
}

coordinate_ticks <- function(length_bp, angle, chrom) {
  macro <- chrom %in% c(as.character(1:9), "Z")
  major_step <- if (chrom %in% as.character(6:9)) 1e7 else if (macro) 5e7 else 1e7
  major <- seq(0, length_bp - 1, by = major_step)
  step <- if (macro) 1e7 else 5e6
  minor <- setdiff(seq(0, length_bp - 1, by = step), major)
  list(major = major, minor = minor,
       labels = ifelse(chrom == "1" & major == 0, "0 MB", format(major / 1e6, trim = TRUE)))
}

draw_coordinate_track_v15 <- function(inputs, layout, cfg) {
  store <- new.env(parent = emptyenv())
  circlize::circos.trackPlotRegion(ylim = c(0, 1),
    track.height = cfg$geometry$coordinate_height, bg.col = "white", bg.border = NA,
    panel.fun = function(x, y) {
      chrom <- circlize::get.cell.meta.data("sector.index")
      z <- coordinate_ticks(inputs$lengths[[chrom]], layout[chrom, "arc_degrees"], chrom)
      store[[chrom]] <- z
      # Natural axis: baseline and thin ticks only; no closed coordinate frames.
      circlize::circos.lines(c(0, inputs$lengths[[chrom]]), c(0.04, 0.04),
        col = "#A9B2BB", lwd = 0.55)
      if (length(z$minor)) circlize::circos.segments(z$minor, 0.04, z$minor, 0.19,
        col = "#B0BAC3", lwd = 0.45)
      circlize::circos.segments(z$major, 0.04, z$major, 0.27,
        col = "#8997A5", lwd = 0.80)
    })
  mget(inputs$chromosomes, envir = store, inherits = FALSE)
}

draw_chromosome_label_track_v15 <- function(inputs, layout, cfg) {
  # Names occupy a separate white annotation band, outside the ideogram.
  circlize::circos.trackPlotRegion(ylim = c(0, 1),
    track.height = cfg$geometry$chromosome_label_height,
    bg.col = "white", bg.border = NA, panel.fun = function(x, y) NULL)
}

annotation_box <- function(chrom, bp, y, track, label, pt, font, kind, cfg, adj_x = 0.5) {
  polar <- circlize::circlize(bp, y, sector.index = chrom, track.index = track)
  theta <- unname(polar[1, "theta"]) * pi / 180
  radius <- unname(polar[1, "rou"])
  w <- strwidth(label, units = "inches", cex = pt / cfg$typography$base_pt, font = font)
  h <- strheight("Mg", units = "inches", cex = pt / cfg$typography$base_pt, font = font)
  tx <- sin(theta); ty <- -cos(theta)
  shift <- (adj_x - 0.5) * w
  data.frame(kind = kind, chromosome = chrom, label = label, position_bp = bp,
    y = y, track_index = track, angle_degrees = theta * 180 / pi,
    font_pt = pt, font = font, adj_x = adj_x,
    center_x_in = grconvertX(radius * cos(theta), "user", "inches") + shift * tx,
    center_y_in = grconvertY(radius * sin(theta), "user", "inches") + shift * ty,
    width_in = w, height_in = h, tangent_x = tx, tangent_y = ty,
    radial_x = cos(theta), radial_y = sin(theta), stringsAsFactors = FALSE)
}

boxes_overlap <- function(a, b, padding_mm) {
  at <- c(a$tangent_x, a$tangent_y); ar <- c(a$radial_x, a$radial_y)
  bt <- c(b$tangent_x, b$tangent_y); br <- c(b$radial_x, b$radial_y)
  delta <- c(b$center_x_in - a$center_x_in, b$center_y_in - a$center_y_in)
  pad <- padding_mm / 25.4 / 2
  aw <- a$width_in / 2 + pad; ah <- a$height_in / 2 + pad
  bw <- b$width_in / 2 + pad; bh <- b$height_in / 2 + pad
  dot <- function(x, y) sum(x * y)
  for (axis in list(at, ar, bt, br)) {
    extent <- aw * abs(dot(at, axis)) + ah * abs(dot(ar, axis)) +
      bw * abs(dot(bt, axis)) + bh * abs(dot(br, axis))
    if (abs(dot(delta, axis)) > extent) return(FALSE)
  }
  TRUE
}

box_on_device <- function(box, cfg) {
  center <- c(box$center_x_in, box$center_y_in)
  corners <- do.call(rbind, lapply(c(-1, 1), function(i) do.call(rbind,
    lapply(c(-1, 1), function(j) center + i * box$width_in / 2 *
      c(box$tangent_x, box$tangent_y) + j * box$height_in / 2 * c(box$radial_x, box$radial_y)))))
  all(corners[, 1] > 0.02 & corners[, 1] < cfg$geometry$width_in - 0.02 &
      corners[, 2] > 0.02 & corners[, 2] < cfg$geometry$height_in - 0.02)
}

plan_outer_annotations <- function(inputs, layout, ticks, cfg) {
  placed <- list(); names <- list(); coordinates <- list(); tick_rows <- list()
  hits <- function(box) length(placed) && any(vapply(placed,
    function(other) boxes_overlap(box, other, cfg$geometry$label_padding_mm), logical(1)))
  for (chrom in inputs$chromosomes) {
    selected <- NULL
    for (y in c(0.05, 0.15, 0.25)) {
      box <- annotation_box(chrom, inputs$lengths[[chrom]] / 2, y, 1, chrom,
        chromosome_font(chrom, cfg), 2,
        "chromosome", cfg)
      if (box_on_device(box, cfg) && !hits(box)) { selected <- box; break }
    }
    if (is.null(selected)) stopf("Chr%s 名称无法在标签带内放下，请增大扇区或画布。", chrom)
    names[[length(names) + 1L]] <- selected
    placed[[length(placed) + 1L]] <- selected
  }
  for (chrom in inputs$chromosomes) {
    z <- ticks[[chrom]]
    for (j in seq_along(z$major)) {
      pt <- if (chrom %in% c(as.character(1:9), "Z")) cfg$typography$tick_macro_pt else cfg$typography$tick_micro_pt
      tick_font <- 2
      candidate <- annotation_box(chrom, z$major[j], 0.47, 2, z$labels[j], pt, tick_font, "coordinate", cfg)
      radius_in <- abs(diff(grconvertX(c(0, 1), "user", "inches")))
      enough <- layout[chrom, "arc_degrees"] * pi / 180 * radius_in > candidate$width_in + 0.7 / 25.4
      status <- if (enough) "collision" else "insufficient_sector_arc"
      if (enough) for (adj in c(0.5, 0, 1)) {
        box <- annotation_box(chrom, z$major[j], 0.47, 2, z$labels[j], pt, tick_font, "coordinate", cfg, adj)
        if (box_on_device(box, cfg) && !hits(box)) {
          candidate <- box; status <- "shown"; placed[[length(placed) + 1L]] <- box; break
        }
      }
      candidate$label_status <- status
      coordinates[[length(coordinates) + 1L]] <- candidate
    }
    for (type in c("major", "minor")) if (length(z[[type]])) {
      tick_rows[[length(tick_rows) + 1L]] <- data.frame(chromosome = chrom,
        position_bp = z[[type]], position_Mb = z[[type]] / 1e6, tick_type = type)
    }
  }
  if (length(placed) > 1) for (i in seq_len(length(placed) - 1)) for (j in (i + 1):length(placed)) {
    if (boxes_overlap(placed[[i]], placed[[j]], cfg$geometry$label_padding_mm)) stopf("文字碰撞检查失败。")
  }
  list(names = do.call(rbind, names), coordinates = do.call(rbind, coordinates),
       ticks = do.call(rbind, tick_rows), collision_pairs = 0L)
}

draw_outer_annotations <- function(plan, cfg) {
  rows <- rbind(plan$names, plan$coordinates[plan$coordinates$label_status == "shown",
    names(plan$names), drop = FALSE])
  for (i in seq_len(nrow(rows))) {
    z <- rows[i, ]
    circlize::circos.text(z$position_bp, z$y, z$label, sector.index = z$chromosome,
      track.index = z$track_index, facing = "outside", niceFacing = TRUE,
      adj = c(z$adj_x, 0.5), cex = z$font_pt / cfg$typography$base_pt,
      font = z$font, col = cfg$text_color)
  }
}

draw_ideogram_v15 <- function(inputs, layout, cfg) {
  markers <- list()
  fills <- ideogram_fill_table(inputs$chromosomes, cfg)
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = cfg$geometry$ideogram_height,
    track.margin = cfg$geometry$ideogram_track_margin,
    bg.col = "white", bg.border = NA,
    panel.fun = function(x, y) {
      chrom <- circlize::get.cell.meta.data("sector.index")
      circlize::circos.rect(0, 0, inputs$lengths[[chrom]], 1,
        col = fills$fill_color[match(chrom, fills$chromosome)], border = NA)
      centers <- inputs$centers
      if (!is.null(centers)) {
        centers <- centers[centers$chrom == chrom, , drop = FALSE]
        if (nrow(centers)) for (i in seq_len(nrow(centers))) {
          z <- centers[i, ]; midpoint <- (z$start + z$end) / 2
          half <- min(1.0, layout[chrom, "arc_degrees"] * 0.40) / layout[chrom, "degrees_per_bp"] / 2
          a <- max(0, midpoint - half); b <- min(inputs$lengths[[chrom]], midpoint + half)
          circlize::circos.rect(a, 0, b, 1,
            col = cfg$centromere_color, border = NA)
          markers[[length(markers) + 1L]] <<- data.frame(chromosome = chrom,
            midpoint_bp = midpoint, source_start_bp = z$start, source_end_bp = z$end,
            marker_start_bp = a, marker_end_bp = b, evidence = z$evidence,
            rendering = "solid red rectangle; full ideogram height",
            fill = cfg$centromere_color, bottom_y = 0, top_y = 1)
        }
      }
      draw_track_boundaries(inputs$lengths[[chrom]], 1, cfg, cfg$geometry$ideogram_lwd)
    })
  if (length(markers)) do.call(rbind, markers) else data.frame()
}

draw_track_boundaries <- function(length_bp, upper, cfg, lwd = NULL) {
  if (is.null(lwd)) lwd <- cfg$geometry$track_lwd
  circlize::circos.rect(0, 0, length_bp, upper,
    col = NA, border = cfg$border_color, lwd = lwd)
}

shape_preserving_profile <- function(x, y, grid) {
  # Weighted harmonic-mean Hermite slopes. Turning points get zero slope;
  # unlike unconstrained splines, each interval stays within its two nodes.
  h <- diff(x); delta <- diff(y) / h
  slope <- numeric(length(x))
  slope[1] <- delta[1]; slope[length(x)] <- tail(delta, 1)
  for (i in 2:(length(x) - 1)) {
    if (delta[i - 1] * delta[i] > 0) {
      w1 <- 2 * h[i] + h[i - 1]
      w2 <- h[i] + 2 * h[i - 1]
      slope[i] <- (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i])
    }
  }
  j <- findInterval(grid, x, all.inside = TRUE)
  t <- (grid - x[j]) / h[j]
  value <- (2 * t^3 - 3 * t^2 + 1) * y[j] +
    (t^3 - 2 * t^2 + t) * h[j] * slope[j] +
    (-2 * t^3 + 3 * t^2) * y[j + 1] +
    (t^3 - t^2) * h[j] * slope[j + 1]
  if (any(value < pmin(y[j], y[j + 1]) - 1e-10 |
          value > pmax(y[j], y[j + 1]) + 1e-10)) stopf("保形插值超出局部节点包络。")
  value
}

make_te_profiles <- function(inputs, cfg) {
  nodes <- list(); curves <- list()
  for (class in cfg$te_classes) for (chrom in inputs$chromosomes) {
    d <- inputs$data[inputs$data$Chromosome == chrom, , drop = FALSE]
    raw <- d[[paste0(class, "_coverage_fraction")]]
    smooth <- raw
    if (length(raw) >= cfg$te_smoothing$minimum_windows) {
      for (i in seq_along(raw)) {
        j <- intersect(seq_along(raw), (i - 1):(i + 1))
        weight <- cfg$te_smoothing$kernel[j - i + 2] * d$effective_window_bp[j]
        smooth[i] <- sum(raw[j] * weight) / sum(weight)
      }
    }
    length_bp <- inputs$lengths[[chrom]]
    x <- c(0, d$Midpoint, length_bp)
    y <- c(smooth[1], smooth, tail(smooth, 1))
    grid <- sort(unique(c(x, seq(0, length_bp,
      length.out = max(30L, nrow(d) * cfg$te_smoothing$points_per_window)))))
    fit <- shape_preserving_profile(x, y, grid)
    if (any(fit < min(y) - 1e-10 | fit > max(y) + 1e-10)) stopf("TE 拟合超出节点包络：%s Chr%s", class, chrom)
    fit <- pmax(0, pmin(1, fit))
    nodes[[length(nodes) + 1L]] <- data.frame(TE_class = class, chromosome = chrom,
      start_bp = d$Start, end_bp = d$End, midpoint_bp = d$Midpoint,
      effective_window_bp = d$effective_window_bp, raw_fraction = raw,
      display_smoothed_fraction = smooth, smoothing_applied = length(raw) >= cfg$te_smoothing$minimum_windows)
    curves[[length(curves) + 1L]] <- data.frame(TE_class = class,
      chromosome = chrom, position_bp = grid, display_fraction = fit)
  }
  list(nodes = do.call(rbind, nodes), curves = do.call(rbind, curves))
}

te_axis_ticks <- function(scale) {
  # Round only the middle tick, then place it at its true percentage coordinate.
  c(0, signif(scale$cap * 50, 2) / 100, scale$cap)
}

draw_te_panel_clipped <- function(d, length_bp, scale, class, chrom, cfg) {
  # A polar polygon must explicitly return along the INNER ARC. Letting a
  # graphics device close it would create a Cartesian chord below the baseline.
  # A real per-panel clip path also contains all outline stroke/antialias edges.
  track <- circlize::get.current.track.index()
  polar <- circlize::circlize(c(0, length_bp), c(0, scale$cap),
    sector.index = chrom, track.index = track)
  theta0 <- polar[1, "theta"]; theta1 <- polar[2, "theta"]
  inner <- circlize::get.cell.meta.data("cell.bottom.radius", chrom, track)
  outer <- circlize::get.cell.meta.data("cell.top.radius", chrom, track)
  angle <- seq(theta0, theta1, length.out = max(32L, ceiling(abs(theta1-theta0)/0.05)+1L))
  radians <- angle * pi / 180
  clip_x <- c(outer*cos(radians), rev(inner*cos(radians)))
  clip_y <- c(outer*sin(radians), rev(inner*sin(radians)))
  normalized <- pmax(0, pmin(1, d$display_fraction / scale$cap))
  curve <- circlize::circlize(d$position_bp, normalized*scale$cap,
    sector.index = chrom, track.index = track)
  curve_angle <- curve[, "theta"] * pi/180
  top_x <- curve[, "rou"]*cos(curve_angle); top_y <- curve[, "rou"]*sin(curve_angle)
  # Reverse sampled baseline: the closing edge is radial, never a long chord.
  fill_x <- c(top_x, rev(inner*cos(radians)))
  fill_y <- c(top_y, rev(inner*sin(radians)))
  if (any(!is.finite(normalized) | normalized<0 | normalized>1) ||
      any(curve[,"rou"] < inner-1e-10 | curve[,"rou"] > outer+1e-10)) stopf("TE 绘图超出轨道。")
  bounds <- data.frame(TE_class=class,chromosome=chrom,track_index=track,
    inner_radius=inner,outer_radius=outer,start_degrees=theta0,end_degrees=theta1,
    normalized_min=min(normalized),normalized_max=max(normalized),
    raw_curve_max=max(d$display_fraction),cap=scale$cap,
    baseline_vertices=length(angle),actual_clip_path=TRUE,baseline_closed_by_arc=TRUE,
    clip_component=chrom,clip_group=class,
    legacy_chord_inward_sagitta=inner*(1-cos(abs(theta1-theta0)*pi/360)))
  list(bounds=bounds,clip_x=clip_x,clip_y=clip_y,fill_x=fill_x,fill_y=fill_y,top_x=top_x,top_y=top_y)

}

draw_te_tracks_v15 <- function(inputs, scales, cfg) {
  profiles <- make_te_profiles(inputs, cfg)
  bounds <- list(); layers <- list()
  for (class in cfg$te_classes) {
    scale <- scales[[class]]
    layers[[class]] <- list()
    circlize::circos.trackPlotRegion(ylim=c(0,scale$cap),
      track.height=cfg$geometry$te_height,track.margin=cfg$geometry$te_track_margin,
      bg.col="white",bg.border=NA,
      panel.fun=function(x,y) {
        chrom <- circlize::get.cell.meta.data("sector.index")
        d <- profiles$curves[profiles$curves$TE_class==class & profiles$curves$chromosome==chrom,]
        geom <- draw_te_panel_clipped(d,inputs$lengths[[chrom]],scale,class,chrom,cfg)
        layers[[class]][[chrom]] <<- geom
        bounds[[length(bounds)+1L]] <<- geom$bounds
      })
  }
  ndcx <- function(x) grconvertX(x,"user","ndc")
  ndcy <- function(y) grconvertY(y,"user","ndc")
  for (i in seq_along(cfg$te_classes)) {
    class <- cfg$te_classes[i]; group <- layers[[class]]
    # Four compound paths, each with 41 disjoint sector components, avoid the
    # limited clipping-path resource pool of R 4.1 Cairo. Gaps stay excluded.
    clip_x <- unlist(lapply(group,`[[`,"clip_x"),use.names=FALSE)
    clip_y <- unlist(lapply(group,`[[`,"clip_y"),use.names=FALSE)
    mask <- grid::pathGrob(x=ndcx(clip_x),y=ndcy(clip_y),
      id.lengths=vapply(group,function(z)length(z$clip_x),integer(1)),
      rule="winding",gp=grid::gpar(fill="black",col=NA))
    grid::pushViewport(grid::viewport(clip=mask,name=paste0("TEclip_",class)))
    for (geom in group) {
      grid::grid.polygon(x=ndcx(geom$fill_x),y=ndcy(geom$fill_y),
        gp=grid::gpar(fill=grDevices::adjustcolor(cfg$colors[[class]],alpha.f=0.58),col=NA))
      grid::grid.lines(x=ndcx(geom$top_x),y=ndcy(geom$top_y),
        gp=grid::gpar(col=cfg$colors[[class]],lwd=0.65,lineend="butt"))
    }
    grid::popViewport()
    for (chrom in inputs$chromosomes) circlize::circos.rect(0,0,inputs$lengths[[chrom]],scales[[class]]$cap,
      sector.index=chrom,track.index=i+3L,col=NA,border=cfg$border_color,lwd=cfg$geometry$track_lwd)
  }
  profiles$render_bounds <- do.call(rbind,bounds)
  profiles
}

draw_te_axis_v15 <- function(scales,cfg) {
  # One set of opening-side axes. TE labels are left-anchored next to their
  # actual axes; variant gradients/units replace the removed page legends.
  tracks<-c(cfg$te_classes,"SNP","SV","TIP")
  rows<-list();strips<-list()
  for(i in seq_along(tracks)) {
    track<-tracks[i];index<-i+3L;scale<-scales[[track]];is_te<-i<=4
    ticks<-if(is_te)te_axis_ticks(scale) else c(0,0.5,1)
    raw<-if(is_te)ticks else if(scale$transform=="log1p")expm1(ticks*log1p(scale$cap)) else ticks*scale$cap
    if(!is_te) {
      edges<-seq(0,1,length.out=65L)
      xl<--circlize::convert_x(0.85,"mm",sector.index="1",track.index=index)
      xr<--circlize::convert_x(0.15,"mm",sector.index="1",track.index=index)
      colors<-density_colors((head(edges,-1)+tail(edges,-1))/2,cfg$colors[[track]])
      circlize::circos.rect(xl,head(edges,-1),xr,tail(edges,-1),sector.index="1",track.index=index,col=colors,border=NA)
      strips[[track]]<-data.frame(track=track,fraction=(head(edges,-1)+tail(edges,-1))/2,color=colors)
    }
    circlize::circos.yaxis(side="left",at=ticks,labels=FALSE,
      sector.index="1",track.index=index,
      tick.length=circlize::convert_x(0.30,"mm",sector.index="1",track.index=index),
      col=cfg$text_color,lwd=0.45)
    polar<-circlize::circlize(rep(0,length(ticks)),ticks,sector.index="1",track.index=index)
    rows[[i]]<-data.frame(track=track,TE_class=if(is_te)track else NA_character_,
      coverage_fraction=if(is_te)raw else NA_real_,raw_density_per_MB=if(is_te)NA_real_ else raw,
      y=ticks,track_index=index,
      label=if(is_te)paste0(vapply(raw*100,format_value,character(1)),"%") else vapply(raw,format_value,character(1)),
      kind=if(is_te)"TE_axis_tick" else "variant_axis_tick",radius=polar[,"rou"],anchor_degrees=polar[,"theta"])
  }
  d<-do.call(rbind,rows);d<-d[order(d$radius),]
  pt<-cfg$typography$coverage_axis_pt
  widths<-strwidth(d$label,units="inches",cex=pt/cfg$typography$base_pt,font=2)
  height<-strheight("Mg",units="inches",cex=pt/cfg$typography$base_pt,font=2)
  unit_in<-abs(diff(grconvertX(c(0,1),"user","inches")))
  offset<-seq(0,nrow(d)-1)*(height+0.50/25.4)
  r_in<-stats::isoreg(d$radius*unit_in-offset)$yf+offset
  is_te<-d$kind=="TE_axis_tick"
  angle<-d$anchor_degrees-atan(ifelse(is_te,0.40,1.05)/25.4/
    (ifelse(is_te,min(d$radius[is_te]),min(d$radius[!is_te]))*unit_in))*180/pi
  theta<-angle*pi/180;radius<-r_in/unit_in
  xx<-radius*cos(theta);yy<-radius*sin(theta)
  for(j in seq_len(nrow(d)))graphics::text(xx[j],yy[j],d$label[j],srt=angle[j]-90,
    cex=pt/cfg$typography$base_pt,font=2,col=cfg$text_color,adj=c(0,0.5))
  result<-data.frame(kind=d$kind,chromosome="1",label=d$label,
    position_bp=0,y=d$y,track_index=d$track_index,
    angle_degrees=angle,font_pt=pt,font=2,adj_x=0,
    center_x_in=grconvertX(xx,"user","inches")+widths/2*sin(theta),
    center_y_in=grconvertY(yy,"user","inches")-widths/2*cos(theta),
    width_in=widths,height_in=height,tangent_x=sin(theta),tangent_y=-cos(theta),
    radial_x=cos(theta),radial_y=sin(theta),track=d$track,TE_class=d$TE_class,
    coverage_fraction=d$coverage_fraction,raw_density_per_MB=d$raw_density_per_MB,
    label_radius=radius,anchor_radius=d$radius,leaders_drawn=FALSE,
    radial_adjustment_mm=(r_in-d$radius*unit_in)*25.4)
  # V15: omit variant unit/transform titles; keep true ticks and gradients.
  for(j in seq_len(nrow(result)))if(!box_on_device(result[j,],cfg))stopf("开口坐标标签超出画布。")
  for(i in seq_len(nrow(result)-1))for(j in (i+1):nrow(result)) {
    if(boxes_overlap(result[i,],result[j,],0.25))stopf("开口坐标标签重叠：%s / %s",result$label[i],result$label[j])
  }
  attr(result,"gradient_strips")<-do.call(rbind,strips)
  result
}

draw_heat_tracks_v15 <- function(inputs, scales, cfg) {
  for (name in c("SNP", "SV", "TIP")) {
    scale <- scales[[name]]
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = cfg$geometry$heat_height,
      bg.col = "white", bg.border = NA,
      panel.fun = function(x, y) {
        chrom <- circlize::get.cell.meta.data("sector.index")
        d <- inputs$data[inputs$data$Chromosome == chrom, , drop = FALSE]
        raw <- d[[scale$column]]
        valid <- !is.na(raw)
        if (any(valid)) {
          circlize::circos.rect(d$Start[valid], 0, d$End[valid], 1,
            col = density_colors(scale_fraction(raw[valid], scale), cfg$colors[[name]]),
            border = NA)
        }
        if (any(!valid)) {
          circlize::circos.rect(d$Start[!valid], 0, d$End[!valid], 1,
            col = cfg$na_fill, border = NA)
          circlize::circos.rect(d$Start[!valid], 0, d$End[!valid], 1,
            density = 22, angle = 45, col = cfg$na_hatch, border = NA, lwd = 0.3)
        }
        draw_track_boundaries(inputs$lengths[[chrom]], 1, cfg)
      }
    )
  }
}

capture_circlize_geometry <- function(inputs, layout, scales, cfg) {
  for (chrom in inputs$chromosomes) {
    transformed <- circlize::circlize(c(0, inputs$lengths[[chrom]] / 2, inputs$lengths[[chrom]]),
      c(0.5, 0.5, 0.5), sector.index = chrom, track.index = 3)
    expected <- c(layout[chrom, "start_degrees"], layout[chrom, "end_degrees"])
    # Circular angles may differ by an integer number of turns.
    error <- (transformed[c(1, 3), "theta"] - expected + 180) %% 360 - 180
    if (any(abs(error) > 1e-7) ||
        abs(transformed[2, "theta"] - mean(transformed[c(1, 3), "theta"])) > 1e-7) {
      stopf("Chr%s 的扇区方向或真实物理坐标线性映射错误。", chrom)
    }
  }
  tracks <- c("Chromosome labels", "Coordinates", "Chromosomes", cfg$te_classes, "SNP", "SV", "TIP")
  geometry <- do.call(rbind, lapply(seq_along(tracks), function(i) data.frame(
    track_index = i, track = tracks[i],
    inner_radius = circlize::get.cell.meta.data("cell.bottom.radius", sector.index = "1", track.index = i),
    outer_radius = circlize::get.cell.meta.data("cell.top.radius", sector.index = "1", track.index = i))))
  axes <- do.call(rbind, lapply(seq_along(cfg$te_classes), function(i) {
    class <- cfg$te_classes[i]; index <- i + 3L; scale <- scales[[class]]
    ticks <- te_axis_ticks(scale)
    polar <- circlize::circlize(rep(0, length(ticks)), ticks, sector.index = "1", track.index = index)
    inner <- geometry$inner_radius[index]; outer <- geometry$outer_radius[index]
    normalized <- (polar[, "rou"] - inner) / (outer - inner)
    if (any(abs(normalized - ticks / scale$cap) > 1e-10)) stopf("独立 TE 纵轴位置错误。")
    data.frame(track = class, percentage = ticks * 100, coverage_fraction = ticks,
      track_cap = scale$cap, angle_degrees = polar[, "theta"], radius = polar[, "rou"],
      fraction_of_track_height = normalized)
  }))
  panels <- expand.grid(chromosome = inputs$chromosomes,
    track = c("Chromosomes", cfg$te_classes, "SNP", "SV", "TIP"), stringsAsFactors = FALSE)
  panels$start_bp <- 0
  panels$end_bp <- unname(inputs$lengths[panels$chromosome])
  panels$start_degrees <- layout[panels$chromosome, "start_degrees"]
  panels$end_degrees <- layout[panels$chromosome, "end_degrees"]
  unit_in <- abs(diff(grconvertX(c(0, 1), "user", "inches")))
  diameter <- 2 * unit_in * geometry$outer_radius[1]
  list(layout = layout, track_geometry = geometry, te_axis_positions = axes, panels = panels,
    circle_diameter_in = diameter, circle_diameter_canvas_fraction = diameter / grDevices::dev.size("in")[1])
}

draw_grid_center <- function(state, cfg) {
  radius <- tail(state$track_geometry$inner_radius, 1)
  center <- c(grconvertX(0, "user", "ndc"), grconvertY(0, "user", "ndc"))
  width <- abs(diff(grconvertX(c(-radius, radius), "user", "ndc")))
  height <- abs(diff(grconvertY(c(-radius, radius), "user", "ndc")))
  grid::pushViewport(grid::viewport(x = center[1], y = center[2], width = width, height = height))
  on.exit(grid::popViewport(), add = TRUE)
  grid::grid.circle(x = 0.5, y = 0.5, r = 0.47,
    gp = grid::gpar(fill = "#F1F3F8", col = NA))
  labels <- c("Chromosomes", "LINE density", "LTR density", "SINE density", "DNA density", "SNV density", "SV density", "TEV density")
  colors <- c(cfg$border_color, unname(cfg$colors[c(cfg$te_classes, "SNP", "SV", "TIP")]))
  label_gp <- grid::gpar(fontfamily = cfg$typography$family,
    fontsize = cfg$typography$legend_pt, fontface = "bold", col = cfg$text_color)
  direction_gp <- grid::gpar(fontfamily = cfg$typography$family,
    fontsize = 8.0, fontface = "bold", col = cfg$text_color)
  measure <- function(label, gp, dim = "width") {
    grob <- grid::textGrob(label, gp = gp)
    if (dim == "width") grid::convertWidth(grid::grobWidth(grob), "inches", valueOnly = TRUE)
    else grid::convertHeight(grid::grobHeight(grob), "inches", valueOnly = TRUE)
  }
  viewport_w <- grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  viewport_h <- grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  label_widths <- vapply(labels, measure, numeric(1), gp = label_gp)
  arrow_width <- max(vapply(c("Outer", "Inner"), measure, numeric(1), gp = direction_gp))
  line_width <- 0.14; arrow_gap <- 0.12; text_gap <- 0.10
  total_width <- arrow_width + arrow_gap + line_width + text_gap + max(label_widths)
  left <- (viewport_w - total_width) / 2
  arrow_x <- (left + arrow_width / 2) / viewport_w
  line_a <- (left + arrow_width + arrow_gap) / viewport_w
  line_b <- line_a + line_width / viewport_w
  text_x <- line_b + text_gap / viewport_w
  yy <- seq(0.82, 0.18, length.out = 8)
  grid::grid.segments(x0 = arrow_x, x1 = arrow_x, y0 = 0.32, y1 = 0.68,
    arrow = grid::arrow(ends = "both", type = "open", length = grid::unit(1.3, "mm")),
    gp = grid::gpar(col = cfg$text_color, lwd = 0.75))
  records <- list()
  record <- function(label, x, y, w, h, role, color) {
    data.frame(role = role, label = label,
      center_x_in = center[1] * cfg$geometry$width_in + (x - 0.5) * viewport_w,
      center_y_in = center[2] * cfg$geometry$height_in + (y - 0.5) * viewport_h,
      width_in = w, height_in = h, color = color,
      circle_center_x_in = center[1] * cfg$geometry$width_in,
      circle_center_y_in = center[2] * cfg$geometry$height_in)
  }
  for (i in 1:2) {
    label <- c("Outer", "Inner")[i]; y <- c(0.75, 0.25)[i]
    grid::grid.text(label, x = arrow_x, y = y, gp = direction_gp)
    records[[length(records)+1L]] <- record(label, arrow_x, y,
      measure(label, direction_gp), measure(label, direction_gp, "height"), "direction", cfg$text_color)
  }
  for (i in seq_along(labels)) {
    grid::grid.segments(x0 = line_a, x1 = line_b, y0 = yy[i], y1 = yy[i],
      gp = grid::gpar(col = colors[i], lwd = if (i == 1) 2.0 else 3.0, lineend = "round"))
    grid::grid.text(labels[i], x = text_x, y = yy[i], just = "left", gp = label_gp)
    records[[length(records)+1L]] <- record(labels[i], text_x + label_widths[i]/viewport_w/2,
      yy[i], label_widths[i], measure(labels[i], label_gp, "height"), "track_label", colors[i])
    records[[length(records)+1L]] <- record(labels[i], (line_a+line_b)/2,
      yy[i], line_width, 3/72, "color_line", colors[i])
  }
  d <- do.call(rbind, records)
  bbox <- c(min(d$center_x_in-d$width_in/2), max(d$center_x_in+d$width_in/2),
    min(d$center_y_in-d$height_in/2), max(d$center_y_in+d$height_in/2))
  if (abs(mean(bbox[1:2])-d$circle_center_x_in[1]) > 1e-8 ||
      abs(mean(bbox[3:4])-d$circle_center_y_in[1]) > 1e-8) stopf("中心图例整体未严格居中。")
  d
}

draw_grid_page <- function(scales, cfg) {
  # V15: no page title, window subtitle or bottom density legends.
  invisible(NULL)
}

draw_complete_figure <- function(inputs, layout, scales, cfg) {
  load_circlize(cfg)
  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  par(mar = c(0, 0, 0, 0), fig = cfg$geometry$figure_region,
    family = cfg$typography$family, ps = cfg$typography$base_pt, cex = 1,
    xaxs = "i", yaxs = "i", xpd = NA)
  g <- cfg$geometry
  gaps <- setNames(c(rep(g$gap_deg, 40), g$opening_deg), inputs$chromosomes)
  first <- g$opening_center_deg + if (g$clockwise) -g$opening_deg / 2 else g$opening_deg / 2
  circlize::circos.par(start.degree = first, gap.after = gaps,
    clock.wise = g$clockwise, xaxis.clock.wise = g$clockwise,
    track.margin = g$track_margin, cell.padding = g$cell_padding,
    canvas.xlim = g$xlim, canvas.ylim = g$ylim,
    unit.circle.segments = g$unit_circle_segments, points.overflow.warning = FALSE)
  limits <- cbind(rep(0, 41), inputs$sizes$length_bp)
  rownames(limits) <- inputs$chromosomes
  circlize::circos.initialize(factors = factor(inputs$chromosomes, levels = inputs$chromosomes),
    xlim = limits, sector.width = layout$displayed_sector_width)
  draw_chromosome_label_track_v15(inputs, layout, cfg)
  ticks <- draw_coordinate_track_v15(inputs, layout, cfg)
  markers <- draw_ideogram_v15(inputs, layout, cfg)
  profiles <- draw_te_tracks_v15(inputs, scales, cfg)
  draw_heat_tracks_v15(inputs, scales, cfg)
  axis_annotations <- draw_te_axis_v15(scales, cfg)
  state <- capture_circlize_geometry(inputs, layout, scales, cfg)
  state$annotations <- plan_outer_annotations(inputs, layout, ticks, cfg)
  state$centromere_markers <- markers
  state$te_profiles <- profiles
  state$te_axis_annotations <- axis_annotations[axis_annotations$kind == "TE_axis_tick", ]
  state$variant_axis_annotations <- axis_annotations[axis_annotations$kind != "TE_axis_tick", ]
  state$axis_gradient_strips <- attr(axis_annotations, "gradient_strips")
  all_boxes <- rbind(state$annotations$names,
    state$annotations$coordinates[state$annotations$coordinates$label_status == "shown",
      names(state$annotations$names), drop = FALSE],
    axis_annotations[, names(state$annotations$names), drop = FALSE])
  if (nrow(all_boxes) > 1) for (i in seq_len(nrow(all_boxes)-1)) for (j in (i+1):nrow(all_boxes)) {
    if (boxes_overlap(all_boxes[i, ], all_boxes[j, ], 0.30)) {
      stopf("外圈/TE 轴标签碰撞：%s / %s", all_boxes$label[i], all_boxes$label[j])
    }
  }
  draw_outer_annotations(state$annotations, cfg)
  state$center_legend <- draw_grid_center(state, cfg)
  state$ideogram_fills <- ideogram_fill_table(inputs$chromosomes, cfg)
  draw_grid_page(scales, cfg)
  state
}

open_output_device <- function(path, format, cfg, preview = FALSE) {
  g <- cfg$geometry
  if (!capabilities("cairo")) stopf("当前 R 未编译 Cairo 支持，无法无界面导出。")
  if (format == "svg") {
    grDevices::svg(path, width = g$width_in, height = g$height_in,
      pointsize = cfg$typography$base_pt, family = cfg$typography$family, bg = "white", onefile = FALSE)
  } else if (format == "pdf") {
    grDevices::cairo_pdf(path, width = g$width_in, height = g$height_in,
      pointsize = cfg$typography$base_pt, family = cfg$typography$family, bg = "white", onefile = FALSE)
  } else if (format == "png") {
    grDevices::png(path, width = g$width_in, height = g$height_in, units = "in",
      res = if (preview) cfg$preview_dpi else cfg$dpi, type = "cairo",
      pointsize = cfg$typography$base_pt, bg = "white")
  } else stopf("不支持的输出格式：%s", format)
}

compact_svg_clip_definitions <- function(path) {
  # Cairo may repeat the same compound clipping path for each paint operation.
  # Deduplicate identical XML bodies and redirect references; no coordinates,
  # paths, colors or text are changed. Base R only, no external dependencies.
  lines <- readLines(path, warn = FALSE)
  starts <- which(grepl("^<clipPath id=", lines))
  ends <- which(grepl("^</clipPath>", lines))
  if (!length(starts)) return(invisible(NULL))
  if (length(starts) != length(ends) || any(ends <= starts)) stopf("SVG 裁剪定义结构异常。")
  ids <- sub('^<clipPath id="([^"]+)".*', '\\1', lines[starts])
  bodies <- vapply(seq_along(starts), function(i)
    paste(lines[seq.int(starts[i]+1L, ends[i]-1L)], collapse="\n"), character(1))
  first <- match(bodies, bodies)
  aliases <- setNames(ids[first], ids)
  remove <- unlist(lapply(which(first != seq_along(first)), function(i)
    seq.int(starts[i], ends[i])), use.names=FALSE)
  if (length(remove)) lines <- lines[-remove]
  matches <- gregexpr("(?<=url\\(#)clip[0-9]+(?=\\))", lines, perl=TRUE)
  old <- regmatches(lines, matches)
  regmatches(lines, matches) <- lapply(old, function(x) unname(aliases[x]))
  writeLines(lines, path, useBytes=TRUE)
  invisible(NULL)
}

save_figure <- function(path, format, inputs, layout, scales, cfg, preview = FALSE, recorded_plot = NULL) {
  open_output_device(path, format, cfg, preview)
  on.exit({
    grDevices::dev.off()
    if (identical(format, "svg")) compact_svg_clip_definitions(path)
  }, add = TRUE)
  if (!is.null(recorded_plot)) {
    grDevices::replayPlot(recorded_plot)
    return(invisible(NULL))
  }
  # 几何只计算一次；同尺寸 PDF/PNG/预览复用显示列表，保证版式一致并减少等待。
  grDevices::dev.control(displaylist = "enable")
  state <- draw_complete_figure(inputs, layout, scales, cfg)
  state$recorded_plot <- grDevices::recordPlot()
  state
}

write_te_profile_qc <- function(path, inputs, scales, profiles, cfg) {
  grDevices::cairo_pdf(path, width = 10, height = 8.0,
    family = cfg$typography$family, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(4, 1), mar = c(2.1, 4.2, 2.0, 0.8), mgp = c(2.3, 0.6, 0),
    las = 1, family = cfg$typography$family, cex = 0.75)
  offsets <- setNames(c(0, head(cumsum(inputs$sizes$length_bp), -1)), inputs$chromosomes)
  for (class in cfg$te_classes) {
    d <- profiles$nodes[profiles$nodes$TE_class == class, ]
    xx <- (offsets[d$chromosome] + d$midpoint_bp) / 1e6
    plot(xx, d$raw_fraction * 100, pch = 16, cex = 0.22,
      col = grDevices::adjustcolor(cfg$colors[[class]], alpha.f = 0.35),
      xlab = "", ylab = "coverage (%)", ylim = c(0, max(d$raw_fraction) * 105), bty = "l",
      main = paste0(class, ": raw windows (points), display smoothing (line); independent main-plot cap ",
        format_value(scales[[class]]$cap * 100), "%"))
    for (chrom in inputs$chromosomes) {
      z <- profiles$curves[profiles$curves$TE_class == class & profiles$curves$chromosome == chrom, ]
      lines((offsets[chrom] + z$position_bp) / 1e6, z$display_fraction * 100,
        col = cfg$colors[[class]], lwd = 0.7)
    }
    abline(h = scales[[class]]$cap * 100, lty = 2, col = "#59656F", lwd = 0.6)
  }
}

write_plot_records <- function(prefix, inputs, layout, scales, cfg, output_paths) {
  scale_rows <- lapply(names(scales), function(name) {
    s <- scales[[name]]
    raw <- inputs$data[[s$column]]
    valid <- raw[is.finite(raw)]
    count_column <- sub("_per_Mb$", "_count", s$column)
    total <- if (grepl("_per_Mb$", s$column)) sum(inputs$data[[count_column]], na.rm = TRUE) else NA_real_
    qv <- function(prob) if (length(valid)) unname(quantile(valid, prob, names = FALSE, type = 7)) else NA_real_
    data.frame(track = name, input_column = s$column, transform = s$transform,
                 raw_upper_limit = s$cap, displayed_upper_limit = s$cap / s$display_divisor,
                 displayed_unit = s$unit, windows = length(raw), missing_windows = sum(is.na(raw)),
                 clipped_windows = sum(raw > s$cap, na.rm = TRUE),
                 minimum_raw = if (length(valid)) min(valid) else NA_real_,
                 mean_raw = if (length(valid)) mean(valid) else NA_real_,
                 median_raw = qv(0.50), q90_raw = qv(0.90), q95_raw = qv(0.95),
                 q99_raw = qv(0.99), maximum_raw = if (length(valid)) max(valid) else NA_real_,
                 class_specific_q99 = if (!is.null(s$class_q99)) s$class_q99 else NA_real_,
                 selected_display_quantile = if (!is.null(s$selected_quantile)) s$selected_quantile else NA_real_,
                 selected_quantile_raw = if (!is.null(s$quantile_value)) s$quantile_value else NA_real_,
                 independent_site_total = total,
                 scale_basis = if (!is.null(s$scale_basis)) s$scale_basis else "unchanged V10 variant density scale",
                 TE_display_smoothing = if (name %in% cfg$te_classes) cfg$te_smoothing$method else "none")
  })
  scales_path <- paste0(prefix, ".track_scales.tsv")
  layout_path <- paste0(prefix, ".chromosome_display_mapping.tsv")
  input_path <- paste0(prefix, ".input_files.tsv")
  settings_path <- paste0(prefix, ".settings.R")
  session_path <- paste0(prefix, ".sessionInfo.txt")
  write_tsv(do.call(rbind, scale_rows), scales_path)
  write_tsv(layout, layout_path)
  source_files <- c(windows = cfg$table, chromosome_sizes = cfg$chrom_sizes)
  if (!identical(cfg$centromeres, "none")) source_files <- c(source_files, centromeres = cfg$centromeres)
  source_files <- c(source_files, R_code = file.path(cfg$project_dir, "06_scripts", "plot_chicken_circos_v15.R"),
                   shell_entrypoint = file.path(cfg$project_dir, "06_scripts", "run_circos_R_v15.sh"))
  sexchr_dir <- file.path(cfg$project_dir, "03_variants", "v11", "sexchr_snps")
  sexchr_sources <- file.path(sexchr_dir, c("assembly51.ZW.discovery.biallelic.sites.vcf.gz",
    "assembly51.ZW.callable_union.bed.gz", "method.json"))
  if (grepl("/v11/", cfg$table, fixed = TRUE) && all(file.exists(sexchr_sources))) {
    source_files <- c(source_files, setNames(sexchr_sources,
      c("ZW_SNP_site_catalog", "ZW_callable_union", "ZW_extraction_method")))
  }
  write_tsv(data.frame(role = names(source_files), path = unname(normalizePath(source_files)),
                        bytes = unname(file.info(source_files)$size),
                        md5 = unname(tools::md5sum(source_files))), input_path)
  dput(list(reference = "GRCg7b", coordinates = "0-based half-open", configuration = cfg,
             centers_drawn = if (is.null(inputs$centers)) 0L else nrow(inputs$centers),
             chromosome_arcs_proportional = all(abs(layout$enlargement_factor - 1) < 1e-7),
             microchromosomes_enlarged = any(layout$enlargement_factor > 1 + 1e-7),
             microchromosome_set = cfg$geometry$micro_chromosomes,
             microchromosome_min_angle_deg = cfg$geometry$min_micro_deg,
             within_chromosome_coordinates_linear = TRUE,
             centromere_marker = "one red (#E75A46), full ideogram height rectangular sectors; evidence retained in QC; illustrative angular length",
             cytobands_displayed = FALSE,
             cytoband_reason = "No validated GRCg7b assembly-coordinate G-band annotation was supplied; synthetic bands were not fabricated.",
             palette = cfg$palette,
             palette_source_repository = "https://github.com/1577377232/Genomics-Figure-Styles",
             palette_source_commit = "13f1fbb7a8d67da1c335934ce8ca839a5c6be6de",
             TE_scale = "independent class scales: LTR q95; LINE/SINE/DNA q99; rounded up to two significant digits; coverage percent",
             center_legend = "measured centered legend; Outer/Inner 8pt; shortened centered double arrow; unchanged background and track colors",
             TE_smoothing = cfg$te_smoothing,
             heat_color_mapping = "monotone Lab ramp: white, light palette tint, palette base, dark palette shade; unchanged density transforms",
             plotting_engine = "circlize + grid",
             circlize_version = as.character(utils::packageVersion("circlize")),
             graphics_tracks = c("Chromosome labels", "Coordinates", "Chromosomes", cfg$te_classes, "SNP", "SV", "TIP"),
             TE_structures_displayed = FALSE,
             reference_TE_tracks = "merged within-class covered bp / actual window bp",
             variant_tracks = "independent loci / actual window bp * 1000000",
             SNP_ZW_supplement = "V11 default table adds real Z/W assembly-discovery SNPs from 51 existing BAMs under original 50kb callable filters; no sex-stratified AF or joint genotypes",
             plotting_code_origin = "V15 refinement of V13; fixed V12 track colors, uniform narrow ideograms, 48 degree opening, left-aligned adjacent TE labels, opening-side variants/MB axes and gradients, no page annotations; retained annular TE clipping and arc baselines; all V11 data and scales unchanged"), settings_path)
  writeLines(capture.output(sessionInfo()), session_path)
  paths <- c(output_paths, scales_path, layout_path, input_path, settings_path, session_path)
  write_tsv(data.frame(path = basename(paths), bytes = unname(file.info(paths)$size),
                        md5 = unname(tools::md5sum(paths))), paste0(prefix, ".output_manifest.tsv"))
}


# ------------------------------------------------------------------------------
# 10. 主执行函数与命令行参数
# ------------------------------------------------------------------------------

validate_configuration <- function(cfg) {
  if (!identical(cfg$palette, "pangenome_spectrum")) {
    stopf("V15 --palette 固定为 pangenome_spectrum。")
  }
  if (length(cfg$window_bp) != 1 || !is.finite(cfg$window_bp) ||
      cfg$window_bp <= 0 || cfg$window_bp != floor(cfg$window_bp)) stopf("--window 必须为正整数 bp。")
  if (length(cfg$dpi) != 1 || !is.finite(cfg$dpi) || cfg$dpi <= 0) stopf("--dpi 必须为正数。")
  if (!length(cfg$formats) || any(!cfg$formats %in% c("svg", "pdf", "png"))) {
    stopf("--formats 只能包含 svg,pdf,png。")
  }
  if (!cfg$scales$snp_display_divisor %in% c(1, 1000)) stopf("SNP 显示单位的除数应为 1 或 1000。")
  invisible(track_columns(cfg))
}

canonical_output_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(getwd(), path)
  missing <- character(0)
  while (!file.exists(path)) {
    parent <- dirname(path)
    if (identical(parent, path)) break
    missing <- c(basename(path), missing)
    path <- parent
  }
  path <- normalizePath(path, mustWork = TRUE)
  if (length(missing)) path <- do.call(file.path, as.list(c(path, missing)))
  path
}

render_circos <- function(cfg = default_config()) {
  cfg <- apply_palette(cfg)
  validate_configuration(cfg)
  cfg <- resolve_paths(cfg)
  candidate <- canonical_output_path(cfg$out_dir)
  old_roots <- file.path(cfg$project_dir, "07_figures", paste0("R_v", 3:14))
  if (any(vapply(old_roots, function(p) identical(candidate, p) || startsWith(candidate, paste0(p, "/")), logical(1)))) {
    stopf("V15 不能写入已有 R_v3/R_v4/R_v5/R_v6/R_v7/R_v8/R_v9/R_v10/R_v11/R_v12/R_v13 目录，请使用 R_v15 或新的输出目录。")
  }
  load_circlize(cfg)
  inputs <- read_and_validate_inputs(cfg)
  layout <- make_layout(inputs, cfg)
  scales <- make_scales(inputs, cfg)
  message(sprintf("输入检查通过：%s bp，%s 个真实窗口，41 条染色体；SNP=%s，SV=%s，TIP=%s。",
                    cfg$window_bp, nrow(inputs$data), cfg$snp_panel, cfg$sv_set, cfg$tip_set))
  if (cfg$validate_only) return(invisible(list(inputs = inputs, layout = layout, scales = scales)))
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(cfg$out_dir, sprintf("chicken_GRCg7b_R_v15_%skb_%s_TIP-%s_SV-%s_palette-%s",
                                          cfg$window_bp / 1000, cfg$snp_panel, cfg$tip_set, cfg$sv_set, cfg$palette))
  output_paths <- character(0)
  state <- NULL
  for (format in unique(cfg$formats)) {
    path <- paste0(prefix, ".", format)
    if (is.null(state)) {
      state <- save_figure(path, format, inputs, layout, scales, cfg)
    } else {
      save_figure(path, format, inputs, layout, scales, cfg, recorded_plot = state$recorded_plot)
    }
    output_paths <- c(output_paths, path)
    message("已输出：", path)
  }
  preview_path <- paste0(prefix, ".preview.png")
  save_figure(preview_path, "png", inputs, layout, scales, cfg, preview = TRUE, recorded_plot = state$recorded_plot)
  layout <- state$layout
  geometry_path <- paste0(prefix, ".track_geometry.tsv")
  write_tsv(state$track_geometry, geometry_path)
  annotation_path <- paste0(prefix, ".chromosome_labels.tsv")
  coordinate_path <- paste0(prefix, ".coordinate_labels.tsv")
  coordinate_ticks_path <- paste0(prefix, ".coordinate_ticks.tsv")
  figure_path <- paste0(prefix, ".figure_geometry.tsv")
  axis_path <- paste0(prefix, ".TE_axis_positions.tsv")
  write_tsv(state$annotations$names, annotation_path)
  write_tsv(state$annotations$coordinates, coordinate_path)
  write_tsv(state$annotations$ticks, coordinate_ticks_path)
  write_tsv(state$te_axis_positions, axis_path)
  write_tsv(state$te_axis_annotations, paste0(prefix, ".TE_axis_labels.tsv"))
  write_tsv(state$variant_axis_annotations, paste0(prefix, ".variant_axis_labels.tsv"))
  write_tsv(state$axis_gradient_strips, paste0(prefix, ".variant_axis_gradients.tsv"))
  write_tsv(state$center_legend, paste0(prefix, ".center_legend.tsv"))
  write_tsv(state$ideogram_fills, paste0(prefix, ".ideogram_fills.tsv"))
  write_tsv(data.frame(track = names(cfg$colors), color = unname(cfg$colors),
    use = c(rep("reference TE coverage area, outline, center legend", 4),
      rep("variant density ramp base, colorbar, center legend", 3))), paste0(prefix, ".color_mapping.tsv"))
  write_tsv(state$panels, paste0(prefix, ".sector_panels.tsv"))
  write_tsv(state$centromere_markers, paste0(prefix, ".centromere_markers.tsv"))
  write_tsv(state$te_profiles$nodes, paste0(prefix, ".TE_profile_nodes.tsv"))
  write_tsv(state$te_profiles$curves, paste0(prefix, ".TE_profile_curves.tsv"))
  write_tsv(state$te_profiles$render_bounds, paste0(prefix, ".TE_render_bounds.tsv"))
  color_rows <- lapply(c("SNP", "SV", "TIP"), function(name) {
    f <- seq(0, 1, length.out = 256); scale <- scales[[name]]
    raw <- if (scale$transform == "log1p") expm1(f * log1p(scale$cap)) else f * scale$cap
    data.frame(track = name, scale_fraction = f, raw_density_per_Mb = raw,
      displayed_density = raw / scale$display_divisor, color = density_colors(f, cfg$colors[[name]]))
  })
  write_tsv(do.call(rbind, color_rows), paste0(prefix, ".heat_color_scale.tsv"))
  write_tsv(data.frame(width_in = cfg$geometry$width_in, height_in = cfg$geometry$height_in,
    circle_diameter_in = state$circle_diameter_in,
    circle_diameter_canvas_fraction = state$circle_diameter_canvas_fraction,
    opening_degrees = cfg$geometry$opening_deg, gap_degrees = cfg$geometry$gap_deg,
    chromosome_names = nrow(state$annotations$names),
    coordinate_labels = sum(state$annotations$coordinates$label_status == "shown"),
    opening_center_degrees = cfg$geometry$opening_center_deg, clockwise = cfg$geometry$clockwise,
    annotation_collisions = state$annotations$collision_pairs,
    separate_coordinate_track = TRUE, separate_chromosome_label_track = TRUE,
    sector_track_frames = TRUE, coordinate_track_frames = FALSE,
    TE_independent_axes = TRUE, centromere_rendering = "short red rectangle"), figure_path)
  output_paths <- c(output_paths, preview_path, geometry_path, annotation_path,
    coordinate_path, coordinate_ticks_path, figure_path, axis_path,
    paste0(prefix, ".sector_panels.tsv"), paste0(prefix, ".centromere_markers.tsv"),
    paste0(prefix, ".TE_profile_nodes.tsv"), paste0(prefix, ".TE_profile_curves.tsv"),
    paste0(prefix, ".TE_render_bounds.tsv"),
    paste0(prefix, ".TE_axis_labels.tsv"),
    paste0(prefix, ".variant_axis_labels.tsv"), paste0(prefix, ".variant_axis_gradients.tsv"),
    paste0(prefix, ".center_legend.tsv"), paste0(prefix, ".ideogram_fills.tsv"),
    paste0(prefix, ".color_mapping.tsv"),
    paste0(prefix, ".heat_color_scale.tsv"))
  if (cfg$window_bp == 1000000) {
    profile_path <- file.path(cfg$out_dir, "TE_smoothing_QC.pdf")
    write_te_profile_qc(profile_path, inputs, scales, state$te_profiles, cfg)
    output_paths <- c(output_paths, profile_path)
  }
  write_plot_records(prefix, inputs, layout, scales, cfg, output_paths)
  # 简短文件名作为当前最后一次绘图的索引，另有每张图独立的完整尺度/映射记录。
  file.copy(paste0(prefix, ".track_scales.tsv"), file.path(cfg$out_dir, "track_scales.tsv"), overwrite = TRUE)
  file.copy(paste0(prefix, ".chromosome_display_mapping.tsv"), file.path(cfg$out_dir, "chromosome_display_mapping.tsv"), overwrite = TRUE)
  invisible(list(prefix = prefix, layout = layout, scales = scales, configuration = cfg))
}

print_usage <- function() {
  cat(paste(c(
    "V15 GRCg7b：circlize + grid，48 度右侧开口；四类 TE 独立百分比纵轴，变异密度单位 /MB。",
    "用法：Rscript plot_chicken_circos_v15.R [选项]",
    "  --project-dir PATH       chicken_7b_TE_circos 项目目录；默认从脚本位置识别",
    "  --window BP              窗口大小，默认 1000000；已有 100000/500000/1000000",
    "  --snp-panel NAME         population95 / discovery50 / adaptive_discovery（默认）",
    "  --tip-set NAME           strict / extended（默认）",
    "  --sv-set NAME            all（默认）/ nonTE_assessed / nonrepeat_assessed",
    "  --opening DEG            右侧开口角度，默认 48",
    "  --min-micro-deg DEG      GGA10-39 最小角度，默认 2.1；W 至少 3 度",
    "  --gap DEG                染色体间隔，默认 1.1",
    "  --palette NAME           pangenome_spectrum（V15 固定映射）",
    "  --dpi INT                PNG 分辨率，默认 600",
    "  --formats LIST           默认 svg,pdf,png",
    "  --out-dir PATH           默认 07_figures/R_v15；禁止写入旧版本目录",
    "  --table PATH             自定义同样字段的窗口表",
    "  --chrom-sizes PATH       自定义主染色体长度文件",
    "  --centromeres PATH|none  自定义近似着丝粒表，或 none 不显示",
    "  --validate-only          仅校验真实输入，不画图",
    "  --all-windows            依次绘制 100 kb、500 kb、1 Mb",
    "  --help                   显示帮助"
  ), collapse = "\n"), "\n")
}

parse_arguments <- function(args, cfg) {
  character_options <- c("project-dir" = "project_dir", "snp-panel" = "snp_panel",
                          "tip-set" = "tip_set", "sv-set" = "sv_set", "palette" = "palette", "out-dir" = "out_dir",
                          "table" = "table", "chrom-sizes" = "chrom_sizes", "centromeres" = "centromeres")
  numeric_options <- c("window" = "window_bp", "dpi" = "dpi")
  geometry_options <- c("opening" = "opening_deg", "min-micro-deg" = "min_micro_deg", "gap" = "gap_deg")
  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--help") { print_usage(); return(NULL) }
    if (arg == "--validate-only") { cfg$validate_only <- TRUE; i <- i + 1L; next }
    if (arg == "--all-windows") { cfg$all_windows <- TRUE; i <- i + 1L; next }
    if (!grepl("^--", arg)) stopf("无法识别参数：%s", arg)
    key <- sub("^--", "", arg)
    allowed <- c(names(character_options), names(numeric_options), names(geometry_options), "formats")
    if (!key %in% allowed) stopf("未知选项：%s；使用 --help 查看说明。", arg)
    if (i == length(args) || grepl("^--", args[i + 1L])) stopf("选项 %s 缺少值。", arg)
    value <- args[i + 1L]
    if (key %in% names(character_options)) cfg[[character_options[[key]]]] <- value
    if (key %in% names(numeric_options)) cfg[[numeric_options[[key]]]] <- suppressWarnings(as.numeric(value))
    if (key %in% names(geometry_options)) cfg$geometry[[geometry_options[[key]]]] <- suppressWarnings(as.numeric(value))
    if (key == "formats") cfg$formats <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
    i <- i + 2L
  }
  cfg
}

main <- function() {
  cfg <- parse_arguments(commandArgs(trailingOnly = TRUE), default_config())
  if (is.null(cfg)) return(invisible(NULL))
  if (cfg$all_windows && !is.null(cfg$table)) stopf("--all-windows 不能与 --table 同时使用。")
  windows <- if (cfg$all_windows) c(100000, 500000, 1000000) else cfg$window_bp
  for (window in windows) {
    one <- cfg
    one$window_bp <- window
    render_circos(one)
  }
}

if (sys.nframe() == 0L) {
  tryCatch(main(), error = function(e) {
    message("错误：", conditionMessage(e))
    quit(status = 1L)
  })
}
