suppressPackageStartupMessages(library(HiCcompare))
suppressPackageStartupMessages(library(strawr))

options(scipen = 999)
options(warn = -1)
options(help_type = "text")

data_dir <- Sys.getenv("LAB_DATA_DIR", path.expand("~/lab-data"))

Args <- commandArgs(trailingOnly = T)

analysis_type <- Args[1]
flag_inter <- as.logical(Args[2])

#################################################################
##                            Intra                            ##
#################################################################
compute_initial_layout <- function(A, xspan = 100, yspan = 100, pad_x = 0, pad_y = 0) {
  d <- min(nrow(A), ncol(A))
  n <- ncol(A)

  X <- xspan - 2 * pad_x
  Y <- yspan - 2 * pad_y

  # Max w allowed by horizontal and vertical constraints
  w_x <- X / (d + 1)
  w_y <- 2 * Y / (n + 1)
  w <- min(w_x, w_y)

  # Actual footprint size with that w
  W <- (d + 1) * w
  H <- ((n + 1) * w) / 2

  # Centering
  left <- (xspan - W) / 2 + pad_x
  base <- (yspan - H) / 2 + pad_y

  start_x <- 2 * (left + w)
  start_y <- 2 * base

  # the y of the triangle's top
  top <- start_y / 2 + ((n + 1) * w) / 2

  list(
    w = w, start_x = start_x, start_y = start_y, top = top,
    bbox = list(left = left, right = left + W, bottom = base, top = base + H)
  )
}

compute_final_layout <- function(A_focus,
                                 w,
                                 n_custom,
                                 xspan = 100,
                                 yspan = 100,
                                 pad_x = 0,
                                 legend_line_units = 1.8, # height per legend row
                                 legend_gap_units = 2.0, # extra gap under the legend
                                 min_top_y = 40, # don't push plot below this
                                 title_reserve = 6 # don't push plot above this
  ) {
  # Basic dimensions of the focus matrix
  d_focus <- min(nrow(A_focus), ncol(A_focus))
  n_focus <- ncol(A_focus)

  # Footprint of the contact map in w units
  W <- (d_focus + 1) * w
  H <- ((n_focus + 1) * w) / 2

  # Center the footprint horizontally
  left <- (xspan - W) / 2 + pad_x
  start_x <- 2 * (left + w)

  # Legend height and vertical placement
  legend_height <- n_custom * legend_line_units + legend_gap_units
  max_top_y <- yspan - title_reserve # e.g. 100 - 12 = 88
  Y_top_target <- yspan - legend_height
  Y_top_target <- min(Y_top_target, max_top_y) # don't go above title area
  Y_top_target <- max(Y_top_target, min_top_y) # don't go below 40

  start_y <- 2 * (Y_top_target - H)
  legend_mid_y <- Y_top_target + 5

  list(
    w            = w,
    start_x      = start_x,
    start_y      = start_y,
    legend_mid_y = legend_mid_y
  )
}

zero_diag <- function(matrix, width) {
  if (nrow(matrix) != ncol(matrix)) {
    stop("zero_diag: matrix must be square")
  }
  if (nrow(matrix) <= width) {
    stop(sprintf(
      "zero_diag: matrix is too small (%d x %d) for requested diagonal width (%d). Try increase the range.",
      nrow(matrix), ncol(matrix), width
    ))
  }

  n <- nrow(matrix)
  for (j in 0:width) {
    for (i in 1:(n - j)) {
      matrix[i + j, i] <- 0
      matrix[i, i + j] <- 0
    }
  }
  return(matrix)
}

draw_scale <- function(
  A, C, w, start_x, start_y, breakList,
  analysis_type = "single_sample",
  gap_in = 0, # start offset along the slope
  chip_len_in = 0.09, # length of each color chip
  thickness_in = 0.06, # bar width
  inset_in = -0.20, # distance from the edge
  label_offset_in = -0.05, # offset for labels
  label_cex = 0.30, label_col = "#888888"
  ) {
  # center and diamond vertices
  get_center <- function(i, j, w, start_x, start_y) {
    c(
      start_x + (j - 1) * w + (i - 1) * w,
      start_y + (i - 1) * (-w) + j * w
    )
  }
  get_vertices <- function(center, w) {
    list(
      left   = c((center[1] - w) / 2, center[2] / 2),
      top    = c(center[1] / 2, (center[2] + w) / 2),
      right  = c((center[1] + w) / 2, center[2] / 2),
      bottom = c(center[1] / 2, (center[2] - w) / 2)
    )
  }

  # Decide which breaks to show
  break_steps <- if (analysis_type == "single_sample") {
    seq(1, length(breakList), by = 5)
  } else {
    seq(1, length(breakList), by = 10)
  }
  breaks_to_show <- breakList[break_steps]

  # Check flag for expanded viewport
  use_expand_viewport <- FALSE
  if (exists("flag_expand_viewport", inherits = TRUE)) {
    use_expand_viewport <- isTRUE(get("flag_expand_viewport", inherits = TRUE))
  }

  old_xpd <- par("xpd")
  on.exit(par(xpd = old_xpd), add = TRUE)
  par(xpd = NA)

  if (!use_expand_viewport) {

    ## Rotated scale bar along plot edge
    n_cols <- ncol(A)
    d <- min(nrow(A), ncol(A))

    # Triangle corners
    TOP_v <- get_vertices(get_center(1, n_cols, w, start_x, start_y), w)$top
    BR_v <- get_vertices(get_center(d, d, w, start_x, start_y), w)$right
    BL_v <- get_vertices(get_center(1, 1, w, start_x, start_y), w)$left

    # Direction along the right edge (TOP -> BR)
    u <- BR_v - TOP_v # vector from TOP corner to BR corner
    u <- u / sqrt(sum(u^2)) # parallel along plot edge
    v <- c(-u[2], u[1]) # perpendicular dir
    if (sum(((TOP_v + BR_v + BL_v) / 3 - (TOP_v + BR_v) / 2) * v) < 0) v <- -v # point inward

    # Convert inches -> plot units
    ux_per_in <- diff(grconvertX(c(0, 1), from = "in", to = "user"))
    uy_per_in <- diff(grconvertY(c(0, 1), from = "in", to = "user"))
    scale_u <- function(L_in) L_in / sqrt((u[1] / ux_per_in)^2 + (u[2] / uy_per_in)^2)
    scale_v <- function(L_in) L_in / sqrt((v[1] / ux_per_in)^2 + (v[2] / uy_per_in)^2)

    step_u <- scale_u(chip_len_in)
    gap_u <- scale_u(gap_in)
    inset_v <- scale_v(inset_in)
    thick_v <- scale_v(thickness_in)
    lab_v <- scale_v(thickness_in + label_offset_in - 0.05)

    start_pt <- TOP_v + gap_u * u + inset_v * v

    # Draw chips
    for (k in seq_along(breaks_to_show)) {
      b <- breaks_to_show[k]
      color <- rownames(C)[order(abs(C[, "breaks"] - b), decreasing = FALSE)[1]]
      p0 <- start_pt + (k - 1) * step_u * u
      p1 <- start_pt + k * step_u * u
      poly <- rbind(p0, p1, p1 + thick_v * v, p0 + thick_v * v)
      polygon(poly[, 1], poly[, 2], col = color, border = NA)
    }

    # Labels at ends, rotated to the slope
    ang <- atan2(u[2], u[1]) * 180 / pi
    end_pt <- start_pt + length(breaks_to_show) * step_u * u
    start_edge <- start_pt
    end_edge <- end_pt

    # Offset perpendicular
    lab_off <- lab_v * v

    # start label at left edge, left-justified
    text(start_edge[1] + lab_off[1], start_edge[2] + lab_off[2],
      labels = format(round(breakList[1], 4)),
      cex = label_cex, col = label_col, srt = ang, adj = c(0, 0.5)
    )

    # end label at right edge, right-justified
    text(end_edge[1] + lab_off[1], end_edge[2] + lab_off[2],
      labels = format(round(tail(breakList, 1), 4), nsmall = 2),
      cex = label_cex, col = label_col, srt = ang, adj = c(1, 0.5)
    )
  } else {

    ## Top-center scale bar for expanded viewport 
    # Convert inches to user units for a horizontal bar
    ux_per_in <- diff(grconvertX(c(0, 1), from = "in", to = "user"))
    uy_per_in <- diff(grconvertY(c(0, 1), from = "in", to = "user"))

    bar_width_in <- chip_len_in * length(breaks_to_show)
    bar_width <- bar_width_in * ux_per_in
    bar_height <- thickness_in * uy_per_in
    lab_offset <- abs(label_offset_in) * uy_per_in

    x_mid <- grconvertX(0.895, from = "ndc", to = "user")
    y_mid <- start_y - 2.5

    x0 <- x_mid - bar_width / 2
    x1 <- x_mid + bar_width / 2
    y0 <- y_mid - bar_height / 2
    y1 <- y_mid + bar_height / 2

    step_x <- bar_width / length(breaks_to_show)

    # Draw chips
    for (k in seq_along(breaks_to_show)) {
      b <- breaks_to_show[k]
      color <- rownames(C)[order(abs(C[, "breaks"] - b), decreasing = FALSE)[1]]

      cx0 <- x0 + (k - 1) * step_x
      cx1 <- x0 + k * step_x

      rect(cx0, y0, cx1, y1, col = color, border = NA)
    }

    # Labels at ends, below the bar
    lab_y <- y0 - lab_offset

    text(x0, lab_y,
      labels = format(round(breakList[1], 4)),
      adj    = c(0, 1),
      cex    = label_cex,
      col    = label_col
    )

    text(x1, lab_y,
      labels = format(round(tail(breakList, 1), 4), nsmall = 2),
      adj    = c(1, 1),
      cex    = label_cex,
      col    = label_col
    )
  }
}

draw_scale_inh <- function(
  A, C, w, start_x, start_y, breakList,
  analysis_type = "single_sample",
  gap_in = 0,
  chip_len_in = 0.04,
  thickness_in = 0.04,
  inset_in = -0.20,
  label_offset_in = -0.05,
  label_cex = 0.30,
  label_col = "#888888",
  steps = 1
  ) {
  get_center <- function(i, j, w, start_x, start_y) {
    c(
      start_x + (j - 1) * w + (i - 1) * w,
      start_y + (i - 1) * (-w) + j * w
    )
  }
  get_vertices <- function(center, w) {
    list(
      left   = c((center[1] - w) / 2, center[2] / 2),
      top    = c(center[1] / 2, (center[2] + w) / 2),
      right  = c((center[1] + w) / 2, center[2] / 2),
      bottom = c(center[1] / 2, (center[2] - w) / 2)
    )
  }

  idx <- seq(1, length(breakList), by = steps)
  breaks_to_show <- breakList[idx]

  # expanded viewport flag
  use_expand_viewport <- FALSE
  if (exists("flag_expand_viewport", inherits = TRUE)) {
    use_expand_viewport <- isTRUE(get("flag_expand_viewport", inherits = TRUE))
  }

  old_xpd <- par("xpd")
  on.exit(par(xpd = old_xpd), add = TRUE)
  par(xpd = NA)

  if (!use_expand_viewport) {

    ## Rotated scale bar along plot edge
    n_cols <- ncol(A)
    d <- min(nrow(A), ncol(A))

    TOP_v <- get_vertices(get_center(1, n_cols, w, start_x, start_y), w)$top
    BR_v <- get_vertices(get_center(d, d, w, start_x, start_y), w)$right
    BL_v <- get_vertices(get_center(1, 1, w, start_x, start_y), w)$left

    # Direction along TOP->BR and point inward
    u <- BR_v - TOP_v
    u <- u / sqrt(sum(u^2))
    v <- c(-u[2], u[1])
    if (sum(((TOP_v + BR_v + BL_v) / 3 - (TOP_v + BR_v) / 2) * v) < 0) v <- -v

    # inches -> user units
    ux_per_in <- diff(grconvertX(c(0, 1), from = "in", to = "user"))
    uy_per_in <- diff(grconvertY(c(0, 1), from = "in", to = "user"))

    scale_u <- function(L_in) L_in / sqrt((u[1] / ux_per_in)^2 + (u[2] / uy_per_in)^2)
    scale_v <- function(L_in) L_in / sqrt((v[1] / ux_per_in)^2 + (v[2] / uy_per_in)^2)

    step_u <- scale_u(chip_len_in)
    gap_u <- scale_u(gap_in)
    inset_v <- scale_v(inset_in)
    thick_v <- scale_v(thickness_in)
    lab_v <- scale_v(thickness_in + label_offset_in - 0.05)

    start_pt <- TOP_v + gap_u * u + inset_v * v

    # Draw chips
    for (k in seq_along(breaks_to_show)) {
      b <- breaks_to_show[k]
      color <- rownames(C)[order(abs(C[, "breaks"] - b), decreasing = FALSE)[1]]
      p0 <- start_pt + (k - 1) * step_u * u
      p1 <- start_pt + k * step_u * u
      poly <- rbind(p0, p1, p1 + thick_v * v, p0 + thick_v * v)
      polygon(poly[, 1], poly[, 2], col = color, border = NA)
    }

    # Labels at ends
    ang <- atan2(u[2], u[1]) * 180 / pi
    end_pt <- start_pt + length(breaks_to_show) * step_u * u
    start_edge <- start_pt
    end_edge <- end_pt

    lab_off <- lab_v * v

    text(start_edge[1] + lab_off[1], start_edge[2] + lab_off[2],
      labels = format(round(breakList[1], 4)),
      cex = label_cex, col = label_col, srt = ang, adj = c(0, 0.5)
    )

    text(end_edge[1] + lab_off[1], end_edge[2] + lab_off[2],
      labels = format(round(tail(breakList, 1), 4)),
      cex = label_cex, col = label_col, srt = ang, adj = c(1, 0.5)
    )

    mid_pt <- (start_edge + end_edge) / 2
    text(mid_pt[1] + lab_off[1], mid_pt[2] + lab_off[2],
      labels = "1",
      cex = label_cex, col = label_col, srt = ang, adj = c(0.5, 0.5)
    )
  } else {

    ## Top-center scale bar for expanded viewport
    ux_per_in <- diff(grconvertX(c(0, 1), from = "in", to = "user"))
    uy_per_in <- diff(grconvertY(c(0, 1), from = "in", to = "user"))

    chip_width_user <- chip_len_in * ux_per_in
    bar_width <- chip_width_user * length(breaks_to_show)
    bar_height <- thickness_in * uy_per_in
    lab_offset <- abs(label_offset_in) * uy_per_in

    x_mid <- grconvertX(0.895, from = "ndc", to = "user")
    y_mid <- start_y - 2.5

    x0 <- x_mid - bar_width / 2
    x1 <- x_mid + bar_width / 2
    y0 <- y_mid - bar_height / 2
    y1 <- y_mid + bar_height / 2

    step_x <- chip_width_user

    # Draw chips
    for (k in seq_along(breaks_to_show)) {
      b <- breaks_to_show[k]
      color <- rownames(C)[order(abs(C[, "breaks"] - b), decreasing = FALSE)[1]]

      cx0 <- x0 + (k - 1) * step_x
      cx1 <- x0 + k * step_x

      rect(cx0, y0, cx1, y1, col = color, border = NA)
    }

    # Labels at ends, below bar
    lab_y <- y0 - lab_offset

    text(x0, lab_y,
      labels = format(round(breakList[1], 4)),
      adj    = c(0, 1),
      cex    = label_cex,
      col    = label_col
    )

    text(x1, lab_y,
      labels = format(round(tail(breakList, 1), 4)),
      adj    = c(1, 1),
      cex    = label_cex,
      col    = label_col
    )

    mid_x <- (x0 + x1) / 2
    text(mid_x, lab_y,
      labels = "1",
      adj    = c(0.5, 1),
      cex    = label_cex,
      col    = label_col
    )
  }
}

is_gene_in_bedpe <- function(gene_start, gene_end, pairs) {
  if (!is.data.frame(pairs) || nrow(pairs) == 0L) {
    return(FALSE)
  }
  ov1 <- (gene_end >= pairs$start1) & (gene_start <= pairs$end1)
  ov2 <- (gene_end >= pairs$start2) & (gene_start <= pairs$end2)
  any(ov1 | ov2, na.rm = TRUE)
}

is_gene_in_any_bedpe <- function(start_bp, end_bp, bedpe_list) {
  if (identical(bedpe_list, "FALSE") || !length(bedpe_list)) {
    return(FALSE)
  }
  stopifnot(is.list(bedpe_list))
  any(vapply(
    bedpe_list,
    function(p) is_gene_in_bedpe(start_bp, end_bp, p),
    logical(1)
  ))
}

draw_default_annotations <- function(bed_file,
                                     depth = 0,
                                     label_genes = TRUE,
                                     merge_exons = TRUE,
                                     show_arrows = TRUE,
                                     pack_lanes = TRUE,
                                     lane_gap = 1.25,
                                     group_gap = 1,
                                     alt_tss_file = NULL,
                                     alt_tss_label = TRUE,
                                     alt_tss_label_cex = 0.26,
                                     label_cex = 0.26,
                                     label_pad_in = 0.03,
                                     page_bottom_y = -Inf) {

  # Constants
  exon_height  = 0.25
  base_lwd     = 0.4
  arrow_dx     = 0.20
  arrow_dy     = 0.14
  arrow_gap    = 0.3
  alt_tss_col  = "#FF5656"
  tss_half_w   = 0.25 * 0.3
  label_family = "DejaVu Sans Mono"

  if (ann_style == "axiotl") {
    lane_gap <- lane_gap + 0.4
  }

  # Vertical layout for a single gene lane
  lane_layout <- function(y_base,
                          gene_label_gap = 0.40,
                          tss_tip_gap = 0.10,
                          tss_label_gap = 0.05) {
    y_gene <- y_base
    y_gene_label <- y_gene - exon_height / 2 - gene_label_gap
    y_tss_tip <- y_gene + exon_height / 2 + tss_tip_gap
    y_tss_base <- y_tss_tip + exon_height
    y_tss_label <- y_tss_base + tss_label_gap

    list(
      y_gene       = y_gene,
      y_gene_label = y_gene_label,
      y_tss_tip    = y_tss_tip,
      y_tss_base   = y_tss_base,
      y_tss_label  = y_tss_label
    )
  }

  # Draw a single alt TSS triangle
  draw_alt_tss_triangle <- function(bp_pos, y_tip_user, col) {
    x_tip <- map_x(bp_pos)
    polygon(
      c(x_tip - tss_half_w, x_tip + tss_half_w, x_tip),
      c(y_tip_user + exon_height, y_tip_user + exon_height, y_tip_user),
      col = col, border = NA, xpd = NA
    )
  }

  # Load data
  min_y_drawn <- Inf

  bed <- load_and_subset_bed(bed_file, interval_chr, interval_start, interval_end)
  if (is.null(bed)) return(invisible(NULL))

  alt_tss_tbl <- load_and_subset_alt_tss(alt_tss_file, interval_chr,
                                          interval_start, interval_end)

  # Coordinate transform
  n_bins_eff <- ncol(A) - 1
  scale_factor <- (n_bins_eff * w) / interval_len
  x0 <- start_x / 2 - 0.5 * w

  map_x <- function(pos_bp) {
    x0 + (pos_bp - interval_start) * scale_factor
  }

  ruler_top <- start_y / 2 + 0.5 * w - 0.001 * w
  base_y0 <- ruler_top - depth

  # Build gene spans and assign lanes
  span_tbl <- build_gene_spans(bed, alt_tss_tbl)
  if (is.null(span_tbl)) return(invisible(NULL))

  genes <- sort(unique(span_tbl$gene))

  lane_info <- assign_all_lanes(span_tbl, label_genes, label_cex, label_pad_in,
                                 scale_factor, lane_gap, group_gap)

  group_offsets  <- lane_info$group_offsets
  lane_by_gene   <- lane_info$lane_by_gene
  group_by_gene  <- lane_info$group_by_gene
  strand_by_gene <- lane_info$strand_by_gene

  # BEDPE-aware dimming
  has_bedpe_for_genes <- FALSE
  bedpe_for_genes <- NULL

  if (exists("bedpe", inherits = TRUE)) {
    bedpe_for_genes <- get("bedpe", inherits = TRUE)
    if (is.list(bedpe_for_genes) &&
      length(bedpe_for_genes) &&
      !identical(bedpe_for_genes, "FALSE")) {
      has_bedpe_for_genes <- TRUE
    }
  }

  # Main drawing loop
  genes_skipped <- 0L

  for (g in genes) {
    df <- bed[bed$gene == g, , drop = FALSE]
    if (!nrow(df)) next

    g_label <- toupper(g)
    bt <- df$biotype[1]
    grp <- group_by_gene[[g]]
    lane_idx <- if (pack_lanes) lane_by_gene[[g]] else 1L

    gx1_bp <- min(df$start)
    gx2_bp <- max(df$end)
    if (!(is.finite(gx1_bp) && is.finite(gx2_bp)) || gx2_bp <= gx1_bp) next

    gcol_base <- map_biotype_col(bt)
    bcol <- gcol_base

    # BEDPE dimming
    dim_this_gene <- FALSE
    if (isTRUE(has_bedpe_for_genes)) {
      in_any_anchor <- is_gene_in_any_bedpe(gx1_bp, gx2_bp, bedpe_for_genes)
      dim_this_gene <- !isTRUE(in_any_anchor)
    }
    if (dim_this_gene) {
      bcol <- grDevices::adjustcolor(gcol_base, alpha.f = 0.40)
    }

    tss_col <- grDevices::adjustcolor(alt_tss_col,
      alpha.f = if (dim_this_gene) 0.40 else 1.00)

    # Vertical layout
    class_offset <- group_offsets[[grp]]
    y_base <- base_y0 - class_offset - (lane_idx - 1L) * lane_gap

    if (is.finite(page_bottom_y) && y_base < page_bottom_y) {
      genes_skipped <- genes_skipped + 1L
      next
    }

    layout <- lane_layout(y_base)

    gx1 <- map_x(gx1_bp)
    gx2 <- map_x(gx2_bp)

    # Gene baseline
    inset <- (gx2 - gx1) * 0.001
    segments(gx1 + inset, layout$y_gene, gx2 - inset, layout$y_gene,
      col = bcol, lwd = base_lwd, xpd = NA, lend = "butt"
    )

    # Exons
    ex_draw <- prepare_exon_draw(df, gx1_bp, gx2_bp, g, strand_by_gene[[g]], bt, merge_exons)

    if (nrow(ex_draw)) {
      y1 <- layout$y_gene - exon_height / 2
      y2 <- layout$y_gene + exon_height / 2

      for (k in seq_len(nrow(ex_draw))) {
        xs <- map_x(ex_draw$start[k])
        xe <- map_x(ex_draw$end[k])
        if (is.finite(xs) && is.finite(xe) && xe > xs) {
          rect(xs, y1, xe, y2,
            col = bcol, border = bcol, lwd = 0.1, xpd = NA
          )
        }
      }
    }

    # Strand arrow
    if (show_arrows && grp != "mirna") {
      dir <- if (strand_by_gene[[g]] == "-") "left" else "right"

      three_prime_bp <- if (nrow(ex_draw)) {
        if (dir == "right") max(ex_draw$end, na.rm = TRUE) else min(ex_draw$start, na.rm = TRUE)
      } else {
        if (dir == "right") gx2_bp else gx1_bp
      }

      sgn <- if (dir == "right") 1 else -1
      tip_x <- map_x(three_prime_bp) + sgn * arrow_gap

      segments(tip_x - sgn * arrow_dx, layout$y_gene - arrow_dy,
        tip_x, layout$y_gene,
        col = bcol, lwd = 0.4, xpd = NA
      )
      segments(tip_x - sgn * arrow_dx, layout$y_gene + arrow_dy,
        tip_x, layout$y_gene,
        col = bcol, lwd = 0.4, xpd = NA
      )
    }

    # Alt TSS cluster
    if (!is.null(alt_tss_tbl)) {
      alt_g <- alt_tss_tbl[alt_tss_tbl$gene == g, , drop = FALSE]
      if (nrow(alt_g)) {
        xs_alt <- map_x(alt_g$start)
        xs_alt <- xs_alt[is.finite(xs_alt)]

        if (length(xs_alt)) {
          segments(
            min(xs_alt) - (tss_half_w * 0.6), layout$y_tss_base,
            max(xs_alt) + (tss_half_w * 0.6), layout$y_tss_base,
            col = tss_col, lwd = base_lwd, xpd = NA
          )

          for (ii in seq_len(nrow(alt_g))) {
            bp_pos <- alt_g$start[ii]
            if (!is.finite(bp_pos)) next
            draw_alt_tss_triangle(bp_pos, layout$y_tss_tip, tss_col)
          }
        }
      }
    }

    # Gene label
    if (label_genes && nzchar(g)) {
      text(
        (gx1 + gx2) / 2, layout$y_gene_label,
        labels = g_label, col = bcol,
        cex = label_cex, xpd = NA, family = label_family
      )
      label_height <- strheight(g_label, cex = label_cex, units = "user")
      min_y_drawn <- min(min_y_drawn, layout$y_gene_label - label_height)
    } else {
      min_y_drawn <- min(min_y_drawn, layout$y_gene - exon_height / 2)
    }
  }

  # Return legend info 
  if (!is.finite(min_y_drawn)) min_y_drawn <- NA_real_
  last_annotation_min_y <<- min_y_drawn

  invisible(build_annotation_legend_info(genes_skipped, min_y_drawn))
}

draw_feature_2 <- function(chr, start_pos, end_pos, label_txt, color,
                           y_offset, flag_text, ruler_top,
                           gene_cex = gene_cex, coord_cex = coord_cex,
                           point_mode = FALSE, point_cex = 0.3) {
  
  # scale (ncol(A) -1 to match strawr boundary indexing (prevents annotation drift)
  n_bins_eff <- ncol(A) - 1
  scale_factor <- (n_bins_eff * w) / interval_len

  # left edge of the first bin on the base
  x0 <- start_x / 2 - 0.5 * w

  # map bp -> user x
  x_start_edge <- x0 + (start_pos - interval_start) * scale_factor
  x_end_edge <- x0 + (end_pos - interval_start) * scale_factor

  x_anchor <- x_start_edge
  y_pos <- ruler_top - y_offset

  # Optional span line
  # lines(c(x_start_edge, x_end_edge), rep(y_pos, 2), col = color, lwd = 0.5)

  if (flag_text) {
    composite_label <- sprintf("%s:%d-%d", chr, start_pos, end_pos)

    text(x_anchor, y_pos,
      labels = composite_label, col = color,
      offset = -0.16, pos = 2, cex = coord_cex, srt = 45, xpd = NA
    )
    text(x_anchor, y_pos,
      labels = label_txt, col = color,
      offset = 0.16, pos = 2, cex = gene_cex, srt = 45, xpd = NA
    )

    segments(x_anchor, y_pos, x_anchor, ruler_top,
      col = "#B2B2B2", lwd = 0.2, lty = "dotted", xpd = NA
    )
  } else {
    if (point_mode) {
      x_mid <- (x_start_edge + x_end_edge) / 2
      points(x_mid, y_pos, pch = 16, col = color, cex = point_cex, xpd = NA)
    } else {
      segments(x_start_edge, y_pos, x_end_edge, y_pos, col = color, lwd = 0.8, xpd = NA)
      segments(x_start_edge, y_pos, x_start_edge, ruler_top, col = "#B2B2B2", lwd = 0.2, lty = "dotted", xpd = NA)
      segments(x_end_edge, y_pos, x_end_edge, ruler_top, col = "#B2B2B2", lwd = 0.2, lty = "dotted", xpd = NA)
    }
  }
}

draw_bed <- function(bed_file, color, depth, flag_text, ruler_top,
                     point_mode = FALSE, point_cex = 0.3) {
  if (!file.exists(bed_file)) {
    stop("draw_bed error: file not found")
  }

  bed <- read.table(bed_file, as.is = TRUE, header = FALSE)
  colnames(bed) <- c("chr", "start", "end")

  # Subset the BED features to visualized genomic interval
  bed <- bed[bed$chr == interval_chr & bed$start < interval_end & bed$end > interval_start, ]
  if (nrow(bed) == 0) {
    warning("draw_bed warning: no features in the specified interval")
    return()
  }

  # Clip the feature boundaries to visualized interval
  bed$start <- pmax(bed$start, interval_start)
  bed$end <- pmin(bed$end, interval_end)

  if (!flag_text) {
    for (i in seq_len(nrow(bed))) {
      b <- bed[i, ]
      draw_feature_2(b$chr, b$start, b$end,
        label_txt = "",
        color = color, y_offset = depth, flag_text = FALSE,
        ruler_top = ruler_top,
        point_mode = point_mode, point_cex = point_cex
      )
    }
    return()
  }

  # Gene annotations
  bin_index <- floor(bed$start / bin_size)
  bed <- bed[order(bin_index), ]
  bin_counts <- table(bin_index)

  for (bin in names(bin_counts)) {
    feature_indices <- which(bin_index == bin)
    for (j in seq_along(feature_indices)) {
      b <- bed[feature_indices[j], ]
      # Adjust vertical offset for multiple genes per bin
      current_offset <- depth - (j - 1) * (-4)

      label_txt <- if (ncol(b) >= 4) b[[4]] else ""

      # define bounds for the overlap test
      gene_start <- b$start
      gene_end <- b$end - 1L

      if (is_gene_in_any_bedpe(gene_start, gene_end, bedpe)) {
        gene_cex <- 0.4
        coord_cex <- 0.2
        colour <- color
      } else {
        gene_cex <- 0.22
        coord_cex <- 0.18
        colour <- "#aaaaaa"
      }

      draw_feature_2(b$chr, b$start, b$end, label_txt,
        color = colour, current_offset,
        flag_text = TRUE, ruler_top = ruler_top,
        gene_cex = gene_cex, coord_cex = coord_cex
      )
    }
  }
}

draw_ruler <- function(start_x, start_y, w, n_bins, bin_coordinates,
                       gap_y = 1.0,
                       col_base = "#B0B0B0",
                       col_ticks = "#afafafff",
                       col_labels = "#6E6E6E",
                       cex_labels = 0.45,
                       y_offset = 0) {
  bin_bp <- as.numeric(bin_size)

  # Horizontal span
  x_left <- start_x / 2 - (0.5 * w)
  x_right <- (start_x + (2 * n_bins - 2) * w) / 2 + (0.5 * w)

  # flat tips
  old_lend <- par("lend")
  on.exit(par(lend = old_lend), add = TRUE)
  par(lend = "butt")

  # Vertical anchoring
  target_top <- start_y / 2 + 0.5 * w - 0.001 * w + y_offset
  H <- gap_y
  len_mb <- H * 0.89
  len_500k <- H * 0.65
  len_100k <- H * 0.65
  len_50k <- H * 0.50
  len_10k <- H * 0.35

  base_y <- target_top - len_mb

  # Background panel
  label_offset_mb <- 0.05
  label_offset_mid <- 0
  cex_mb <- 0.35
  cex_mid <- 0.25

  bg_pad_x <- 0.40 * w
  bg_pad_y <- 2 * w
  rect(
    xleft = x_left - bg_pad_x,
    ybottom = base_y - max(label_offset_mb, label_offset_mid) - bg_pad_y,
    xright = x_right + bg_pad_x,
    ytop = target_top,
    col = "#eeeeee", border = NA, xpd = NA
  )

  # Baseline
  segments(x_left, base_y, x_right, base_y, col = col_base, lwd = 0.2, xpd = NA)

  # Which tick families are relevant for this resolution
  show_10kb <- (bin_bp <= 5000)
  show_50kb <- (bin_bp %in% c(1000, 10000))
  show_100kb <- (bin_bp %in% c(1000, 5000, 10000, 25000)) ||
                  (bin_bp == 50000  && w >= 1.5) ||
                  (bin_bp == 100000 && w >= 4)
  show_500kb <- (bin_bp %in% c(50000, 100000)) ||
                  (bin_bp == 1000  && w < 0.1)  ||
                  (bin_bp == 5000  && w < 0.2)  ||
                  (bin_bp == 10000 && w < 0.30) ||
                  (bin_bp == 25000 && w >= 0.15 && w < 0.75)
  major_step_bp <- if (bin_bp >= 1e6) bin_bp else 1e6

  # Label suppression rules for clutter
  suppress_label_100kb <-
      (bin_bp == 1000  && w < 0.05)  ||
      (bin_bp == 5000  && w < 0.2)  ||
      (bin_bp == 10000 && w < 0.4) ||
      (bin_bp == 25000 && w < 0.9) ||
      (bin_bp == 50000 && w < 2)
  suppress_label_50kb  <- (bin_bp == 5000 && w <= 0.8) ||  (bin_bp == 1000 && w < 0.1)
  suppress_label_10kb  <- (bin_bp == 5000 && w <= 2) ||  (bin_bp == 1000 && w < 0.5)


  # First bin coordinate (bp)
  first_bp <- suppressWarnings(as.numeric(bin_coordinates[1]))
  if (is.na(first_bp)) stop("ruler error: bin_coordinates[1] must be numeric base-pair coordinate")

  # Draw ticks
  for (k in 1:n_bins) {
    xk <- ((start_x + (2 * k - 2) * w) / 2) - (0.5 * w)
    o <- (k - 1)
    coord_bp <- first_bp + o * bin_bp

    at_mb <- (coord_bp %% major_step_bp == 0L)
    at_500kb <- show_500kb && (coord_bp %% 500000 == 0L)
    at_100kb <- show_100kb && (coord_bp %% 100000 == 0L)
    at_50kb <- show_50kb && (coord_bp %% 50000 == 0L)
    at_10kb <- show_10kb && (coord_bp %% 10000 == 0L)

    if (at_mb) {
      segments(xk, base_y, xk, base_y + len_mb, col = col_ticks, lwd = 0.25, xpd = NA)
      mb_lab <- if (major_step_bp == 1e6) {
        paste0(formatC(coord_bp / 1e6, digits = 0, format = "f"), " MB")
      } else {
        paste0(formatC(coord_bp / 1e6, digits = 1, format = "f"), " MB")
      }
      text(xk, base_y - label_offset_mb,
        labels = mb_lab,
        cex = cex_mb, col = col_labels, adj = c(0.5, 1), xpd = NA
      )
    } else if (at_500kb) {
      segments(xk, base_y, xk, base_y + len_500k, col = col_ticks, lwd = 0.18, xpd = NA)
      lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
      text(xk, base_y - label_offset_mid,
        labels = lab,
        cex = cex_mid, col = col_labels, adj = c(0.5, 1), xpd = NA
      )
    } else if (at_100kb) {
      segments(xk, base_y, xk, base_y + len_100k, col = col_ticks, lwd = 0.18, xpd = NA)
      if (!suppress_label_100kb) {
        lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
        text(xk, base_y - label_offset_mid,
          labels = lab,
          cex = cex_mid, col = col_labels, adj = c(0.5, 1), xpd = NA
        )
      }
    } else if (at_50kb) {
      segments(xk, base_y, xk, base_y + len_50k, col = col_ticks, lwd = 0.16, xpd = NA)
      # Label 50 kb only at 1 kb resolution; suppress if cramped
      if (bin_bp %in% c(1000, 5000) && !suppress_label_50kb) {
        lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
        text(xk, base_y - label_offset_mid,
          labels = lab,
          cex = cex_mid, col = col_labels, adj = c(0.5, 1), xpd = NA
        )
      }
    } else if (at_10kb) {
      segments(xk, base_y, xk, base_y + len_10k, col = col_ticks, lwd = 0.14, xpd = NA)
      if (!suppress_label_10kb) {                                  # NEW
        lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
        text(xk, base_y - label_offset_mid,
          labels = lab,
          cex = cex_mid, col = col_labels, adj = c(0.5, 1), xpd = NA
        )
      }
    }
  }
  invisible(list(target_top = target_top, base_y = base_y))
}

# Clip bedpe polygons that cross below diagonal
# Sutherland–Hodgman algorithm
clip_poly_ymin <- function(x, y, ymin) {
  n <- length(x)
  if (n < 3) {
    return(list(x = numeric(0), y = numeric(0)))
  }

  outx <- c()
  outy <- c()
  for (k in seq_len(n)) {
    i <- k
    j <- if (k == n) 1 else (k + 1)
    x1 <- x[i]
    y1 <- y[i]
    x2 <- x[j]
    y2 <- y[j]
    in1 <- (y1 >= ymin)
    in2 <- (y2 >= ymin)

    # intersection with horizontal line y = ymin
    inter <- function(x1, y1, x2, y2, ymin) {
      if (y2 == y1) {
        return(c(NA_real_, NA_real_))
      }
      t <- (ymin - y1) / (y2 - y1)
      c(x1 + t * (x2 - x1), ymin)
    }

    if (in1 && in2) {
      # both inside: keep end
      outx <- c(outx, x2)
      outy <- c(outy, y2)
    } else if (in1 && !in2) {
      # leaving: add intersection
      p <- inter(x1, y1, x2, y2, ymin)
      outx <- c(outx, p[1])
      outy <- c(outy, p[2])
    } else if (!in1 && in2) {
      # entering: add intersection THEN end
      p <- inter(x1, y1, x2, y2, ymin)
      outx <- c(outx, p[1], x2)
      outy <- c(outy, p[2], y2)
    } else {
      # both outside: add nothing
    }
  }
  list(x = outx, y = outy)
}

draw_contacts <- function(A,
                          C,
                          left_side,
                          bedpe,
                          focus_start_idx = 1L,
                          focus_len = min(nrow(A), ncol(A)),
                          band_height = focus_len - 1L,
                          include_expand_viewport = TRUE) {
  # Helpers
  get_center <- function(i, j, w) {
    cx <- start_x + (j - 1) * w + (i - 1) * w
    cy <- start_y + (i - 1) * (-w) + j * w
    c(cx, cy)
  }

  get_vertices <- function(center, w) {
    list(
      left   = c((center[1] - w) / 2, center[2] / 2),
      top    = c(center[1] / 2, (center[2] + w) / 2),
      right  = c((center[1] + w) / 2, center[2] / 2),
      bottom = c(center[1] / 2, (center[2] - w) / 2)
    )
  }

  draw_diamond <- function(center, w, color, highlight = FALSE) {
    verts <- get_vertices(center, w)
    x_coords <- c(verts$left[1], verts$top[1], verts$right[1], verts$bottom[1])
    y_coords <- c(verts$left[2], verts$top[2], verts$right[2], verts$bottom[2])
    polygon(x_coords, y_coords, col = color, border = NA)
    if (highlight) {
      polygon(x_coords, y_coords, col = color, border = bedpe_color, lwd = 0.4)
    }
  }

  n_rows <- nrow(A)
  n_cols <- ncol(A)
  legend_out <- NULL

  # ---------------- TRIANGLE MODE -----------------------------------
  if (!include_expand_viewport) {
    # Footprint outline
    d <- min(nrow(A), ncol(A))

    BL <- c(start_x / 2 - w, start_y / 2)
    TOP <- c((start_x + (ncol(A) - 1) * w) / 2, (start_y + (ncol(A) + 1) * w) / 2)
    BR <- c((start_x + (2 * d - 2) * w) / 2 + w, start_y / 2)

    polygon(
      x = c(BL[1], TOP[1], BR[1]),
      y = c(BL[2], TOP[2], BR[2]),
      col = "white", border = "#eeeeee", lwd = 0.08
    )

    # Draw diamonds
    for (i in seq_len(n_rows)) {
      for (j in i:n_cols) {
        val <- A[i, j]
        if (is.na(val) || val == 0) next
        color <- rownames(C[order(abs(C[, "breaks"] - val), decreasing = FALSE), ])[1]
        center <- get_center(i, j, w)
        draw_diamond(center, w, color, highlight = FALSE)
      }
    }

    # Ruler uses all columns
    ruler <- draw_ruler(start_x, start_y, w, n_cols, colnames(A))

    # BEDPE highlighting
    if (!identical(bedpe, "FALSE")) {
      prep <- prep_bedpe_layers(bedpe, bedpe_color_param = bedpe_color)
      pairs_list <- prep$pairs_list
      bedpe_colors <- prep$color_map

      base_y <- start_y / 2 + 0.5 * w

      for (nm in names(pairs_list)) {
        pairs <- pairs_list[[nm]]
        if (is.null(pairs) || nrow(pairs) == 0L) next
        bedpe_col <- unname(bedpe_colors[[nm]])

        for (p in seq_len(nrow(pairs))) {
          i_start <- pairs$i[p]
          i_end <- pairs$i_end[p]
          j_start <- pairs$j[p]
          j_end <- pairs$j_end[p]

          if (i_start == i_end && j_start == j_end) {
            center <- get_center(i_start, j_start, w)
            v <- get_vertices(center, w)
            px <- c(v$left[1], v$top[1], v$right[1], v$bottom[1])
            py <- c(v$left[2], v$top[2], v$right[2], v$bottom[2])
          } else {
            vL <- get_vertices(get_center(i_start, j_start, w), w)$left
            vT <- get_vertices(get_center(i_start, j_end, w), w)$top
            vR <- get_vertices(get_center(i_end, j_end, w), w)$right
            vB <- get_vertices(get_center(i_end, j_start, w), w)$bottom
            px <- c(vL[1], vT[1], vR[1], vB[1])
            py <- c(vL[2], vT[2], vR[2], vB[2])
          }

          clp <- clip_poly_ymin(px, py, base_y)
          if (length(clp$x) >= 3) {
            polygon(clp$x, clp$y, border = bedpe_col, col = NA, lwd = 0.8, xpd = NA)
          }
        }
      }

      labs <- names(bedpe_colors)
      legend_out <- list(labels = labs, colors = unname(bedpe_colors))
    }

    return(invisible(list(legend = legend_out, ruler = ruler)))
  }

  # ---------------- expand_viewport ---------------------

  max_offset <- as.integer(band_height) # max allowed j - i in bins
  offset <- as.integer(focus_start_idx) - 1L # number of bins before the focus block
  focus_end_idx <- offset + focus_len # 1-based index of last focus bin in A_big
  focus_idx <- (offset + 1L):focus_end_idx # 1-based indices of all focus bins in A_big

  # First pass: bounding box
  xmin <- Inf
  xmax <- -Inf
  ymin <- Inf
  ymax <- -Inf

  for (i_big in seq_len(n_rows)) {
    for (j_big in seq_len(n_cols)) {
      if (j_big < i_big) next

      i_plot <- i_big - offset
      j_plot <- j_big - offset

      if (j_plot < i_plot) next # skip lower triangle
      if ((j_plot - i_plot) > max_offset) next # too far from diagonal

      # Do not let expand_viewport extend past the NxN footprint
      if (i_plot < 1L || j_plot < 1L || i_plot > focus_len || j_plot > focus_len) next # outside NxN footprint

      center <- get_center(i_plot, j_plot, w)
      verts <- get_vertices(center, w)

      xs <- c(verts$left[1], verts$top[1], verts$right[1], verts$bottom[1])
      ys <- c(verts$left[2], verts$top[2], verts$right[2], verts$bottom[2])

      xmin <- min(xmin, xs)
      xmax <- max(xmax, xs)
      ymin <- min(ymin, ys)
      ymax <- max(ymax, ys)
    }
  }

  if (is.finite(xmin)) {
    polygon(
      x = c(xmin, xmax, xmax, xmin),
      y = c(ymin, ymin, ymax, ymax),
      col = "white", border = "#eeeeee", lwd = 0.08
    )
  }

  # Save for clipping in the draw pass
  bbox_xmin <- xmin
  bbox_xmax <- xmax
  bbox_ymin <- ymin
  bbox_ymax <- ymax

  # Second pass: draw diamonds (skipping zeros)
  # Set clip so expand_viewport get drawn but can't exceed the square background
  usr0 <- par("usr")
  graphics::clip(bbox_xmin, bbox_xmax, bbox_ymin, bbox_ymax)
  on.exit(graphics::clip(usr0[1], usr0[2], usr0[3], usr0[4]), add = TRUE)

  for (i_big in seq_len(n_rows)) {
    for (j_big in seq_len(n_cols)) {
      if (j_big < i_big) next

      i_plot <- i_big - offset
      j_plot <- j_big - offset

      if (j_plot < i_plot) next
      if ((j_plot - i_plot) > max_offset) next

      val <- A[i_big, j_big]
      if (is.na(val) || val == 0) next

      color <- rownames(C[order(abs(C[, "breaks"] - val), decreasing = FALSE), ])[1]
      center <- get_center(i_plot, j_plot, w)
      draw_diamond(center, w, color, highlight = FALSE)
    }
  }

  # Ruler only over the focus interval
  all_bins <- colnames(A)
  focus_labs <- all_bins[focus_idx]
  ruler <- draw_ruler(start_x, start_y, w, focus_len, focus_labs)

  # BEDPE
  if (!identical(bedpe, "FALSE")) {
    prep <- prep_bedpe_layers(bedpe, bedpe_color_param = bedpe_color)
    pairs_list <- prep$pairs_list
    bedpe_colors <- prep$color_map

    base_y <- start_y / 2 + 0.5 * w

    for (nm in names(pairs_list)) {
      pairs <- pairs_list[[nm]]
      if (is.null(pairs) || nrow(pairs) == 0L) next
      bedpe_col <- unname(bedpe_colors[[nm]])

      for (p in seq_len(nrow(pairs))) {
        i_start_big <- pairs$i[p]
        i_end_big <- pairs$i_end[p]
        j_start_big <- pairs$j[p]
        j_end_big <- pairs$j_end[p]

        i_start <- i_start_big
        i_end <- i_end_big
        j_start <- j_start_big
        j_end <- j_end_big

        # skip anything outside the focus footprint
        if (i_start < 1L || j_start < 1L || i_end > focus_len || j_end > focus_len) next

        if (i_start == i_end && j_start == j_end) {
          center <- get_center(i_start, j_start, w)
          v <- get_vertices(center, w)
          px <- c(v$left[1], v$top[1], v$right[1], v$bottom[1])
          py <- c(v$left[2], v$top[2], v$right[2], v$bottom[2])
        } else {
          vL <- get_vertices(get_center(i_start, j_start, w), w)$left
          vT <- get_vertices(get_center(i_start, j_end, w), w)$top
          vR <- get_vertices(get_center(i_end, j_end, w), w)$right
          vB <- get_vertices(get_center(i_end, j_start, w), w)$bottom
          px <- c(vL[1], vT[1], vR[1], vB[1])
          py <- c(vL[2], vT[2], vR[2], vB[2])
        }

        clp <- clip_poly_ymin(px, py, base_y)
        if (length(clp$x) >= 3) {
          polygon(clp$x, clp$y, border = bedpe_col, col = NA, lwd = 0.8, xpd = NA)
        }
      }
    }

    labs <- names(bedpe_colors)
    legend_out <- list(labels = labs, colors = unname(bedpe_colors))
  }
  return(invisible(list(legend = legend_out, ruler = ruler)))
}

draw_title <- function(
  interval_chr,
  interval_start,
  interval_end
  ) {
  label_family <- "DejaVu Sans Mono"

  if (analysis_type == "single_sample") {
    sample_name <- basename(sample_dir)
    title <- sprintf(
      "Sample: %s\nRange: %s:%d-%d \nGenome Build: %s\nResolution: %s",
      sample_name,
      interval_chr,
      interval_start,
      interval_end,
      genome_build,
      bin_size
    )
    text(0, top, labels = title, col = "#888888", pos = 4, cex = 0.5, family = label_family)
  } else if (analysis_type == "two_sample") {
    sample_nameA <- basename(sample_dirA)
    sample_nameB <- basename(sample_dirB)
    title <- sprintf(
      "Delta = \n%s - %s\nRange: %s:%d-%d\nGenome Build: %s\nResolution: %s",
      sample_nameB, sample_nameA,
      interval_chr,
      interval_start,
      interval_end,
      genome_build,
      bin_size
    )
    text(0, top - 1, labels = title, col = "#888888", pos = 4, cex = 0.5, family = label_family)
  }
}

estimate_crop <- function(min_y_user,
                          page_height_in = 8.5,
                          min_crop_pts = 0,
                          max_crop_pts = 220,
                          padding_pts = 15) {
  if (!is.finite(min_y_user)) {
    return(0L)
  }

  # page height in points
  page_height_pts <- page_height_in * 72

  # portion of page below the lowest annotation 
  blank_pts <- (min_y_user / 100) * page_height_pts

  # keep a bit of padding below the labels
  crop_raw <- blank_pts - padding_pts

  cut <- max(min_crop_pts, min(round(crop_raw), max_crop_pts))
  if (cut < 0) cut <- 0L

  cut
}

crop_pdf <- function(path, cut_bottom = 20) {
  # Only run for .pdf outputs
  if (!grepl("\\.pdf$", path, ignore.case = TRUE)) {
    return(invisible(FALSE))
  }

  # Need both pdfcrop and pdfinfo
  if (nzchar(Sys.which("pdfcrop")) == FALSE) {
    message("crop_pdf: pdfcrop not found on PATH; skipping ", path)
    return(invisible(FALSE))
  }
  if (nzchar(Sys.which("pdfinfo")) == FALSE) {
    message("crop_pdf: pdfinfo not found on PATH; skipping ", path)
    return(invisible(FALSE))
  }

  # Get page size from pdfinfo (in points)
  info <- system2("pdfinfo", shQuote(path), stdout = TRUE)
  size_line <- grep("^Page size:", info, value = TRUE)
  if (length(size_line) == 0) {
    return(invisible(FALSE))
  }

  nums <- as.numeric(unlist(regmatches(size_line, gregexpr("[0-9.]+", size_line))))
  width <- nums[1]
  height <- nums[2]

  # Define bounding box: left, bottom, right, top
  # Keep left = 0 and right/top = full size, only crop the bottom
  lly <- cut_bottom
  if (lly >= height) {
    message("crop_pdf_bottom_only: cut_bottom >= page height; skipping ", path)
    return(invisible(FALSE))
  }

  bbox <- sprintf("%f %f %f %f", 0, lly, width, height)

  cmd <- sprintf(
    "pdfcrop --bbox '%s' --margins '0 0 0 0' %s %s",
    bbox,
    shQuote(path),
    shQuote(path) # overwrite in place
  )

  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  invisible(TRUE)
}

preflight_custom_beds <- function(ann_custom, interval_chr, interval_start, interval_end) {
  if (identical(ann_custom, "NONE") || !nzchar(ann_custom)) {
    return(character(0))
  }

  files <- strsplit(ann_custom, "[[:space:]]+")[[1]]
  files <- files[nzchar(files)]
  if (!length(files)) {
    return(character(0))
  }

  valid <- character(0)

  for (f in files) {
    if (!file.exists(f)) {
      cat(sprintf("Warning: custom annotation file not found: %s\n", f))
      next
    }

    x <- tryCatch(read.table(f, header = FALSE, as.is = TRUE), error = function(e) NULL)
    if (is.null(x) || ncol(x) < 3) {
      cat(sprintf("Warning: could not read (or <3 cols): %s\n", basename(f)))
      next
    }

    x <- x[, 1:3, drop = FALSE]
    colnames(x) <- c("chr", "start", "end")

    in_range <- x$chr == interval_chr & x$start < interval_end & x$end > interval_start
    if (any(in_range, na.rm = TRUE)) {
      valid <- c(valid, f)
    } else {
      cat(sprintf(
        "Warning: %s has no features in interval %s:%d-%d\n",
        basename(f), interval_chr, interval_start, interval_end
      ))
    }
  }

  valid
}

render_triangle_plot <- function(output_name,
                                 A_mat,
                                 C_tbl,
                                 breaks_vec,
                                 bedpe,
                                 scale_fun = draw_scale,
                                 device = c("svg", "pdf"),
                                 ann_style = "standard",
                                 focus_start_idx = 1L,
                                 focus_len = min(nrow(A_mat), ncol(A_mat)),
                                 band_height = focus_len - 1L,
                                 left_side = 0,
                                 legend_mid_y = NULL) {
  # Constants
  legend_line_lwd <- 1.0
  custom_start <- 2.5
  custom_step <- 0.25
  gap_after_custom <- 1.2
  legend_max_chars <- 36
  bbox_gap <- 1.0
  bbox_pad_x <- 0.8
  bbox_pad_y <- 0.1

  # Base palette for custom annotations
  triangle_base_colors <- c("#E69F00", "#377EB8", "#4DAF4A", "#C77CFF", "#984EA3")

  legend_add <- function(label, color, type,
                         lwd = NA, ptbg = NA, ptlwd = NA, ptcex = NA) {
    legend_labels <<- c(legend_labels, label)
    legend_colors <<- c(legend_colors, color)
    legend_types <<- c(legend_types, type)
    legend_lwd <<- c(legend_lwd, lwd)
    legend_ptbg <<- c(legend_ptbg, ptbg)
    legend_ptlwd <<- c(legend_ptlwd, ptlwd)
    legend_ptcex <<- c(legend_ptcex, ptcex)
  }

  legend_add_many <- function(labels, colors, types, lwds = NA, ptbgs = NA, ptlwds = NA, ptcexs = NA) {
    stopifnot(length(labels) == length(colors), length(colors) == length(types))
    for (i in seq_along(labels)) {
      legend_add(
        labels[i], colors[i], types[i],
        if (length(lwds) == 1) lwds else lwds[i],
        if (length(ptbgs) == 1) ptbgs else ptbgs[i],
        if (length(ptlwds) == 1) ptlwds else ptlwds[i],
        if (length(ptcexs) == 1) ptcexs else ptcexs[i]
      )
    }
  }

  # Setup device and plot surface
  test_con <- tryCatch(file(output_name, open = "ab"), error = function(e) NULL)
  if (is.null(test_con)) {
    stop(sprintf(
      "Cannot write to '%s' -- the file may be open in another application. Please close it and try again.",
      output_name
    ))
  }
  close(test_con)

  device <- match.arg(device)
  if (device == "svg") {
    svg(output_name, width = 8.5, height = 8.5)
  } else {
    grDevices::cairo_pdf(output_name, width = 8.5, height = 8.5)
  }
  on.exit(dev.off(), add = TRUE)

  par(family = "DejaVu Sans Mono", omi = rep(0, 4), mai = rep(0, 4), bg = "#eeeeee")
  plot(NULL,
    xlim = c(0, 100), ylim = c(0, 100),
    xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", bty = "n", asp = 1
  )

  # Compute page bottom in user coordinates
  page_bottom_y <- grconvertY(0.02, from = "ndc", to = "user")

  # Initialize custom legend
  legend_labels <- c()
  legend_colors <- c()
  legend_types <- c()
  legend_lwd <- c()
  legend_ptbg <- c()
  legend_ptlwd <- c()
  legend_ptcex <- c()

  # Initialize default legend
  default_labels <- c()
  default_colors <- c()
  default_types <- c()
  default_lwd <- c()
  default_ptbg <- c()
  default_ptlwd <- c()
  default_ptcex <- c()

  # Draw title
  top <<- grconvertY(1 - 0.06, from = "ndc", to = "user")
  draw_title(interval_chr, interval_start, interval_end)

  scale_y <- if (isTRUE(flag_expand_viewport) && !is.null(legend_mid_y)) {
    legend_mid_y
  } else {
    start_y
  }

  # Draw scale bar
  scale_fun(A_mat, C_tbl,
    w = w, start_x = start_x, start_y = scale_y,
    breakList = breaks_vec, analysis_type = analysis_type,
    gap_in = 1, chip_len_in = 0.05, thickness_in = 0.05
  )

  # Draw contacts and bedpe highlights
  legend_info <- draw_contacts(
    A = A_mat,
    C = C_tbl,
    left_side = left_side,
    bedpe = bedpe,
    focus_start_idx = focus_start_idx,
    focus_len = focus_len,
    band_height = band_height,
    include_expand_viewport = flag_expand_viewport
  )

  # Compute annotation region horizontal extent (match ruler width)
  x_left <- start_x / 2 - 0.5 * w
  x_right <- (start_x + (2 * focus_len - 2) * w) / 2 + 0.5 * w

  # Track the number of custom annotation files actually drawn
  n_bed <- 0
  custom_bbox_bottom_actual <- NULL

  if (!identical(ann_custom, "NONE") && nzchar(ann_custom)) {
    annotation_files <- strsplit(ann_custom, "[[:space:]]+")[[1]]
    annotation_files <- annotation_files[nzchar(annotation_files)]

    n_files <- length(annotation_files)
    n_bed <- n_files

    track_depths <- seq(from = custom_start, by = custom_step, length.out = n_files)

    final_colors <- resolve_custom_colors(
      ann_custom_colors, n_files, base_colors = triangle_base_colors
    )

    # Draw each track
    for (i in seq_along(annotation_files)) {
      if (!file.exists(annotation_files[i])) next

      draw_bed(annotation_files[i], final_colors[i], track_depths[i], FALSE,
        ruler_top = legend_info$ruler$target_top
      )

      legend_add(basename(annotation_files[i]), final_colors[i], "line", lwd = 1.4)
    }

    # Custom annotation bounding box
    if (n_bed > 0) {
      ruler_top_val <- legend_info$ruler$target_top
      ruler_base_val <- legend_info$ruler$base_y

      custom_bbox_top <- ruler_base_val - bbox_gap
      custom_bbox_bottom <- ruler_top_val - track_depths[n_bed] - 0.5

      custom_bbox_bottom_actual <- custom_bbox_bottom - bbox_pad_y

      rect(
        xleft = x_left - bbox_pad_x,
        ybottom = custom_bbox_bottom - bbox_pad_y,
        xright = x_right + bbox_pad_x,
        ytop = custom_bbox_top + bbox_pad_y,
        col = "white", border = NA, xpd = NA
      )

      # Redraw custom tracks on top of the bounding box
      for (i in seq_along(annotation_files)) {
        if (!file.exists(annotation_files[i])) next
        draw_bed(annotation_files[i], final_colors[i], track_depths[i], FALSE,
          ruler_top = legend_info$ruler$target_top
        )
      }
    }
  }

  # Add bedpe info to legend
  if (!is.null(legend_info$legend) && length(legend_info$legend$labels)) {
    legend_add_many(legend_info$legend$labels, legend_info$legend$colors,
      types = rep("bedpe", length(legend_info$legend$labels))
    )
  }

  # Draw default annotations
  exon_depth <- if (n_bed > 0) {
    last_custom_depth <- custom_start + custom_step * (n_bed - 1)
    last_custom_depth + gap_after_custom
  } else {
    2.5
  }

  # Let genes draw to page bottom
  # Legend and overflow notice overlay on top if needed
  legend_floor_y <- page_bottom_y

  # Draw annotations based on style
  legend_info_exon <- NULL
  genes_skipped <- 0L

  if (identical(ann_style, "standard")) {
    exon_bed_file <- switch(genome_build,
      "hg38" = file.path(data_dir, genome_build, "reference", "gencode.hg38.annotation.bed"),
      "hg19" = file.path(data_dir, genome_build, "reference", "gencode.hg19.annotation.bed"),
      "mm10" = file.path(data_dir, genome_build, "reference", "gencode.mm10.annotation.bed"),
      NULL
    )

    if (is.null(exon_bed_file) || !file.exists(exon_bed_file)) {
      warning(
        "No exon BED file configured for genome_build = ", genome_build,
        "; skipping standard exon annotations."
      )
    } else {
      legend_info_exon <- draw_default_annotations(
        bed_file       = exon_bed_file,
        depth          = exon_depth,
        page_bottom_y  = legend_floor_y
      )
      genes_skipped <- legend_info_exon$genes_skipped
    }
  } else if (identical(ann_style, "axiotl")) {
    alt_tss_file <- switch(genome_build,
      "hg38" = file.path(data_dir, genome_build, "reference", "gencode.hg38.transcript_tss.bed"),
      "hg19" = file.path(data_dir, genome_build, "reference", "gencode.hg19.transcript_tss.bed"),
      "mm10" = file.path(data_dir, genome_build, "reference", "gencode.mm10.transcript_tss.bed"),
      NULL
    )

    exon_bed_file <- switch(genome_build,
      "hg38" = file.path(data_dir, genome_build, "reference", "gencode.hg38.annotation.bed"),
      "hg19" = file.path(data_dir, genome_build, "reference", "gencode.hg19.annotation.bed"),
      "mm10" = file.path(data_dir, genome_build, "reference", "gencode.mm10.annotation.bed"),
      NULL
    )

    missing_exon <- is.null(exon_bed_file) || !file.exists(exon_bed_file)
    missing_tss <- is.null(alt_tss_file) || !file.exists(alt_tss_file)

    if (missing_exon || missing_tss) {
      warning(
        "Missing annotation file(s) for genome_build = ", genome_build,
        ": exon_bed_file=", exon_bed_file, ", alt_tss_file=", alt_tss_file,
        "; skipping exon / TSS annotations."
      )
    } else {
      legend_info_exon <- draw_default_annotations(
        bed_file       = exon_bed_file,
        depth          = exon_depth,
        alt_tss_file   = alt_tss_file,
        page_bottom_y  = legend_floor_y
      )
      genes_skipped <- legend_info_exon$genes_skipped
    }

    # Add TSS info to legend
    if (!is.null(alt_tss_file)) {
      default_labels <- c(default_labels, "TSSs")
      default_colors <- c(default_colors, "#FF5656")
      default_types <- c(default_types, "tss")
      default_lwd <- c(default_lwd, NA)
      default_ptbg <- c(default_ptbg, "#FF5656")
      default_ptlwd <- c(default_ptlwd, NA)
      default_ptcex <- c(default_ptcex, 0.5)
    }
  } else if (identical(ann_style, "none")) {
    # Draw nothing
  }

  # Default annotations bounding box and bracket
  if (!is.null(legend_info_exon) && is.finite(last_annotation_min_y)) {
    ruler_top_val <- legend_info$ruler$target_top
    x_left_ann <- start_x / 2 - 0.5 * w

    y_top <- ruler_top_val - exon_depth
    y_bottom <- last_annotation_min_y - 1.5
    ref_bbox_top <- y_top + bbox_pad_y

    rect(
      xleft = x_left - bbox_pad_x,
      ybottom = y_bottom - bbox_pad_y,
      xright = x_right + bbox_pad_x,
      ytop = ref_bbox_top,
      col = "white", border = NA, xpd = NA
    )

    # Redraw default annotations on top of the bounding box
    if (identical(ann_style, "standard")) {
      if (!is.null(exon_bed_file) && file.exists(exon_bed_file)) {
        legend_info_exon <- draw_default_annotations(
          bed_file       = exon_bed_file,
          depth          = exon_depth,
          page_bottom_y  = legend_floor_y
        )
        genes_skipped <- legend_info_exon$genes_skipped
      }
    } else if (identical(ann_style, "axiotl")) {
      if (!missing_exon && !missing_tss) {
        legend_info_exon <- draw_default_annotations(
          bed_file       = exon_bed_file,
          depth          = exon_depth,
          alt_tss_file   = alt_tss_file,
          page_bottom_y  = legend_floor_y
        )
        genes_skipped <- legend_info_exon$genes_skipped
      }
    }

    # Draw pivots
    if (identical(genome_build, "hg38") && isTRUE(inherent) && identical(analysis_type, "single_sample")) {
      bed_enh <- file.path(data_dir, genome_build, "reference", "encode-3_pivots_hg38_enh_off.bed")
      bed_pro <- file.path(data_dir, genome_build, "reference", "encode-3_pivots_hg38_pro_on.bed")

      if (file.exists(bed_enh)) {
        draw_bed(bed_enh, "#40E0D0", 0.30, FALSE,
          ruler_top = legend_info$ruler$target_top,
          point_mode = TRUE, point_cex = 0.2
        )
        default_labels <- c(default_labels, "inherent off")
        default_colors <- c(default_colors, "#40E0D0")
        default_types <- c(default_types, "dot")
        default_ptcex <- c(default_ptcex, 0.5)
      }

      if (file.exists(bed_pro)) {
        draw_bed(bed_pro, "#EE5C42", 0.30, FALSE,
          ruler_top = legend_info$ruler$target_top,
          point_mode = TRUE, point_cex = 0.2
        )
        default_labels <- c(default_labels, "inherent on")
        default_colors <- c(default_colors, "#EE5C42")
        default_types <- c(default_types, "dot")
        default_ptcex <- c(default_ptcex, 0.5)
      }
    }
  }

  # Draw default legend
  if (!is.null(legend_info_exon) && is.finite(last_annotation_min_y)) {
    src <- list(
      legend_info_exon,
      list(labels = default_labels, colors = default_colors, types = default_types)
    )
    src <- Filter(function(s) length(s$labels) > 0, src)

    if (length(src) > 0) {
      bottom_labels <- truncate_labels(unlist(lapply(src, `[[`, "labels")))
      bottom_colors <- unlist(lapply(src, `[[`, "colors"))
      bottom_types <- unlist(lapply(src, `[[`, "types"))

      m <- compute_legend_markers(bottom_types, lwd_default = legend_line_lwd)

      legend(
        x = x_left - bbox_pad_x, y = y_bottom - bbox_pad_y,
        legend = bottom_labels, col = bottom_colors,
        lty = m$lty, lwd = m$lwd, pch = m$pch,
        pt.bg = m$pt_bg, pt.cex = m$pt_cex, pt.lwd = m$pt_lwd,
        bty = "o",
        bg = adjustcolor("white", alpha.f = 0.75),
        box.col = NA,
        box.lwd = 0.5,
        cex = 0.30,
        text.col = "#888888", xjust = 0, yjust = 0, xpd = NA
      )
    }
  }

  # Draw overflow notice if genes were skipped
  if (genes_skipped > 0L) {
    notice_x <- grconvertX(0.96, from = "ndc", to = "user")
    notice_y <- page_bottom_y + 1
    notice_w <- 22
    notice_h <- 4

    rect(
      xleft = notice_x - notice_w,
      ybottom = notice_y - notice_h / 2,
      xright = notice_x,
      ytop = notice_y + notice_h / 2,
      col = adjustcolor("white", alpha.f = 0.75),
      border = adjustcolor("#bbbbbb", alpha.f = 0.9),
      lwd = 0.5, xpd = NA
    )

    notice_text <- sprintf("%d gene(s) omitted\n(insufficient page space)", genes_skipped)
    text(notice_x - notice_w / 2, notice_y,
      labels = notice_text, cex = 0.30, col = "#888888", xpd = NA
    )
  }

  # Draw custom legend (top right)
  if (length(legend_labels) > 0) {
    legend_labels <- truncate_labels(legend_labels)

    m <- compute_legend_markers(legend_types, lwd_default = legend_line_lwd)

    leg_x <- grconvertX(0.96, from = "ndc", to = "user")
    leg_y <- grconvertY(0.97, from = "ndc", to = "user")

    legend(
      x = leg_x, y = leg_y,
      legend = legend_labels, col = legend_colors,
      lty = m$lty, lwd = m$lwd, pch = m$pch,
      pt.bg = m$pt_bg, pt.cex = m$pt_cex, pt.lwd = m$pt_lwd,
      bty = "n", cex = 0.30, text.col = "#888888",
      xjust = 1, yjust = 1, xpd = NA, y.intersp = 1.15
    )
  }
}

#################################################################
##                            Inter                            ##
#################################################################

calculate_w_inter <- function(plot_width, plot_height, A) {
  n_rows <- nrow(A)
  n_cols <- ncol(A)

  # Choose the smaller bin size to ensure cells fit
  w <- min(plot_width / n_cols, plot_height / n_rows)

  return(w)
}

calculate_plot_dimensions <- function(interval_start1, interval_end1,
                                      interval_start2, interval_end2,
                                      bin_size,
                                      max_width  = 100,
                                      max_height = 100) {
  bins1 <- ceiling((interval_end1 - interval_start1) / bin_size)
  bins2 <- ceiling((interval_end2 - interval_start2) / bin_size)

  aspect_ratio <- bins2 / bins1

  # Start by filling the width
  plot_width  <- max_width
  plot_height <- max_width * aspect_ratio

  # If that overflows the height, scale down to fit height instead
  if (plot_height > max_height) {
    plot_height <- max_height
    plot_width  <- max_height / aspect_ratio
  }

  return(list(width = plot_width, height = plot_height))
}

draw_ruler_inter <- function(left_side, top, w, A,
                             axis = "x",
                             gap = 1.0,
                             col_base = "#B0B0B0",
                             col_ticks = "#afafafff",
                             col_labels = "#6E6E6E",
                             cex_mb = 0.25,
                             cex_mid = 0.18) {
  bin_bp <- as.numeric(bin_size)
  n_rows <- nrow(A)
  n_cols <- ncol(A)

  # flat tips
  old_lend <- par("lend")
  on.exit(par(lend = old_lend), add = TRUE)
  par(lend = "butt")

  # Hierarchical tick lengths
  H <- gap
  len_mb <- H * 0.89
  len_500k <- H * 0.65
  len_100k <- H * 0.65
  len_50k <- H * 0.50
  len_10k <- H * 0.35

    # --- Tick visibility ---
  show_10kb  <- (bin_bp <= 5000)
  show_50kb  <- (bin_bp %in% c(1000, 5000, 10000))
  show_100kb <- (bin_bp %in% c(1000, 5000, 10000, 25000)) ||
                (bin_bp == 50000  && w >= 0.30) ||
                (bin_bp == 100000 && w >= 1.50)
  show_500kb <- (bin_bp %in% c(50000, 100000)) || (bin_bp == 5000 && w < 1)
  major_step_bp <- if (bin_bp >= 1e6) bin_bp else 1e6

  # --- Label suppression ---
  suppress_label_100kb <- (bin_bp == 5000  && w < 0.08) ||
                          (bin_bp == 25000 && w < 0.25)
  suppress_label_50kb  <- (bin_bp %in% c(1000, 5000) && w < 0.3) ||
                          (bin_bp == 10000 && w < 0.20)
  suppress_label_10kb  <- (bin_bp <= 5000 && w < 2)

  # Label offset from baseline
  label_offset_mb <- 0.15
  label_offset_mid <- 0.5

  # X-axis ruler
  if (axis == "x") {
    first_bp <- suppressWarnings(as.numeric(colnames(A)[1]))
    n_bins <- n_cols

    # Span of ruler
    x_left <- left_side
    x_right <- left_side + n_cols * w
    base_y <- top - n_rows * w # bottom edge of the contact rectangle

    # Baseline along bottom edge
    segments(x_left, base_y, x_right, base_y, col = col_base, lwd = 0.2, xpd = NA)

    # Draw ticks
    for (k in 1:n_bins) {
      xk <- left_side + (k - 1) * w
      coord_bp <- first_bp + (k - 1) * bin_bp

      at_mb <- (coord_bp %% major_step_bp == 0L)
      at_500kb <- show_500kb && (coord_bp %% 500000 == 0L)
      at_100kb <- show_100kb && (coord_bp %% 100000 == 0L)
      at_50kb <- show_50kb && (coord_bp %% 50000 == 0L)
      at_10kb <- show_10kb && (coord_bp %% 10000 == 0L)

      if (at_mb) {
        segments(xk, base_y, xk, base_y - len_mb, col = col_ticks, lwd = 0.25, xpd = NA)
        mb_lab <- if (major_step_bp == 1e6) {
          paste0(formatC(coord_bp / 1e6, digits = 0, format = "f"), " MB")
        } else {
          paste0(formatC(coord_bp / 1e6, digits = 1, format = "f"), " MB")
        }
        text(xk, base_y - len_mb - label_offset_mb,
          labels = mb_lab,
          cex = cex_mb, col = col_labels, adj = c(0.5, 1), xpd = NA
        )
      } else if (at_500kb) {
        segments(xk, base_y, xk, base_y - len_500k, col = col_ticks, lwd = 0.18, xpd = NA)
        lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
        text(xk, base_y - len_500k - label_offset_mid,
          labels = lab,
          cex = cex_mid, col = col_labels, adj = c(0.5, 1), xpd = NA
        )
      } else if (at_100kb) {
        segments(xk, base_y, xk, base_y - len_100k, col = col_ticks, lwd = 0.18, xpd = NA)
        if (!suppress_label_100kb) {
          lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
          text(xk, base_y - len_100k - label_offset_mid,
            labels = lab,
            cex = cex_mid, col = col_labels, adj = c(0.5, 1), xpd = NA
          )
        }
      } else if (at_50kb) {
        segments(xk, base_y, xk, base_y - len_50k, col = col_ticks, lwd = 0.16, xpd = NA)
        if (bin_bp %in% c(1000, 5000) && !suppress_label_50kb) {
          lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
          text(xk, base_y - len_50k - label_offset_mid,
            labels = lab,
            cex = cex_mid, col = col_labels, adj = c(0.5, 1), xpd = NA
          )
        }
      } else if (at_10kb) {
        segments(xk, base_y, xk, base_y - len_10k, col = col_ticks, lwd = 0.14, xpd = NA)
        if (bin_bp <= 5000 && !suppress_label_10kb) {
          lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
          text(xk, base_y - len_10k - label_offset_mid,
            labels = lab,
            cex = cex_mid, col = col_labels, adj = c(0.5, 1), xpd = NA
          )
        }
      }
    }

    # Y-axis ruler (vertical, left of the rectangle)
  } else if (axis == "y") {
    first_bp <- suppressWarnings(as.numeric(rownames(A)[1]))
    n_bins <- n_rows

    # Span of ruler
    y_top <- top
    y_bottom <- top - n_rows * w
    base_x <- left_side # left edge of the contact rectangle

    # Baseline along left edge
    segments(base_x, y_bottom, base_x, y_top, col = col_base, lwd = 0.2, xpd = NA)

    # Draw ticks (leftward from baseline)
    for (k in 1:n_bins) {
      yk <- top - (k - 1) * w
      coord_bp <- first_bp + (k - 1) * bin_bp

      at_mb <- (coord_bp %% major_step_bp == 0L)
      at_500kb <- show_500kb && (coord_bp %% 500000 == 0L)
      at_100kb <- show_100kb && (coord_bp %% 100000 == 0L)
      at_50kb <- show_50kb && (coord_bp %% 50000 == 0L)
      at_10kb <- show_10kb && (coord_bp %% 10000 == 0L)

      if (at_mb) {
        segments(base_x, yk, base_x - len_mb, yk, col = col_ticks, lwd = 0.25, xpd = NA)
        mb_lab <- if (major_step_bp == 1e6) {
          paste0(formatC(coord_bp / 1e6, digits = 0, format = "f"), " MB")
        } else {
          paste0(formatC(coord_bp / 1e6, digits = 1, format = "f"), " MB")
        }
        text(base_x - len_mb - label_offset_mb, yk,
          labels = mb_lab,
          cex = cex_mb, col = col_labels, adj = c(1, 0.5), xpd = NA
        )
      } else if (at_500kb) {
        segments(base_x, yk, base_x - len_500k, yk, col = col_ticks, lwd = 0.18, xpd = NA)
        lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
        text(base_x - len_500k - label_offset_mid, yk,
          labels = lab,
          cex = cex_mid, col = col_labels, adj = c(1, 0.5), xpd = NA
        )
      } else if (at_100kb) {
        segments(base_x, yk, base_x - len_100k, yk, col = col_ticks, lwd = 0.18, xpd = NA)
        if (!suppress_label_100kb) {
          lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
          text(base_x - len_100k - label_offset_mid, yk,
            labels = lab,
            cex = cex_mid, col = col_labels, adj = c(1, 0.5), xpd = NA
          )
        }
      } else if (at_50kb) {
        segments(base_x, yk, base_x - len_50k, yk, col = col_ticks, lwd = 0.16, xpd = NA)
        if (bin_bp %in% c(1000, 5000) && !suppress_label_50kb) {
          lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
          text(base_x - len_50k - label_offset_mid, yk,
            labels = lab,
            cex = cex_mid, col = col_labels, adj = c(1, 0.5), xpd = NA
          )
        }
      } else if (at_10kb) {
        segments(base_x, yk, base_x - len_10k, yk, col = col_ticks, lwd = 0.14, xpd = NA)
        if (bin_bp <= 5000 && !suppress_label_10kb) {
          lab <- prettyNum(coord_bp, big.mark = ",", scientific = FALSE, preserve.width = "none")
          text(base_x - len_10k - label_offset_mid, yk,
            labels = lab,
            cex = cex_mid, col = col_labels, adj = c(1, 0.5), xpd = NA
          )
        }
      }
    }
  }
    
  invisible(NULL)
}

draw_scale_inter <- function() {
  if (analysis_type == "single_sample") {
    break_steps <- seq(1, length(breakList), by = 5)
  } else {
    break_steps <- seq(1, length(breakList), by = 10)
  }

  wid <- 0.5
  n_boxes <- length(break_steps)

  # Position scale bar in top-right area using NDC coordinates
  x_right <- grconvertX(0.96, from = "ndc", to = "user")
  x_offset <- x_right - n_boxes * wid

  # Place vertically near the top of the page
  y_top_bar <- grconvertY(0.96, from = "ndc", to = "user")
  y_bottom_bar <- y_top_bar - wid

  for (j in seq_len(n_boxes)) {
    color <- rownames(C[order(abs(C[, "breaks"] - breakList[break_steps[j]])), ])[1]
    rect(x_offset + (j - 1) * wid, y_bottom_bar,
      x_offset + j * wid, y_top_bar,
      col = color, border = NA, xpd = NA
    )
  }

  # Labels above the bar
  label_y <- y_top_bar + wid

  text(x_offset, label_y,
    labels = format(round(breakList[1], 4)),
    cex = 0.4, col = "#888888", pos = 4, offset = 0, xpd = NA
  )
  text(x_offset + n_boxes * wid, label_y,
    labels = format(round(breakList[length(breakList)], 4), nsmall = 2),
    cex = 0.4, col = "#888888", pos = 2, offset = 0, xpd = NA
  )

  # Return the bounding box info
  return(list(
    x_left   = x_offset,
    x_right  = x_offset + n_boxes * wid,
    y_bottom = y_bottom_bar,
    y_top    = label_y
  ))
}

draw_feature_inter <- function(
  chr, up, do,
  gene_name,
  col,
  flag_text,
  w, top, axis,
  interval_start, interval_len,
  plot_width, plot_height,
  depth,
  gap = 1.0,
  gene_cex = 0.2,
  coord_cex = 0.1,
  left_side = 0
  ) {
  # Clip the start and end positions to the interval boundaries
  clipped_start <- max(up, interval_start)
  clipped_end <- min(do, interval_start + interval_len)

  genomic_region_label <- sprintf("%s:%d-%d", chr, clipped_start, clipped_end)

  # Convert genomic positions to plot coordinates, offset by left_side
  calculate_position <- function(pos) {
    return(left_side + (pos - interval_start) / bin_size * w)
  }

  str <- calculate_position(clipped_start)
  end <- calculate_position(clipped_end)

  # Ruler geometry
  n_rows <- nrow(A)
  n_cols <- ncol(A)
  H <- gap
  len_mb <- H * 0.89
  label_offset_mb <- 0.05

  base_x <- left_side # left edge of contact rectangle
  base_y <- top - n_rows * w # bottom edge of contact rectangle

  # Where the ruler's longest ticks end
  ruler_tick_end_x <- base_x - len_mb
  ruler_tick_end_y <- base_y - len_mb

  if (axis == "x") {
    annotation_pad <- 3
    if (flag_text && gene_name != "") {
      label <- sprintf("%s (%s)", gene_name, genomic_region_label)
      label_x <- (str + end) / 2
      label_y <- ruler_tick_end_y - label_offset_mb - annotation_pad

      # Dotted guide line from baseline down to annotation
      segments(label_x, base_y, label_x, label_y,
        col = "#d7d7d7", lwd = 0.4, lty = "dotted", xpd = NA
      )
      # Label below the ruler, rotated 45 degrees
      text(label_x, label_y, label,
        col = col, cex = gene_cex,
        srt = 45, adj = c(1, 1), xpd = NA
      )
    } else {
      bar_y <- base_y - depth
      lines(c(str, end), c(bar_y, bar_y), col = col, lwd = 0.5, xpd = NA)
    }
  } else if (axis == "y") {
    annotation_pad <- 5
    y_start <- top - (clipped_start - interval_start) / bin_size * w
    y_end <- top - (clipped_end - interval_start) / bin_size * w

    if (flag_text && gene_name != "") {
      label <- sprintf("%s (%s)", gene_name, genomic_region_label)
      label_x <- ruler_tick_end_x - label_offset_mb - annotation_pad
      label_y <- (y_start + y_end) / 2

      # Dotted guide line from baseline leftward to annotation
      segments(base_x, label_y, label_x - 2, label_y,
        col = "#d7d7d7", lwd = 0.4, lty = "dotted", xpd = NA
      )
      # Label to the left of the ruler
      text(label_x, label_y, label,
        col = col, cex = gene_cex,
        pos = 2, xpd = NA
      )
    } else {
      x_position <- base_x - depth
      lines(c(x_position, x_position), c(y_start, y_end), col = col, lwd = 0.5, xpd = NA)
    }
  }
}

draw_bed_inter <- function(
  bed_file,
  color,
  depth,
  flag_text,
  w, top, axis,
  interval_start1, interval_len1, interval_start2, interval_len2,
  plot_width, plot_height, left_side = 0
  ) {
  if (!file.exists(bed_file)) {
    cat("bed file not found\n")
    stop()
  }

  bed <- read.table(bed_file, as.is = TRUE, header = FALSE)

  # Assign column names based on the number of columns
  if (ncol(bed) == 3) {
    colnames(bed) <- c("chr", "start", "end")
  } else if (ncol(bed) >= 4) {
    colnames(bed) <- c("chr", "start", "end", "name")
  } else {
    cat("bed file must have at least 3 columns\n")
    stop()
  }

  if (nrow(bed) == 0) {
    cat("draw_bed_inter error: empty bed file\n")
    return()
  }

  # Check if there's a fourth column for text
  if (ncol(bed) == 3 && flag_text) {
    flag_text <- FALSE
  }

  # Determine which axis to process
  if (axis == "x") {
    bed_data <- bed[bed$chr == new_interval_chr1, ]
    axis_start <- interval_start1
    axis_len <- interval_len1
  } else if (axis == "y") {
    bed_data <- bed[bed$chr == new_interval_chr2, ]
    axis_start <- interval_start2
    axis_len <- interval_len2
  } else {
    stop("Invalid axis value: must be 'x' or 'y'")
  }

  process_axis <- function(bed_data, axis, interval_start, interval_len) {
    bed_data <- bed_data[bed_data$start < interval_start + interval_len & bed_data$end > interval_start, ]
    bed_data$bin <- floor(bed_data$start / bin_size)
    bin_counts <- table(bed_data$bin)

    for (bin in names(bin_counts)) {
      bin_genes <- bed_data[bed_data$bin == as.numeric(bin), ]

      for (i in seq_len(nrow(bin_genes))) {
        b <- bin_genes[i, ]

        # Non-gene annotations
        if (!flag_text) {
          draw_feature_inter(
            b$chr, b$start, b$end, "",
            col = color, flag_text = FALSE,
            w = w, top = top, axis = axis,
            interval_start = interval_start, interval_len = interval_len,
            plot_width = plot_width, plot_height = plot_height,
            depth = depth, left_side = left_side
          )
          next
        }

        # Gene annotations: evaluate BEDPE overlap
        gene_name <- if (ncol(bed) >= 4) b$name else ""

        if (exists("pairs") && nrow(pairs) > 0 &&
          is_gene_in_bedpe_inter(b$chr, b$start, b$end, axis, pairs, flip_axes, gene_name)) {
          gene_cex <- 0.3
          coord_cex <- 0.2
          colour <- "#ef0000"
        } else {
          gene_cex <- 0.3
          coord_cex <- 0.2
          colour <- "#aaaaaa"
        }

        draw_feature_inter(
          b$chr, b$start, b$end, gene_name,
          col = colour, flag_text = TRUE,
          w = w, top = top, axis = axis,
          interval_start = interval_start, interval_len = interval_len,
          plot_width = plot_width, plot_height = plot_height,
          depth = depth,
          gene_cex = gene_cex, coord_cex = coord_cex, left_side = left_side
        )
      }
    }
  }

  process_axis(bed_data, axis, axis_start, axis_len)
}

draw_contacts_inter <- function(A, C, left_side, bedpe, top, w) {
  n_rows <- nrow(A)
  n_cols <- ncol(A)
 
  # Draw contact matrix
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      color <- rownames(C[order(abs(C[, "breaks"] - A[i, j]), decreasing = FALSE), ])[1]
      rect(
        left_side + (j - 1) * w,
        top - (i - 1) * w,
        left_side + j * w,
        top - i * w,
        col = color,
        border = NA
      )
    }
  }
 
  # Rulers
  draw_ruler_inter(left_side, top, w, A, axis = "x")
  draw_ruler_inter(left_side, top, w, A, axis = "y")
 
  # Draw BEDPE highlights
  legend_out <- NULL
 
  if (!identical(bedpe, "FALSE")) {
    prep <- prep_bedpe_layers(bedpe, bedpe_color_param = bedpe_color)
    pairs_list  <- prep$pairs_list
    bedpe_colors <- prep$color_map
 
    for (nm in names(pairs_list)) {
      pairs <- pairs_list[[nm]]
      if (is.null(pairs) || nrow(pairs) == 0L) next
      bedpe_col <- unname(bedpe_colors[[nm]])
 
      for (p in seq_len(nrow(pairs))) {
        x1 <- left_side + (pairs$i[p] - 1) * w
        y1 <- top - (pairs$j[p] - 1) * w
        x2 <- left_side + (pairs$i_end[p]) * w
        y2 <- top - (pairs$j_end[p]) * w
 
        rect(
          x1, y1, x2, y2,
          border = bedpe_col,
          col = NA,
          lwd = 0.8
        )
      }
    }
 
    labs <- names(bedpe_colors)
    legend_out <- list(labels = labs, colors = unname(bedpe_colors))
  }
 
  return(invisible(list(legend = legend_out)))
}
 
is_gene_in_bedpe_inter <- function(gene_chr, gene_start, gene_end, axis, bedpe, flip_axes, gene_name = NULL) {
  # Handle list of data.frames (multi-file) or single data.frame
  if (is.data.frame(bedpe)) {
    pairs_list <- list(bedpe)
  } else if (is.list(bedpe) && !identical(bedpe, "FALSE")) {
    pairs_list <- bedpe
  } else {
    return(FALSE)
  }
 
  for (pairs in pairs_list) {
    if (!is.data.frame(pairs) || nrow(pairs) == 0) next
 
    for (p in seq_len(nrow(pairs))) {
      if (!flip_axes) {
        if (axis == "x" && gene_chr == new_interval_chr1) {
          region_start <- pairs$start1[p]
          region_end <- pairs$end1[p]
        } else if (axis == "y" && gene_chr == new_interval_chr2) {
          region_start <- pairs$start2[p]
          region_end <- pairs$end2[p]
        } else {
          next
        }
      } else {
        if (axis == "x" && gene_chr == new_interval_chr1) {
          region_start <- pairs$start2[p]
          region_end <- pairs$end2[p]
        } else if (axis == "y" && gene_chr == new_interval_chr2) {
          region_start <- pairs$start1[p]
          region_end <- pairs$end1[p]
        } else {
          next
        }
      }
 
      if (gene_start < region_end && gene_end > region_start) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

draw_title_inter <- function(
  interval_chr1, interval_start1, interval_end1,
  interval_chr2, interval_start2, interval_end2
  ) {
  label_family <- "DejaVu Sans Mono"

  # Position title using NDC -> user conversion, matching intra draw_title
  title_y <- grconvertY(1 - 0.06, from = "ndc", to = "user")
  title_x <- grconvertX(0.02, from = "ndc", to = "user")

  if (analysis_type == "single_sample") {
    sample_name <- basename(sample_dir)
    title <- sprintf(
      "Sample: %s\nRegions: %s:%d:%d - %s:%d:%d\nGenome Build: %s | Resolution: %s",
      sample_name,
      interval_chr1,
      interval_start1,
      interval_end1,
      interval_chr2,
      interval_start2,
      interval_end2,
      genome_build,
      bin_size
    )
  } else if (analysis_type == "two_sample") {
    sample_nameA <- basename(sample_dirA)
    sample_nameB <- basename(sample_dirB)
    title <- sprintf(
      "Samples: %s \n         and %s\nRegions: %s:%d:%d - %s:%d:%d\nGenome Build: %s | Resolution: %s",
      sample_nameA, sample_nameB,
      interval_chr1,
      interval_start1,
      interval_end1,
      interval_chr2,
      interval_start2,
      interval_end2,
      genome_build,
      bin_size
    )
  }

  # Draw title in user coordinates with xpd=NA to allow drawing outside plot region
  text(title_x, title_y, labels = title, col = "#888888",
    pos = 4, cex = 0.5, family = label_family, xpd = NA
  )
}

draw_default_annotations_inter <- function(
    bed_file,
    axis,                        # "x" or "y"
    A,                           # contact matrix
    w,                           # bin width in user coords
    top,                         # top-edge y of the contact rectangle
    left_side,                   # left-edge x of the contact rectangle
    interval_chr,                # chromosome for this axis
    interval_start,              # start bp for this axis
    interval_end,                # end bp for this axis
    bin_size,                    # resolution
    bedpe          = "FALSE",    # bedpe list/data.frame or "FALSE"
    pairs          = NULL,       # parsed pairs data.frame
    flip_axes      = FALSE,      # whether axes were flipped
    new_interval_chr1 = NULL,    # chr for axis 1 (needed for bedpe check)
    new_interval_chr2 = NULL,    # chr for axis 2 (needed for bedpe check)
    depth          = 2.5,        # gap between ruler baseline and first gene lane
    label_genes    = TRUE,
    merge_exons    = TRUE,
    show_arrows    = TRUE,
    pack_lanes     = TRUE,
    lane_gap       = 1.25,
    group_gap      = 1,
    alt_tss_file   = NULL,
    alt_tss_label  = TRUE,
    alt_tss_label_cex = 0.26,
    label_cex      = 0.26,
    label_pad_in   = 0.03,
    ann_style      = "standard",
    page_boundary  = -Inf        # min y for x-axis, min x for y-axis
  ) {

  # Constants

  exon_height  = 0.25
  base_lwd     = 0.4
  arrow_dx     = 0.20
  arrow_dy     = 0.14
  arrow_gap    = 0.3
  alt_tss_col  = "#FF5656"
  tss_half_w   = 0.25 * 0.3
  label_family = "DejaVu Sans Mono"

  interval_len <- interval_end - interval_start

  if (identical(ann_style, "axiotl")) {
    lane_gap <- lane_gap + 0.4
  }

  n_rows <- nrow(A)
  n_cols <- ncol(A)

  # --- Geometry anchors ---
  rect_bottom <- top - n_rows * w
  rect_left   <- left_side

  if (axis == "x") {
    ruler_baseline <- rect_bottom
  } else {
    ruler_baseline <- rect_left
  }

  # --- Coordinate mappers: bp -> user coords ---
  if (axis == "x") {
    map_pos <- function(bp) {
      left_side + (bp - interval_start) / bin_size * w
    }
    scale_factor <- (n_cols * w) / interval_len
  } else {
    map_pos <- function(bp) {
      top - (bp - interval_start) / bin_size * w
    }
    scale_factor <- (n_rows * w) / interval_len
  }

  # ---- BEDPE dimming helper ----
  is_gene_in_bedpe_for_axis <- function(gene_start, gene_end) {
    if (identical(bedpe, "FALSE")) return(FALSE)

    all_pairs <- if (is.data.frame(bedpe)) list(bedpe)
                 else if (is.list(bedpe)) bedpe
                 else return(FALSE)

    if (!is.null(pairs) && is.data.frame(pairs)) {
      all_pairs <- c(all_pairs, list(pairs))
    }

    for (p_df in all_pairs) {
      if (!is.data.frame(p_df) || nrow(p_df) == 0L) next
      for (p in seq_len(nrow(p_df))) {
        if (!flip_axes) {
          if (axis == "x" && !is.null(new_interval_chr1) && interval_chr == new_interval_chr1) {
            rs <- p_df$start1[p]; re <- p_df$end1[p]
          } else if (axis == "y" && !is.null(new_interval_chr2) && interval_chr == new_interval_chr2) {
            rs <- p_df$start2[p]; re <- p_df$end2[p]
          } else { next }
        } else {
          if (axis == "x" && !is.null(new_interval_chr1) && interval_chr == new_interval_chr1) {
            rs <- p_df$start2[p]; re <- p_df$end2[p]
          } else if (axis == "y" && !is.null(new_interval_chr2) && interval_chr == new_interval_chr2) {
            rs <- p_df$start1[p]; re <- p_df$end1[p]
          } else { next }
        }
        if (gene_start < re && gene_end > rs) return(TRUE)
      }
    }
    FALSE
  }

  # Load data
  bed <- load_and_subset_bed(bed_file, interval_chr, interval_start, interval_end)
  if (is.null(bed)) return(invisible(NULL))

  alt_tss_tbl <- load_and_subset_alt_tss(alt_tss_file, interval_chr,
                                          interval_start, interval_end)

  # Build gene spans and assign lanes
  span_tbl <- build_gene_spans(bed, alt_tss_tbl)
  if (is.null(span_tbl)) return(invisible(NULL))

  genes <- sort(unique(span_tbl$gene))

  lane_info <- assign_all_lanes(span_tbl, label_genes, label_cex, label_pad_in,
                                 scale_factor, lane_gap, group_gap)

  group_offsets  <- lane_info$group_offsets
  lane_by_gene   <- lane_info$lane_by_gene
  group_by_gene  <- lane_info$group_by_gene
  strand_by_gene <- lane_info$strand_by_gene

  # Drawing loop
  old_xpd <- par("xpd"); on.exit(par(xpd = old_xpd), add = TRUE); par(xpd = NA)

  min_extent_drawn <- Inf
  genes_skipped <- 0L

  for (g in genes) {
    df <- bed[bed$gene == g, , drop = FALSE]
    if (!nrow(df)) next

    g_label <- toupper(g)
    bt      <- df$biotype[1]
    grp     <- group_by_gene[[g]]
    lane_idx <- if (pack_lanes) lane_by_gene[[g]] else 1L

    gx1_bp <- min(df$start); gx2_bp <- max(df$end)
    if (!(is.finite(gx1_bp) && is.finite(gx2_bp)) || gx2_bp <= gx1_bp) next

    gcol_base <- map_biotype_col(bt)
    bcol <- gcol_base

    # BEDPE dimming
    dim_this <- is_gene_in_bedpe_for_axis(gx1_bp, gx2_bp)
    has_any_bedpe <- !identical(bedpe, "FALSE") && !is.null(pairs) && is.data.frame(pairs) && nrow(pairs) > 0L
    if (has_any_bedpe && !dim_this) {
      bcol <- grDevices::adjustcolor(gcol_base, alpha.f = 0.40)
    }

    tss_col <- grDevices::adjustcolor(alt_tss_col,
      alpha.f = if (has_any_bedpe && !dim_this) 0.40 else 1.00)

    # Lane offset
    class_offset <- group_offsets[[grp]]
    lane_offset  <- depth + class_offset + (lane_idx - 1L) * lane_gap

    # Primary coordinates
    p1 <- map_pos(gx1_bp)
    p2 <- map_pos(gx2_bp)

    # Prepare exon draw data using shared helper
    ex_draw <- prepare_exon_draw(df, gx1_bp, gx2_bp, g, strand_by_gene[[g]], bt, merge_exons)

    # Axis-specific drawing
    if (axis == "x") {
      y_gene       <- ruler_baseline - lane_offset
      y_label      <- y_gene - exon_height / 2 - 0.40
      y_tss_tip    <- y_gene + exon_height / 2 + 0.10
      y_tss_base   <- y_tss_tip + exon_height

      if (is.finite(page_boundary) && y_label < page_boundary) {
        genes_skipped <- genes_skipped + 1L
        next
      }

      # Gene baseline
      inset <- (p2 - p1) * 0.001
      segments(p1 + inset, y_gene, p2 - inset, y_gene,
               col = bcol, lwd = base_lwd, xpd = NA, lend = "butt")

      # Exons
      for (k in seq_len(nrow(ex_draw))) {
        xs <- map_pos(ex_draw$start[k]); xe <- map_pos(ex_draw$end[k])
        if (is.finite(xs) && is.finite(xe) && xe > xs) {
          rect(xs, y_gene - exon_height / 2, xe, y_gene + exon_height / 2,
               col = bcol, border = bcol, lwd = 0.1, xpd = NA)
        }
      }

      # Strand arrow
      if (show_arrows && grp != "mirna") {
        dir <- if (strand_by_gene[[g]] == "-") "left" else "right"
        tip_bp <- if (dir == "right") max(ex_draw$end) else min(ex_draw$start)
        sgn <- if (dir == "right") 1 else -1
        tip_x <- map_pos(tip_bp) + sgn * arrow_gap
        segments(tip_x - sgn * arrow_dx, y_gene - arrow_dy,
                 tip_x, y_gene, col = bcol, lwd = 0.4, xpd = NA)
        segments(tip_x - sgn * arrow_dx, y_gene + arrow_dy,
                 tip_x, y_gene, col = bcol, lwd = 0.4, xpd = NA)
      }

      # Alt TSS triangles
      if (!is.null(alt_tss_tbl)) {
        alt_g <- alt_tss_tbl[alt_tss_tbl$gene == g, , drop = FALSE]
        if (nrow(alt_g)) {
          xs_alt <- map_pos(alt_g$start); xs_alt <- xs_alt[is.finite(xs_alt)]
          if (length(xs_alt)) {
            segments(min(xs_alt) - tss_half_w * 0.6, y_tss_base,
                     max(xs_alt) + tss_half_w * 0.6, y_tss_base,
                     col = tss_col, lwd = base_lwd, xpd = NA)
            for (ii in seq_len(nrow(alt_g))) {
              bp_pos <- alt_g$start[ii]; if (!is.finite(bp_pos)) next
              x_tip <- map_pos(bp_pos)
              polygon(c(x_tip - tss_half_w, x_tip + tss_half_w, x_tip),
                      c(y_tss_base, y_tss_base, y_tss_tip),
                      col = tss_col, border = NA, xpd = NA)
            }
          }
        }
      }

      # Label
      if (label_genes && nzchar(g)) {
        text((p1 + p2) / 2, y_label, labels = g_label, col = bcol,
             cex = label_cex, xpd = NA, family = label_family)
        label_height <- strheight(g_label, cex = label_cex, units = "user")
        min_extent_drawn <- min(min_extent_drawn, y_label - label_height)
      } else {
        min_extent_drawn <- min(min_extent_drawn, y_gene - exon_height / 2)
      }

    } else {
      # -- Y-AXIS: genes left of the rectangle, vertical --
      x_gene       <- ruler_baseline - lane_offset
      x_label      <- x_gene - exon_height / 2 - 0.40
      x_tss_tip    <- x_gene + exon_height / 2 + 0.10
      x_tss_base   <- x_tss_tip + exon_height

      if (is.finite(page_boundary) && x_label < page_boundary) {
        genes_skipped <- genes_skipped + 1L
        next
      }

      y1 <- p1; y2 <- p2

      # Gene baseline (vertical)
      inset <- abs(y1 - y2) * 0.001
      segments(x_gene, min(y1, y2) + inset, x_gene, max(y1, y2) - inset,
               col = bcol, lwd = base_lwd, xpd = NA, lend = "butt")

      # Exons
      for (k in seq_len(nrow(ex_draw))) {
        ys <- map_pos(ex_draw$start[k]); ye <- map_pos(ex_draw$end[k])
        if (is.finite(ys) && is.finite(ye) && ys != ye) {
          rect(x_gene - exon_height / 2, min(ys, ye),
               x_gene + exon_height / 2, max(ys, ye),
               col = bcol, border = bcol, lwd = 0.1, xpd = NA)
        }
      }

      # Strand arrow (vertical)
      if (show_arrows && grp != "mirna") {
        dir <- if (strand_by_gene[[g]] == "-") "up" else "down"
        if (dir == "down") {
          tip_y <- map_pos(max(ex_draw$end)) - arrow_gap
          segments(x_gene - arrow_dy, tip_y + arrow_dx,
                   x_gene, tip_y, col = bcol, lwd = 0.4, xpd = NA)
          segments(x_gene + arrow_dy, tip_y + arrow_dx,
                   x_gene, tip_y, col = bcol, lwd = 0.4, xpd = NA)
        } else {
          tip_y <- map_pos(min(ex_draw$start)) + arrow_gap
          segments(x_gene - arrow_dy, tip_y - arrow_dx,
                   x_gene, tip_y, col = bcol, lwd = 0.4, xpd = NA)
          segments(x_gene + arrow_dy, tip_y - arrow_dx,
                   x_gene, tip_y, col = bcol, lwd = 0.4, xpd = NA)
        }
      }

      # Alt TSS triangles (point rightward)
      if (!is.null(alt_tss_tbl)) {
        alt_g <- alt_tss_tbl[alt_tss_tbl$gene == g, , drop = FALSE]
        if (nrow(alt_g)) {
          ys_alt <- map_pos(alt_g$start); ys_alt <- ys_alt[is.finite(ys_alt)]
          if (length(ys_alt)) {
            segments(x_tss_base, min(ys_alt) - tss_half_w * 0.6,
                     x_tss_base, max(ys_alt) + tss_half_w * 0.6,
                     col = tss_col, lwd = base_lwd, xpd = NA)
            for (ii in seq_len(nrow(alt_g))) {
              bp_pos <- alt_g$start[ii]; if (!is.finite(bp_pos)) next
              y_tip <- map_pos(bp_pos)
              polygon(c(x_tss_base, x_tss_base, x_tss_tip),
                      c(y_tip - tss_half_w, y_tip + tss_half_w, y_tip),
                      col = tss_col, border = NA, xpd = NA)
            }
          }
        }
      }

      # Label (rotated)
      if (label_genes && nzchar(g)) {
        text(x_label, (min(y1, y2) + max(y1, y2)) / 2, labels = g_label,
             col = bcol, cex = label_cex, srt = -90, xpd = NA, family = label_family)
        label_width <- strheight(g_label, cex = label_cex, units = "user")
        min_extent_drawn <- min(min_extent_drawn, x_label - label_width)
      } else {
        min_extent_drawn <- min(min_extent_drawn, x_gene - exon_height / 2)
      }
    }
  }

  # Return legend info using shared helper
  if (!is.finite(min_extent_drawn)) min_extent_drawn <- NA_real_

  invisible(build_annotation_legend_info(genes_skipped, min_extent_drawn))
}

render_inter_plot <- function(
    output_name,
    A_mat,
    C_tbl,
    breaks_vec,
    bedpe,
    bedpe_color,
    w,
    ann_style         = "standard",
    ann_custom        = "NONE",
    ann_custom_colors = "NONE",
    genome_build      = "hg38",
    bin_size,
    new_interval_chr1,
    interval_start1,
    interval_end1,
    interval_len1,
    new_interval_chr2,
    interval_start2,
    interval_end2,
    interval_len2,
    flip_axes         = FALSE,
    pairs             = NULL,
    flag_profiles     = FALSE,
    data_dir          = path.expand("~/lab-data")
  ) {

  # Constants
  legend_max_chars    <- 36
  bbox_gap            <- 1.0
  bbox_pad_x          <- 0.8
  bbox_pad_y          <- 0.4
  gap_after_custom    <- 1.2
  base_depth_custom_x <- 3.0
  base_depth_custom_y <- 4.0
  step_depth          <- 0.3
  inter_base_colors <- c("#E69F00", "#00829d", "#8622b7", "#4C8219")
  left_side_contact <- ann_reserve_left

  # Open device
  test_con <- tryCatch(file(output_name, open = "ab"), error = function(e) NULL)
  if (is.null(test_con)) {
    stop(sprintf(
      "Cannot write to '%s' -- the file may be open in another application. Please close it and try again.",
      output_name
    ))
  }
  close(test_con)

  grDevices::cairo_pdf(output_name, width = 8.5, height = 8.5)
  on.exit(dev.off(), add = TRUE)

  par(family = "DejaVu Sans Mono",
      omi = rep(0, 4), mai = rep(0, 4), bg = "#eeeeee")

  plot(NULL,
    xlim = c(0, 100), ylim = c(0, 100),
    xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", bty = "n", asp = 1
  )

  # Title
  draw_title_inter(
    new_interval_chr1, interval_start1, interval_end1,
    new_interval_chr2, interval_start2, interval_end2
  )

  # Scale bar
  scale_pos <- draw_scale_inter()

  # Draw contacts
  # Count expected custom legend entries to nudge top down if needed
  n_legend <- 0L
  if (!identical(ann_custom, "NONE") && nzchar(ann_custom)) {
    n_legend <- n_legend + length(strsplit(ann_custom, "[[:space:]]+")[[1]])
  }
  if (!identical(bedpe, "FALSE")) {
    if (is.list(bedpe)) n_legend <- n_legend + length(bedpe)
    else                n_legend <- n_legend + 1L
  }

  # If more than 6 custom legend items, nudge top down to make room
  legend_overflow <- max(0L, n_legend - 6L)
  top <- 90 - legend_overflow * 1.5

  draw_contacts_inter(A_mat, C_tbl, left_side = left_side_contact, bedpe, top, w)

  # Geometry anchors
  n_rows <- nrow(A_mat)
  n_cols <- ncol(A_mat)

  rect_left   <- left_side_contact
  rect_right  <- left_side_contact + n_cols * w
  rect_bottom <- top - n_rows * w
  rect_top    <- top

  page_boundary_x <- grconvertY(0.02, from = "ndc", to = "user")
  page_boundary_y <- grconvertX(0.02, from = "ndc", to = "user")

  # Legend accumulators
  # Default legend (bottom-left, anchored to default annotation bbox)
  default_legend_labels <- c()
  default_legend_colors <- c()
  default_legend_types  <- c()

  # Custom legend (upper-right, below scale bar)
  custom_legend_labels <- c()
  custom_legend_colors <- c()
  custom_legend_types  <- c()

  legend_info_x     <- NULL
  legend_info_y     <- NULL
  x_def_bbox_bottom <- NULL
  y_def_bbox_left   <- NULL

  # Custom annotations
  max_custom_x_depth <- 0
  max_custom_y_depth <- 0

  custom_x_depths         <- c()
  custom_y_depths         <- c()
  custom_annotation_files <- c()
  custom_final_colors     <- c()

  if (!identical(ann_custom, "NONE") && nzchar(ann_custom)) {
    cat("\nDrawing custom annotations for", genome_build, "\n")
    custom_annotation_files <- strsplit(ann_custom, "[[:space:]]+")[[1]]
    custom_annotation_files <- custom_annotation_files[nzchar(custom_annotation_files)]
    n_files <- length(custom_annotation_files)

    custom_final_colors <- resolve_custom_colors(
      ann_custom_colors, n_files, base_colors = inter_base_colors
    )

    # First pass: draw custom tracks
    x_depth <- base_depth_custom_x
    y_depth <- base_depth_custom_y

    for (i in seq_along(custom_annotation_files)) {
      current_file <- custom_annotation_files[i]
      this_color   <- custom_final_colors[i]

      bed_tmp <- read.table(current_file, as.is = TRUE, header = FALSE)
      has_x <- any(bed_tmp[, 1] == new_interval_chr1 &
                    bed_tmp[, 2] < interval_end1 &
                    bed_tmp[, 3] > interval_start1)
      has_y <- any(bed_tmp[, 1] == new_interval_chr2 &
                    bed_tmp[, 2] < interval_end2 &
                    bed_tmp[, 3] > interval_start2)

      if (has_x) {
        custom_x_depths <- c(custom_x_depths, x_depth)
        draw_bed_inter(
          current_file, this_color, x_depth, FALSE, w, top, "x",
          interval_start1, interval_len1, interval_start2, interval_len2,
          plot_width = NULL, plot_height = NULL, left_side_contact
        )
        x_depth <- x_depth + step_depth
      }

      if (has_y) {
        custom_y_depths <- c(custom_y_depths, y_depth)
        draw_bed_inter(
          current_file, this_color, y_depth, FALSE, w, top, "y",
          interval_start1, interval_len1, interval_start2, interval_len2,
          plot_width = NULL, plot_height = NULL, left_side_contact
        )
        y_depth <- y_depth + step_depth
      }
    }

    if (length(custom_x_depths) > 0) max_custom_x_depth <- max(custom_x_depths)
    if (length(custom_y_depths) > 0) max_custom_y_depth <- max(custom_y_depths)

    # White-box + redraw: X-axis custom
    if (length(custom_x_depths) > 0) {
      cust_x_bbox_top    <- rect_bottom - min(custom_x_depths) + 0.5
      cust_x_bbox_bottom <- rect_bottom - max(custom_x_depths) - 0.5

      rect(xleft = rect_left - bbox_pad_x, ybottom = cust_x_bbox_bottom,
           xright = rect_right + bbox_pad_x, ytop = cust_x_bbox_top,
           col = "white", border = NA, xpd = NA)

      x_depth_redraw <- base_depth_custom_x
      for (i in seq_along(custom_annotation_files)) {
        bed_tmp <- read.table(custom_annotation_files[i], as.is = TRUE, header = FALSE)
        if (any(bed_tmp[, 1] == new_interval_chr1 &
                bed_tmp[, 2] < interval_end1 &
                bed_tmp[, 3] > interval_start1)) {
          draw_bed_inter(
            custom_annotation_files[i], custom_final_colors[i],
            x_depth_redraw, FALSE, w, top, "x",
            interval_start1, interval_len1, interval_start2, interval_len2,
            plot_width = NULL, plot_height = NULL, left_side_contact
          )
          x_depth_redraw <- x_depth_redraw + step_depth
        }
      }
    }

    # White-box + redraw: Y-axis custom
    if (length(custom_y_depths) > 0) {
      cust_y_bbox_right <- rect_left - min(custom_y_depths) + 0.5
      cust_y_bbox_left  <- rect_left - max(custom_y_depths) - 0.5

      rect(xleft = cust_y_bbox_left, ybottom = rect_bottom - bbox_pad_y,
           xright = cust_y_bbox_right, ytop = rect_top + bbox_pad_y,
           col = "white", border = NA, xpd = NA)

      y_depth_redraw <- base_depth_custom_y
      for (i in seq_along(custom_annotation_files)) {
        bed_tmp <- read.table(custom_annotation_files[i], as.is = TRUE, header = FALSE)
        if (any(bed_tmp[, 1] == new_interval_chr2 &
                bed_tmp[, 2] < interval_end2 &
                bed_tmp[, 3] > interval_start2)) {
          draw_bed_inter(
            custom_annotation_files[i], custom_final_colors[i],
            y_depth_redraw, FALSE, w, top, "y",
            interval_start1, interval_len1, interval_start2, interval_len2,
            plot_width = NULL, plot_height = NULL, left_side_contact
          )
          y_depth_redraw <- y_depth_redraw + step_depth
        }
      }
    }

    # Add custom tracks to legend
    for (i in seq_along(custom_annotation_files)) {
      custom_legend_labels <- c(custom_legend_labels, basename(custom_annotation_files[i]))
      custom_legend_colors <- c(custom_legend_colors, custom_final_colors[i])
      custom_legend_types  <- c(custom_legend_types, "line")
    }
  }

  # Default annotations
  if (ann_style != "none") {
    exon_bed_file <- switch(genome_build,
      "hg38" = file.path(data_dir, genome_build, "reference", "gencode.hg38.annotation.bed"),
      "hg19" = file.path(data_dir, genome_build, "reference", "gencode.hg19.annotation.bed"),
      "mm10" = file.path(data_dir, genome_build, "reference", "gencode.mm10.annotation.bed"),
      NULL
    )

    alt_tss_file <- if (identical(ann_style, "axiotl")) {
      switch(genome_build,
        "hg38" = file.path(data_dir, genome_build, "reference", "gencode.hg38.transcript_tss.bed"),
        "hg19" = file.path(data_dir, genome_build, "reference", "gencode.hg19.transcript_tss.bed"),
        "mm10" = file.path(data_dir, genome_build, "reference", "gencode.mm10.transcript_tss.bed"),
        NULL
      )
    } else {
      NULL
    }

    if (!is.null(exon_bed_file) && file.exists(exon_bed_file)) {
      cat("\nDrawing default annotations for", genome_build, "\n")

      # Compute default annotation depths
      def_depth_x <- if (max_custom_x_depth > 0) {
        max_custom_x_depth + gap_after_custom
      } else {
        3
      }
      def_depth_y <- if (max_custom_y_depth > 0) {
        max_custom_y_depth + gap_after_custom
      } else {
        base_depth_custom_y
      }

      left_side_ann <- rect_left - 1

      # Shared args for draw_default_annotations_inter
      common_args <- list(
        bed_file          = exon_bed_file,
        A                 = A_mat,
        w                 = w,
        top               = top,
        bin_size          = bin_size,
        bedpe             = bedpe,
        pairs             = pairs,
        flip_axes         = flip_axes,
        new_interval_chr1 = new_interval_chr1,
        new_interval_chr2 = new_interval_chr2,
        ann_style         = ann_style,
        alt_tss_file      = alt_tss_file
      )

      x_args <- c(common_args, list(
        axis           = "x",
        left_side      = left_side_contact,
        interval_chr   = new_interval_chr1,
        interval_start = interval_start1,
        interval_end   = interval_end1,
        depth          = def_depth_x,
        page_boundary  = page_boundary_x
      ))

      y_args <- c(common_args, list(
        axis           = "y",
        left_side      = left_side_ann,
        interval_chr   = new_interval_chr2,
        interval_start = interval_start2,
        interval_end   = interval_end2,
        depth          = def_depth_y,
        page_boundary  = page_boundary_y
      ))

      # First pass
      legend_info_x <- do.call(draw_default_annotations_inter, x_args)
      legend_info_y <- do.call(draw_default_annotations_inter, y_args)

      # ---- X-axis default white box + redraw ----
      if (max_custom_x_depth > 0) {
        x_def_bbox_top <- rect_bottom - max_custom_x_depth - gap_after_custom
      } else {
        x_def_bbox_top <- rect_bottom - def_depth_x + 0.5
      }

      if (!is.null(legend_info_x$min_extent) && is.finite(legend_info_x$min_extent)) {
        x_def_bbox_bottom <- legend_info_x$min_extent - 1.5
      } else {
        x_def_bbox_bottom <- x_def_bbox_top - 8
      }

      rect(xleft = rect_left - bbox_pad_x, ybottom = x_def_bbox_bottom,
           xright = rect_right + bbox_pad_x, ytop = x_def_bbox_top,
           col = "white", border = NA, xpd = NA)

      legend_info_x <- do.call(draw_default_annotations_inter, x_args)

      # ---- Y-axis default white box + redraw ----
      if (max_custom_y_depth > 0) {
        y_def_bbox_right <- rect_left - max_custom_y_depth - gap_after_custom
      } else {
        y_def_bbox_right <- left_side_ann - def_depth_y + 0.5
      }

      if (!is.null(legend_info_y$min_extent) && is.finite(legend_info_y$min_extent)) {
        y_def_bbox_left <- legend_info_y$min_extent - 1.5
      } else {
        y_def_bbox_left <- y_def_bbox_right - 8
      }

      rect(xleft = y_def_bbox_left, ybottom = rect_bottom - bbox_pad_y,
           xright = y_def_bbox_right, ytop = rect_top + bbox_pad_y,
           col = "white", border = NA, xpd = NA)

      legend_info_y <- do.call(draw_default_annotations_inter, y_args)

      # Default legend entries (types match build_annotation_legend_info vocabulary)
      default_legend_labels <- c(default_legend_labels, "protein coding", "lncRNA", "miRNA")
      default_legend_colors <- c(default_legend_colors, "#4DAF4A", "#C77CFF", "#377EB8")
      default_legend_types  <- c(default_legend_types, "exon", "lncrna", "mirna")

      if (identical(ann_style, "axiotl") && !is.null(alt_tss_file)) {
        default_legend_labels <- c(default_legend_labels, "TSSs")
        default_legend_colors <- c(default_legend_colors, "#FF5656")
        default_legend_types  <- c(default_legend_types, "tss")
      }

    } else {
      warning("Exon BED not found for genome_build = ", genome_build,
              "; skipping default annotations.")
    }
  }

  # BEDPE legend entry
  if (!identical(bedpe, "FALSE")) {
    prep <- prep_bedpe_layers(bedpe, bedpe_color_param = bedpe_color)

    for (nm in names(prep$color_map)) {
      custom_legend_labels <- c(custom_legend_labels, nm)
      custom_legend_colors <- c(custom_legend_colors, unname(prep$color_map[[nm]]))
      custom_legend_types  <- c(custom_legend_types, "bedpe")
    }
  }

  # Default legend (bottom-left, anchored to annotation bbox)
  if (length(default_legend_labels) > 0) {
    bl_x <- rect_left
    bl_y <- if (!is.null(x_def_bbox_bottom)) {
      x_def_bbox_bottom
    } else {
      grconvertY(0.04, from = "ndc", to = "user")
    }

    m <- compute_legend_markers(default_legend_types)

    legend(
      x = bl_x, y = bl_y,
      legend = truncate_labels(default_legend_labels),
      col = default_legend_colors,
      lty = m$lty, lwd = m$lwd, pch = m$pch,
      pt.bg = m$pt_bg, pt.cex = m$pt_cex, pt.lwd = m$pt_lwd,
      bty = "o",
      bg = adjustcolor("white", alpha.f = 0.75),
      box.col = NA, box.lwd = 0.5,
      cex = 0.30,
      text.col = "#888888",
      xjust = 0, yjust = 0, xpd = NA
    )
  }

  # Custom legend (upper-right, below scale bar)
  if (length(custom_legend_labels) > 0) {
    leg_x <- grconvertX(0.96, from = "ndc", to = "user")
    leg_y <- scale_pos$y_bottom - 2

    m <- compute_legend_markers(custom_legend_types)

    legend(
      x = leg_x, y = leg_y,
      legend = truncate_labels(custom_legend_labels),
      col = custom_legend_colors,
      lty = m$lty, lwd = m$lwd, pch = m$pch,
      pt.bg = m$pt_bg, pt.cex = m$pt_cex, pt.lwd = m$pt_lwd,
      bty = "n",
      cex = 0.30,
      text.col = "#888888",
      xjust = 1, yjust = 1, xpd = NA, y.intersp = 1.15
    )
  }

  # Overflow notices
  x_genes_skipped <- if (!is.null(legend_info_x)) legend_info_x$genes_skipped else 0L
  y_genes_skipped <- if (!is.null(legend_info_y)) legend_info_y$genes_skipped else 0L

  if (x_genes_skipped > 0L && !is.null(x_def_bbox_bottom)) {
    notice_w <- 22; notice_h <- 4
    notice_x_right  <- rect_right + bbox_pad_x
    notice_y_center <- x_def_bbox_bottom + notice_h / 2 + 0.5

    rect(xleft = notice_x_right - notice_w,
         ybottom = notice_y_center - notice_h / 2,
         xright = notice_x_right,
         ytop = notice_y_center + notice_h / 2,
         col = adjustcolor("white", alpha.f = 0.75),
         border = adjustcolor("#bbbbbb", alpha.f = 0.9),
         lwd = 0.5, xpd = NA)

    text(notice_x_right - notice_w / 2, notice_y_center,
      labels = sprintf("%d gene(s) omitted on x-axis\n(insufficient page space)",
                       x_genes_skipped),
      cex = 0.30, col = "#888888", xpd = NA)
  }

  if (y_genes_skipped > 0L && !is.null(y_def_bbox_left)) {
    notice_w <- 22; notice_h <- 4
    notice_x_center <- y_def_bbox_left + notice_h / 2 + 0.5
    notice_y_bottom <- rect_bottom - bbox_pad_y

    rect(xleft = notice_x_center - notice_h / 2,
         ybottom = notice_y_bottom,
         xright = notice_x_center + notice_h / 2,
         ytop = notice_y_bottom + notice_w,
         col = adjustcolor("white", alpha.f = 0.75),
         border = adjustcolor("#bbbbbb", alpha.f = 0.9),
         lwd = 0.5, xpd = NA)

    text(notice_x_center, notice_y_bottom + notice_w / 2,
      labels = sprintf("%d gene(s) omitted on y-axis\n(insufficient page space)",
                       y_genes_skipped),
      cex = 0.30, col = "#888888", srt = -90, xpd = NA)
  }

  invisible(NULL)
}

#################################################################
##                            Shared                           ##
#################################################################

calculate_cap <- function(
  denMat,
  max_cap,
  quant_cut,
  flag_inter = "FALSE"
  ) {
  calculate_skewness <- function(values) {
    n <- length(values)
    if (n < 2) {
      return(0)
    }
    mean_val <- mean(values)
    sd_val <- sd(values)
    if (sd_val == 0) {
      return(0)
    }
    return((n / ((n - 1) * (n - 2))) * sum(((values - mean_val) / sd_val)^3))
  }

  if (all(denMat == 0)) {
    cat(sprintf("\nNo contact values found for this interval at resolution %d\n", bin_size))
    quit(save = "no", status = 1)
  }

  # Check for 'max_cap' and 'quant_cut' conditions
  if (max_cap != "none") {
    numeric_max_cap <- as.numeric(max_cap)

    if (numeric_max_cap < 0) {
      cat("\nMax_cap values cannot be negative. Switching to its absolute value.\n\n")
      return(abs(numeric_max_cap))
    } else {
      return(numeric_max_cap)
    }
  } else if (quant_cut != 1.00) {
    non_zero_values <- denMat[denMat != 0]
    quantile_value <- quantile(non_zero_values, probs = c(quant_cut), na.rm = TRUE)
    calculated_cap <- max(quantile_value)

    if (calculated_cap < 0) {
      cat("\nThe calculated cap value based on quant_cut is negative. Switching to its absolute value.\n\n")
      return(abs(calculated_cap))
    } else {
      return(calculated_cap)
    }
  } else {
    # Check for interchromosomal plot for no zero diag
    if (isFALSE(flag_inter)) {
      modified_mat <- zero_diag(denMat, 4)
      non_zero_values <- modified_mat[modified_mat != 0]
    } else {
      non_zero_values <- denMat[denMat != 0]
    }

    if (length(non_zero_values) == 0) {
      cat("\nNot enough data to plot. Please decrease the resolution or increase the range\n")
      quit(save = "no", status = 1)
    }

    # Calculate skewness
    skewness <- calculate_skewness(non_zero_values)

    # Define thresholds for skewness to apply scaling factors
    skewness_threshold1 <- 1 # For moderate skewness
    skewness_threshold2 <- 2 # For high skewness

    # Determine the scaling factor based on skewness
    scaling_factor <- ifelse(abs(skewness) > skewness_threshold2, 3,
      ifelse(abs(skewness) > skewness_threshold1, 2, 1)
    )

    # Calculate the 90th percentile of the absolute values
    high_percentile <- quantile(abs(non_zero_values), probs = 0.90, na.rm = TRUE)

    # Apply scaling factor to high percentile
    cap_value <- high_percentile * scaling_factor

    return(abs(cap_value))
  }
}

truncate_labels <- function(x, max_chars = 36L, use_unicode_ellipsis = TRUE) {
  if (is.null(x) || !length(x)) return(x)
  ell <- if (use_unicode_ellipsis) "\u2026" else "..."
  nc_widths <- nchar(x, type = "width", allowNA = FALSE, keepNA = FALSE)
  idx <- which(nc_widths > max_chars)
  if (length(idx)) {
    x[idx] <- paste0(substr(x[idx], 1L, max_chars - nchar(ell)), ell)
  }
  x
}

prep_bedpe_layers <- function(bedpe, bedpe_color_param = "#000000") {
  # Case 1: No highlights
  if (identical(bedpe, "FALSE")) {
    return(list(pairs_list = list(), color_map = character(0)))
  }
  # Case 2: highlights
  if (is.list(bedpe)) {
    pairs_list <- bedpe
    if (is.null(names(pairs_list)) || any(!nzchar(names(pairs_list)))) {
      names(pairs_list) <- paste0("bedpe_", seq_along(pairs_list))
    }
  } else {
    stop("`bedpe` must be 'FALSE' or a list of data.frames.")
  }

  n <- length(pairs_list)
  nm <- names(pairs_list)

  # Parse user color string
  toks <- unlist(strsplit(bedpe_color_param, "[[:space:]]+"))
  toks <- toks[nzchar(toks)]
  toks <- ifelse(substr(toks, 1, 1) == "#", toupper(toks), toupper(paste0("#", toks)))

  # Defaults: 1st black, 2nd pink; rest automatic
  first_two <- c("#000000", "#FF69B4")
  if (length(toks) == 0L || (length(toks) == 1L && toks[1] == "#000000")) {
    # Only default black provided, use default sequence
    if (n <= 2L) {
      cols <- first_two[seq_len(n)]
    } else {
      fill_n <- n - 2L
      fill <- grDevices::hsv(seq(0, 1, length.out = fill_n + 1)[-1], s = 0.65, v = 0.9)
      cols <- c(first_two, fill)
    }
  } else {
    # User provided colors first, fill remainder
    if (length(toks) >= n) {
      cols <- toks[seq_len(n)]
    } else {
      need <- n - length(toks)
      seed <- c("#FF69B4", "#377EB8", "#E69F00", "#4DAF4A", "#C77CFF", "#984EA3")
      fill <- if (need <= length(seed)) seed[seq_len(need)] else grDevices::colorRampPalette(seed)(need)
      cols <- c(toks, fill)
    }
  }

  names(cols) <- nm
  list(pairs_list = pairs_list, color_map = cols)
}

resolve_custom_colors <- function(ann_custom_colors, n_files,
                                  base_colors = c("#E69F00", "#377EB8",
                                                  "#4DAF4A", "#C77CFF",
                                                  "#984EA3")) {
  if (identical(ann_custom_colors, "NONE")) {
    return(colorRampPalette(base_colors)(n_files))
  }
  cols <- strsplit(ann_custom_colors, "[[:space:]]+")[[1]]
  cols <- cols[nzchar(cols)]
  cols <- sapply(cols, function(x) {
    if (substr(x, 1, 1) == "#") x else paste0("#", x)
  }, USE.NAMES = FALSE)
  if (length(cols) >= n_files) {
    cols[seq_len(n_files)]
  } else {
    fill <- colorRampPalette(base_colors)(n_files - length(cols))
    c(cols, fill)
  }
}

compute_legend_markers <- function(types, lwd_default = 0.8) {
  is_line <- types %in% c("coding", "lncrna", "mirna", "line", "exon")
  list(
    pch    = ifelse(types == "bedpe", 0,
              ifelse(types == "dot", 16,
                ifelse(types == "tss", 25, NA))),
    pt_bg  = ifelse(types == "dot", NA,
              ifelse(types == "tss", "#FF5656", NA)),
    pt_cex = ifelse(types == "bedpe", 0.6,
              ifelse(types == "dot", 0.5,
                ifelse(types == "tss", 0.5, NA))),
    pt_lwd = ifelse(types == "bedpe", 0.8,
              ifelse(types == "tss", NA, NA)),
    lty    = ifelse(is_line, 1, NA),
    lwd    = ifelse(is_line, lwd_default, NA)
  )
}

validate_intervals <- function(df) {
  df[!is.na(df$start) & !is.na(df$end) &
    is.finite(df$start) & is.finite(df$end) &
    (df$end > df$start), , drop = FALSE]
}

map_biotype_col <- function(bt) {
  key <- tolower(gsub("_", ".", bt, fixed = TRUE))
  pal <- c(
    "protein.coding" = "#4DAF4A",
    "lncrna"         = "#C77CFF",
    "mirna"          = "#377EB8"
  )
  if (!is.na(pal[key])) pal[key] else "#888888"
}

merge_intervals <- function(df) {
  df <- df[, c("chr", "start", "end", "gene", "strand", "biotype")]
  df$start <- suppressWarnings(as.numeric(df$start))
  df$end   <- suppressWarnings(as.numeric(df$end))
  df <- validate_intervals(df)
  if (!nrow(df)) return(df[0, ])

  df <- df[order(df$start, df$end), , drop = FALSE]
  s <- df$start[1]; e <- df$end[1]
  out <- vector("list", max(1L, nrow(df))); j <- 0L

  if (nrow(df) > 1L) {
    for (i in 2:nrow(df)) {
      si <- df$start[i]; ei <- df$end[i]
      if (si <= e) { e <- max(e, ei) }
      else { j <- j + 1L; out[[j]] <- c(s, e); s <- si; e <- ei }
    }
  }
  j <- j + 1L; out[[j]] <- c(s, e)
  m <- do.call(rbind, out[seq_len(j)])
  colnames(m) <- c("start", "end")

  data.frame(
    chr = df$chr[1], start = m[, 1], end = m[, 2],
    gene = df$gene[1], strand = df$strand[1], biotype = df$biotype[1],
    stringsAsFactors = FALSE
  )
}

assign_lanes <- function(df, label_cex = NULL, label_pad_in = NULL,
                         scale_factor = NULL) {
  if (!nrow(df)) return(df)

  use_labels <- !is.null(label_cex) && !is.null(label_pad_in) &&
                !is.null(scale_factor)

  if (use_labels) {
    # Widen packing envelope to include label footprint
    mid_bp <- (df$start + df$end) / 2
    w_in   <- strwidth(df$gene, cex = label_cex, units = "inches")

    inch_to_bp <- function(dx_in) {
      dx_user <- diff(grconvertX(c(0, dx_in), from = "in", to = "user"))
      dx_user / scale_factor
    }

    half_bp <- vapply(w_in / 2 + label_pad_in, inch_to_bp, numeric(1))
    lab_s   <- mid_bp - half_bp
    lab_e   <- mid_bp + half_bp
    pack_s  <- pmin(df$start, lab_s)
    pack_e  <- pmax(df$end,   lab_e)
  } else {
    pack_s <- df$start
    pack_e <- df$end
  }

  ord <- order(pack_s, pack_e)
  lane <- integer(nrow(df)); last_end <- numeric(0)

  for (jj in seq_along(ord)) {
    i <- ord[jj]; placed <- FALSE
    for (L in seq_along(last_end)) {
      if (pack_s[i] > last_end[L]) {
        lane[i] <- L; last_end[L] <- pack_e[i]; placed <- TRUE; break
      }
    }
    if (!placed) { last_end <- c(last_end, pack_e[i]); lane[i] <- length(last_end) }
  }
  df$lane <- lane
  df
}

load_and_subset_bed <- function(bed_file, interval_chr, interval_start, interval_end) {
  bed <- read.table(
    bed_file,
    header = FALSE, sep = "\t", quote = "", comment.char = "",
    stringsAsFactors = FALSE,
    colClasses = c("character", "numeric", "numeric", "character", "character", "character")
  )
  if (ncol(bed) != 6L) {
    stop(sprintf("Expected 6 BED columns; got %d", ncol(bed)))
  }
  colnames(bed) <- c("chr", "start", "end", "strand", "gene", "biotype")

  in_view <- (bed$chr == interval_chr) &
    (bed$start < interval_end) &
    (bed$end > interval_start)
  bed <- bed[in_view, , drop = FALSE]
  if (!nrow(bed)) return(NULL)

  bed$start <- pmax(bed$start, interval_start)
  bed$end   <- pmin(bed$end, interval_end)
  bed <- validate_intervals(bed)
  if (!nrow(bed)) return(NULL)
  unique(bed)
}

load_and_subset_alt_tss <- function(alt_tss_file, interval_chr,
                                     interval_start, interval_end) {
  if (is.null(alt_tss_file) || !file.exists(alt_tss_file)) return(NULL)

  alt_raw <- read.table(
    alt_tss_file, header = FALSE, sep = "\t",
    quote = "", comment.char = "", stringsAsFactors = FALSE
  )
  if (ncol(alt_raw) < 7L) {
    warning("alt_tss_file must have >= 7 columns; skipping")
    return(NULL)
  }
  colnames(alt_raw)[1:7] <- c("chr", "start", "end", "gene", "biotype", "strand", "level")
  alt_raw$level <- suppressWarnings(as.integer(alt_raw$level))

  alt_tss_tbl <- alt_raw[
    alt_raw$chr == interval_chr &
      alt_raw$start < interval_end &
      alt_raw$end > interval_start, , drop = FALSE
  ]
  if (!nrow(alt_tss_tbl)) return(NULL)

  alt_tss_tbl$start <- pmax(alt_tss_tbl$start, interval_start)
  alt_tss_tbl$end   <- pmin(alt_tss_tbl$end, interval_end)
  alt_tss_tbl <- validate_intervals(alt_tss_tbl)
  if (!nrow(alt_tss_tbl)) NULL else alt_tss_tbl
}

build_gene_spans <- function(bed, alt_tss_tbl = NULL) {
  genes <- sort(unique(bed$gene))
  if (!length(genes)) return(NULL)

  span_tbl <- do.call(rbind, lapply(genes, function(g) {
    df <- bed[bed$gene == g, , drop = FALSE]
    s <- min(df$start); e <- max(df$end)
    if (!is.null(alt_tss_tbl)) {
      alt_g <- alt_tss_tbl[alt_tss_tbl$gene == g, , drop = FALSE]
      if (nrow(alt_g)) {
        s <- min(s, min(alt_g$start, na.rm = TRUE))
        e <- max(e, max(alt_g$end, na.rm = TRUE))
      }
    }
    if (!is.finite(s) || !is.finite(e) || e <= s) return(NULL)
    bt  <- df$biotype[1]
    st  <- if (all(df$strand == df$strand[1])) df$strand[1] else "+"

    # Inline classify_group
    k <- tolower(gsub("_", ".", bt, fixed = TRUE))
    grp <- if (k == "mirna") "mirna"
           else if (k == "lncrna") "lncrna"
           else if (k == "protein.coding") "coding"
           else NA_character_

    data.frame(gene = g, start = s, end = e, strand = st,
               biotype = bt, group = grp, stringsAsFactors = FALSE)
  }))

  if (is.null(span_tbl) || !nrow(span_tbl)) return(NULL)
  span_tbl <- span_tbl[!is.na(span_tbl$group), , drop = FALSE]
  if (!nrow(span_tbl)) NULL else span_tbl
}

assign_all_lanes <- function(span_tbl, label_genes, label_cex, label_pad_in,
                              scale_factor, lane_gap, group_gap) {
  S_mir <- assign_lanes(subset(span_tbl, group == "mirna",  drop = FALSE))
  S_lnc <- assign_lanes(subset(span_tbl, group == "lncrna", drop = FALSE))
  S_cod <- assign_lanes(subset(span_tbl, group == "coding", drop = FALSE))

  if (label_genes) {
    for (grp_name in c("mir", "lnc", "cod")) {
      S <- get(paste0("S_", grp_name))
      if (nrow(S)) {
        assign(paste0("S_", grp_name),
               assign_lanes(S, label_cex, label_pad_in, scale_factor))
      }
    }
  }

  n_mir <- if (nrow(S_mir)) max(S_mir$lane) else 0L
  n_lnc <- if (nrow(S_lnc)) max(S_lnc$lane) else 0L
  n_cod <- if (nrow(S_cod)) max(S_cod$lane) else 0L

  # Inline compute_group_offsets
  offset_cod <- 1
  offset_lnc <- if (n_cod > 0L) n_cod * lane_gap + group_gap else 1
  offset_mir <- offset_lnc + if (n_lnc > 0L) n_lnc * (lane_gap * 0.98) else 1
  group_offsets <- c(mirna = offset_mir, lncrna = offset_lnc, coding = offset_cod)

  span_lanes <- rbind(S_mir, S_lnc, S_cod)

  list(
    span_lanes    = span_lanes,
    group_offsets = group_offsets,
    lane_by_gene  = setNames(span_lanes$lane, span_lanes$gene),
    group_by_gene = setNames(span_lanes$group, span_lanes$gene),
    strand_by_gene = setNames(span_lanes$strand, span_lanes$gene)
  )
}

prepare_exon_draw <- function(df, gx1_bp, gx2_bp, g, strand, bt, merge_exons) {
  df_ex <- unique(df[, c("chr", "start", "end", "strand", "gene", "biotype")])
  df_ex <- df_ex[!(df_ex$start == gx1_bp & df_ex$end == gx2_bp), , drop = FALSE]
  df_ex <- validate_intervals(df_ex)

  if (nrow(df_ex)) {
    if (merge_exons) merge_intervals(df_ex)
    else df_ex[order(df_ex$start, df_ex$end), , drop = FALSE]
  } else {
    data.frame(
      chr = df$chr[1], start = gx1_bp, end = gx2_bp,
      gene = g, strand = strand, biotype = bt,
      stringsAsFactors = FALSE
    )
  }
}

build_annotation_legend_info <- function(genes_skipped, min_extent) {
  list(
    labels = c("protein coding", "lncRNA", "miRNA"),
    colors = c(
      map_biotype_col("protein_coding"),
      map_biotype_col("lncRNA"),
      map_biotype_col("miRNA")
    ),
    types         = c("exon", "lncrna", "mirna"),
    min_extent    = min_extent,
    genes_skipped = genes_skipped
  )
}

n_bedpe <- 0L

###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if (isFALSE(flag_inter)) {
  if (analysis_type == "single_sample") {

    path_hic            <- Args[3]
    path_mgs            <- Args[4]
    norm                <- Args[5]
    unit                <- Args[6]
    bin_size            <- as.numeric(Args[7])
    prefix              <- Args[8]
    range_text          <- Args[9]
    out_dir             <- Args[10]
    interval_chr        <- Args[11]
    interval_start      <- as.numeric(Args[12])
    interval_end        <- as.numeric(Args[13])
    genome_build        <- Args[14]
    flag_profiles       <- as.logical(Args[15])
    contact_color       <- sprintf("#%s", Args[16])
    ann_style           <- Args[17]
    ann_custom          <- Args[18]
    ann_custom_colors   <- Args[19]
    max_cap             <- Args[20]
    flag_norm           <- Args[21]
    quant_cut           <- as.numeric(Args[22])
    output_name         <- Args[23]
    bedpe               <- Args[24]
    bedpe_color         <- sprintf("#%s", Args[25])
    inherent            <- as.logical(Args[26])
    sample_dir          <- Args[27]
    flag_w              <- Args[28]
    inh_col_floor       <- sprintf("#%s", Args[29])
    inh_col_off         <- sprintf("#%s", Args[30])
    inh_col_on          <- sprintf("#%s", Args[31])
    inh_col_ceil        <- sprintf("#%s", Args[32])
    flag_matrix         <- as.logical(Args[33])
    flag_plexus         <- as.logical(Args[34])
    flag_expand_viewport <- as.logical(Args[35])

    output_name <- basename(path.expand(output_name))
    out_path <- file.path(out_dir, output_name)

    interval_len <- interval_end - interval_start

    if (flag_expand_viewport) {
      shoulder_bins <- ceiling((interval_len / 2) / bin_size)
      shoulder_bp <- shoulder_bins * bin_size
      plot_start <- interval_start - shoulder_bp
      plot_end <- interval_end + shoulder_bp
    } else {
      plot_start <- interval_start
      plot_end <- interval_end
    }

    if (any(grepl("chr", interval_chr, ignore.case = T)) == FALSE) {
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save = "no", status = 1)
    }

    if (isFALSE(flag_matrix)) {
      if (bedpe != "FALSE") {
        # Accept a single space-joined string or a character vector
        bedpe_files <- if (length(bedpe) == 1L) strsplit(bedpe, "[[:space:]]+")[[1]] else bedpe
        bedpe_files <- bedpe_files[nzchar(bedpe_files)]

        # Pre-allocate list (names are basenames of files)
        pairs_by_file <- setNames(vector("list", length(bedpe_files)), basename(bedpe_files))

        n_bedpe <- length(bedpe_files) # original number of bedpe files provided
        k <- 0L # number of bedpe files that have >=1 pair fully in-range

        for (bedpe_path in bedpe_files) {
          pairs <- read.table(bedpe_path, as.is = TRUE)[, 1:6]
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

          if (sum(grepl("chr", pairs$chr1)) == 0) {
            pairs$chr1 <- paste0("chr", pairs$chr1)
          }
          if (sum(grepl("chr", pairs$chr2)) == 0) {
            pairs$chr2 <- paste0("chr", pairs$chr2)
          }

          # BEDPE is typically half-open [start, end).
          # Convert to an inclusive end coordinate [start, end-1]
          pairs$end1 <- ifelse(pairs$end1 > pairs$start1, pairs$end1 - 1, pairs$end1)
          pairs$end2 <- ifelse(pairs$end2 > pairs$start2, pairs$end2 - 1, pairs$end2)

          # STRICT: require both anchors fully contained in the plot interval
          pairs <- pairs[pairs[, "chr1"] == interval_chr, ]
          pairs <- pairs[
            pairs[, "end2"] <= interval_end &
              pairs[, "start1"] >= interval_start &
              pairs[, "end1"] <= interval_end &
              pairs[, "start2"] >= interval_start,
          ]

          if (nrow(pairs) == 0) {
            next # no rows from this file in range, skip
          }

          pairs <- pairs[pairs$chr1 %in% c(paste0(rep("chr", 22), 1:22), "chrX", "chrY"), ]

          # convert coordinates to bins
          pairs$start1_bin <- as.integer(floor(pairs$start1 / bin_size) * bin_size)
          pairs$end1_bin <- as.integer(floor(pairs$end1 / bin_size) * bin_size)
          pairs$start2_bin <- as.integer(floor(pairs$start2 / bin_size) * bin_size)
          pairs$end2_bin <- as.integer(floor(pairs$end2 / bin_size) * bin_size)

          for (i in 1:nrow(pairs)) {
            if (pairs[i, "start2_bin"] < pairs[i, "start1_bin"]) {
              a <- pairs[i, "start1_bin"]
              b <- pairs[i, "end1_bin"]
              c <- pairs[i, "start2_bin"]
              d <- pairs[i, "end2_bin"]
              pairs[i, "start1_bin"] <- c
              pairs[i, "end1_bin"] <- d
              pairs[i, "start2_bin"] <- a
              pairs[i, "end2_bin"] <- b
            }
            if (pairs[i, "chr1"] != pairs[i, "chr2"]) {
              cat("Please provide intra-chromosomal bedpes\n")
              quit(save = "no", status = 1)
            }
          }

          # Convert bp positions to bin indices
          # Add +1 to make first bin row/column (0 to 1 indexing to match straw)
          pairs$i <- (pairs$start1_bin - interval_start) / bin_size + 1
          pairs$j <- (pairs$start2_bin - interval_start) / bin_size + 1
          pairs$i_end <- (pairs$end1_bin - interval_start) / bin_size + 1
          pairs$j_end <- (pairs$end2_bin - interval_start) / bin_size + 1

          pairs_by_file[[basename(bedpe_path)]] <- pairs
          k <- k + 1L
        }

        # Keep only files that actually contributed pairs
        pairs_by_file <- pairs_by_file[vapply(pairs_by_file, function(x) !is.null(x) && nrow(x) > 0, logical(1))]

        # Update bedpe count to reflect only files in range
        n_bedpe <- length(pairs_by_file)

        if (n_bedpe == 0L) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          bedpe <- pairs_by_file
        }
      }
    }

    cat("\n  interval\n")
    cat(sprintf("  interval_chr:   %s\n", interval_chr))
    cat(sprintf("  interval_start: %d\n", interval_start))
    cat(sprintf("  interval_end:   %d\n", interval_end))

    ## 'chr' prefix handler:

    hic_A_chroms <- readHicChroms(path_hic)$name

    if (sum(grepl("chr", hic_A_chroms)) > 0) {
      flag_strip_chr <- FALSE
    }
    if (sum(grepl("chr", hic_A_chroms)) == 0) {
      flag_strip_chr <- TRUE
    }

    if (flag_strip_chr) {
      new_interval_chr <- sub("chr", "", interval_chr)
    } else {
      new_interval_chr <- interval_chr
    }

    ## CPM or AQuA factors

    mergeStats_A <- read.table(path_mgs, as.is = T)
    spikeVar <- ncol(mergeStats_A)
    hg_total1 <- as.numeric(mergeStats_A["valid_interaction_rmdup", 1])
    mm_total1 <- as.numeric(mergeStats_A["valid_interaction_rmdup", 2])
    total1 <- sum(hg_total1, mm_total1)
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1

    # set normalization factor for spike-in and non-spike-in data.
    # non-spike-in data defaults to cpm.
    # spike in data defaults to aqua.

    if (spikeVar == 1) {
      if (flag_norm == "blank") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        flag_norm <- "cpm"
      } else if (flag_norm == "none") {
        norm_factor1 <- 1
        aqua_factor1 <- 1
      } else if (flag_norm == "cpm") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      } else if (flag_norm == "aqua") {
        cat("\n\n--norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        flag_norm <- "cpm"
      }
    } else if (spikeVar == 2) {
      if (flag_norm == "blank") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
        flag_norm <- "aqua"
      } else if (flag_norm == "none") {
        norm_factor1 <- 1
        aqua_factor1 <- 1
      } else if (flag_norm == "cpm") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      } else if (flag_norm == "aqua") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
      }
    }

    if (isTRUE(inherent)) {
      flag_norm <- "inherent"
    }

    cat("\n  factors\n")
    cat(sprintf("  norm_factor: %f\n", norm_factor1))
    cat(sprintf("  aqua_factor: %f\n", aqua_factor1))

    ## Process matrix:

    cat("\n  parameters\n")
    cat(sprintf("       norm: %s\n", flag_norm))
    cat(sprintf("        hic: %s\n", path_hic))
    cat(sprintf("   interval: %s\n", paste(new_interval_chr, interval_start, interval_end, sep = ":")))
    cat(sprintf("       unit: %s\n", unit))
    cat(sprintf("   bin_size: %s\n", bin_size))

    sparMat1 <- straw(
      norm,
      path_hic,
      paste(new_interval_chr, plot_start, plot_end, sep = ":"),
      paste(new_interval_chr, plot_start, plot_end, sep = ":"),
      unit,
      bin_size
    )

    if (nrow(sparMat1) == 0) {
      cat("\n No contact values found for this interval\n")
      quit(save = "no", status = 1)
    }
    # Bin coords for expanded plot if expand_viewport = true
    bin_coords <- seq(plot_start, plot_end, by = bin_size)

    # Focus bins = interval bins if expand_viewport = false
    focus_bins <- which(bin_coords >= interval_start & bin_coords <= interval_end)
    if (!length(focus_bins)) {
      cat("\n No bins found overlapping the focus interval; check coordinates.\n")
      quit(save = "no", status = 1)
    }

    focus_start_idx <- focus_bins[1L]
    focus_end_idx <- focus_bins[length(focus_bins)]
    focus_len <- length(focus_bins)

    if (isFALSE(inherent)) {
      sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1

      denMat1 <- matrix(
        data = 0,
        nrow = length(bin_coords),
        ncol = length(bin_coords)
      )

      rownames(denMat1) <- bin_coords
      colnames(denMat1) <- bin_coords

      # Populate denMat1
      for (i in 1:nrow(sparMat1)) {
        x_index <- as.character(sparMat1[i, "x"])
        y_index <- as.character(sparMat1[i, "y"])

        if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
          denMat1[x_index, y_index] <- sparMat1[i, "counts"]
        }
      }

      if (flag_matrix) {
        # Subset to the focus interval
        denMat1_focus <- denMat1[
          focus_start_idx:focus_end_idx,
          focus_start_idx:focus_end_idx,
          drop = FALSE
        ]

        # Build bin labels
        focus_coords <- bin_coords[focus_start_idx:focus_end_idx]

        bin_labels <- sprintf(
          "%s:%d:%d",
          interval_chr,
          focus_coords,
          focus_coords + bin_size - 1L
        )

        if (length(bin_labels) == nrow(denMat1_focus)) {
          rownames(denMat1_focus) <- bin_labels
          colnames(denMat1_focus) <- bin_labels
        } else {
          warning("bin_labels length does not match denMat1 dimensions; keeping original row/col names")
        }
        # Use the output name as passed from the shell
        file_name <- basename(output_name)
        out_file <- file.path(out_dir, file_name)
        write.table(denMat1_focus, file = out_file, quote = FALSE, sep = "\t")
        quit(save = "no", status = 0)
      }

      # Get cap value
      cap <- calculate_cap(denMat1, max_cap, quant_cut)
      cap <- round(cap, 3)

      # Print the cap
      cat(sprintf("        cap: %s\n", cap))

      # Cap the values in the original matrix
      denMat1[denMat1 > cap] <- cap
      denMat1max <- denMat1

      color_ramp <- colorRampPalette(c("white", contact_color))(101)

      breakList <- seq(0, cap, by = cap / 100)

      C <- matrix(ncol = 2, nrow = length(breakList), data = 0)
      rownames(C) <- color_ramp
      colnames(C) <- c("rank", "breaks")
      C[, 1] <- 1:nrow(C)
      C[, 2] <- breakList


      A_big <- denMat1max
      A_focus <- A_big[focus_start_idx:focus_end_idx,
        focus_start_idx:focus_end_idx,
        drop = FALSE
      ]
      A <- A_focus

      initial_layout <- compute_initial_layout(
        A          = A_focus,
        xspan      = 100,
        yspan      = 100,
        pad_x      = 0,
        pad_y      = 0
      )
      max_w <- initial_layout$w

      # Choose w
      if (identical(flag_w, "blank")) {
        w <- max_w
      } else {
        w_req <- suppressWarnings(as.numeric(flag_w))
        if (!is.finite(w_req) || w_req <= 0) {
          w <- max_w
        } else if (w_req > max_w) {
          cat(sprintf("    Requested -w %.3f > max allowable w %.3f; using max w...\n", w_req, max_w))
          w <- max_w
        } else {
          w <- w_req
        }
      }
      cat(sprintf("   plotting bin width: %s\n\n", round(w, 3)))

      # check annotation files
      valid_custom_files <- preflight_custom_beds(ann_custom, interval_chr, interval_start, interval_end)
      n_bed <- length(valid_custom_files)
      n_custom <- n_bedpe + n_bed
      ann_custom <- if (n_bed > 0) paste(valid_custom_files, collapse = " ") else "NONE"

      final_layout <- compute_final_layout(A_focus, w, n_custom)
      w <- final_layout$w
      start_x <- final_layout$start_x
      start_y <- final_layout$start_y
      legend_mid_y <- final_layout$legend_mid_y

      fmt <- if (flag_plexus) "svg" else "pdf"

      render_triangle_plot(
        output_name     = out_path,
        A_mat           = A_big,
        C_tbl           = C,
        breaks_vec      = breakList,
        bedpe           = bedpe,
        scale_fun       = draw_scale,
        device          = fmt,
        ann_style       = ann_style,
        focus_start_idx = focus_start_idx,
        focus_len       = focus_len,
        band_height     = focus_len - 1L,
        legend_mid_y    = legend_mid_y
      )
      if (identical(fmt, "pdf")) {
        # Read lowest annotation y from draw_default_annotations
        min_y <- get0("last_annotation_min_y", ifnotfound = NA_real_)
        min_y <- min_y - 2

        if (is.finite(min_y)) {
          cut_bottom <- estimate_crop(min_y)
        } else {
          # when annotations_default is none
          cut_bottom <- 200
        }
        crop_pdf(out_path, cut_bottom = cut_bottom)
      }
    }

    if (isTRUE(inherent)) {
      cat(sprintf("   inh_col_floor: %s\n", inh_col_floor))
      cat(sprintf("   inh_col_off:   %s\n", inh_col_off))
      cat(sprintf("   inh_col_on:    %s\n", inh_col_on))
      cat(sprintf("   inh_col_ceil:  %s\n", inh_col_ceil))

      denMat1 <- matrix(
        data = 0,
        nrow = length(bin_coords),
        ncol = length(bin_coords)
      )

      rownames(denMat1) <- bin_coords
      colnames(denMat1) <- bin_coords

      # Populate denMat1
      for (i in 1:nrow(sparMat1)) {
        x_index <- as.character(sparMat1[i, "x"])
        y_index <- as.character(sparMat1[i, "y"])

        if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
          denMat1[x_index, y_index] <- sparMat1[i, "counts"]
        }
      }

      power_law_path <- list.files(
        sample_dir,
        pattern = paste0("inherentStats", ".txt"),
        full.names = T
      )

      if (length(power_law_path) == 0) {
        cat("Power laws unavailable for this sample! \n")
        q(save = "no")
      }

      power_laws <- read.table(
        power_law_path,
        as.is = T,
        skip = 1
      )

      colnames(power_laws) <- c("off", "on", "res")
      power_laws <- power_laws[!power_laws$off %in% "off", ]

      power_laws$off <- as.numeric(power_laws$off)
      power_laws$on <- as.numeric(power_laws$on)
      power_laws$res <- as.numeric(power_laws$res)
      power_laws <- power_laws[power_laws$res == bin_size, ]

      if (nrow(power_laws) == 0) {
        cat("Power laws unavailable at this resolution! \n")
        q(save = "no")
      }

      if (nrow(power_laws) < ncol(denMat1)) {
        diff <- ncol(denMat1) - nrow(power_laws)

        power_laws <- rbind(
          power_laws,
          do.call(
            "rbind",
            replicate(diff, power_laws[nrow(power_laws), ], simplify = FALSE)
          )
        )
      }

      background <- power_laws[, "off"]
      foreground <- power_laws[, "on"]

      # Remove background
      tad_matrix_bg_removed <- denMat1

      for (j in 0:(ncol(denMat1) - 1)) {
        for (i in 1:(nrow(denMat1) - j)) {
          tad_matrix_bg_removed[i, i + j] <-
            denMat1[i, i + j] - background[j + 1]
        }
      }

      # Standardize
      tad_matrix_sd <- tad_matrix_bg_removed

      for (j in 0:(ncol(tad_matrix_bg_removed) - 1)) {
        for (i in 1:(nrow(tad_matrix_bg_removed) - j)) {
          tad_matrix_sd[i, i + j] <-
            tad_matrix_bg_removed[i, i + j] / (foreground[j + 1] - background[j + 1])
        }
      }

      A <- tad_matrix_sd
      for (i in 2:nrow(A)) {
        for (j in 1:(i - 1)) {
          A[i, j] <- 0
        }
      }

      if (flag_matrix) {
        # Subset to the focus interval
        A_focus <- A[
          focus_start_idx:focus_end_idx,
          focus_start_idx:focus_end_idx,
          drop = FALSE
        ]

        # Build bin labels
        focus_coords <- bin_coords[focus_start_idx:focus_end_idx]

        bin_labels <- sprintf(
          "%s:%d:%d",
          interval_chr,
          focus_coords,
          focus_coords + bin_size - 1L
        )

        if (length(bin_labels) == nrow(A_focus)) {
          rownames(A_focus) <- bin_labels
          colnames(A_focus) <- bin_labels
        } else {
          warning("bin_labels length does not match A dimensions; keeping original row/col names")
        }

        # Use the output name as passed from the shell
        file_name <- basename(output_name)
        out_file <- file.path(out_dir, file_name)
        write.table(A_focus, file = out_file, quote = FALSE, sep = "\t")
        quit(save = "no", status = 0)
      }

      if (max_cap != "none") {
        cat("\nParameter --max_cap is not applicable when inherent = TRUE\nContinuing without --max_cap...\n\n")
      }
      if (quant_cut != 1) {
        cat("\nParameter --quant_cut is not applicable when inherent = TRUE\nContinuing without --quant_cut...\n\n")
      }

      colors_n <- 10
      colors_0_1 <- colorRampPalette(c(inh_col_off, inh_col_on))(colors_n)
      colors_1_2 <- colorRampPalette(c(inh_col_on, inh_col_ceil))(colors_n)
      colors <- c(inh_col_floor, colors_0_1, colors_1_2)

      breaks <- c((1:colors_n) / colors_n, (1:colors_n) / colors_n + 1)
      breaks <- c(0, breaks) # 0 … 2
      C <- matrix(ncol = 2, nrow = length(breaks), data = 0)
      rownames(C) <- colors
      colnames(C) <- c("rank", "breaks")
      C[, 1] <- 1:nrow(C)
      C[, 2] <- breaks

      breakList <- breaks

      A_big <- A
      A_focus <- A_big[focus_start_idx:focus_end_idx,
        focus_start_idx:focus_end_idx,
        drop = FALSE
      ]
      A <- A_focus

      # Compute layout
      initial_layout <- compute_initial_layout(
        A          = A_focus,
        xspan      = 100,
        yspan      = 100,
        pad_x      = 0,
        pad_y      = 0
      )
      max_w <- initial_layout$w

      # Choose w
      if (identical(flag_w, "blank")) {
        w <- max_w
      } else {
        w_req <- suppressWarnings(as.numeric(flag_w))
        if (!is.finite(w_req) || w_req <= 0) {
          w <- max_w
        } else if (w_req > max_w) {
          cat(sprintf("    Requested -w %.3f > max allowable w %.3f; using max w...\n", w_req, max_w))
          w <- max_w
        } else {
          w <- w_req
        }
      }

      cat(sprintf("   plotting bin width: %s\n", round(w, 3)))

      # check annotation files
      valid_custom_files <- preflight_custom_beds(ann_custom, interval_chr, interval_start, interval_end)
      n_bed <- length(valid_custom_files)
      n_custom <- n_bedpe + n_bed
      ann_custom <- if (n_bed > 0) paste(valid_custom_files, collapse = " ") else "NONE"

      final_layout <- compute_final_layout(A_focus, w, n_custom)

      w <- final_layout$w
      start_x <- final_layout$start_x
      start_y <- final_layout$start_y
      legend_mid_y <- final_layout$legend_mid_y

      fmt <- if (flag_plexus) "svg" else "pdf"

      render_triangle_plot(
        output_name     = out_path,
        A_mat           = A_big,
        C_tbl           = C,
        breaks_vec      = breakList,
        bedpe           = bedpe,
        scale_fun       = draw_scale_inh,
        device          = fmt,
        ann_style       = ann_style,
        focus_start_idx = focus_start_idx,
        focus_len       = focus_len,
        band_height     = focus_len - 1L,
        legend_mid_y    = legend_mid_y
      )

      if (identical(fmt, "pdf")) {
        # Read lowest annotation y from draw_default_annotations
        min_y <- get0("last_annotation_min_y", ifnotfound = NA_real_)
        min_y <- min_y - 3

        if (is.finite(min_y)) {
          cut_bottom <- estimate_crop(min_y)
        } else {
          # when annotations_default is none
          cut_bottom <- 200
        }
        crop_pdf(out_path, cut_bottom = cut_bottom)
      }
    }
  }
}

###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if (isFALSE(flag_inter)) {
  if (analysis_type == "two_sample") {

    path_hic_A           <- Args[3]
    path_hic_B           <- Args[4]
    path_mgs_A           <- Args[5]
    path_mgs_B           <- Args[6]
    norm                 <- Args[7]
    unit                 <- Args[8]
    bin_size             <- as.numeric(Args[9])
    prefix               <- Args[10]
    range_text           <- Args[11]
    out_dir              <- Args[12]
    interval_chr         <- Args[13]
    interval_start       <- as.numeric(Args[14])
    interval_end         <- as.numeric(Args[15])
    genome_build         <- Args[16]
    flag_profiles        <- as.logical(Args[17])
    contact_color        <- Args[18]
    ann_style            <- Args[19]
    ann_custom           <- Args[20]
    ann_custom_colors    <- Args[21]
    max_cap              <- Args[22]
    flag_norm            <- Args[23]
    quant_cut            <- as.numeric(Args[24])
    output_name          <- Args[25]
    bedpe                <- Args[26]
    bedpe_color          <- sprintf("#%s", Args[27])
    flag_w               <- Args[28]
    sample_dirA          <- Args[29]
    sample_dirB          <- Args[30]
    flag_matrix          <- as.logical(Args[31])
    inh_col_floor        <- sprintf("#%s", Args[32])
    inh_col_off          <- sprintf("#%s", Args[33])
    inh_col_on           <- sprintf("#%s", Args[34])
    inh_col_ceil         <- sprintf("#%s", Args[35])
    inherent             <- as.logical(Args[36])
    flag_plexus          <- as.logical(Args[37])
    flag_expand_viewport <- as.logical(Args[38])

    output_name <- basename(path.expand(output_name))
    out_path <- file.path(out_dir, output_name)

    interval_len <- interval_end - interval_start

    if (flag_expand_viewport) {
      shoulder_bins <- ceiling((interval_len / 2) / bin_size)
      shoulder_bp <- shoulder_bins * bin_size
      plot_start <- interval_start - shoulder_bp
      plot_end <- interval_end + shoulder_bp
    } else {
      plot_start <- interval_start
      plot_end <- interval_end
    }

    if (any(grepl("chr", interval_chr, ignore.case = T)) == FALSE) {
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save = "no", status = 1)
    }

    if (isFALSE(flag_matrix)) {
      if (bedpe != "FALSE") {
        # Accept a single space-joined string or a character vector
        bedpe_files <- if (length(bedpe) == 1L) strsplit(bedpe, "[[:space:]]+")[[1]] else bedpe
        bedpe_files <- bedpe_files[nzchar(bedpe_files)]

        # Pre-allocate list (names are basenames of files)
        pairs_by_file <- setNames(vector("list", length(bedpe_files)), basename(bedpe_files))

        n_bedpe <- length(bedpe_files) # original number of bedpe files provided
        k <- 0L # number of bedpe files that have >=1 pair fully in-range

        for (bedpe_path in bedpe_files) {
          pairs <- read.table(bedpe_path, as.is = TRUE)[, 1:6]
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

          if (sum(grepl("chr", pairs$chr1)) == 0) {
            pairs$chr1 <- paste0("chr", pairs$chr1)
          }
          if (sum(grepl("chr", pairs$chr2)) == 0) {
            pairs$chr2 <- paste0("chr", pairs$chr2)
          }

          # BEDPE is typically half-open [start, end).
          # Convert to an inclusive end coordinate [start, end-1]
          pairs$end1 <- ifelse(pairs$end1 > pairs$start1, pairs$end1 - 1, pairs$end1)
          pairs$end2 <- ifelse(pairs$end2 > pairs$start2, pairs$end2 - 1, pairs$end2)

          # STRICT: require both anchors fully contained in the plot interval
          pairs <- pairs[pairs[, "chr1"] == interval_chr, ]
          pairs <- pairs[
            pairs[, "end2"] <= interval_end &
              pairs[, "start1"] >= interval_start &
              pairs[, "end1"] <= interval_end &
              pairs[, "start2"] >= interval_start,
          ]

          if (nrow(pairs) == 0) {
            next # no rows from this file in range, skip
          }

          pairs <- pairs[pairs$chr1 %in% c(paste0(rep("chr", 22), 1:22), "chrX", "chrY"), ]

          # convert coordinates to bins
          pairs$start1_bin <- as.integer(floor(pairs$start1 / bin_size) * bin_size)
          pairs$end1_bin <- as.integer(floor(pairs$end1 / bin_size) * bin_size)
          pairs$start2_bin <- as.integer(floor(pairs$start2 / bin_size) * bin_size)
          pairs$end2_bin <- as.integer(floor(pairs$end2 / bin_size) * bin_size)

          for (i in 1:nrow(pairs)) {
            if (pairs[i, "start2_bin"] < pairs[i, "start1_bin"]) {
              a <- pairs[i, "start1_bin"]
              b <- pairs[i, "end1_bin"]
              c <- pairs[i, "start2_bin"]
              d <- pairs[i, "end2_bin"]
              pairs[i, "start1_bin"] <- c
              pairs[i, "end1_bin"] <- d
              pairs[i, "start2_bin"] <- a
              pairs[i, "end2_bin"] <- b
            }
            if (pairs[i, "chr1"] != pairs[i, "chr2"]) {
              cat("Please provide intra-chromosomal bedpes\n")
              quit(save = "no", status = 1)
            }
          }

          # Convert bp positions to bin indices
          # Add +1 to make first bin row/column (0 to 1 indexing to match straw)
          pairs$i <- (pairs$start1_bin - interval_start) / bin_size + 1
          pairs$j <- (pairs$start2_bin - interval_start) / bin_size + 1
          pairs$i_end <- (pairs$end1_bin - interval_start) / bin_size + 1
          pairs$j_end <- (pairs$end2_bin - interval_start) / bin_size + 1

          pairs_by_file[[basename(bedpe_path)]] <- pairs
          k <- k + 1L
        }

        # Keep only files that actually contributed pairs
        pairs_by_file <- pairs_by_file[vapply(pairs_by_file, function(x) !is.null(x) && nrow(x) > 0, logical(1))]

        # Update bedpe count to reflect only files in range
        n_bedpe <- length(pairs_by_file)

        if (n_bedpe == 0L) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          bedpe <- pairs_by_file
        }
      }
    }

    cat("\n  interval\n")
    cat(sprintf("  interval_chr:   %s\n", interval_chr))
    cat(sprintf("  interval_start: %d\n", interval_start))
    cat(sprintf("  interval_end:   %d\n", interval_end))

    ## 'chr' prefix handler:

    hic_A_chroms <- readHicChroms(path_hic_A)$name
    hic_B_chroms <- readHicChroms(path_hic_B)$name

    if (sum(grepl("chr", hic_A_chroms)) > 0 && sum(grepl("chr", hic_B_chroms)) > 0) {
      flag_strip_chr <- "no"
    }
    if (sum(grepl("chr", hic_A_chroms)) > 0 && sum(grepl("chr", hic_B_chroms)) == 0) {
      flag_strip_chr <- "A"
    }
    if (sum(grepl("chr", hic_B_chroms)) > 0 && sum(grepl("chr", hic_A_chroms)) == 0) {
      flag_strip_chr <- "B"
    }
    if (sum(grepl("chr", hic_A_chroms)) == 0 && sum(grepl("chr", hic_B_chroms)) == 0) {
      flag_strip_chr <- "yes"
    }

    if (flag_strip_chr == "yes") {
      interval_chr_A <- sub("chr", "", interval_chr)
      interval_chr_B <- sub("chr", "", interval_chr)
    } else if (flag_strip_chr == "B") {
      interval_chr_A <- interval_chr
      interval_chr_B <- sub("chr", "", interval_chr)
    } else if (flag_strip_chr == "A") {
      interval_chr_A <- sub("chr", "", interval_chr)
      interval_chr_B <- interval_chr
    } else if (flag_strip_chr == "no") {
      interval_chr_A <- interval_chr
      interval_chr_B <- interval_chr
    }

    ## CPM or AQuA factors
    mergeStats_A <- read.table(path_mgs_A, as.is = T)
    spikeVar_A <- ncol(mergeStats_A)
    hg_total1 <- as.numeric(mergeStats_A["valid_interaction_rmdup", 1])
    mm_total1 <- as.numeric(mergeStats_A["valid_interaction_rmdup", 2])
    total1 <- sum(hg_total1, mm_total1)
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1

    mergeStats_B <- read.table(path_mgs_B, as.is = T)
    spikeVar_B <- ncol(mergeStats_B)
    hg_total2 <- as.numeric(mergeStats_B["valid_interaction_rmdup", 1])
    mm_total2 <- as.numeric(mergeStats_B["valid_interaction_rmdup", 2])
    total2 <- sum(hg_total2, mm_total2)
    norm_factor2 <- 1000000 / total2
    aqua_factor2 <- hg_total2 / mm_total2

    # set normalization factor for spike-in and non-spike-in data.
    # non-spike-in data defaults to cpm
    # spike-in data defaults to aqua

    if (spikeVar_A == 1 || spikeVar_B == 1) {
      if (flag_norm == "blank") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- 1
        flag_norm <- "cpm"
      } else if (flag_norm == "none") {
        norm_factor1 <- 1
        aqua_factor1 <- 1
        norm_factor2 <- 1
        aqua_factor2 <- 1
      } else if (flag_norm == "cpm") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- 1
      } else if (flag_norm == "aqua") {
        cat("\n\n--norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- 1
        flag_norm <- "cpm"
      }
    } else if (spikeVar_A == 2 && spikeVar_B == 2) {
      if (flag_norm == "blank") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- aqua_factor2
        flag_norm <- "aqua"
      } else if (flag_norm == "none") {
        norm_factor1 <- 1
        aqua_factor1 <- 1
        norm_factor2 <- 1
        aqua_factor2 <- 1
      } else if (flag_norm == "cpm") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- 1
      } else if (flag_norm == "aqua") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- aqua_factor2
      }
    }

    cat("\n  factors\n")
    cat(sprintf("  norm_factor1: %f\n", norm_factor1))
    cat(sprintf("  norm_factor2: %f\n", norm_factor2))
    cat(sprintf("  aqua_factor1: %f\n", aqua_factor1))
    cat(sprintf("  aqua_factor2: %f\n", aqua_factor2))

    ## Process matrix:
    cat("\n  parameters\n")
    cat(sprintf("       norm: %s\n", flag_norm))
    cat(sprintf("      hic_A: %s\n", path_hic_A))
    cat(sprintf("      hic_B: %s\n", path_hic_B))
    cat(sprintf("   interval: %s\n", paste(interval_chr, interval_start, interval_end, sep = ":")))
    cat(sprintf("       unit: %s\n", unit))
    cat(sprintf("   bin_size: %s\n", bin_size))

    sparMat1 <- straw(
      norm,
      path_hic_A,
      paste(interval_chr_A, plot_start, plot_end, sep = ":"),
      paste(interval_chr_A, plot_start, plot_end, sep = ":"),
      unit,
      bin_size
    )

    if (nrow(sparMat1) == 0) {
      cat("\n No contact values found for this interval in sample A\n")
      quit(save = "no", status = 1)
    }

    sparMat2 <- straw(
      norm,
      path_hic_B,
      paste(interval_chr_B, plot_start, plot_end, sep = ":"),
      paste(interval_chr_B, plot_start, plot_end, sep = ":"),
      unit,
      bin_size
    )

    if (nrow(sparMat2) == 0) {
      cat("\n No contact values found for this interval in sample B\n")
      quit(save = "no", status = 1)
    }

    bin_coords <- seq(plot_start, plot_end, by = bin_size)

    focus_bins <- which(bin_coords >= interval_start & bin_coords <= interval_end)
    if (!length(focus_bins)) {
      cat("\n No bins found overlapping the focus interval; check coordinates.\n")
      quit(save = "no", status = 1)
    }

    focus_start_idx <- focus_bins[1L]
    focus_end_idx <- focus_bins[length(focus_bins)]
    focus_len <- length(focus_bins)

    valid_custom_files <- preflight_custom_beds(ann_custom, interval_chr, interval_start, interval_end)

    n_bed <- length(valid_custom_files)
    n_custom <- n_bedpe + n_bed

    ann_custom <- if (n_bed > 0) paste(valid_custom_files, collapse = " ") else "NONE"

    if (isFALSE(inherent)) {
      sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1
      sparMat2$counts <- sparMat2$counts * norm_factor2 * aqua_factor2

      denMat1 <- matrix(
        data = 0,
        nrow = length(bin_coords),
        ncol = length(bin_coords)
      )

      rownames(denMat1) <- bin_coords
      colnames(denMat1) <- bin_coords

      denMat2 <- matrix(
        data = 0,
        nrow = length(bin_coords),
        ncol = length(bin_coords)
      )

      rownames(denMat2) <- bin_coords
      colnames(denMat2) <- bin_coords

      # Populate denMat1
      for (i in 1:nrow(sparMat1)) {
        x_index <- as.character(sparMat1[i, "x"])
        y_index <- as.character(sparMat1[i, "y"])

        if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
          denMat1[x_index, y_index] <- sparMat1[i, "counts"]
        }
      }

      # Populate denMat2
      for (i in 1:nrow(sparMat2)) {
        x_index <- as.character(sparMat2[i, "x"])
        y_index <- as.character(sparMat2[i, "y"])

        if (x_index %in% rownames(denMat2) && y_index %in% colnames(denMat2)) {
          denMat2[x_index, y_index] <- sparMat2[i, "counts"]
        }
      }

      if (nrow(denMat1) != nrow(denMat2)) {
        intersect <- intersect(rownames(denMat1), rownames(denMat2))
        denMat1 <- denMat1[intersect, intersect]
        denMat2 <- denMat2[intersect, intersect]
      }

      denDelta <- denMat2 - denMat1

      if (flag_matrix) {
        # Subset to the focus interval
        denDelta_focus <- denDelta[
          focus_start_idx:focus_end_idx,
          focus_start_idx:focus_end_idx,
          drop = FALSE
        ]

        # Build bin labels
        focus_coords <- bin_coords[focus_start_idx:focus_end_idx]

        bin_labels <- sprintf(
          "%s:%d:%d",
          interval_chr,
          focus_coords,
          focus_coords + bin_size - 1L
        )

        if (length(bin_labels) == nrow(denDelta_focus)) {
          rownames(denDelta_focus) <- bin_labels
          colnames(denDelta_focus) <- bin_labels
        } else {
          warning("bin_labels length does not match denMat1 dimensions; keeping original row/col names")
        }
        # Use the output name as passed from the shell
        file_name <- basename(output_name)
        out_file <- file.path(out_dir, file_name)
        write.table(denDelta_focus, file = out_file, quote = FALSE, sep = "\t")
        quit(save = "no", status = 0)
      }

      # Max_cap and quant_cut handling
      cap <- calculate_cap(denDelta, max_cap, quant_cut)
      cap <- round(cap, 3)

      denDelta[denDelta > cap] <- cap
      denDelta[denDelta < -cap] <- -cap

      cat(sprintf("        cap: %s\n", cap))

      if (isFALSE(grepl("-", contact_color))) {
        cat("\nPlease separate the two delta colors with '-' only! \n")
        quit(save = "no", status = 1)
      }

      color_neg <- sprintf("#%s", unlist(strsplit(contact_color, "-"))[1])
      color_pos <- sprintf("#%s", unlist(strsplit(contact_color, "-"))[2])

      breakList <- seq(
        -cap,
        cap,
        by = cap / 100
      )

      color_ramp <- colorRampPalette(c(color_neg, "white", color_pos))(length(breakList))

      C <- matrix(ncol = 2, nrow = length(breakList), data = 0)
      rownames(C) <- color_ramp
      colnames(C) <- c("rank", "breaks")
      C[, 1] <- 1:nrow(C)
      C[, 2] <- breakList

      A_big <- denDelta
      A_focus <- A_big[focus_start_idx:focus_end_idx,
        focus_start_idx:focus_end_idx,
        drop = FALSE
      ]
      A <- A_focus

      n_bed <- if (!identical(ann_custom, "NONE")) length(strsplit(ann_custom, " ")[[1]]) else 0
      n_custom <- n_bedpe + n_bed

      initial_layout <- compute_initial_layout(
        A          = A_focus,
        xspan      = 100,
        yspan      = 100,
        pad_x      = 0,
        pad_y      = 0
      )
      max_w <- initial_layout$w

      # Choose w
      if (identical(flag_w, "blank")) {
        w <- max_w
      } else {
        w_req <- suppressWarnings(as.numeric(flag_w))
        if (!is.finite(w_req) || w_req <= 0) {
          w <- max_w
        } else if (w_req > max_w) {
          cat(sprintf("    Requested -w %.3f > max allowable w %.3f; using max w...\n", w_req, max_w))
          w <- max_w
        } else {
          w <- w_req
        }
      }
      cat(sprintf("   plotting bin width: %s\n", round(w, 3)))

      final_layout <- compute_final_layout(A_focus, w, n_custom)
      w <- final_layout$w
      start_x <- final_layout$start_x
      start_y <- final_layout$start_y
      legend_mid_y <- final_layout$legend_mid_y

      fmt <- if (flag_plexus) "svg" else "pdf"

      render_triangle_plot(
        output_name     = out_path,
        A_mat           = A_big,
        C_tbl           = C,
        breaks_vec      = breakList,
        bedpe           = bedpe,
        scale_fun       = draw_scale,
        device          = fmt,
        ann_style       = ann_style,
        focus_start_idx = focus_start_idx,
        focus_len       = focus_len,
        band_height     = focus_len - 1L,
        legend_mid_y    = legend_mid_y
      )
      if (identical(fmt, "pdf")) {
        # Read lowest annotation y from draw_default_annotations
        min_y <- get0("last_annotation_min_y", ifnotfound = NA_real_)
        min_y <- min_y - 2

        if (is.finite(min_y)) {
          cut_bottom <- estimate_crop(min_y)
        } else {
          # when annotations_default is none
          cut_bottom <- 200
        }
        crop_pdf(out_path, cut_bottom = cut_bottom)
      }
    }

    if (isTRUE(inherent)) {
      cat("\n--inherent TRUE is currently only available for one sample tests only.\n")
    }
  }
}

###########################################################################
###########################################################################
###                                                                     ###
###                   INTER-CHROMOSOMAL ONE-SAMPLE                      ###
###                                                                     ###
###########################################################################
###########################################################################

if (isTRUE(flag_inter)) {
  if (analysis_type == "single_sample") {
    path_hic             <- Args[3]
    path_mgs             <- Args[4]
    norm                 <- Args[5]
    unit                 <- Args[6]
    bin_size             <- as.numeric(Args[7])
    prefix               <- Args[8]
    range_text           <- Args[9]
    out_dir              <- Args[10]
    interval_chr         <- Args[11]
    interval_start       <- as.numeric(Args[12])
    interval_end         <- as.numeric(Args[13])
    genome_build         <- Args[14]
    flag_profiles        <- as.logical(Args[15])
    contact_color        <- sprintf("#%s", Args[16])
    ann_style            <- Args[17]
    ann_custom           <- Args[18]
    ann_custom_colors    <- Args[19]
    max_cap              <- Args[20]
    flag_norm            <- Args[21]
    quant_cut            <- as.numeric(Args[22])
    output_name          <- Args[23]
    bedpe                <- Args[24]
    bedpe_color          <- sprintf("#%s", Args[25])
    inherent             <- as.logical(Args[26])
    sample_dir           <- Args[27]
    flag_w               <- Args[28]
    inh_col_floor        <- sprintf("#%s", Args[29])
    inh_col_off          <- sprintf("#%s", Args[30])
    inh_col_on           <- sprintf("#%s", Args[31])
    inh_col_ceil         <- sprintf("#%s", Args[32])
    flag_matrix          <- as.logical(Args[33])
    flag_plexus          <- as.logical(Args[34])
    flag_expand_viewport <- as.logical(Args[35])
    orig_chr1            <- Args[36]
    orig_chr2            <- Args[37]

    output_name <- basename(path.expand(output_name))
    out_path <- file.path(out_dir, output_name)

    flag_expand_viewport <- FALSE

    # Split the range_text at underscores
    range_components <- unlist(strsplit(range_text, "_"))

    # Extract the individual components
    interval_chr1 <- range_components[1]
    interval_start1 <- as.numeric(range_components[2])
    interval_end1 <- as.numeric(range_components[3])
    interval_chr2 <- range_components[4]
    interval_start2 <- as.numeric(range_components[5])
    interval_end2 <- as.numeric(range_components[6])

    interval_len1 <- interval_end1 - interval_start1
    interval_len2 <- interval_end2 - interval_start2

    if (ann_style != "none"){
      ann_reserve_bottom <- 35
      ann_reserve_left   <- 25
    } else {
      ann_reserve_bottom <- 16
      ann_reserve_left   <- 8
    }

    available_width  <- 100 - ann_reserve_left
    available_height <- 100 - ann_reserve_bottom

    dimensions <- calculate_plot_dimensions(
      interval_start1, interval_end1,
      interval_start2, interval_end2,
      bin_size,
      max_width  = available_width,
      max_height = available_height
    )

    plot_width <- dimensions$width
    plot_height <- dimensions$height

    flip_axes <- (orig_chr1 != interval_chr1)

    if (flip_axes) {
      temp_val <- plot_width
      plot_width <- plot_height
      plot_height <- temp_val
    }

    if (any(grepl("chr", interval_chr1, ignore.case = T)) == FALSE) {
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save = "no", status = 1)
    }
    if (any(grepl("chr", interval_chr2, ignore.case = T)) == FALSE) {
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save = "no", status = 1)
    }


    ## 'chr' prefix handler:
    hic_A_chroms <- readHicChroms(path_hic)$name

    if (sum(grepl("chr", hic_A_chroms)) > 0) {
      flag_strip_chr <- FALSE
    }
    if (sum(grepl("chr", hic_A_chroms)) == 0) {
      flag_strip_chr <- TRUE
    }

    if (flag_strip_chr) {
      new_interval_chr1 <- sub("chr", "", interval_chr1)
      new_interval_chr2 <- sub("chr", "", interval_chr2)
    } else {
      new_interval_chr1 <- interval_chr1
      new_interval_chr2 <- interval_chr2
    }

    if (isFALSE(flag_matrix)) {
      if (bedpe != "FALSE") {
        # Accept a single space-joined string or a character vector
        bedpe_files <- if (length(bedpe) == 1L) strsplit(bedpe, "[[:space:]]+")[[1]] else bedpe
        bedpe_files <- bedpe_files[nzchar(bedpe_files)]

        # Pre-allocate list (names are basenames of files)
        pairs_by_file <- setNames(vector("list", length(bedpe_files)), basename(bedpe_files))

        n_bedpe <- length(bedpe_files)
        k <- 0L

        for (bedpe_path in bedpe_files) {
          pairs <- read.table(bedpe_path, as.is = TRUE, header = FALSE)[, 1:6]
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

          # Make sure chromosome names have 'chr' prefix
          pairs$chr1 <- ifelse(grepl("chr", pairs$chr1), pairs$chr1, paste0("chr", pairs$chr1))
          pairs$chr2 <- ifelse(grepl("chr", pairs$chr2), pairs$chr2, paste0("chr", pairs$chr2))

          # Filter: both anchors must be fully contained in the two intervals
          # (allow either orientation in the file)
          pairs <- pairs[
            (
              (pairs$chr1 == new_interval_chr1 &
              pairs$start1 >= interval_start1 & pairs$end1 <= interval_end1 &
              pairs$chr2 == new_interval_chr2 &
              pairs$start2 >= interval_start2 & pairs$end2 <= interval_end2) |
              (pairs$chr1 == new_interval_chr2 &
              pairs$start1 >= interval_start2 & pairs$end1 <= interval_end2 &
              pairs$chr2 == new_interval_chr1 &
              pairs$start2 >= interval_start1 & pairs$end2 <= interval_end1)
            ),
          ]

          if (nrow(pairs) == 0) {
            next
          }

          # Reorder columns so chr1 always corresponds to interval 1
          for (idx in 1:nrow(pairs)) {
            if (pairs[idx, "chr1"] != new_interval_chr1) {
              pairs[idx, c("chr1", "start1", "end1", "chr2", "start2", "end2")] <-
                pairs[idx, c("chr2", "start2", "end2", "chr1", "start1", "end1")]
            }
          }
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

          # Adjust end coordinates for half-open binning
          pairs$end1 <- ifelse(pairs$end1 > pairs$start1, pairs$end1 - 1, pairs$end1)
          pairs$end2 <- ifelse(pairs$end2 > pairs$start2, pairs$end2 - 1, pairs$end2)

          # Convert coordinates to bins
          pairs$start1_bin <- as.integer(floor(pairs$start1 / bin_size) * bin_size)
          pairs$end1_bin   <- as.integer(floor(pairs$end1 / bin_size) * bin_size)
          pairs$start2_bin <- as.integer(floor(pairs$start2 / bin_size) * bin_size)
          pairs$end2_bin   <- as.integer(floor(pairs$end2 / bin_size) * bin_size)

          # Calculate plotting coordinates
          pairs$i     <- (pairs$start1_bin - interval_start1) / bin_size + 1
          pairs$j     <- (pairs$start2_bin - interval_start2) / bin_size + 1
          pairs$i_end <- (pairs$end1_bin - interval_start1) / bin_size + 1
          pairs$j_end <- (pairs$end2_bin - interval_start2) / bin_size + 1

          if (flip_axes) {
            tmp_i     <- pairs$i
            tmp_i_end <- pairs$i_end
            pairs$i     <- pairs$j
            pairs$i_end <- pairs$j_end
            pairs$j     <- tmp_i
            pairs$j_end <- tmp_i_end
          }

          pairs_by_file[[basename(bedpe_path)]] <- pairs
          k <- k + 1L
        }

        # Keep only files that actually contributed pairs
        pairs_by_file <- pairs_by_file[vapply(pairs_by_file, function(x) !is.null(x) && nrow(x) > 0, logical(1))]

        n_bedpe <- length(pairs_by_file)

        if (n_bedpe == 0L) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          bedpe <- pairs_by_file
        }
      }
    }

    ## CPM or AQuA factors

    mergeStats_A <- read.table(path_mgs, as.is = T)
    spikeVar <- ncol(mergeStats_A)
    hg_total1 <- as.numeric(mergeStats_A["valid_interaction_rmdup", 1])
    mm_total1 <- as.numeric(mergeStats_A["valid_interaction_rmdup", 2])
    total1 <- sum(hg_total1, mm_total1)
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1

    # set normalization factor for spike-in and non-spike-in data.
    # non-spike-in data defaults to cpm.
    # spike in data defaults to aqua.

    if (spikeVar == 1) {
      if (flag_norm == "blank") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        flag_norm <- "cpm"
      } else if (flag_norm == "none") {
        norm_factor1 <- 1
        aqua_factor1 <- 1
      } else if (flag_norm == "cpm") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      } else if (flag_norm == "aqua") {
        cat("\n\n--norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        flag_norm <- "cpm"
      }
    } else if (spikeVar == 2) {
      if (flag_norm == "blank") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
        flag_norm <- "aqua"
      } else if (flag_norm == "none") {
        norm_factor1 <- 1
        aqua_factor1 <- 1
      } else if (flag_norm == "cpm") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      } else if (flag_norm == "aqua") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
      }
    }

    cat("\n  factors\n")
    cat(sprintf("  norm_factor1: %f\n", norm_factor1))
    cat(sprintf("  aqua_factor1: %f\n", aqua_factor1))

    ## Process matrix:

    cat("\n  parameters\n")
    cat(sprintf("       norm: %s\n", flag_norm))
    cat(sprintf("        hic: %s\n", path_hic))
    cat(sprintf("   interval: %s\n", paste(new_interval_chr1, interval_start1, interval_end1, sep = ":")))
    cat(sprintf("   interval: %s\n", paste(new_interval_chr2, interval_start2, interval_end2, sep = ":")))
    cat(sprintf("       unit: %s\n", unit))
    cat(sprintf("   bin_size: %s\n", bin_size))

    sparMat1 <- straw(
      norm,
      path_hic,
      paste(new_interval_chr1, interval_start1, interval_end1, sep = ":"),
      paste(new_interval_chr2, interval_start2, interval_end2, sep = ":"),
      unit,
      bin_size
    )

    if (nrow(sparMat1) == 0) {
      cat("\n No contact values found for this interval\n")
      quit(save = "no", status = 1)
    }

    sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1

    denMat1 <- matrix(
      data = 0,
      nrow = length(seq(interval_start1, interval_end1, bin_size)),
      ncol = length(seq(interval_start2, interval_end2, bin_size))
    )

    rownames(denMat1) <- seq(interval_start1, interval_end1, bin_size)
    colnames(denMat1) <- seq(interval_start2, interval_end2, bin_size)

    # Populate denMat1
    for (i in 1:nrow(sparMat1)) {
      x_index <- as.character(sparMat1[i, "x"])
      y_index <- as.character(sparMat1[i, "y"])

      if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
        denMat1[x_index, y_index] <- sparMat1[i, "counts"]
      }
    }

    denMat1 <- t(denMat1)

    if (flip_axes) {
      denMat1 <- t(denMat1) # transpose back

      tmp_chr1 <- new_interval_chr1
      tmp_start1 <- interval_start1
      tmp_end1 <- interval_end1
      tmp_len1 <- interval_len1

      new_interval_chr1 <- new_interval_chr2
      interval_start1 <- interval_start2
      interval_end1 <- interval_end2
      interval_len1 <- interval_len2

      new_interval_chr2 <- tmp_chr1
      interval_start2 <- tmp_start1
      interval_end2 <- tmp_end1
      interval_len2 <- tmp_len1
    }

    if (flag_matrix) {

      # first region given is x, second is y
      row_coords <- seq(interval_start2, interval_end2, by = bin_size)
      col_coords <- seq(interval_start1, interval_end1, by = bin_size)

      row_labels <- sprintf("%s:%d:%d", new_interval_chr2, row_coords, row_coords + bin_size - 1L)
      col_labels <- sprintf("%s:%d:%d", new_interval_chr1, col_coords, col_coords + bin_size - 1L)

      if (length(row_labels) == nrow(denMat1)) rownames(denMat1) <- row_labels
      if (length(col_labels) == ncol(denMat1)) colnames(denMat1) <- col_labels

      file_name <- basename(output_name)
      out_file <- file.path(out_dir, file_name)
      write.table(denMat1, file = out_file, quote = FALSE, sep = "\t")
      quit(save = "no", status = 0)
    }

    cap <- calculate_cap(denMat1, max_cap, quant_cut, flag_inter)
    cap <- round(cap, 3)

    # Print the cap
    cat(sprintf("        cap: %s\n", cap))

    # Cap the values in the original matrix
    denMat1[denMat1 > cap] <- cap
    denMat1max <- denMat1

    color_ramp <- colorRampPalette(c("white", contact_color))(101)

    breakList <- seq(0, cap, by = cap / 100)

    C <- matrix(ncol = 2, nrow = length(breakList), data = 0)
    rownames(C) <- color_ramp
    colnames(C) <- c("rank", "breaks")
    C[, 1] <- 1:nrow(C)
    C[, 2] <- breakList

    A <- denMat1max

    if (flag_w == "blank") {
      w <- calculate_w_inter(plot_width, plot_height, A)
    }

    if (flag_w != "blank") {
      w <- as.numeric(flag_w)
      max_w <- calculate_w_inter(plot_width, plot_height, A)
      if (w > max_w) {
        cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
        w <- max_w
      }
    }

    cat(sprintf("plotting bin width: %s\n", round(w, 3)))

    render_inter_plot(
      output_name       = out_path,
      A_mat             = A,
      C_tbl             = C,
      breaks_vec        = breakList,
      bedpe             = bedpe,
      bedpe_color       = bedpe_color,
      w                 = w,
      ann_style         = ann_style,
      ann_custom        = ann_custom,
      ann_custom_colors = ann_custom_colors,
      genome_build      = genome_build,
      bin_size          = bin_size,
      new_interval_chr1 = new_interval_chr1,
      interval_start1   = interval_start1,
      interval_end1     = interval_end1,
      interval_len1     = interval_len1,
      new_interval_chr2 = new_interval_chr2,
      interval_start2   = interval_start2,
      interval_end2     = interval_end2,
      interval_len2     = interval_len2,
      flip_axes         = flip_axes,
      pairs             = if (exists("pairs") && is.data.frame(pairs)) pairs else NULL,
      flag_profiles     = flag_profiles
    )
  }
}

###########################################################################
###########################################################################
###                                                                     ###
###                   INTER-CHROMOSOMAL TWO-SAMPLE                      ###
###                                                                     ###
###########################################################################
###########################################################################

if (isTRUE(flag_inter)) {
  if (analysis_type == "two_sample") {
    path_hic_A           <- Args[3]
    path_hic_B           <- Args[4]
    path_mgs_A           <- Args[5]
    path_mgs_B           <- Args[6]
    norm                 <- Args[7]
    unit                 <- Args[8]
    bin_size             <- as.numeric(Args[9])
    prefix               <- Args[10]
    range_text           <- Args[11]
    out_dir              <- Args[12]
    interval_chr         <- Args[13]
    interval_start       <- as.numeric(Args[14])
    interval_end         <- as.numeric(Args[15])
    genome_build         <- Args[16]
    flag_profiles        <- as.logical(Args[17])
    contact_color        <- Args[18]
    ann_style            <- Args[19]
    ann_custom           <- Args[20]
    ann_custom_colors    <- Args[21]
    max_cap              <- Args[22]
    flag_norm            <- Args[23]
    quant_cut            <- as.numeric(Args[24])
    output_name          <- Args[25]
    bedpe                <- Args[26]
    bedpe_color          <- sprintf("#%s", Args[27])
    flag_w               <- Args[28]
    sample_dirA          <- Args[29]
    sample_dirB          <- Args[30]
    flag_matrix          <- as.logical(Args[31])
    inh_col_floor        <- sprintf("#%s", Args[32])
    inh_col_off          <- sprintf("#%s", Args[33])
    inh_col_on           <- sprintf("#%s", Args[34])
    inh_col_ceil         <- sprintf("#%s", Args[35])
    inherent             <- as.logical(Args[36])
    flag_plexus          <- as.logical(Args[37])
    flag_expand_viewport <- as.logical(Args[38])
    orig_chr1            <- Args[39]
    orig_chr2            <- Args[40]

    output_name <- basename(path.expand(output_name))
    out_path <- file.path(out_dir, output_name)

    flag_expand_viewport <- FALSE
  
    # Split the range_text at underscores
    range_components <- unlist(strsplit(range_text, "_"))

    # Extract the individual components
    interval_chr1 <- range_components[1]
    interval_start1 <- as.numeric(range_components[2])
    interval_end1 <- as.numeric(range_components[3])
    interval_chr2 <- range_components[4]
    interval_start2 <- as.numeric(range_components[5])
    interval_end2 <- as.numeric(range_components[6])


    if (any(grepl("chr", interval_chr1, ignore.case = T)) == FALSE) {
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save = "no", status = 1)
    }
    if (any(grepl("chr", interval_chr2, ignore.case = T)) == FALSE) {
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save = "no", status = 1)
    }

    interval_len1 <- interval_end1 - interval_start1
    interval_len2 <- interval_end2 - interval_start2

    if (ann_style != "none"){
      ann_reserve_bottom <- 35
      ann_reserve_left   <- 25
    } else {
      ann_reserve_bottom <- 16
      ann_reserve_left   <- 8
    }

    available_width  <- 100 - ann_reserve_left
    available_height <- 100 - ann_reserve_bottom

    dimensions <- calculate_plot_dimensions(
      interval_start1, interval_end1,
      interval_start2, interval_end2,
      bin_size,
      max_width  = available_width,
      max_height = available_height
    )

    plot_width <- dimensions$width
    plot_height <- dimensions$height

    flip_axes <- (orig_chr1 != interval_chr1)

    if (flip_axes) {
      temp_val <- plot_width
      plot_width <- plot_height
      plot_height <- temp_val
    }

    ## 'chr' prefix handler:

    hic_A_chroms <- readHicChroms(path_hic_A)$name
    hic_B_chroms <- readHicChroms(path_hic_B)$name

    if (sum(grepl("chr", hic_A_chroms)) > 0 && sum(grepl("chr", hic_B_chroms)) > 0) {
      flag_strip_chr <- "no"
    }
    if (sum(grepl("chr", hic_A_chroms)) > 0 && sum(grepl("chr", hic_B_chroms)) == 0) {
      flag_strip_chr <- "A"
    }
    if (sum(grepl("chr", hic_B_chroms)) > 0 && sum(grepl("chr", hic_A_chroms)) == 0) {
      flag_strip_chr <- "B"
    }
    if (sum(grepl("chr", hic_A_chroms)) == 0 && sum(grepl("chr", hic_B_chroms)) == 0) {
      flag_strip_chr <- "yes"
    }


    if (flag_strip_chr == "yes") {
      new_interval_chr1 <- sub("chr", "", interval_chr1)
      new_interval_chr2 <- sub("chr", "", interval_chr2)
    } else if (flag_strip_chr == "B") {
      new_interval_chr1 <- interval_chr1
      new_interval_chr2 <- sub("chr", "", interval_chr2)
    } else if (flag_strip_chr == "A") {
      new_interval_chr1 <- sub("chr", "", interval_chr1)
      new_interval_chr2 <- interval_chr2
    } else if (flag_strip_chr == "no") {
      new_interval_chr1 <- interval_chr1
      new_interval_chr2 <- interval_chr2
    }


    if (isFALSE(flag_matrix)) {
      if (bedpe != "FALSE") {
        # Accept a single space-joined string or a character vector
        bedpe_files <- if (length(bedpe) == 1L) strsplit(bedpe, "[[:space:]]+")[[1]] else bedpe
        bedpe_files <- bedpe_files[nzchar(bedpe_files)]

        # Pre-allocate list (names are basenames of files)
        pairs_by_file <- setNames(vector("list", length(bedpe_files)), basename(bedpe_files))

        n_bedpe <- length(bedpe_files)
        k <- 0L

        for (bedpe_path in bedpe_files) {
          pairs <- read.table(bedpe_path, as.is = TRUE, header = FALSE)[, 1:6]
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

          # Make sure chromosome names have 'chr' prefix
          pairs$chr1 <- ifelse(grepl("chr", pairs$chr1), pairs$chr1, paste0("chr", pairs$chr1))
          pairs$chr2 <- ifelse(grepl("chr", pairs$chr2), pairs$chr2, paste0("chr", pairs$chr2))

          # Filter: both anchors must be fully contained in the two intervals
          # (allow either orientation in the file)
          pairs <- pairs[
            (
              (pairs$chr1 == new_interval_chr1 &
              pairs$start1 >= interval_start1 & pairs$end1 <= interval_end1 &
              pairs$chr2 == new_interval_chr2 &
              pairs$start2 >= interval_start2 & pairs$end2 <= interval_end2) |
              (pairs$chr1 == new_interval_chr2 &
              pairs$start1 >= interval_start2 & pairs$end1 <= interval_end2 &
              pairs$chr2 == new_interval_chr1 &
              pairs$start2 >= interval_start1 & pairs$end2 <= interval_end1)
            ),
          ]

          if (nrow(pairs) == 0) {
            next
          }

          # Reorder columns so chr1 always corresponds to interval 1
          for (idx in 1:nrow(pairs)) {
            if (pairs[idx, "chr1"] != new_interval_chr1) {
              pairs[idx, c("chr1", "start1", "end1", "chr2", "start2", "end2")] <-
                pairs[idx, c("chr2", "start2", "end2", "chr1", "start1", "end1")]
            }
          }
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

          # Adjust end coordinates for half-open binning
          pairs$end1 <- ifelse(pairs$end1 > pairs$start1, pairs$end1 - 1, pairs$end1)
          pairs$end2 <- ifelse(pairs$end2 > pairs$start2, pairs$end2 - 1, pairs$end2)

          # Convert coordinates to bins
          pairs$start1_bin <- as.integer(floor(pairs$start1 / bin_size) * bin_size)
          pairs$end1_bin   <- as.integer(floor(pairs$end1 / bin_size) * bin_size)
          pairs$start2_bin <- as.integer(floor(pairs$start2 / bin_size) * bin_size)
          pairs$end2_bin   <- as.integer(floor(pairs$end2 / bin_size) * bin_size)

          # Calculate plotting coordinates
          pairs$i     <- (pairs$start1_bin - interval_start1) / bin_size + 1
          pairs$j     <- (pairs$start2_bin - interval_start2) / bin_size + 1
          pairs$i_end <- (pairs$end1_bin - interval_start1) / bin_size + 1
          pairs$j_end <- (pairs$end2_bin - interval_start2) / bin_size + 1

          if (flip_axes) {
            tmp_i     <- pairs$i
            tmp_i_end <- pairs$i_end
            pairs$i     <- pairs$j
            pairs$i_end <- pairs$j_end
            pairs$j     <- tmp_i
            pairs$j_end <- tmp_i_end
          }

          pairs_by_file[[basename(bedpe_path)]] <- pairs
          k <- k + 1L
        }

        # Keep only files that actually contributed pairs
        pairs_by_file <- pairs_by_file[vapply(pairs_by_file, function(x) !is.null(x) && nrow(x) > 0, logical(1))]

        n_bedpe <- length(pairs_by_file)

        if (n_bedpe == 0L) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          bedpe <- pairs_by_file
        }
      }
    }

    ## CPM or AQuA factors

    mergeStats_A <- read.table(path_mgs_A, as.is = T)
    spikeVar_A <- ncol(mergeStats_A)
    hg_total1 <- as.numeric(mergeStats_A["valid_interaction_rmdup", 1])
    mm_total1 <- as.numeric(mergeStats_A["valid_interaction_rmdup", 2])
    total1 <- sum(hg_total1, mm_total1)
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1

    mergeStats_B <- read.table(path_mgs_B, as.is = T)
    spikeVar_B <- ncol(mergeStats_B)
    hg_total2 <- as.numeric(mergeStats_B["valid_interaction_rmdup", 1])
    mm_total2 <- as.numeric(mergeStats_B["valid_interaction_rmdup", 2])
    total2 <- sum(hg_total2, mm_total2)
    norm_factor2 <- 1000000 / total2
    aqua_factor2 <- hg_total2 / mm_total2

    # set normalization factor for spike-in and non-spike-in data.
    # non-spike-in data defaults to cpm
    # spike-in data defaults to aqua

    if (spikeVar_A == 1 || spikeVar_B == 1) {
      if (flag_norm == "blank") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- 1
        flag_norm <- "cpm"
      } else if (flag_norm == "none") {
        norm_factor1 <- 1
        aqua_factor1 <- 1
        norm_factor2 <- 1
        aqua_factor2 <- 1
      } else if (flag_norm == "cpm") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- 1
      } else if (flag_norm == "aqua") {
        cat("\n\n--norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- 1
        flag_norm <- "cpm"
      }
    } else if (spikeVar_A == 2 || spikeVar_B == 2) {
      if (flag_norm == "blank") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- aqua_factor2
        flag_norm <- "aqua"
      } else if (flag_norm == "none") {
        norm_factor1 <- 1
        aqua_factor1 <- 1
        norm_factor2 <- 1
        aqua_factor2 <- 1
      } else if (flag_norm == "cpm") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- 1
      } else if (flag_norm == "aqua") {
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- aqua_factor2
      }
    }

    cat("\n  factors\n")
    cat(sprintf("  norm_factor1: %f\n", norm_factor1))
    cat(sprintf("  norm_factor2: %f\n", norm_factor2))
    cat(sprintf("  aqua_factor1: %f\n", aqua_factor1))
    cat(sprintf("  aqua_factor2: %f\n", aqua_factor2))

    ## Process matrix:

    cat("\n  parameters\n")
    cat(sprintf("       norm: %s\n", flag_norm))
    cat(sprintf("      hic_A: %s\n", path_hic_A))
    cat(sprintf("      hic_B: %s\n", path_hic_B))
    cat(sprintf("   interval: %s\n", paste(new_interval_chr1, interval_start1, interval_end1, sep = ":")))
    cat(sprintf("   interval: %s\n", paste(new_interval_chr2, interval_start2, interval_end2, sep = ":")))
    cat(sprintf("       unit: %s\n", unit))
    cat(sprintf("   bin_size: %s\n", bin_size))


    sparMat1 <- straw(
      norm,
      path_hic_A,
      paste(new_interval_chr1, interval_start1, interval_end1, sep = ":"),
      paste(new_interval_chr2, interval_start2, interval_end2, sep = ":"),
      unit,
      bin_size
    )

    if (nrow(sparMat1) == 0) {
      cat("\n No contact values found for this interval in sample A\n")
      quit(save = "no", status = 1)
    }

    sparMat2 <- straw(
      norm,
      path_hic_B,
      paste(new_interval_chr1, interval_start1, interval_end1, sep = ":"),
      paste(new_interval_chr2, interval_start2, interval_end2, sep = ":"),
      unit,
      bin_size
    )

    if (nrow(sparMat2) == 0) {
      cat("\n No contact values found for this interval in sample B\n")
      quit(save = "no", status = 1)
    }

    sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1
    sparMat2$counts <- sparMat2$counts * norm_factor2 * aqua_factor2

    denMat1 <- matrix(
      data = 0,
      nrow = length(seq(interval_start1, interval_end1, bin_size)),
      ncol = length(seq(interval_start2, interval_end2, bin_size))
    )
    rownames(denMat1) <- seq(interval_start1, interval_end1, bin_size)
    colnames(denMat1) <- seq(interval_start2, interval_end2, bin_size)

    denMat2 <- matrix(
      data = 0,
      nrow = length(seq(interval_start1, interval_end1, bin_size)),
      ncol = length(seq(interval_start2, interval_end2, bin_size))
    )
    rownames(denMat2) <- seq(interval_start1, interval_end1, bin_size)
    colnames(denMat2) <- seq(interval_start2, interval_end2, bin_size)

    # Populate denMat1
    for (i in 1:nrow(sparMat1)) {
      x_index <- as.character(sparMat1[i, "x"])
      y_index <- as.character(sparMat1[i, "y"])

      if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
        denMat1[x_index, y_index] <- sparMat1[i, "counts"]
      }
    }

    # Populate denMat2
    for (i in 1:nrow(sparMat2)) {
      x_index <- as.character(sparMat2[i, "x"])
      y_index <- as.character(sparMat2[i, "y"])

      if (x_index %in% rownames(denMat2) && y_index %in% colnames(denMat2)) {
        denMat2[x_index, y_index] <- sparMat2[i, "counts"]
      }
    }

    denMat1 <- t(denMat1)
    denMat2 <- t(denMat2)

    if (nrow(denMat1) != nrow(denMat2)) {
      intersect <- intersect(rownames(denMat1), rownames(denMat2))
      denMat1 <- denMat1[intersect, intersect]
      denMat2 <- denMat2[intersect, intersect]
    }

    denDelta <- denMat2 - denMat1

    if (flip_axes) {
      denDelta <- t(denDelta) # transpose back

      tmp_chr1 <- new_interval_chr1
      tmp_start1 <- interval_start1
      tmp_end1 <- interval_end1
      tmp_len1 <- interval_len1

      new_interval_chr1 <- new_interval_chr2
      interval_start1 <- interval_start2
      interval_end1 <- interval_end2
      interval_len1 <- interval_len2

      new_interval_chr2 <- tmp_chr1
      interval_start2 <- tmp_start1
      interval_end2 <- tmp_end1
      interval_len2 <- tmp_len1
    }

    if (flag_matrix) {

      # first region given is x, second is y
      row_coords <- seq(interval_start2, interval_end2, by = bin_size)
      col_coords <- seq(interval_start1, interval_end1, by = bin_size)

      row_labels <- sprintf("%s:%d:%d", new_interval_chr2, row_coords, row_coords + bin_size - 1L)
      col_labels <- sprintf("%s:%d:%d", new_interval_chr1, col_coords, col_coords + bin_size - 1L)

      if (length(row_labels) == nrow(denMat1)) rownames(denMat1) <- row_labels
      if (length(col_labels) == ncol(denMat1)) colnames(denMat1) <- col_labels

      file_name <- basename(output_name)
      out_file <- file.path(out_dir, file_name)
      write.table(denMat1, file = out_file, quote = FALSE, sep = "\t")
      quit(save = "no", status = 0)
    }

    # Max_cap and quant_cut handling
    cap <- calculate_cap(denDelta, max_cap, quant_cut, flag_inter)
    cap <- round(cap, 3)

    denDelta[denDelta > cap] <- cap
    denDelta[denDelta < -cap] <- -cap

    cat(sprintf("        cap: %s\n", cap))

    if (isFALSE(grepl("-", contact_color))) {
      cat("\nPlease separate the two delta colors with '-' only! \n")
      quit(save = "no", status = 1)
    }

    color_neg <- sprintf("#%s", unlist(strsplit(contact_color, "-"))[1])
    color_pos <- sprintf("#%s", unlist(strsplit(contact_color, "-"))[2])

    # Create the break list
    breakList <- seq(
      from = -cap,
      to = cap,
      by = cap / 100
    )

    color_ramp <- colorRampPalette(c(color_neg, "white", color_pos))(length(breakList))

    C <- matrix(ncol = 2, nrow = length(breakList), data = 0)
    rownames(C) <- color_ramp
    colnames(C) <- c("rank", "breaks")
    C[, 1] <- 1:nrow(C)
    C[, 2] <- breakList

    A <- denDelta

    if (flag_w == "blank") {
      w <- calculate_w_inter(plot_width, plot_height, A)
    }

    if (flag_w != "blank") {
      w <- as.numeric(flag_w)
      max_w <- calculate_w_inter(plot_width, plot_height, A)
      if (w > max_w) {
        cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
        w <- max_w
      }
    }

    cat(sprintf("plotting bin width: %s\n", round(w, 3)))

    render_inter_plot(
      output_name       = out_path,
      A_mat             = A,
      C_tbl             = C,
      breaks_vec        = breakList,
      bedpe             = bedpe,
      bedpe_color       = bedpe_color,
      w                 = w,
      ann_style         = ann_style,
      ann_custom        = ann_custom,
      ann_custom_colors = ann_custom_colors,
      genome_build      = genome_build,
      bin_size          = bin_size,
      new_interval_chr1 = new_interval_chr1,
      interval_start1   = interval_start1,
      interval_end1     = interval_end1,
      interval_len1     = interval_len1,
      new_interval_chr2 = new_interval_chr2,
      interval_start2   = interval_start2,
      interval_end2     = interval_end2,
      interval_len2     = interval_len2,
      flip_axes         = flip_axes,
      pairs             = if (exists("pairs") && is.data.frame(pairs)) pairs else NULL,
      flag_profiles     = flag_profiles
    )
  }
}
