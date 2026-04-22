suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(strawr)
})
options(scipen = 999)
options(help_type = "text")
Args <- commandArgs(trailingOnly = TRUE)

# ==========================================================================
#                              FUNCTIONS
# ==========================================================================

# Resolve chromosome name between bedpe and .hic conventions
match_chr_name <- function(chr_bedpe, hic_chroms) {
  if (chr_bedpe %in% hic_chroms) return(chr_bedpe)
  stripped <- sub("^chr", "", chr_bedpe)
  if (stripped %in% hic_chroms) return(stripped)
  prefixed <- paste0("chr", chr_bedpe)
  if (prefixed %in% hic_chroms) return(prefixed)
  return(chr_bedpe)
}

pre_check_power_law <- function(sample_dir, bin_size) {
  power_law_path <- list.files(
    sample_dir,
    pattern = paste0("inherentStats", ".txt"),
    full.names = TRUE
  )
  if (length(power_law_path) == 0) {
    cat("\n--inherent unavailable for", sample_dir, "\n")
    quit(save="no")
  }
  power_laws <- read.table(power_law_path, as.is = TRUE, skip = 1)
  colnames(power_laws) <- c("off", "on", "res")
  power_laws$res <- as.numeric(power_laws$res)
  
  if (!bin_size %in% power_laws$res) {
    cat("\nNo power law data available for resolution", bin_size, "at", sample_dir,"\n")
    quit(save="no")
  }
}

# Read power law data once and cache
read_power_law_data <- function(sample_dir, bin_size) {
  power_law_path <- list.files(
    sample_dir,
    pattern = paste0("inherentStats", ".txt"),
    full.names = TRUE
  )
  
  power_laws <- read.table(power_law_path, as.is = TRUE, skip = 1)
  colnames(power_laws) <- c("off", "on", "res")
  power_laws <- power_laws[! power_laws$off %in% "off",]
  power_laws$off <- as.numeric(power_laws$off)
  power_laws$on <- as.numeric(power_laws$on)
  power_laws$res <- as.numeric(power_laws$res)
  power_laws <- power_laws[power_laws$res == bin_size,]
  
  return(power_laws)
}

inherent_normalization <- function(records, power_laws, resolution) {
  bin_index <- (records$y - records$x) / resolution + 1
  bin_index <- pmin(bin_index, nrow(power_laws))
  bin_index <- pmax(bin_index, 1L)
  off_vals <- power_laws$off[bin_index]
  on_vals  <- power_laws$on[bin_index]
  denom <- on_vals - off_vals
  # If on == off, set normalized value to 0 instead of Inf/NaN
  safe <- denom != 0
  records$counts <- ifelse(safe,
                           (records$counts - off_vals) / denom,
                           0)
  return(records)
}

build_inherent_colors <- function(min_cap_value, cap,
                                  col_floor = "#ffffff",
                                  col_off   = "#d4e4fb",
                                  col_on    = "#ff0000",
                                  col_ceil  = "#ff8d4a") {
  colors_n <- 10
  colors_0_1 <- colorRampPalette(c(col_off, col_on))(colors_n)
  colors_1_2 <- colorRampPalette(c(col_on, col_ceil))(colors_n)
  color_ramp <- c(col_floor, colors_0_1, colors_1_2)

  breaks <- c(0, (1:colors_n) / colors_n, (1:colors_n) / colors_n + 1)

  # Rescale breaks to [0,1] relative to the plot's cap range
  values <- (breaks - min_cap_value) / (cap - min_cap_value)
  values <- pmin(pmax(values, 0), 1)

  list(colors = color_ramp, values = values, breaks = breaks)
}

# Sort key that extracts the numeric part of the chromosome name
chr_sort_key <- function(chr_name) {
  stripped <- sub("^chr", "", chr_name)
  if (stripped == "X")  return(23)
  if (stripped == "Y")  return(24)
  if (stripped == "M")  return(25)
  if (stripped == "MT") return(25)
  num <- suppressWarnings(as.numeric(stripped))
  if (is.na(num)) return(999)
  return(num)
}

# Build a hash-bucket environment from straw records sorted by x-bin
build_bucket_env <- function(rec_xb, rec_yb, rec_cts) {
  ord     <- order(rec_xb)
  rec_xb  <- rec_xb[ord]
  rec_yb  <- rec_yb[ord]
  rec_cts <- rec_cts[ord]

  r <- rle(rec_xb)
  bucket_end   <- cumsum(r$lengths)
  bucket_start <- bucket_end - r$lengths + 1L

  bucket_env <- new.env(hash = TRUE, size = length(r$values))
  for (k in seq_along(r$values)) {
    s <- bucket_start[k]; e <- bucket_end[k]
    assign(as.character(r$values[k]),
           list(yb = rec_yb[s:e], xb = rec_xb[s:e], cts = rec_cts[s:e]),
           envir = bucket_env)
  }
  bucket_env
}

# Accumulate contact values for one loop into the local APA matrix
# via hash-bucket lookup on x-bins, then filtering y-bins
accumulate_loop <- function(i_rs, i_re, i_cs, i_ce, bucket_env,
                            mat_size, swap = FALSE) {
  xbins_needed <- seq.int(i_rs, i_re)

  all_yb  <- vector("list", mat_size)
  all_xb  <- vector("list", mat_size)
  all_cts <- vector("list", mat_size)
  count <- 0L
  for (xb in xbins_needed) {
    key <- as.character(xb)
    if (exists(key, envir = bucket_env, inherits = FALSE)) {
      count <- count + 1L
      b <- get(key, envir = bucket_env, inherits = FALSE)
      all_yb[[count]]  <- b$yb
      all_xb[[count]]  <- b$xb
      all_cts[[count]] <- b$cts
    }
  }
  if (count == 0L) return(NULL)

  sub_yb  <- unlist(all_yb[1:count],  use.names = FALSE)
  sub_xb  <- unlist(all_xb[1:count],  use.names = FALSE)
  sub_cts <- unlist(all_cts[1:count], use.names = FALSE)

  y_ok <- sub_yb >= i_cs & sub_yb <= i_ce
  if (!any(y_ok)) return(NULL)

  ri  <- sub_xb[y_ok] - i_rs + 1L
  ci  <- sub_yb[y_ok] - i_cs + 1L
  cts <- sub_cts[y_ok]

  if (swap) {
    tmp <- ri; ri <- ci; ci <- tmp
  }

  valid <- ri >= 1L & ri <= mat_size & ci >= 1L & ci <= mat_size
  if (!any(valid)) return(NULL)

  list(idx = cbind(ri[valid], ci[valid]), cts = cts[valid])
}

process_cis_group <- function(indices, chr1_q_vec, rs_bin, re_bin, cs_bin, ce_bin,
                              hic_path, resolution, mat_size,
                              flag_inherent, power_laws) {
  local_agg     <- matrix(0, nrow = mat_size, ncol = mat_size)
  local_counted <- 0L
  local_skipped <- 0L

  chr_q <- chr1_q_vec[indices[1]]

  grp_rs <- min(rs_bin[indices]) * resolution
  grp_re <- max(re_bin[indices]) * resolution
  grp_cs <- min(cs_bin[indices]) * resolution
  grp_ce <- max(ce_bin[indices]) * resolution

  region1 <- paste0(chr_q, ":", max(0L, grp_rs), ":", grp_re)
  region2 <- paste0(chr_q, ":", max(0L, grp_cs), ":", grp_ce)

  records <- tryCatch(
    straw("NONE", hic_path, region1, region2, "BP", resolution),
    error = function(e) NULL
  )
  if (is.null(records) || nrow(records) == 0) {
    return(list(agg = local_agg, counted = 0L, skipped = length(indices)))
  }

  if (flag_inherent) {
    records <- inherent_normalization(records, power_laws, resolution)
  }

  bucket_env <- build_bucket_env(
    as.integer(records$x / resolution),
    as.integer(records$y / resolution),
    records$counts
  )

  for (j in seq_along(indices)) {
    i <- indices[j]
    result <- accumulate_loop(rs_bin[i], re_bin[i], cs_bin[i], ce_bin[i],
                              bucket_env, mat_size, swap = FALSE)
    if (is.null(result)) {
      local_skipped <- local_skipped + 1L
    } else {
      local_agg[result$idx] <- local_agg[result$idx] + result$cts
      local_counted <- local_counted + 1L
    }
  }

  list(agg = local_agg, counted = local_counted, skipped = local_skipped)
}

process_trans_group <- function(indices, chr1_q_vec, chr2_q_vec,
                                rs_bin, re_bin, cs_bin, ce_bin,
                                hic_path, resolution, mat_size) {
  local_agg     <- matrix(0, nrow = mat_size, ncol = mat_size)
  local_counted <- 0L
  local_skipped <- 0L

  chr1_q <- chr1_q_vec[indices[1]]
  chr2_q <- chr2_q_vec[indices[1]]

  # Determine if we need to swap anchors to match straw's axis convention
  swap_anchors <- (chr_sort_key(chr1_q) > chr_sort_key(chr2_q))

  if (swap_anchors) {
    loc_rs_bin <- cs_bin[indices]    # anchor2 bins -> straw x (rows)
    loc_re_bin <- ce_bin[indices]
    loc_cs_bin <- rs_bin[indices]    # anchor1 bins -> straw y (cols)
    loc_ce_bin <- re_bin[indices]
    query_chr1 <- chr2_q             # numerically lower chr
    query_chr2 <- chr1_q
  } else {
    loc_rs_bin <- rs_bin[indices]
    loc_re_bin <- re_bin[indices]
    loc_cs_bin <- cs_bin[indices]
    loc_ce_bin <- ce_bin[indices]
    query_chr1 <- chr1_q
    query_chr2 <- chr2_q
  }

  grp_rs <- min(loc_rs_bin) * resolution
  grp_re <- max(loc_re_bin) * resolution
  grp_cs <- min(loc_cs_bin) * resolution
  grp_ce <- max(loc_ce_bin) * resolution

  region1 <- paste0(query_chr1, ":", max(0L, grp_rs), ":", grp_re)
  region2 <- paste0(query_chr2, ":", max(0L, grp_cs), ":", grp_ce)

  records <- tryCatch(
    straw("NONE", hic_path, region1, region2, "BP", resolution),
    error = function(e) NULL
  )
  if (is.null(records) || nrow(records) == 0) {
    return(list(agg = local_agg, counted = 0L, skipped = length(indices)))
  }

  bucket_env <- build_bucket_env(
    as.integer(records$x / resolution),
    as.integer(records$y / resolution),
    records$counts
  )

  for (j in seq_along(indices)) {
    result <- accumulate_loop(loc_rs_bin[j], loc_re_bin[j],
                              loc_cs_bin[j], loc_ce_bin[j],
                              bucket_env, mat_size, swap = swap_anchors)
    if (is.null(result)) {
      local_skipped <- local_skipped + 1L
    } else {
      local_agg[result$idx] <- local_agg[result$idx] + result$cts
      local_counted <- local_counted + 1L
    }
  }

  list(agg = local_agg, counted = local_counted, skipped = local_skipped)
}

run_apa <- function(hic_path, bedpe, resolution, window,
                    flag_inherent = FALSE, power_laws = NULL) {
  mat_size  <- 2L * window + 1L
  aggregate <- matrix(0, nrow = mat_size, ncol = mat_size)

  hic_chrom_vec <- readHicChroms(hic_path)$name

  # Vectorize bedpe columns
  chr1_vec <- as.character(bedpe[, 1])
  s1_vec   <- as.numeric(bedpe[, 2])
  e1_vec   <- as.numeric(bedpe[, 3])
  chr2_vec <- as.character(bedpe[, 4])
  s2_vec   <- as.numeric(bedpe[, 5])
  e2_vec   <- as.numeric(bedpe[, 6])

  # Resolve chromosome names
  chr1_q_vec <- vapply(chr1_vec, match_chr_name, character(1),
                       hic_chroms = hic_chrom_vec, USE.NAMES = FALSE)
  chr2_q_vec <- vapply(chr2_vec, match_chr_name, character(1),
                       hic_chroms = hic_chrom_vec, USE.NAMES = FALSE)

  # Deduplicate by exact interval identity (matches Juicer HashSet)
  dedup_key <- paste(chr1_q_vec, s1_vec, e1_vec, chr2_q_vec, s2_vec, e2_vec, sep = "_")
  dup_mask  <- duplicated(dedup_key)

  if (any(dup_mask)) {
    keep <- !dup_mask
    chr1_q_vec <- chr1_q_vec[keep]; chr2_q_vec <- chr2_q_vec[keep]
    s1_vec <- s1_vec[keep]; e1_vec <- e1_vec[keep]
    s2_vec <- s2_vec[keep]; e2_vec <- e2_vec[keep]
  }

  # Snap midpoints to resolution grid
  mid1_vec <- floor(((s1_vec + e1_vec) / 2) / resolution) * resolution
  mid2_vec <- floor(((s2_vec + e2_vec) / 2) / resolution) * resolution

  # Window boundaries in bin coordinates
  win_res <- as.integer(window)
  rs_bin  <- as.integer(mid1_vec / resolution) - win_res
  re_bin  <- as.integer(mid1_vec / resolution) + win_res
  cs_bin  <- as.integer(mid2_vec / resolution) - win_res
  ce_bin  <- as.integer(mid2_vec / resolution) + win_res

  # Split into cis and trans groups
  is_cis    <- (chr1_q_vec == chr2_q_vec)
  chr_pair  <- paste0(chr1_q_vec, "_", chr2_q_vec)

  cis_indices  <- which(is_cis)
  trans_indices <- which(!is_cis)

  # Group by chromosome pair within cis and trans separately
  results <- list()

  if (length(cis_indices) > 0) {
    cis_groups <- split(cis_indices, chr_pair[cis_indices])
    results <- c(results, lapply(cis_groups, process_cis_group,
      chr1_q_vec = chr1_q_vec, rs_bin = rs_bin, re_bin = re_bin,
      cs_bin = cs_bin, ce_bin = ce_bin,
      hic_path = hic_path, resolution = resolution, mat_size = mat_size,
      flag_inherent = flag_inherent, power_laws = power_laws))
  }

  if (length(trans_indices) > 0) {
    trans_groups <- split(trans_indices, chr_pair[trans_indices])
    results <- c(results, lapply(trans_groups, process_trans_group,
      chr1_q_vec = chr1_q_vec, chr2_q_vec = chr2_q_vec,
      rs_bin = rs_bin, re_bin = re_bin,
      cs_bin = cs_bin, ce_bin = ce_bin,
      hic_path = hic_path, resolution = resolution, mat_size = mat_size))
  }

  # Reduce
  counted <- 0L; skipped <- 0L
  for (res_item in results) {
    aggregate <- aggregate + res_item$agg
    counted   <- counted + res_item$counted
    skipped   <- skipped + res_item$skipped
  }
  return(list(aggregate = aggregate, counted = counted, skipped = skipped))
}

# Compute intrachromosomal distance for each loop
compute_distance <- function(i, pairs) {
  if (pairs[i, 1] == pairs[i, 4]) return(abs(pairs[i, 5] - pairs[i, 2]))
  return(NA)
}

# Compute peak value and peak-to-corner ratios
compute_peak_scores <- function(apa, res) {
  apa_m <- as.matrix(apa)
  n   <- nrow(apa_m)
  mid <- (n + 1) %/% 2
  q   <- if (res <= 5000) 3 else 6

  peak <- as.numeric(apa_m[mid, mid])
  corners <- list(
    p2ll = apa_m[(n - q + 1):n,         1:q],
    p2ul = apa_m[1:q,                    1:q],
    p2ur = apa_m[1:q,          (n - q + 1):n],
    p2lr = apa_m[(n - q + 1):n, (n - q + 1):n]
  )
  scores <- lapply(corners, function(region) peak / mean(region, na.rm = TRUE))
  scores$peak <- peak
  return(scores)
}

# Build annotation data.frames for APA score overlay on heatmaps
build_box_annotations <- function(matrix_norm, res, peak_scores,
                                  peak_x_shift = 3, peak_y_shift = 1.2) {
  n   <- nrow(matrix_norm)
  mid <- (n + 1) %/% 2
  q   <- if (res <= 5000) 3 else 6

  region_names <- c("P2LL", "P2UL", "P2UR", "P2LR", "peak")
  scores <- c(peak_scores$p2ll, peak_scores$p2ul,
              peak_scores$p2ur, peak_scores$p2lr, peak_scores$peak)

  box_data <- data.frame(
    region     = region_names,
    xmin       = c(1, 1, n - q + 1, n - q + 1, mid) - 0.5,
    xmax       = c(q, q, n, n, mid) + 0.5,
    ymin       = n + 1 - c(n, q, q, n, mid) - 0.5,
    ymax       = n + 1 - c(n - q + 1, 1, 1, n - q + 1, mid) + 0.5,
    label_head = paste0(region_names, ":"),
    label_val  = sprintf("%.3f", scores)
  )
  box_data$x_label <- (box_data$xmin + box_data$xmax) / 2
  box_data$y_label <- (box_data$ymin + box_data$ymax) / 2
  box_data$y_head  <- box_data$y_label + 0.35
  box_data$y_val   <- box_data$y_label - 0.35

  is_peak <- box_data$region == "peak"
  box_data$label_one_line <- NA_character_
  box_data$label_one_line[is_peak] <- sprintf("peak: %s", box_data$label_val[is_peak])
  box_data$x_label[is_peak] <- box_data$x_label[is_peak] + peak_x_shift
  box_data$y_label[is_peak] <- box_data$y_label[is_peak] + peak_y_shift
  box_data
}

# Choose label format based on scale magnitude
get_label_fmt <- function(cap) {
  if (cap < 1)         scales::label_number(accuracy = 0.001, big.mark = "")
  else if (cap < 1000) scales::label_number(accuracy = 0.1,   big.mark = "")
  else                 scales::label_number(accuracy = 1,     big.mark = "")
}

compute_norm_factors <- function(mergeStats, flag_norm) {
  hg_total   <- as.numeric(mergeStats["valid_interaction_rmdup", 1])

  if (flag_norm == "none") {
    norm_factor <- 1
    aqua_factor <- 1
  } else if (flag_norm == "cpm") {
    if (ncol(mergeStats) >= 2) {
      mm_total    <- as.numeric(mergeStats["valid_interaction_rmdup", 2])
      norm_factor <- 1000000 / (hg_total + mm_total)
    } else {
      norm_factor <- 1000000 / hg_total
    }
    aqua_factor <- 1
  } else if (flag_norm == "aqua") {
    mm_total    <- as.numeric(mergeStats["valid_interaction_rmdup", 2])
    total       <- hg_total + mm_total
    norm_factor <- 1000000 / total
    aqua_factor <- hg_total / mm_total
  }

  list(norm_factor = norm_factor, aqua_factor = aqua_factor)
}

# Build APA heatmap ggplot
build_apa_plot <- function(matrix_long, color_ramp, breakList, cap, min_cap_value,
                           box_data = NULL, apa_scores = FALSE,
                           flag_inherent = FALSE) {

  if (flag_inherent) {
    inh <- build_inherent_colors(min_cap_value, cap)
    fill_scale <- scale_fill_gradientn(
      colors = inh$colors,
      values = inh$values,
      limits = c(min_cap_value, cap),
      oob    = scales::squish,
      labels = get_label_fmt(cap),
      name   = NULL
    )
  } else {
    fill_scale <- scale_fill_gradientn(
      colors = color_ramp,
      breaks = breakList,
      limits = c(min_cap_value, cap),
      oob    = scales::squish,
      labels = get_label_fmt(cap),
      name   = NULL
    )
  }

  p <- ggplot(matrix_long, aes(x = col, y = row, fill = value)) +
    geom_tile(color = NA) +
    fill_scale +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(matrix_long$row)), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y       = element_blank(),
      axis.title        = element_blank(),
      axis.ticks        = element_blank(),
      panel.grid        = element_blank(),
      plot.title        = element_text(hjust = 0.5),
      legend.key.height = unit(1.2, "cm")
    )

  if (apa_scores && !is.null(box_data)) {
    box_corners <- subset(box_data, region != "peak")
    box_peak    <- subset(box_data, region == "peak")
    p <- p +
      geom_rect(data = box_data,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE, color = "black", fill = NA, linewidth = 0.6) +
      geom_text(data = box_corners, aes(x = x_label, y = y_head, label = label_head),
        inherit.aes = FALSE, color = "black", size = 3.0) +
      geom_text(data = box_corners, aes(x = x_label, y = y_val, label = label_val),
        inherit.aes = FALSE, color = "black", size = 2.5) +
      geom_text(data = box_peak, aes(x = x_label, y = y_label, label = label_one_line),
        inherit.aes = FALSE, color = "black", size = 3.0)
  }
  p
}

# Convert matrix to long format for ggplot
matrix_to_long <- function(mat, cols) {
  colnames(mat) <- cols
  as.data.frame(mat) %>%
    mutate(row = seq_len(nrow(mat))) %>%
    pivot_longer(-row, names_to = "col", values_to = "value") %>%
    mutate(col = factor(col, levels = cols),
           row = factor(row, levels = seq_len(nrow(mat))))
}

build_axis_labels <- function(res, win_size) {
  res_kb <- res / 1000
  c(paste0("-", rev(seq(res_kb, res_kb * win_size, res_kb)), "kb"),
    "0",
    paste0(seq(res_kb, res_kb * win_size, res_kb), "kb"))
}

# Plot bedpe loop-size distribution histogram
plot_size_distribution <- function(pairs, sample_label, out_path,
                                   width = 7, height = 6) {
  distances <- sapply(seq_len(nrow(pairs)), compute_distance, pairs = pairs)
  distances <- distances[!is.na(distances)]

  # If all pairs are trans (no cis distances), skip the histogram
  if (length(distances) == 0) {
    cat("\nAll pairs are inter-chromosomal; skipping size distribution histogram.\n")
    return(invisible(NULL))
  }

  distribution <- data.frame(distance = round(distances))

  n <- nrow(distribution)
  bins <- if (n < 50) 10 else if (n <= 500) 30 else min(ceiling(sqrt(n)), 200)

  test_con <- tryCatch(file(out_path, open = "ab"), error = function(e) NULL)
  if (is.null(test_con)) {
    stop(sprintf(
      "Cannot write to '%s' -- the file may be open in another application. Please close it and try again.",
      out_path
    ))
  }
  close(test_con)

  pdf(out_path, width = width, height = height, useDingbats = FALSE)
  p <- ggplot(distribution, aes(x = distance)) +
    geom_histogram(bins = bins, color = "#757575") +
    xlab("Distance between bedpe feet") + ylab("Count") +
    xlim(0, ceiling(quantile(distribution$distance)[4])) +
    ggtitle(sample_label) + theme_classic()
  suppressWarnings(suppressMessages(print(p)))
  dev.off()
}

# Print APA scores to console
print_apa_scores <- function(label, scores) {
  cat(sprintf("\nAPA scores for %s\n", label))
  cat(sprintf("peak           <- %.3f\n", scores$peak))
  cat(sprintf("P2LL           <- %.3f\n", scores$p2ll))
  cat(sprintf("P2UL           <- %.3f\n", scores$p2ul))
  cat(sprintf("P2UR           <- %.3f\n", scores$p2ur))
  cat(sprintf("P2LR           <- %.3f\n", scores$p2lr))
}

# ==========================================================================
#                          ONE-SAMPLE/COHORT ANALYSIS
# ==========================================================================
if (length(Args) == 17) {

  hic_paths      <- strsplit(Args[1], ",")[[1]]
  sample_labels  <- strsplit(Args[2], ",")[[1]]
  max_cap        <- Args[3]
  max_cap_delta  <- Args[4]
  out_file       <- Args[5]
  win_size       <- as.numeric(Args[6])
  res            <- as.numeric(Args[7])
  pairs_path     <- Args[8]
  flag_norm      <- Args[9]
  stat_paths     <- strsplit(Args[10], ",")[[1]]
  flag_loop_norm <- Args[11]
  num_loops      <- as.numeric(Args[12])
  out_dir        <- Args[13]
  min_cap        <- Args[14]
  apa_scores     <- as.logical(Args[15])
  flag_inherent  <- as.logical(Args[16])
  original_bedpe_path <- Args[17]

  n_samples <- length(hic_paths)

  if (!flag_norm %in% c("blank", "none", "cpm", "aqua")) {
    cat("Norm should be one of: none, cpm, aqua\n"); q(save = "no")
  }

  pairs <- read.table(pairs_path, as.is = TRUE)[, 1:6]

  # Read the original bedpe for the size distribution histogram
  original_pairs <- read.table(original_bedpe_path, as.is = TRUE)[, 1:6]

  # Read mergeStats once per sample and cache
  merge_stats_list <- lapply(stat_paths, read.table, as.is = TRUE)

  # Resolve normalization across all samples
  if (isTRUE(flag_inherent)) {
    flag_norm <- "none"
    flag_loop_norm <- "TRUE"
    norm_label <- "inherent"
  } else {
    all_have_spike <- all(sapply(merge_stats_list, function(x) ncol(x) == 2))

    if (flag_norm == "aqua" && !all_have_spike) {
      cat("\n--norm cannot be aqua: not all samples have spike-in data. Using cpm.\n\n")
      flag_norm <- "cpm"
    } else if (flag_norm == "blank") {
      flag_norm <- if (all_have_spike) "aqua" else "cpm"
    }
    norm_label <- flag_norm
  }

  if (isTRUE(flag_inherent) && max_cap != "no_cap") {
    cat("\nParameter --max_cap is not applicable when inherent = TRUE\nContinuing with inherent...\n\n")
  }

  max_label_width <- max(nchar(sample_labels))

  # Header
  type_header <- if (n_samples > 1) "Cohort" else "Single Sample"
  type_label <- if (n_samples > 1) "Cohort" else "Sample"

  cat(sprintf("\nAPA %s Analysis\n", type_header))
  cat(sprintf("============================================\n"))
  cat(sprintf("%-16s<- %s\n", type_label, paste(sample_labels, collapse = ", ")))
  cat(sprintf("Resolution      <- %d\n", res))
  cat(sprintf("Normalization   <- %s\n", norm_label))
  cat(sprintf("Loop norm       <- %s\n", flag_loop_norm))
  cat(sprintf("Output APA Plot <- %s\n\n", out_file))

  # Per-sample APA + normalization, then average
  matrix_list <- vector("list", n_samples)

  for (k in seq_len(n_samples)) {
    hic_path   <- hic_paths[k]
    sample_dir <- dirname(hic_path)
    cached_power_laws <- NULL
    if (flag_inherent == TRUE) {
      pre_check_power_law(sample_dir, res)
      cached_power_laws <- read_power_law_data(sample_dir, res)
    }

    nf <- compute_norm_factors(merge_stats_list[[k]], flag_norm)

    cat(sprintf("  [%d/%d] %-*s  (norm=%.6f, aqua=%.6f)\n",
                k, n_samples, max_label_width, sample_labels[k],
                nf$norm_factor, nf$aqua_factor))

    apa_result <- run_apa(hic_path, pairs, resolution = res, window = win_size,
                          flag_inherent = flag_inherent, power_laws = cached_power_laws)
    matrix_raw <- as.data.frame(apa_result$aggregate)
    loops_counted <- apa_result$counted

    if (loops_counted == 0) {
      cat(sprintf("Warning: No valid loops contributed contacts for sample %s.\n",
                  sample_labels[k]))
      matrix_norm <- matrix_raw * 0
    } else {
      matrix_norm <- matrix_raw * nf$norm_factor * nf$aqua_factor
      if (flag_loop_norm == "TRUE" || isTRUE(flag_inherent)) {
        matrix_norm <- matrix_norm / loops_counted
      }
    }
    matrix_list[[k]] <- matrix_norm
  }

  # Average across samples
  if (n_samples > 1) {
    cat(sprintf("\nAveraging across %d samples...", n_samples))
  }

  matrix_norm <- Reduce("+", matrix_list) / n_samples

  # Check for both all-zero and all-NaN
  if (all(is.na(matrix_norm) | matrix_norm == 0)) {
    cat("APA matrix is all zeros or NaN - no valid contacts detected.\n")
    q(save = "no", status = 1)
  }

  write.table(matrix_norm, file = file.path(out_dir, "matrix.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # Caps
  if (flag_inherent) {
    cap           <- 2
    min_cap_value <- 0
  } else {
    cap           <- if (max_cap == "no_cap") max(matrix_norm) else as.numeric(max_cap)
    min_cap_value <- if (min_cap == "no_cap") 0 else as.numeric(min_cap)
  }

  cat(sprintf("\nmax_cap         <- %.3f\n", cap))
  cat(sprintf("min_cap         <- %.3f\n", min_cap_value))

  # Plot
  breakList      <- seq(min_cap_value, cap, length.out = 6)
  color_ramp_red <- colorRampPalette(c("white", "red"))(100)
  cols           <- build_axis_labels(res, win_size)
  matrix_long    <- matrix_to_long(matrix_norm, cols)

  plot_label <- if (n_samples == 1) sample_labels[1] else paste0("cohort (n.samples=", n_samples, ")")

  box_data <- NULL
  if (apa_scores) {
    peak_scores <- compute_peak_scores(matrix_norm, res)
    box_data    <- build_box_annotations(matrix_norm, res, peak_scores)
    print_apa_scores(plot_label, peak_scores)
  }

  p1 <- build_apa_plot(matrix_long, color_ramp_red, breakList, cap,
                      min_cap_value, box_data, apa_scores, flag_inherent)

  wrap_label <- function(text, width = 75) {
    paste(strwrap(text, width = width), collapse = "\n")
  }

  if (n_samples > 1) {
    plot_title <- wrap_label(paste0("Cohort: ", paste(sample_labels, collapse = ", ")))
    plot_title <- paste0(plot_title, "\n(n.loop=", num_loops, ")")
  } else {
    plot_title <- paste0(sample_labels[1], " (n.loop=", num_loops, ")")
  }

  title_fontsize <- if (n_samples > 1) 8 else 10

  test_con <- tryCatch(file(out_file, open = "ab"), error = function(e) NULL)
  if (is.null(test_con)) {
    stop(sprintf(
      "Cannot write to '%s' -- the file may be open in another application. Please close it and try again.",
      out_file
    ))
  }
  close(test_con)

  pdf(out_file, width = 4.5, height = 4.5, useDingbats = FALSE)
  grid.draw(arrangeGrob(p1, nrow = 1,
    top = textGrob(plot_title, x = 0.01, hjust = 0, gp = gpar(fontsize = title_fontsize))))
  dev.off()

  plot_size_distribution(original_pairs, plot_label,
    file.path(dirname(out_file), "bedpe_size-distribution.pdf"))

  cat("\nOutput folder created at", out_dir, "\n\n")
}


# ==========================================================================
#                          TWO-SAMPLE ANALYSIS
# ==========================================================================
if (length(Args) == 20) {

  hic_paths_A     <- strsplit(Args[1], ",")[[1]]
  hic_paths_B     <- strsplit(Args[2], ",")[[1]]
  sample_labels_A <- strsplit(Args[3], ",")[[1]]
  sample_labels_B <- strsplit(Args[4], ",")[[1]]
  max_cap         <- Args[5]
  max_cap_delta   <- Args[6]
  out_file        <- Args[7]
  win_size        <- as.numeric(Args[8])
  res             <- as.numeric(Args[9])
  pairs_path      <- Args[10]
  flag_norm       <- Args[11]
  stat_paths_A    <- strsplit(Args[12], ",")[[1]]
  stat_paths_B    <- strsplit(Args[13], ",")[[1]]
  flag_loop_norm  <- Args[14]
  num_loops       <- as.numeric(Args[15])
  out_dir         <- Args[16]
  min_cap         <- Args[17]
  apa_scores      <- as.logical(Args[18])
  flag_inherent   <- as.logical(Args[19])
  original_bedpe_path <- Args[20]

  n_samples_A <- length(hic_paths_A)
  n_samples_B <- length(hic_paths_B)

  if (!flag_norm %in% c("blank", "none", "cpm", "aqua")) {
    cat("Norm should be one of: none, cpm, aqua\n"); q(save = "no")
  }

  pairs <- read.table(pairs_path, as.is = TRUE)[, 1:6]

  # Read the original bedpe for the size distribution histogram
  original_pairs <- read.table(original_bedpe_path, as.is = TRUE)[, 1:6]

  merge_stats_list_A <- lapply(stat_paths_A, read.table, as.is = TRUE)
  merge_stats_list_B <- lapply(stat_paths_B, read.table, as.is = TRUE)

  # Resolve normalization across all samples 
  if (isTRUE(flag_inherent)) {
    flag_norm <- "none"  # so aqua/cpm factors = 1
    flag_loop_norm <- "TRUE"
    norm_label <- "inherent"
  } else {
    all_merge_stats <- c(merge_stats_list_A, merge_stats_list_B)
    all_have_spike <- all(sapply(all_merge_stats, function(x) ncol(x) == 2))

    if (flag_norm == "aqua" && !all_have_spike) {
      cat("\n--norm cannot be aqua: not all samples have spike-in data. Using cpm.\n\n")
      flag_norm <- "cpm"
    } else if (flag_norm == "blank") {
      flag_norm <- if (all_have_spike) "aqua" else "cpm"
    }
    norm_label <- flag_norm
  }

  if (isTRUE(flag_inherent) && max_cap != "no_cap") {
    cat("\nParameter --max_cap is not applicable when inherent = TRUE\nContinuing with inherent...\n\n")
  }

  all_labels <- c(sample_labels_A, sample_labels_B)
  max_label_width <- max(nchar(all_labels))

  # Header
  typeA <- if (n_samples_A > 1) "Cohort A" else "Sample A"
  typeB <- if (n_samples_B > 1) "Cohort B" else "Sample B"
  type_header <- if (n_samples_A > 1 || n_samples_B > 1) "Cohort" else "Two Sample"

  cat(sprintf("\nAPA %s Analysis\n", type_header))
  cat(sprintf("============================================\n"))
  cat(sprintf("%-16s<- %s\n", typeA, paste(sample_labels_A, collapse = ", ")))
  cat(sprintf("%-16s<- %s\n", typeB, paste(sample_labels_B, collapse = ", ")))
  cat(sprintf("Resolution      <- %d\n", res))
  cat(sprintf("Normalization   <- %s\n", norm_label))
  cat(sprintf("Loop norm       <- %s\n", flag_loop_norm))
  cat(sprintf("Output          <- %s\n\n", out_file))

  matrix_list_A <- vector("list", n_samples_A)
  matrix_list_B <- vector("list", n_samples_B)

  for (k in seq_len(n_samples_A)) {
    hic_path   <- hic_paths_A[k]
    sample_dir <- dirname(hic_path)
    cached_power_laws <- NULL
    if(flag_inherent == TRUE){ 
      pre_check_power_law(sample_dir, res)
      cached_power_laws <- read_power_law_data(sample_dir, res)
    }

    nf_A <- compute_norm_factors(merge_stats_list_A[[k]], flag_norm)

    cat(sprintf("  [A %d/%d] %-*s  (norm=%.6f, aqua=%.6f)\n",
              k, n_samples_A, max_label_width, sample_labels_A[k],
              nf_A$norm_factor, nf_A$aqua_factor))
    apa_result_A <- run_apa(hic_path, pairs, resolution = res, window = win_size,
                            flag_inherent = flag_inherent, power_laws = cached_power_laws)
    matrix_raw_A <- as.data.frame(apa_result_A$aggregate)
    loops_counted <- apa_result_A$counted

    if (loops_counted == 0) {
      cat(sprintf("Warning: No valid loops contributed contacts for sample %s.\n",
                  sample_labels_A[k]))
      matrix_norm_A <- matrix_raw_A * 0
    } else {
      matrix_norm_A <- matrix_raw_A * nf_A$norm_factor * nf_A$aqua_factor
      if (flag_loop_norm == "TRUE" || isTRUE(flag_inherent)) {
        matrix_norm_A <- matrix_norm_A / loops_counted
      }
    }
    matrix_list_A[[k]] <- matrix_norm_A
  }

  for (k in seq_len(n_samples_B)) {
    hic_path   <- hic_paths_B[k]
    sample_dir <- dirname(hic_path)
    cached_power_laws <- NULL
    if(flag_inherent == TRUE){ 
      pre_check_power_law(sample_dir, res)
      cached_power_laws <- read_power_law_data(sample_dir, res)
    }

    nf_B <- compute_norm_factors(merge_stats_list_B[[k]], flag_norm)

    cat(sprintf("  [B %d/%d] %-*s  (norm=%.6f, aqua=%.6f)\n",
              k, n_samples_B, max_label_width, sample_labels_B[k],
              nf_B$norm_factor, nf_B$aqua_factor))

    apa_result_B <- run_apa(hic_path, pairs, resolution = res, window = win_size,
                            flag_inherent = flag_inherent, power_laws = cached_power_laws)
    matrix_raw_B <- as.data.frame(apa_result_B$aggregate)
    loops_counted <- apa_result_B$counted

    if (loops_counted == 0) {
      cat(sprintf("Warning: No valid loops contributed contacts for sample %s.\n",
                  sample_labels_B[k]))
      matrix_norm_B <- matrix_raw_B * 0
    } else {
      matrix_norm_B <- matrix_raw_B * nf_B$norm_factor * nf_B$aqua_factor
      if (flag_loop_norm == "TRUE" || isTRUE(flag_inherent)) {
        matrix_norm_B <- matrix_norm_B / loops_counted
      }
    }
    matrix_list_B[[k]] <- matrix_norm_B
  }

  if (n_samples_A > 1 || n_samples_B > 1) {
    cat(sprintf("\nAveraging: %d sample(s) (A), %d sample(s) (B)", n_samples_A, n_samples_B))
  }

  matrix_norm_A <- Reduce("+", matrix_list_A) / n_samples_A
  matrix_norm_B <- Reduce("+", matrix_list_B) / n_samples_B

  if (all(is.na(matrix_norm_A) | matrix_norm_A == 0)) {
    cat("APA matrix for sample/cohort A is all zeros or NaN - no valid contacts detected.\n")
    q(save = "no", status = 1)
  }
  if (all(is.na(matrix_norm_B) | matrix_norm_B == 0)) {
    cat("APA matrix for sample/cohort B is all zeros or NaN - no valid contacts detected.\n")
    q(save = "no", status = 1)
  }

  delta <- matrix_norm_B - matrix_norm_A

  write.table(delta, file = file.path(out_dir, "matrix.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)

  # Caps
  if (flag_inherent) {
    cap           <- 2
    min_cap_value <- 0
  } else {
    cap           <- if (max_cap == "no_cap") max(matrix_norm_A, matrix_norm_B) else as.numeric(max_cap)
    min_cap_value <- if (min_cap == "no_cap") 0 else as.numeric(min_cap)
  }
  delta_cap <- if (max_cap_delta == "no_cap") max(abs(delta)) else as.numeric(max_cap_delta)

  breakList       <- seq(min_cap_value, cap, length.out = 6)
  breakList_delta <- seq(-delta_cap, delta_cap, length.out = 7)
  cat(sprintf("\nmax_cap         <- %.3f\n", cap))
  cat(sprintf("min_cap         <- %.3f\n", min_cap_value))
  cat(sprintf("max_cap_delta   <- %.3f\n", delta_cap))

  color_ramp_red <- colorRampPalette(c("white", "red"))(100)
  color_ramp_vio <- colorRampPalette(c("dodgerblue", "white", "mediumvioletred"))(100)
  cols <- build_axis_labels(res, win_size)

  matrix_long_A <- matrix_to_long(matrix_norm_A, cols)
  matrix_long_B <- matrix_to_long(matrix_norm_B, cols)
  matrix_long_D <- matrix_to_long(delta, cols)

  # Labels
  labelA <- paste(sample_labels_A, collapse = ", ")
  labelB <- paste(sample_labels_B, collapse = ", ")

  box_data_A <- NULL; box_data_B <- NULL
  if (apa_scores) {
    peak_scores_A <- compute_peak_scores(matrix_norm_A, res)
    peak_scores_B <- compute_peak_scores(matrix_norm_B, res)
    box_data_A <- build_box_annotations(matrix_norm_A, res, peak_scores_A)
    box_data_B <- build_box_annotations(matrix_norm_B, res, peak_scores_B)
    print_apa_scores(labelA, peak_scores_A)
    print_apa_scores(labelB, peak_scores_B)
  }

  p1 <- build_apa_plot(matrix_long_A, color_ramp_red, breakList, cap,
                       min_cap_value, box_data_A, apa_scores, flag_inherent)
  p2 <- build_apa_plot(matrix_long_B, color_ramp_red, breakList, cap,
                       min_cap_value, box_data_B, apa_scores, flag_inherent)

  p3 <- ggplot(matrix_long_D, aes(x = col, y = row, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(
      colors = color_ramp_vio, breaks = breakList_delta,
      limits = c(-delta_cap, delta_cap), oob = scales::squish,
      labels = get_label_fmt(delta_cap), name = NULL
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(matrix_long_D$row)), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y       = element_blank(),
      axis.title        = element_blank(),
      axis.ticks        = element_blank(),
      panel.grid        = element_blank(),
      plot.title        = element_text(hjust = 0.5),
      legend.key.height = unit(1.2, "cm")
    )

  if (apa_scores) {
    n_d <- nrow(delta); mid <- (n_d + 1) %/% 2
    delta_peak_box <- data.frame(
      xmin = mid - 0.5, xmax = mid + 0.5,
      ymin = mid - 0.5, ymax = mid + 0.5,
      x_label = mid + 3, y_label = mid + 1.2,
      label = sprintf("peak: %.3f", delta[mid, mid])
    )
    p3 <- p3 +
      geom_rect(data = delta_peak_box,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE, color = "black", fill = NA, linewidth = 0.6) +
      geom_text(data = delta_peak_box,
        aes(x = x_label, y = y_label, label = label),
        inherit.aes = FALSE, color = "black", size = 3.0)
  }

 if (n_samples_A > 1) {
    plot_labelA <- paste0("Cohort A: ", paste(sample_labels_A, collapse = ", "))
  } else {
    plot_labelA <- sample_labels_A[1]
  }

  if (n_samples_B > 1) {
    plot_labelB <- paste0("Cohort B: ", paste(sample_labels_B, collapse = ", "))
  } else {
    plot_labelB <- sample_labels_B[1]
  }

  test_con <- tryCatch(file(out_file, open = "ab"), error = function(e) NULL)
  if (is.null(test_con)) {
    stop(sprintf(
      "Cannot write to '%s' -- the file may be open in another application. Please close it and try again.",
      out_file
    ))
  }
  close(test_con)

  if (n_samples_A > 1 || n_samples_B > 1) {
    plot_title <- paste0(plot_labelA, "\n", plot_labelB, "\nB - A delta (n.loop=", num_loops, ")")
    pdf(out_file, width = 14, height = 4.5, useDingbats = FALSE)
    grid.draw(arrangeGrob(p1, p2, p3, nrow = 1,
      top = textGrob(plot_title, gp = gpar(fontsize = 10))))
    dev.off()
  } else {
    p1 <- p1 + ggtitle(plot_labelA)
    p2 <- p2 + ggtitle(plot_labelB)
    p3 <- p3 + ggtitle(paste0("delta (n.loop=", num_loops, ")"))

    pdf(out_file, width = 14, height = 4.5, useDingbats = FALSE)
    grid.draw(arrangeGrob(p1, p2, p3, nrow = 1))
    dev.off()
  }

  plot_size_distribution(original_pairs, paste(labelA, labelB, sep = " / "),
    file.path(dirname(out_file), "bedpe_size-distribution.pdf"), width = 10, height = 10)

  cat("\nOutput folder created at", out_dir, "\n\n")
}

if (!length(Args) %in% c(17, 20)) {
  cat(sprintf("\nERROR: Unexpected number of arguments: %d (expected 17 for one-sample or 20 for two-sample)\n",
              length(Args)))
  cat("Arguments received:\n")
  for (i in seq_along(Args)) {
    cat(sprintf("  [%d] %s\n", i, Args[i]))
  }
  q(save = "no", status = 1)
}
