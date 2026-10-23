#################################################################
##        Self-contained driver for subprocessing chunks       ##
#################################################################

`%||%` <- function(x, y) if (is.null(x)) y else x

ca_trailing <- commandArgs(trailingOnly = TRUE)

if (length(ca_trailing) > 0 && identical(ca_trailing[1], "--driver")) {
  args <- ca_trailing[-1]
  if (!(length(args) == 17 || length(args) == 20)) {
    stop(sprintf("driver: expected 17 (one-sample) or 20 (two-sample) args; got %d\n%s",
                 length(args), paste(args, collapse = "\n")))
  }

  # Original positional args
  path_pairs <- args[1]
  bin_size   <- suppressWarnings(as.numeric(args[2]))

  # num_loops (7th for one-sample, 9th for two-sample)
  num_loops_all <- suppressWarnings(as.numeric(if (length(args) == 17) args[7] else args[9]))
  if (is.na(num_loops_all) || num_loops_all <= 0) {
    num_loops_all <- as.integer(system2(
      "bash", c("-lc", sprintf("wc -l < %s", shQuote(path_pairs))), stdout = TRUE
    ))
  }

  chunk_size <- if (bin_size == 1000) {
    min(num_loops_all, 2000L)
  } else if (num_loops_all <= 10000) {
    num_loops_all
  } else {
    10000L
  }

  # Resolve path to this script to call itself
  ca_full <- commandArgs(FALSE)
  m <- grep("^--file=", ca_full)
  self_script <- if (length(m)) {
    normalizePath(sub("^--file=", "", ca_full[m][1]))
  } else {
    # Fallback
    file.path(Sys.getenv("AQUA_TOOLS_DIR", file.path(Sys.getenv("HOME"), "aqua_tools")), "query_bedpe.r")
  }

  con <- file(description = path_pairs, open = "r")
  on.exit(close(con), add = TRUE)

  i <- 0L
  repeat {
    chunk <- tryCatch(read.table(con, nrows = chunk_size, as.is = TRUE),
                      error = function(e) data.frame())
    if (!nrow(chunk)) break

    i <- i + 1L
    chunk_path <- tempfile(pattern = sprintf("pairs_chunk_%03d_", i), fileext = ".tsv")
    write.table(chunk, file = chunk_path, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)

    # Build args for chunk mode:
    # - replace original path_pairs with chunk_path
    # - num_loops to nrow(chunk)
    args_tail <- args[-1]
    pos_num_loops_tail <- if (length(args) == 17) 6 else 8  # (7|9) - 1
    args_tail[pos_num_loops_tail] <- as.character(nrow(chunk))

    status <- system2("Rscript", c(self_script, chunk_path, args_tail),
                      stdout = "", stderr = "")
    if (status != 0) {
      stop(sprintf("R failed on chunk %d (temp: %s) with exit status %d",
                   i, chunk_path, status))
    }

    unlink(chunk_path)
  }

  quit(save = "no")
}

suppressPackageStartupMessages( library(strawr  ) )
suppressPackageStartupMessages( library(parallel) )

options(scipen = 999)

source(file.path(Sys.getenv("AQUA_TOOLS_DIR", file.path(Sys.getenv("HOME"), "aqua_tools")), "aqua_functions.r"))

#################################################################
##                          Functions                          ##
#################################################################

parse_coords <- function(txt) {
  parts <- unlist(strsplit(txt, ":"))
  if (length(parts) != 3) return(NULL)
  parts_numeric <- suppressWarnings(as.numeric(parts[-1]))
  if (any(is.na(parts_numeric))) return(NULL)
  list(chr = parts[1], start = parts_numeric[1], end = parts_numeric[2])
}

prefix_filter_bedpe <- function( bedpe ){

  chr_A <- bedpe[1]
  chr_B <- bedpe[4]

  start1_bin <- as.numeric(bedpe[7])
  end1_bin   <- as.numeric(bedpe[8])

  start2_bin <- as.numeric(bedpe[9 ])
  end2_bin   <- as.numeric(bedpe[10])

  midpoint_1_bin <- as.numeric(bedpe[11])
  midpoint_2_bin <- as.numeric(bedpe[12])

  if( end1_bin == start1_bin ){ midpoint_1_bin <- start1_bin }
  if( end2_bin == start2_bin ){ midpoint_2_bin <- start2_bin }

  if(flag_strip_chr == "yes"){

    chr_L <- sub( "chr", "", chr_A )
    chr_R <- sub( "chr", "", chr_B )

    txt_L <- sprintf( "%s:%d:%d", chr_L, start1_bin, end1_bin )
    txt_R <- sprintf( "%s:%d:%d", chr_R, start2_bin, end2_bin )

    obj <- list( txt_L = txt_L,
                 txt_R = txt_R,
                 midpoint_1_bin = midpoint_1_bin,
                 midpoint_2_bin = midpoint_2_bin )
  }

  if(flag_strip_chr == "no"){

    chr_L <- chr_A
    chr_R <- chr_B

    txt_L <- sprintf( "%s:%d:%d", chr_L, start1_bin, end1_bin )
    txt_R <- sprintf( "%s:%d:%d", chr_R, start2_bin, end2_bin )

    obj <- list( txt_L = txt_L,
                 txt_R = txt_R,
                 midpoint_1_bin = midpoint_1_bin,
                 midpoint_2_bin = midpoint_2_bin )
  }

  if(flag_strip_chr == "A"){

    sample_A_chr_L <- chr_A
    sample_A_chr_R <- chr_B

    sample_A_txt_L <- sprintf( "%s:%d:%d", sample_A_chr_L, start1_bin, end1_bin )
    sample_A_txt_R <- sprintf( "%s:%d:%d", sample_A_chr_R, start2_bin, end2_bin )

    sample_B_chr_L <- sub("chr", "", chr_A )
    sample_B_chr_R <- sub("chr", "", chr_B )

    sample_B_txt_L <- sprintf( "%s:%d:%d", sample_B_chr_L, start1_bin, end1_bin )
    sample_B_txt_R <- sprintf( "%s:%d:%d", sample_B_chr_R, start2_bin, end2_bin )

    obj <- list( sample_A_txt_L = sample_A_txt_L,
                 sample_A_txt_R = sample_A_txt_R,
                 sample_B_txt_L = sample_B_txt_L,
                 sample_B_txt_R = sample_B_txt_R,
                 midpoint_1_bin = midpoint_1_bin,
                 midpoint_2_bin = midpoint_2_bin )

  }

  if(flag_strip_chr == "B"){

    sample_B_chr_L <- chr_A
    sample_B_chr_R <- chr_B

    sample_B_txt_L <- sprintf( "%s:%d:%d", sample_B_chr_L, start1_bin, end1_bin )
    sample_B_txt_R <- sprintf( "%s:%d:%d", sample_B_chr_R, start2_bin, end2_bin )

    sample_A_chr_L <- sub("chr","", chr_A)
    sample_A_chr_R <- sub("chr","", chr_B)

    sample_A_txt_L <- sprintf( "%s:%d:%d", sample_A_chr_L, start1_bin, end1_bin )
    sample_A_txt_R <- sprintf( "%s:%d:%d", sample_A_chr_R, start2_bin, end2_bin )

    obj <- list( sample_A_txt_L = sample_A_txt_L,
                 sample_A_txt_R = sample_A_txt_R,
                 sample_B_txt_L = sample_B_txt_L,
                 sample_B_txt_R = sample_B_txt_R,
                 midpoint_1_bin = midpoint_1_bin,
                 midpoint_2_bin = midpoint_2_bin )

  }

  return( obj )

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

inherent_normalization <- function(mat, power_laws, bin_size) {

  if(flag_inherent == TRUE && nrow(mat) > 0){

    for (i in 1:nrow(mat)){
      bin_index <- (mat[i, "y"] - mat[i, "x"]) / bin_size + 1
      if (bin_index >= nrow(power_laws)){
        bin_index <- max(nrow(power_laws))
      }
      raw_val <- mat[i, "counts"]
      inh_val <- round((raw_val - power_laws[bin_index, 1]) / (power_laws[bin_index, 2] - power_laws[bin_index, 1]), 3)
      mat[i, "counts"] <- inh_val

    }
  }
  return(mat)
}

generate_sequences <- function(txt, bin_size) {
  parts <- unlist(strsplit(txt, ":"))

  # Check input splits into exactly three parts (chr, start, end)
  if (length(parts) != 3) {
    return(NA)
  }

  parts_numeric <- suppressWarnings(as.numeric(parts[-1]))

  if (any(is.na(parts_numeric))) {
    return(NA)
  }

  seq(from = parts_numeric[1], to = parts_numeric[2], by = bin_size)
}

initialize_matrix <- function(dim_1, dim_2) {
  sparse_mat <- matrix(0, ncol = length(dim_1), nrow = length(dim_2))
  colnames(sparse_mat) <- dim_1
  rownames(sparse_mat) <- dim_2
  return(sparse_mat)
}

# Power laws for a sample when inherent units are on, otherwise NULL
inherent_laws <- function(sample) {
  if (!isTRUE(flag_inherent)) return(NULL)
  if (isTRUE(two_sample_analysis)) get(paste0("cached_power_laws_", sample)) else cached_power_laws
}

# Fill the window with straw contacts (rows = foot 2 bins, cols = foot 1 bins)
populate_matrix <- function(mat, data, txt_L = NULL, txt_R = NULL, power_laws = NULL) {

  chr_L <- if (is.null(txt_L)) NA_character_ else sub(":.*", "", txt_L)
  chr_R <- if (is.null(txt_R)) NA_character_ else sub(":.*", "", txt_R)
  trans <- !is.na(chr_L) && !is.na(chr_R) && chr_L != chr_R
  attr(mat, "trans") <- trans

  # if inherent: init every box at the inherent value of 0
  if (!is.null(power_laws) && !trans) {
    d_idx <- round(abs(outer(as.numeric(rownames(mat)), as.numeric(colnames(mat)), "-")) / bin_size) + 1
    d_idx[d_idx > nrow(power_laws)] <- nrow(power_laws)
    off <- power_laws[[1]][d_idx]
    on  <- power_laws[[2]][d_idx]
    mat[] <- round((0 - off) / (on - off), 3)
  }

  # Track which boxes straw filled
  filled <- matrix(FALSE, nrow(mat), ncol(mat))
  attr(mat, "filled") <- filled

  if (nrow(data) == 0) {
    mat[is.na(mat)] <- 0
    return(mat)
  }

  row_keys   <- as.character(data[["x"]])
  col_keys   <- as.character(data[["y"]])
  cell_vals  <- data[["counts"]]

  row_names <- rownames(mat)
  col_names <- colnames(mat)

  if (trans) {
    # straw puts the chromosome that comes first in the .hic in x
    hic_order  <- sub("^chr", "", hic_A_chroms)
    idx_L      <- match(sub("^chr", "", chr_L), hic_order)
    idx_R      <- match(sub("^chr", "", chr_R), hic_order)
    x_on_foot1 <- idx_L < idx_R
    if (isTRUE(x_on_foot1)) {
      target_rows <- match(col_keys, row_names)
      target_cols <- match(row_keys, col_names)
    } else {
      target_rows <- match(row_keys, row_names)
      target_cols <- match(col_keys, col_names)
    }
    keep        <- !is.na(target_rows) & !is.na(target_cols)
    target_rows <- target_rows[keep]
    target_cols <- target_cols[keep]
    target_vals <- cell_vals[keep]
  } else {
    # Primary orientation (upstream foot first)
    row_idx_primary <- match(row_keys, row_names)
    col_idx_primary <- match(col_keys, col_names)
    keep_primary    <- !is.na(row_idx_primary) & !is.na(col_idx_primary)

    # Swapped orientation
    row_idx_swapped <- match(col_keys, row_names)
    col_idx_swapped <- match(row_keys, col_names)
    keep_swapped    <- !keep_primary & !is.na(row_idx_swapped) & !is.na(col_idx_swapped)

    target_rows <- c(row_idx_primary[keep_primary], row_idx_swapped[keep_swapped])
    target_cols <- c(col_idx_primary[keep_primary], col_idx_swapped[keep_swapped])
    target_vals <- c(cell_vals      [keep_primary], cell_vals       [keep_swapped])
  }

  if (length(target_rows)) {
    mat[cbind(target_rows, target_cols)] <- target_vals
    filled[cbind(target_rows, target_cols)] <- TRUE
  }

  mat[is.na(mat)] <- 0
  attr(mat, "filled") <- filled
  mat
}

# handles bins in both feet, keep only upper triangle and set
# mirror to NA, or copy the value into it if keep_mirrors. Skips trans
dedupe_mirror_cells <- function(mat, keep_mirrors = FALSE) {
  if (isTRUE(attr(mat, "trans"))) return(mat)
  shared <- intersect(rownames(mat), colnames(mat))
  if (length(shared) > 1) {
    ri  <- match(shared, rownames(mat))
    ci  <- match(shared, colnames(mat))
    blk <- mat[ri, ci, drop = FALSE]
    filled <- attr(mat, "filled")
    real <- if (is.null(filled)) blk != 0 else filled[ri, ci, drop = FALSE]
    sym <- ifelse(real, blk, t(blk))
    diag(sym) <- diag(blk)
    if (!keep_mirrors) sym[lower.tri(sym)] <- NA
    mat[ri, ci] <- sym
  }
  mat
}

process_formula <- function(sparse_mat, formula, midpoint1, midpoint2, bin_size, norm_factor, aqua_factor) {
  result <- NULL
  # center reads a single box, so it keeps both mirror copies; every other
  # formula counts each pair of bins once (upper triangle)
  sparse_mat <- dedupe_mirror_cells(sparse_mat, keep_mirrors = formula == "center")
  switch(
    formula,
    "center" = {
      if(as.character(midpoint2) %in% rownames(sparse_mat) &&
         as.character(midpoint1) %in% colnames(sparse_mat)) {
        result <- sparse_mat[as.character(midpoint2), as.character(midpoint1)]
      } else {
        result <- sparse_mat[ceiling(nrow(sparse_mat)/2), ceiling(ncol(sparse_mat)/2)]
      }
    },
    "max" = {
      result <- max(sparse_mat, na.rm = TRUE)
    },
    "sum" = {
      result <- sum(sparse_mat, na.rm = TRUE)
    },
    "average" = {
      result <- mean(sparse_mat, na.rm = TRUE)
    },
    "median" = {
      result <- median(sparse_mat, na.rm = TRUE)
    }
  )
  if (!is.null(result)) {
    result <- result * norm_factor * aqua_factor
    result <- round(result, 3)
  }
  return(result)
}

merge_tied_bins <- function(coords, connectivity = 4) {
  k <- nrow(coords)
  if (k == 1) return(list(1L))

  offs <- if (connectivity == 8) {
    expand.grid(dr = -1:1, dc = -1:1)
  } else {
    data.frame(dr = c(-1, 1, 0, 0), dc = c(0, 0, -1, 1))
  }
  offs <- offs[!(offs$dr == 0 & offs$dc == 0), , drop = FALSE]

  idx    <- setNames(seq_len(k), paste(coords[, 1], coords[, 2], sep = "-"))
  seen   <- logical(k)
  groups <- list()

  for (s in seq_len(k)) {
    if (seen[s]) next
    seen[s] <- TRUE
    stack   <- s
    comp    <- integer(0)
    while (length(stack)) {
      cur   <- stack[length(stack)]
      stack <- stack[-length(stack)]
      comp  <- c(comp, cur)
      nb  <- paste(coords[cur, 1] + offs$dr, coords[cur, 2] + offs$dc, sep = "-")
      hit <- idx[nb]
      hit <- hit[!is.na(hit)]
      hit <- hit[!seen[hit]]
      if (length(hit)) { seen[hit] <- TRUE; stack <- c(stack, hit) }
    }
    groups[[length(groups) + 1L]] <- unname(comp)
  }
  groups
}

calculate_new_feet <- function(sparse_mat, formula, midpoint_1, midpoint_2, bin_size, coord_L, coord_R) {
  switch(
    formula,
    "center" = {
      s1 <- as.numeric(midpoint_1)
      s2 <- as.numeric(midpoint_2)
    },
    "max" = {
      sparse_mat <- dedupe_mirror_cells(sparse_mat)
      max_val <- max(sparse_mat, na.rm = TRUE)
      if (max_val <= 0) {
        # if no meaningful max to narrow to,
        # so keep the original anchor span
        return(data.frame(
          start1 = coord_L$start, end1 = coord_L$end + bin_size,
          start2 = coord_R$start, end2 = coord_R$end + bin_size))
      }

      coords <- which(sparse_mat == max_val, arr.ind = TRUE)
      groups <- merge_tied_bins(coords, connectivity = 4)

      cols <- as.numeric(colnames(sparse_mat))
      rows <- as.numeric(rownames(sparse_mat))

      return(do.call(rbind, lapply(groups, function(g) {
        cc <- coords[g, , drop = FALSE]
        data.frame(
          start1 = min(cols[cc[, 2]]), end1 = max(cols[cc[, 2]]) + bin_size,
          start2 = min(rows[cc[, 1]]), end2 = max(rows[cc[, 1]]) + bin_size)
      })))
    },
    "sum" = , "average" = , "median" = {
      s1 <- coord_L$start
      s2 <- coord_R$start
    }
  )
  data.frame(start1 = s1, end1 = s1 + bin_size,
             start2 = s2, end2 = s2 + bin_size)
}

obtain_one_sample_counts <- function( mat_1, txt_L, txt_R, midpoint_1, midpoint_2, coord_L = NULL, coord_R = NULL ){

  # Parse coordinates once if not provided
  if (is.null(coord_L)) coord_L <- parse_coords(txt_L)
  if (is.null(coord_R)) coord_R <- parse_coords(txt_R)

  # Run the formula for inherent even when the window is empty
  n_1 <- if (isTRUE(flag_inherent)) max(nrow(mat_1), 1) else nrow(mat_1)

  if( isTRUE(flag_fix)  && shrink_wrap == "FALSE" && split == "FALSE" ){
    if( n_1 >= 1 ){    # build the window and run the formula
      dim_1 <- generate_sequences(txt_L, bin_size)
      dim_2 <- generate_sequences(txt_R, bin_size)

      sparse_mat_1 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_1 <- populate_matrix(sparse_mat_1, mat_1, txt_L, txt_R, inherent_laws("A"))
      count_1 <- process_formula(sparse_mat_1, formula, midpoint_1, midpoint_2, bin_size, norm_factor1, aqua_factor1)

    } else if ( n_1 == 0 ){ # shortcut: report 0
      count_1 <- 0
    } else {
      count_1 <- "*"
    }
    return(count_1)
  }

  if( isFALSE(flag_fix) && shrink_wrap == "FALSE" && split == "FALSE" ){
    if( n_1 >= 1 ){
      dim_1 <- generate_sequences(txt_L, bin_size)
      dim_2 <- generate_sequences(txt_R, bin_size)
      sparse_mat_1 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_1 <- populate_matrix(sparse_mat_1, mat_1, txt_L, txt_R, inherent_laws("A"))
      count_1 <- process_formula(sparse_mat_1, formula, midpoint_1, midpoint_2, bin_size, norm_factor1, aqua_factor1)
      new_feet   <- calculate_new_feet(sparse_mat_1, formula, midpoint_1, midpoint_2, bin_size, coord_L, coord_R)
      return_obj <- list(feet = new_feet, count1 = count_1)
    } else if ( n_1 == 0 ){
      count_1 <- 0
      return_obj <- list(
        feet = data.frame(start1 = coord_L$start, end1 = coord_L$end + bin_size,
                          start2 = coord_R$start, end2 = coord_R$end + bin_size),
        count1 = count_1)
    } else {
      count_1 <- "*"
      return_obj <- list(
        feet = data.frame(start1 = coord_L$start, end1 = coord_L$end + bin_size,
                          start2 = coord_R$start, end2 = coord_R$end + bin_size),
        count1 = count_1)
    }
    return(return_obj)
  }
}

call_one_sample_straw_bedpe <- function( bedpe_list, flag_inherent, flag_fix, formula, norm_factor1, aqua_factor1, shrink_wrap, split ){

  txt_L <- as.character( bedpe_list[1]$txt_L )
  txt_R <- as.character( bedpe_list[2]$txt_R )
  midpoint_1 <- as.character( bedpe_list[3]$midpoint_1_bin )
  midpoint_2 <- as.character( bedpe_list[4]$midpoint_2_bin )

  mat_1 <- straw(norm, hic_A, txt_L, txt_R, "BP", bin_size)

  if(flag_inherent) {
    mat_1 <- inherent_normalization(mat_1, cached_power_laws, bin_size)
  }

  coord_L <- parse_coords(txt_L)
  coord_R <- parse_coords(txt_R)

  res <- obtain_one_sample_counts(mat_1, txt_L, txt_R, midpoint_1, midpoint_2, coord_L, coord_R)
  return(res)
}

obtain_two_sample_counts <- function( mat_1, mat_2, txt_L, txt_R, midpoint_1, midpoint_2, coord_L = NULL, coord_R = NULL ){

  # Parse coordinates once if not provided
  if (is.null(coord_L)) coord_L <- parse_coords(txt_L)
  if (is.null(coord_R)) coord_R <- parse_coords(txt_R)

  # Run the formula for inherent even when the window is empty
  n_1 <- if (isTRUE(flag_inherent)) max(nrow(mat_1), 1) else nrow(mat_1)
  n_2 <- if (isTRUE(flag_inherent)) max(nrow(mat_2), 1) else nrow(mat_2)

  if( isTRUE(flag_fix) ){

    if ( n_1 >= 1 && n_2 == 0 ){
      dim_1 <- generate_sequences(txt_L, bin_size)
      dim_2 <- generate_sequences(txt_R, bin_size)
      sparse_mat_1 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_1 <- populate_matrix(sparse_mat_1, mat_1, txt_L, txt_R, inherent_laws("A"))
      count_1 <- process_formula(sparse_mat_1, formula, midpoint_1, midpoint_2, bin_size, norm_factor1, aqua_factor1)
      count_2 <- 0
    }
    else if ( n_1 == 0 && n_2 >= 1 ){
      count_1 <- 0
      dim_1 <- generate_sequences(txt_L, bin_size)
      dim_2 <- generate_sequences(txt_R, bin_size)
      sparse_mat_2 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_2 <- populate_matrix(sparse_mat_2, mat_2, txt_L, txt_R, inherent_laws("B"))
      count_2 <- process_formula(sparse_mat_2, formula, midpoint_1, midpoint_2, bin_size, norm_factor2, aqua_factor2)
    }
    else if ( n_1 >= 1  && n_2 >=  1 ){
      dim_1 <- generate_sequences(txt_L, bin_size)
      dim_2 <- generate_sequences(txt_R, bin_size)
      sparse_mat_1 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_1 <- populate_matrix(sparse_mat_1, mat_1, txt_L, txt_R, inherent_laws("A"))
      count_1 <- process_formula(sparse_mat_1, formula, midpoint_1, midpoint_2, bin_size, norm_factor1, aqua_factor1)
      sparse_mat_2 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_2 <- populate_matrix(sparse_mat_2, mat_2, txt_L, txt_R, inherent_laws("B"))
      count_2 <- process_formula(sparse_mat_2, formula, midpoint_1, midpoint_2, bin_size, norm_factor2, aqua_factor2)
    }
    else if ( n_1 == 0 && n_2 == 0 ){
      count_1 <- 0
      count_2 <- 0
    } else {
      count_1 <- "*"
      count_2 <- "*"
    }
    obj <- list(count_1 = count_1, count_2 = count_2 )
  }

  if( isFALSE(flag_fix) ){

    anchor_feet <- data.frame(
      start1 = coord_L$start, end1 = coord_L$end + bin_size,
      start2 = coord_R$start, end2 = coord_R$end + bin_size)

    if ( n_1 >= 1 && n_2 == 0 ){
      dim_1 <- generate_sequences(txt_L, bin_size)
      dim_2 <- generate_sequences(txt_R, bin_size)
      sparse_mat_1 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_1 <- populate_matrix(sparse_mat_1, mat_1, txt_L, txt_R, inherent_laws("A"))
      count_1 <- process_formula(sparse_mat_1, formula, midpoint_1, midpoint_2, bin_size, norm_factor1, aqua_factor1)
      feet    <- calculate_new_feet(sparse_mat_1, formula, midpoint_1, midpoint_2, bin_size, coord_L, coord_R)
      count_2 <- 0
    }
    else if ( n_1 == 0 && n_2 >= 1 ){
      dim_1 <- generate_sequences(txt_L, bin_size)
      dim_2 <- generate_sequences(txt_R, bin_size)
      sparse_mat_2 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_2 <- populate_matrix(sparse_mat_2, mat_2, txt_L, txt_R, inherent_laws("B"))
      count_2 <- process_formula(sparse_mat_2, formula, midpoint_1, midpoint_2, bin_size, norm_factor2, aqua_factor2)
      feet    <- calculate_new_feet(sparse_mat_2, formula, midpoint_1, midpoint_2, bin_size, coord_L, coord_R)
      count_1 <- 0
    }
    else if ( n_1 >= 1  && n_2 >=  1 ){
      dim_1 <- generate_sequences(txt_L, bin_size)
      dim_2 <- generate_sequences(txt_R, bin_size)
      sparse_mat_1 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_1 <- populate_matrix(sparse_mat_1, mat_1, txt_L, txt_R, inherent_laws("A"))
      count_1 <- process_formula(sparse_mat_1, formula, midpoint_1, midpoint_2, bin_size, norm_factor1, aqua_factor1)
      sparse_mat_2 <- initialize_matrix(dim_1, dim_2)
      sparse_mat_2 <- populate_matrix(sparse_mat_2, mat_2, txt_L, txt_R, inherent_laws("B"))
      count_2 <- process_formula(sparse_mat_2, formula, midpoint_1, midpoint_2, bin_size, norm_factor2, aqua_factor2)

      # narrow to whichever sample scored higher
      feet <- if (count_2 > count_1) {
        calculate_new_feet(sparse_mat_2, formula, midpoint_1, midpoint_2, bin_size, coord_L, coord_R)
      } else {
        calculate_new_feet(sparse_mat_1, formula, midpoint_1, midpoint_2, bin_size, coord_L, coord_R)
      }
    }
    else if ( n_1 == 0 && n_2 == 0 ){
      count_1 <- 0
      count_2 <- 0
      feet    <- anchor_feet
    }
    else {
      count_1 <- "*"
      count_2 <- "*"
      feet    <- anchor_feet
    }

    obj <- list(feet = feet, count1 = count_1, count2 = count_2 )
  }
  return(obj)
}

call_two_sample_straw_bedpe <- function( bedpe_list, flag_inherent, flag_fix, formula, norm_factor1, aqua_factor1, norm_factor2, aqua_factor2 ){

  if( length( unlist( bedpe_list ) ) == 4 ){
    txt_L <- as.character( bedpe_list[1]$txt_L )
    txt_R <- as.character( bedpe_list[2]$txt_R )
    midpoint_1 <- as.character( bedpe_list[3]$midpoint_1_bin )
    midpoint_2 <- as.character( bedpe_list[4]$midpoint_2_bin )
    mat_1 <- straw( norm, hic_A, txt_L, txt_R, "BP", bin_size )
    mat_2 <- straw( norm, hic_B, txt_L, txt_R, "BP", bin_size )
    if(flag_inherent) {
      mat_1 <- inherent_normalization(mat_1, cached_power_laws_A, bin_size)
      mat_2 <- inherent_normalization(mat_2, cached_power_laws_B, bin_size)
    }
    # Parse coordinates once
    coord_L <- parse_coords(txt_L)
    coord_R <- parse_coords(txt_R)
    ret <- obtain_two_sample_counts( mat_1, mat_2, txt_L, txt_R, midpoint_1, midpoint_2, coord_L, coord_R )
  }

  if( length( unlist( bedpe_list ) ) == 6 ){
    sample_A_txt_L <- as.character( bedpe_list[1]$sample_A_txt_L )
    sample_A_txt_R <- as.character( bedpe_list[2]$sample_A_txt_R )
    sample_B_txt_L <- as.character( bedpe_list[3]$sample_B_txt_L )
    sample_B_txt_R <- as.character( bedpe_list[4]$sample_B_txt_R )
    midpoint_1 <- as.character( bedpe_list[5]$midpoint_1_bin )
    midpoint_2 <- as.character( bedpe_list[6]$midpoint_2_bin )
    mat_1 <- straw( norm, hic_A, sample_A_txt_L, sample_A_txt_R, "BP", bin_size )
    mat_2 <- straw( norm, hic_B, sample_B_txt_L, sample_B_txt_R, "BP", bin_size )
    txt_L <- sample_A_txt_L
    txt_R <- sample_A_txt_R
    if(flag_inherent) {
      mat_1 <- inherent_normalization(mat_1, cached_power_laws_A, bin_size)
      mat_2 <- inherent_normalization(mat_2, cached_power_laws_B, bin_size)
    }
    # Parse coordinates once
    coord_L <- parse_coords(txt_L)
    coord_R <- parse_coords(txt_R)
    ret <- obtain_two_sample_counts( mat_1, mat_2, txt_L, txt_R, midpoint_1, midpoint_2, coord_L, coord_R )
  }
  return(ret)
}

get_c <- function(row, chr_col_name, tss_col_name, enh_col_name, bsize, hic_path){

  norm  <- "NONE"
  unit  <- "BP"

  chr       <- as.character(row[chr_col_name])
  tss_start <-   as.numeric(row[tss_col_name])
  enh_start <-   as.numeric(row[enh_col_name])

  # 1. Obtain contacts w.r.t. TSS ± 5 Mb
  tss_start  <- floor(tss_start/bsize)*bsize

  tss_window_upstream   <- tss_start - 5000000
  tss_window_downstream <- tss_start + 5000000

  chr1loc <- paste(chr,tss_window_upstream,tss_window_downstream,sep=":")
  chr2loc <- paste(chr,tss_window_upstream,tss_window_downstream,sep=":")

  matrix <- straw(
    norm,
    hic_path,
    chr1loc,
    chr2loc,
    unit,
    bsize )

  tss_interactions <- unique(
    rbind(
      matrix[matrix[,"x"]==tss_start,],
      matrix[matrix[,"y"]==tss_start,]
    )
  )

  enh_start  <- floor(enh_start/bsize)*bsize

  if( enh_start %in% unique(c(tss_interactions[,"x"],tss_interactions[,"y"])) ){
    TSS_vec      <- as.numeric(tss_interactions[,"counts"])
    TSS_vec_norm <- TSS_vec/sum(TSS_vec)
    TSS_vec_norm <- TSS_vec_norm/max(TSS_vec_norm)
    tss_interactions[,"counts_norm"] <- TSS_vec_norm
    if(tss_start < enh_start){
      id <- which(tss_interactions[,"y"] == enh_start)
    } else if(tss_start > enh_start) {
      id <- which(tss_interactions[,"x"] == enh_start)
    } else if(tss_start == enh_start) {
      id <- which(tss_interactions[,"x"] == enh_start & tss_interactions[,"y"] == enh_start)
    }
    if (length(id) == 0) return(0) else return(tss_interactions[id, "counts_norm"])
  } else {
    return(0)
  }
}

# auxiliary functions for --split :-

get_neighbours <- function( cell_row, cell_col, radius ){
  neighbours <- matrix( data = NA, nrow = radius*8, ncol = 2 )
  vec_1 <- paste( seq( -radius, radius, 1 ), -radius, sep = "," )
  vec_2 <- paste( seq( -radius, radius, 1 ),  radius, sep = "," )
  vec_3 <- paste( -radius, seq( -radius, radius, 1 ), sep = "," )
  vec_4 <- paste(  radius, seq( -radius, radius, 1 ), sep = "," )
  vec <- unique( c(vec_1, vec_2, vec_3, vec_4) )
  for( i in 1:length(vec) ){
    a <- cell_row + as.numeric(unlist(strsplit(vec[i],","))[1])
    b <- cell_col + as.numeric(unlist(strsplit(vec[i],","))[2])
    neighbours[i,1] <- a
    neighbours[i,2] <- b
  }
  return(neighbours)
}

get_clusters <- function( sparse_mat_1, split, padding, chr, bin_size, formula, norm_factor1, aqua_factor1 ){

  radius <- padding

  split_indices     <- which( sparse_mat_1 >= split, arr.ind = T )
  split_indices     <- split_indices[order(split_indices[,2], split_indices[,1]),]
  split_indices_ids <- paste( split_indices[,1], split_indices[,2], sep = "-")

  boolean_matrix <- matrix( data = NA, nrow = nrow(sparse_mat_1), ncol = ncol(sparse_mat_1) )
  boolean_matrix[split_indices] <- 1

  rownames(boolean_matrix) <- 1:nrow(boolean_matrix)
  colnames(boolean_matrix) <- 1:ncol(boolean_matrix)

  # 1. Get the neighboring indices of all cells greater than split threshold

  cell_clusters <- list()

  for(i in 1:nrow(split_indices)){
    cell_row <- split_indices[i,"row"]
    cell_col <- split_indices[i,"col"]
    cell_clusters[[i]] <- paste( cell_row, cell_col, sep = "-")
    for( j in 1:radius){
      neighbours <- get_neighbours( cell_row, cell_col, j )
      for( k in 1:nrow(neighbours)){
        neighbour_row <- neighbours[ k , 1 ]
        neighbour_col <- neighbours[ k , 2 ]
        if( paste( neighbour_row, neighbour_col, sep = "-") %in% split_indices_ids ){
          cell_clusters[[i]] <- c( cell_clusters[[i]], paste( neighbour_row, neighbour_col, sep = "-") )
        }
      }
    }
  }

  # 2. Construct shared neighbor matrix

  names(cell_clusters) <- split_indices_ids

  shared_neighbours <- matrix( data = NA, nrow = length(split_indices_ids), ncol = length(split_indices_ids) )
  colnames(shared_neighbours) <- split_indices_ids
  rownames(shared_neighbours) <- split_indices_ids

  for( i in 1:length(cell_clusters) ){
    for( j in 1:length(cell_clusters[[i]])){
      row <- names(cell_clusters[i])
      col <- cell_clusters[[i]][j]
      shared_neighbours[row,col] <- 1
    }
  }

  # 3. Obtain cluster indices

  raw_clusters <- list()
  for( i in 1:ncol(shared_neighbours) ){
    vec <- !(is.na(shared_neighbours[,i]))
    vec <- which(vec == TRUE)
    raw_clusters[[i]] <- as.numeric(vec)
  }

  for( i in 1:length(raw_clusters) ){
    members <- c()
    vec_1 <- raw_clusters[[i]]
    for( j in (i+1):length(raw_clusters) ){
      if(j > length(split_indices_ids)){break}
      vec_2 <- raw_clusters[[j]]
      if( any( vec_2 %in% vec_1 )){
        vec_1 <- unique( c(vec_1, vec_2) )
        members <- c( members, j)
      }
    }
    raw_clusters[[i]] <- unique( c(raw_clusters[[i]] , as.numeric(unlist(raw_clusters[members] ))) )
    raw_clusters[members] <- NA
  }

  final_clusters <- Filter(Negate(anyNA),raw_clusters)

  for( i in 1:length(final_clusters) ){
    cluster <- split_indices_ids[ final_clusters[[i]] ]
    cluster <- as.numeric(unlist(strsplit(cluster,"-")))
    row_min <- min(cluster[seq(1,length(cluster),2)])
    row_max <- max(cluster[seq(1,length(cluster),2)])
    col_min <- min(cluster[seq(2,length(cluster),2)])
    col_max <- max(cluster[seq(2,length(cluster),2)])
    mat <- as.matrix(sparse_mat_1[ row_min:row_max , col_min:col_max ])
    if( formula == "center"){
      count_1    <- mat[ ceiling(nrow(mat)/2) , ceiling(ncol(mat)/2) ]
      count_1    <- count_1 * norm_factor1 * aqua_factor1
      count_1    <- round( count_1 , 3 )
    } else if( formula == "max"){
      count_1    <- max( mat, na.rm = T )
      count_1    <- count_1 * norm_factor1 * aqua_factor1
      count_1    <- round( count_1 , 3 )
    } else if( formula == "sum"){
      count_1    <- sum( mat, na.rm = T )
      count_1    <- count_1 * norm_factor1 * aqua_factor1
      count_1    <- round( count_1 , 3 )
    } else if( formula == "average"){
      count_1    <- mean( mat, na.rm = T )
      count_1    <- count_1 * norm_factor1 * aqua_factor1
      count_1    <- round( count_1 , 3 )
    }
    final_clusters[[i]]         <- list()
    final_clusters[[i]]$chr     <- chr
    final_clusters[[i]]$feet1   <- c( as.numeric(rownames(sparse_mat_1)[row_min]) - (radius*bin_size) , as.numeric(rownames(sparse_mat_1)[row_max]) + (radius*bin_size) + bin_size )
    final_clusters[[i]]$feet2   <- c( as.numeric(colnames(sparse_mat_1)[col_min]) - (radius*bin_size) , as.numeric(colnames(sparse_mat_1)[col_max]) + (radius*bin_size) + bin_size )
    final_clusters[[i]]$count_1 <- count_1
    rm(mat)
  }
  return(final_clusters)
}

#################################################################
##                          Arguments                          ##
#################################################################

args        <- commandArgs( trailingOnly = TRUE )
path_pairs  <- args[1]
bin_size    <- as.numeric(args[2])
norm        <- "NONE"
flag_debug  <- FALSE

if( length(args) == 17 ){

  one_sample_analysis <-  TRUE
  two_sample_analysis <- FALSE

  hic_A              <-  args[3]
  path_mergeStats_A  <-  args[4]
  flag_norm          <-  args[5]
  genome             <-  args[6]
  num_loops          <-  as.numeric(args[7])
  formula            <-  args[8]
  flag_fix           <-  as.logical(args[9])
  split              <-  args[10]
  shrink_wrap        <-  args[11]
  padding            <-  as.numeric(args[12])
  expand             <-  as.numeric(args[13])
  flag_inherent      <-  as.logical(args[14])
  sample_dir         <-  args[15]
  flag_meta_string   <-  args[16]
  num_cores          <-  args[17]

  available_cores <- max(detectCores() - 1, 1)
  if (num_cores == "blank") {
    num_cores <- available_cores
  } else {
    if (as.numeric(num_cores) > available_cores){
      num_cores <- available_cores
    } else {
      num_cores <- as.integer(num_cores)
    }
  }

  hic_A_chroms      <-  readHicChroms(hic_A)$name

  if(formula == "mean"){ formula <- "average" }

  if(flag_inherent == TRUE){
    flag_norm <- "none"
    pre_check_power_law(sample_dir, bin_size)
    cached_power_laws <- read_power_law_data(sample_dir, bin_size)
  }

  flag_meta <- flag_meta_string == "TRUE"

  if( split != "FALSE"){
    cat("--split currently deprecated, contact Axiotl if needed. \n")
    q( save = "no" )
  }

  if( shrink_wrap != "FALSE"){
    cat("--shrink_wrap currently deprecated, contact Axiotl if needed. \n")
    q( save = "no" )
  }
}

if( length(args) == 20 ){

  one_sample_analysis <- FALSE
  two_sample_analysis <-  TRUE

  hic_A             <-  args[3]
  path_mergeStats_A <-  args[4]
  hic_B             <-  args[5]
  path_mergeStats_B <-  args[6]
  flag_norm         <-  args[7]
  genome            <-  args[8]
  num_loops         <-  as.numeric(args[9])
  formula           <-  args[10]
  flag_fix          <-  as.logical(args[11])
  split             <-  args[12]
  shrink_wrap       <-  args[13]
  padding           <-  as.numeric(args[14])
  expand            <-  as.numeric(args[15])
  sample_dir_A      <-  args[16]
  sample_dir_B      <-  args[17]
  flag_meta_string  <-  args[18]
  flag_inherent     <-  as.logical(args[19])
  num_cores         <-  args[20]

  available_cores <- max(detectCores() - 1, 1)
  if (num_cores == "blank") {
    num_cores <- available_cores
  } else {
    if (as.numeric(num_cores) > available_cores){
      num_cores <- available_cores
    } else {
      num_cores <- as.integer(num_cores)
    }
  }

  hic_A_chroms      <- readHicChroms(hic_A)$name
  hic_B_chroms      <- readHicChroms(hic_B)$name

  if(formula == "mean"){ formula <- "average" }

  if(flag_inherent == TRUE){
    flag_norm <- "none"
    pre_check_power_law(sample_dir_A, bin_size)
    pre_check_power_law(sample_dir_B, bin_size)
    cached_power_laws_A <- read_power_law_data(sample_dir_A, bin_size)
    cached_power_laws_B <- read_power_law_data(sample_dir_B, bin_size)
  }

  if(shrink_wrap != "FALSE"){ cat("--shrink_wrap currently deprecated, contact Axiotl if needed. \n") ; q(save="no") }
  if(split       != "FALSE"){ cat("--split currently deprecated, contact Axiotl if needed. \n")       ; q(save="no") }

  flag_meta <- flag_meta_string == "TRUE"
}

#################################################################
##                         Code Block                          ##
#################################################################


if( ! flag_norm %in% c( "blank", "none", "cpm", "aqua", "abc" ) ){
  cat("Norm should strictly be none, cpm, aqua, or abc in lower case \n")
  q( save = "no" )
}

if( num_loops <= 10000){
  chunk_size <- num_loops
} else if( num_loops > 10000){
  chunk_size <- 10000
}
# smaller chunks at 1kb
if (bin_size == 1000) chunk_size <- min(num_loops, 2000L)

num_chunks <- ceiling( num_loops / chunk_size )

if( shrink_wrap != "FALSE" || split != "FALSE" ){
  flag_fix <- FALSE
}

if( expand < 0 ){
  cat("Invalid use of --expand \n") ; q(save="no")
}

if( shrink_wrap != "FALSE" && split != "FALSE" ){
  cat("Please use --split or --shrink_wrap one at a time \n") ; q(save="no")
}

if (one_sample_analysis) {
  mergeStats_A <- read.table( path_mergeStats_A, as.is = TRUE)
  spikeVar     <- ncol(mergeStats_A)
  hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
  mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
  total1       <- sum( hg_total1 , mm_total1 )
  norm_factor1 <- 1000000 / total1
  aqua_factor1 <- hg_total1 / mm_total1

  if( sum( grepl("chr", hic_A_chroms ) ) > 0 ){
    flag_strip_chr <- "no"
  }
  if( sum( grepl("chr", hic_A_chroms ) ) == 0 ){
    flag_strip_chr <- "yes"
  }

  if(spikeVar == 1){
    if(flag_norm == "blank"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
    } else if(flag_norm == "none"){
      norm_factor1 <- 1 ; aqua_factor1 <- 1
    } else if(flag_norm == "cpm"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
    } else if(flag_norm == "aqua"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
    }
  }else if(spikeVar == 2){
    if(flag_norm == "blank"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- aqua_factor1
    } else if(flag_norm == "none"){
      norm_factor1 <- 1 ; aqua_factor1 <- 1
    } else if(flag_norm == "cpm"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
    }else if(flag_norm == "aqua"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- aqua_factor1
    }
  }
}

if (two_sample_analysis) {
  mergeStats_A <- read.table( path_mergeStats_A, as.is = TRUE)
  spikeVar_A   <- ncol(mergeStats_A)
  hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
  mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
  total1       <- sum( hg_total1 , mm_total1 )
  norm_factor1 <- 1000000 / total1
  aqua_factor1 <- hg_total1 / mm_total1

  mergeStats_B <- read.table( path_mergeStats_B, as.is = TRUE)
  spikeVar_B   <- ncol(mergeStats_B)
  hg_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 1 ])
  mm_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 2 ])
  total2       <- sum( hg_total2 , mm_total2 )
  norm_factor2 <- 1000000 / total2
  aqua_factor2 <- hg_total2 / mm_total2

  if( sum( grepl("chr", hic_A_chroms ) ) > 0 && sum( grepl("chr", hic_B_chroms ) ) > 0 ){
    flag_strip_chr <- "no"
  }
  if( sum( grepl("chr", hic_A_chroms ) ) > 0 && sum( grepl("chr", hic_B_chroms ) ) == 0 ){
    flag_strip_chr <- "A"
  }
  if( sum( grepl("chr", hic_B_chroms ) ) > 0 && sum( grepl("chr", hic_A_chroms ) ) == 0 ){
    flag_strip_chr <- "B"
  }
  if( sum( grepl("chr", hic_A_chroms ) ) == 0 && sum( grepl("chr", hic_B_chroms ) ) == 0 ){
    flag_strip_chr <- "yes"
  }

  if( spikeVar_A == 1 || spikeVar_B == 1){
    if (flag_norm == "blank"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
    } else if (flag_norm == "none"){
      norm_factor1 <- 1 ; aqua_factor1 <- 1
      norm_factor2 <- 1 ; aqua_factor2 <- 1
    } else if (flag_norm == "cpm"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
    } else if (flag_norm == "aqua"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
    }
  } else if (spikeVar_A == 2 || spikeVar_B == 2){
    if (flag_norm == "blank"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- aqua_factor1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- aqua_factor2
    } else if (flag_norm == "none"){
      norm_factor1 <- 1 ; aqua_factor1 <- 1
      norm_factor2 <- 1 ; aqua_factor2 <- 1
    } else if (flag_norm == "cpm"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
    } else if (flag_norm == "aqua"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- aqua_factor1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- aqua_factor2
    }
  }
}

#################################################################
##                         Chunk loop                          ##
#################################################################

con <- file( description = path_pairs, open = "r" )
pairs <- read.table(con, nrows=chunk_size, as.is = TRUE)

index <- 0
repeat {

  index <- index + 1

  if (ncol(pairs) < 6) {
    cat("Please use a minimum of 6 column bedpe, exiting....\n")
    q(save = "no")
  }

  colnames(pairs) <- c("chr1","start1","end1","chr2","start2","end2")

  # Apply distance restriction if flag_norm is "abc"
  if (flag_norm == "abc") {
    pairs <- pairs[abs(pairs$start1 - pairs$start2) <= 5000000, ]
  }

  if (ncol(pairs) > 6) {
    # Store columns beyond 6 in meta_cols
    meta_cols <- pairs[, 7:ncol(pairs), drop = FALSE]
    meta_cols <- data.frame(meta_cols, stringsAsFactors = FALSE)
    num_meta_cols <- ncol(meta_cols)
    pairs <- pairs[, 1:6]
  }

  # Ensure chr prefix
  if (sum(grepl("chr", pairs$chr1)) == 0) pairs$chr1 <- paste0("chr", pairs$chr1)
  if (sum(grepl("chr", pairs$chr2)) == 0) pairs$chr2 <- paste0("chr", pairs$chr2)

  # Keep only standard chromosomes
  # meta_cols must be subset with the same index
  keep_std_chr <- pairs$chr1 %in% c(paste0(rep("chr",22), 1:22), "chrX", "chrY")
  pairs <- pairs[keep_std_chr, ]
  if (exists("meta_cols")) {
    meta_cols <- meta_cols[keep_std_chr, , drop = FALSE]
  }

  # Bin coordinates
  b1 <- interval_to_bins(pairs$start1, pairs$end1, bin_size)
  b2 <- interval_to_bins(pairs$start2, pairs$end2, bin_size)

  pairs$start1_bin <- as.integer(b1$first * bin_size)
  pairs$end1_bin   <- as.integer(b1$last  * bin_size)
  pairs$start2_bin <- as.integer(b2$first * bin_size)
  pairs$end2_bin   <- as.integer(b2$last  * bin_size)

  if( expand > 0 ){
    expand_size <- expand*bin_size
  } else if(expand == 0){
    expand_size <- 0
  }

  for( i in 1:nrow(pairs) ){

    # expand feet
    pairs[i,"start1_bin"] <- pairs[i,"start1_bin"] - expand_size
    pairs[i,"end1_bin"  ] <- pairs[i,"end1_bin"  ] + expand_size
    pairs[i,"start2_bin"] <- pairs[i,"start2_bin"] - expand_size
    pairs[i,"end2_bin"  ] <- pairs[i,"end2_bin"  ] + expand_size

    if( pairs[i,"chr1"] == pairs[i,"chr2"] ){
          up_end   <- if( pairs[i,"start1_bin"] <= pairs[i,"start2_bin"] ) pairs[i,"end1_bin"] else pairs[i,"end2_bin"]
          dn_start <- max( pairs[i,"start1_bin"], pairs[i,"start2_bin"] )
          if( up_end > dn_start ){
            pairs[i,"start1_bin"] <- pairs[i,"start1_bin"] + expand_size
            pairs[i,"end1_bin"  ] <- pairs[i,"end1_bin"  ] - expand_size
            pairs[i,"start2_bin"] <- pairs[i,"start2_bin"] + expand_size
            pairs[i,"end2_bin"  ] <- pairs[i,"end2_bin"  ] - expand_size
          }
    }
  }

  # expanded feet can't start before the chromosome
  pairs$start1_bin <- pmax(pairs$start1_bin, 0)
  pairs$start2_bin <- pmax(pairs$start2_bin, 0)

  pairs$midpoint_1_bin <- as.numeric( floor( (pairs$start1 + pairs$end1)/2   / bin_size ) * bin_size )
  pairs$midpoint_2_bin <- as.numeric( floor( (pairs$start2 + pairs$end2)/2   / bin_size ) * bin_size )

  if( one_sample_analysis ){

    if (flag_norm == "abc"){
      get_c_results <- mclapply(seq_len(nrow(pairs)), function(i) {
        result <- get_c(
          row = pairs[i, ],
          chr_col_name = "chr1",
          tss_col_name = "start2",
          enh_col_name = "start1",
          bsize = bin_size,
          hic_path = hic_A
        )
        round(result, 4)
      }, mc.cores = num_cores)

      get_c_results_vector <- unlist(get_c_results)
      get_c_results_df <- data.frame(result = get_c_results_vector)
      pairs$c <- get_c_results_df$result

      if (flag_meta && exists("meta_cols")) {
        last_col <- unlist(pairs[, ncol(pairs)])
        C <- cbind(pairs[, 1:6], meta_cols, last_col)
        C[ , "chr1" ] <- as.character( C[ , "chr1" ] )
        C[ , "chr2" ] <- as.character( C[ , "chr2" ] )
        for (i in 1:nrow(C)) {
          try(cat(paste(as.character(C[i, ]), collapse="\t"), "\n", sep=""), silent=TRUE)
        }
      } else {
        for (i in 1:nrow(pairs)) {
          try(cat(paste(c(pairs[i, 1:6], pairs[i, ncol(pairs)]), collapse="\t"), "\n", sep=""), silent=TRUE)
        }
      }

    } else {

      A  <-  apply( pairs, 1, prefix_filter_bedpe )
      B  <-  mclapply( A, call_one_sample_straw_bedpe, mc.cores = num_cores, flag_inherent, flag_fix, formula, norm_factor1, aqua_factor1, shrink_wrap, split )

      rep_idx <- NULL

      if( isTRUE(flag_fix) ){
        C <- cbind( pairs[,1:6], data.frame(count = unlist(B), stringsAsFactors = FALSE) )
      } else if( isFALSE(flag_fix) ){
        if( split == "FALSE" ){
          rep_idx <- rep(seq_along(B), times = vapply(B, function(x) nrow(x$feet), integer(1)))
          result_list <- lapply(1:length(B), function(i) {
            f <- B[[i]]$feet
            data.frame(
              chr1   = pairs[i, "chr1"],
              start1 = f$start1,
              end1 = f$end1,
              chr2   = pairs[i, "chr2"],
              start2 = f$start2,
              end2 = f$end2,
              count  = B[[i]]$count1,
              stringsAsFactors = FALSE
            )
          })
          C <- do.call(rbind, result_list)
        } else {
          # Split mode
          result_list <- list()
          for( i in 1:length(B)){
            length <- as.numeric(lengths(B[i]))
            data_list <- lapply(1:length, function(j) {
              data.frame(
                chr1 = B[[as.character(i)]][[j]]$chr,
                start1 = B[[as.character(i)]][[j]]$feet1[1],
                end1 = B[[as.character(i)]][[j]]$feet1[2],
                chr2 = B[[as.character(i)]][[j]]$chr,
                start2 = B[[as.character(i)]][[j]]$feet2[1],
                end2 = B[[as.character(i)]][[j]]$feet2[2],
                count = B[[as.character(i)]][[j]]$count_1,
                stringsAsFactors = FALSE
              )
            })
            result_list <- c(result_list, data_list)
          }
          C <- do.call(rbind, result_list)
        }
      }

      if (flag_meta && exists("meta_cols")) {
        m <- if (is.null(rep_idx)) meta_cols else meta_cols[rep_idx, , drop = FALSE]
        rownames(m) <- NULL
        rownames(C) <- NULL
        C <- cbind(C[, 1:6], m, C[, 7])
      }

      C[ , "chr1" ] <- as.character( C[ , "chr1" ] )
      C[ , "chr2" ] <- as.character( C[ , "chr2" ] )

      for (i in 1:nrow(C)) {
        try(cat(paste(C[i,], collapse="\t"), "\n", sep=""), silent=TRUE)
      }
    }
  }

  if( two_sample_analysis ){

    if (flag_norm == "abc"){
      get_c_results <- mclapply(seq_len(nrow(pairs)), function(i) {
        result_sample_A <- get_c(
          row = pairs[i, ],
          chr_col_name = "chr1",
          tss_col_name = "start2",
          enh_col_name = "start1",
          bsize = bin_size,
          hic_path = hic_A
        )
        result_sample_B <- get_c(
          row = pairs[i, ],
          chr_col_name = "chr1",
          tss_col_name = "start2",
          enh_col_name = "start1",
          bsize = bin_size,
          hic_path = hic_B
        )
        delta <- result_sample_B - result_sample_A
        c(sample_A = round(result_sample_A, 4),
          sample_B = round(result_sample_B, 4),
          delta    = round(delta, 4))
      }, mc.cores = num_cores)

      get_c_results_df <- do.call(rbind, get_c_results)
      colnames(get_c_results_df) <- c("sample_A", "sample_B", "delta")
      pairs <- cbind(pairs, get_c_results_df)

      if (flag_meta && exists("meta_cols")) {
        C <- cbind(pairs[, 1:6], meta_cols, pairs[, (ncol(pairs)-2):ncol(pairs)])
        C[ , "chr1" ] <- as.character( C[ , "chr1" ] )
        C[ , "chr2" ] <- as.character( C[ , "chr2" ] )
        for (i in 1:nrow(C)) {
          try(cat(paste(as.character(C[i, ]), collapse="\t"), "\n", sep=""), silent=TRUE )
        }
      } else {
        for (i in 1:nrow(pairs)) {
          try(cat(paste(c(pairs[i, 1:6], pairs[i, (ncol(pairs)-2):ncol(pairs)]), collapse="\t"), "\n", sep=""), silent=TRUE )
        }
      }

     } else {
      A <- apply( pairs, 1, prefix_filter_bedpe )
      B <- mclapply( A, call_two_sample_straw_bedpe, mc.cores = num_cores, flag_inherent, flag_fix, formula, norm_factor1, aqua_factor1, norm_factor2, aqua_factor2 )

      rep_idx <- NULL

      if( isTRUE(flag_fix) ){
        result_matrix <- do.call(rbind, lapply(B, function(x) c(x$count_1, x$count_2)))
        C <- cbind( pairs[,1:6], data.frame(result_matrix, stringsAsFactors = FALSE) )
        colnames(C)[7:8] <- c("count1", "count2")
      } else if( isFALSE(flag_fix) ){
        rep_idx <- rep(seq_along(B), times = vapply(B, function(x) nrow(x$feet), integer(1)))
        result_list <- lapply(1:length(B), function(i) {
          f <- B[[i]]$feet
          data.frame(
            chr1   = pairs[i, "chr1"],
            start1 = f$start1,
            end1   = f$end1,
            chr2   = pairs[i, "chr2"],
            start2 = f$start2,
            end2   = f$end2,
            count1 = B[[i]]$count1,
            count2 = B[[i]]$count2,
            stringsAsFactors = FALSE
          )
        })
        C <- do.call(rbind, result_list)
      }

      C[ , "chr1" ] <- as.character( C[ , "chr1" ] )
      C[ , "chr2" ] <- as.character( C[ , "chr2" ] )
      C$delta <- round(C[,8] - C[,7], 3)

      if (isTRUE(flag_meta) && exists("meta_cols")) {
        m <- if (is.null(rep_idx)) meta_cols else meta_cols[rep_idx, , drop = FALSE]
        rownames(m) <- NULL
        rownames(C) <- NULL
        C_meta <- cbind(C[, 1:6], m, C[, 7:9])
        for( i in 1:nrow(C_meta) ){
          try(cat( paste( C_meta[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE)
        }
      } else {
        for( i in 1:nrow(C) ){
          try(cat( paste( C[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE)
        }
      }
    }
  }

  if(index == num_chunks) break
  pairs <- read.table(con, nrows=chunk_size, as.is = T)
}

close(con)
