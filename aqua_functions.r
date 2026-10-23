# aqua_functions.r - shared helpers for aqua_tools


# ============================================================
#  Binning
# ------------------------------------------------------------
#  The single bp <-> bin conversion for aqua_tools.
#  Half-open bp in, closed bin indices out. Bins are 0-indexed:
#  bin k covers [k*R, (k+1)*R)
#
# ============================================================

interval_to_bins <- function(s, e, R) {

  if (length(R) != 1L || !is.finite(R) || R <= 0) {
    stop("interval_to_bins: R must be a single positive number")
  }

  if (length(s) == 0L || length(e) == 0L) {
    return(list(first = numeric(0), last = numeric(0), n = numeric(0)))
  }
  n <- max(length(s), length(e))
  if (n %% length(s) != 0L || n %% length(e) != 0L) {
    stop("interval_to_bins: lengths of s and e are not compatible")
  }
  s <- rep_len(s, n)
  e <- rep_len(e, n)

  k1 <- floor(s / R)
  k2 <- pmax(floor((e - 1) / R), k1)

  list(first = k1, last = k2, n = k2 - k1 + 1)
}


# straw-style coordinate string, "chr:first_bin_start:last_bin_start"
# Returns one string per interval.
bins_to_query <- function(chr, b, R) {

  if (length(chr) == 0L || length(b$first) == 0L) {
    return(character(0))
  }
  n <- max(length(chr), length(b$first))
  if (n %% length(chr) != 0L || n %% length(b$first) != 0L) {
    stop("bins_to_query: lengths of chr and b are not compatible")
  }

  sprintf(
    "%s:%.0f:%.0f",
    rep_len(chr, n),
    rep_len(b$first, n) * R,
    rep_len(b$last,  n) * R
  )
}


# Bin-start coordinates, one per bin in the interval
bins_to_grid <- function(b, R) {
  g <- Map(function(k1, k2) (k1:k2) * R, b$first, b$last)
  if (length(g) == 1L) g[[1L]] else g
}


# Half-open bp span of the binned interval, [start, end)
bins_to_bp <- function(b, R) {

  cbind(
    start = b$first * R,
    end   = (b$last + 1) * R
  )
}


# Binify a bedpe: expand each pair into one row per anchor1-bin x anchor2-bin
# combination, using interval_to_bins for each anchor
binify_bedpe <- function(bedpe, R) {

  if (nrow(bedpe) == 0L) return(bedpe)

  b1 <- interval_to_bins(as.numeric(bedpe[[2]]), as.numeric(bedpe[[3]]), R)
  b2 <- interval_to_bins(as.numeric(bedpe[[5]]), as.numeric(bedpe[[6]]), R)

  n_out <- b1$n * b2$n
  row   <- rep(seq_len(nrow(bedpe)), n_out)
  k     <- sequence(n_out) - 1
  n2    <- b2$n[row]

  k1 <- b1$first[row] + k %/% n2
  k2 <- b2$first[row] + k %%  n2

  out <- bedpe[row, , drop = FALSE]
  out[[2]] <- k1 * R
  out[[3]] <- (k1 + 1) * R
  out[[5]] <- k2 * R
  out[[6]] <- (k2 + 1) * R
  rownames(out) <- NULL

  out
}


# ============================================================
#  Matrices
# ============================================================

# Zero the main diagonal and the 'width' diagonals either side of it
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
