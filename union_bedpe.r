suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(InteractionSet))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(S4Vectors))
options(stringsAsFactors = FALSE)
options(scipen = 999)
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)

# Internal function to find BEDPE format
.findBedpeFormat <- function(df) {
  df_cols <- colnames(df)[seq(1, 10)]
  cols <- c("seqnames1", "start1", "end1", "width1", "strand1",
            "seqnames2", "start2", "end2", "width2", "strand2")
  
  if (ncol(df) > 9) {
    if (identical(df_cols, cols) ||
        identical(df_cols, gsub("seqnames", "chr", cols)) ||
        identical(df_cols, gsub("seqnames", "chrom", cols))) {
      return("bedpe10")
    }
  }
  return("bedpe6")
}

.splitAnchors <- function(df, a1Coords = seq(1, 3), a2Coords = seq(4, 6), starts.in.df.are.0based) {
  a1 <- df[a1Coords]
  colnames(a1) <- c("seqnames", "start", "end")
  a2 <- df[a2Coords]
  colnames(a2) <- c("seqnames", "start", "end")
  
  a1 <- makeGRangesFromDataFrame(df = a1, starts.in.df.are.0based = starts.in.df.are.0based)
  a2 <- makeGRangesFromDataFrame(df = a2, starts.in.df.are.0based = starts.in.df.are.0based)
  
  suppressWarnings(GInteractions(a1, a2))
}

# Main function to convert DataFrame to GInteractions
.as_ginteractions <- function(df, keep.extra.columns = FALSE, starts.in.df.are.0based = FALSE) {
  df <- DataFrame(df)
  
  if (ncol(df) < 6) {
    stop("Improper dimensions in `df`. There must be at least 6 columns.")
  }
  
  bedpeFormat <- .findBedpeFormat(df)
  
  if (identical(bedpeFormat, "bedpe6")) {
    gi <- .splitAnchors(df, a1Coords = seq(1, 3), a2Coords = seq(4, 6), starts.in.df.are.0based = starts.in.df.are.0based)
    if (keep.extra.columns & ncol(df) > 6) {
      mcols(gi) <- df[seq(7, ncol(df))]
    }
  }
  
  if (identical(bedpeFormat, "bedpe10")) {
    gi <- .splitAnchors(df, a1Coords = seq(1, 3), a2Coords = seq(6, 8), starts.in.df.are.0based = starts.in.df.are.0based)
    if (keep.extra.columns & ncol(df) > 10) {
      mcols(gi) <- df[seq(11, ncol(df))]
    }
  }
  
  gi
}

# Alias function
makeGInteractionsFromDataFrame <- function(df, keep.extra.columns = FALSE, starts.in.df.are.0based = FALSE) {
  .as_ginteractions(df, keep.extra.columns, starts.in.df.are.0based)
}

# Single-file merging
.merge_single_pass <- function(bedpe) {
  
  if (nrow(bedpe) <= 1) {
    colnames(bedpe) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
    return(bedpe[, 1:6])
  }
  
  colnames(bedpe)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
  
  # Convert to GI object
  gi <- makeGInteractionsFromDataFrame(bedpe)

  overlaps <- as.data.frame(suppressWarnings(findOverlaps(
    gi, gi, type = "any",
    ignore.strand = TRUE,
    minoverlap = 2L
  )))
  overlaps <- overlaps[overlaps$queryHits != overlaps$subjectHits, ]
  
  # If no overlaps remain, return unique rows
  if (nrow(overlaps) == 0) {
    final_bedpe <- unique(bedpe[, 1:6])
    final_bedpe <- final_bedpe[order(final_bedpe$chr1, final_bedpe$start1), ]
    return(final_bedpe)
  }
  
  # Build undirected graph
  edges <- data.frame(
    from = as.character(overlaps$queryHits),
    to = as.character(overlaps$subjectHits)
  )
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Find connected components
  clusters        <- igraph::components(g)
  cluster_members <- split(V(g)$name, clusters$membership)
  
  # Merge each cluster
  merged_list <- vector("list", length(cluster_members))
  
  for (j in seq_along(cluster_members)) {
    indices <- as.numeric(cluster_members[[j]])
    comb    <- bedpe[indices, 1:6]
    chr1_val <- comb[1, 1]
    chr2_val <- comb[1, 4]
    merged_list[[j]] <- data.frame(
      chr1   = chr1_val,
      start1 = min(comb[, 2]),
      end1   = max(comb[, 3]),
      chr2   = chr2_val,
      start2 = min(comb[, 5]),
      end2   = max(comb[, 6]),
      stringsAsFactors = FALSE
    )
  }
  merged_bedpe <- do.call(rbind, merged_list)
  
  # Non-overlaps
  all_overlapping   <- unique(c(overlaps$queryHits, overlaps$subjectHits))
  non_overlaps_idx  <- setdiff(1:nrow(bedpe), all_overlapping)
  
  if (length(non_overlaps_idx) > 0) {
    non_overlaps <- bedpe[non_overlaps_idx, 1:6]
    colnames(non_overlaps) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
  } else {
    non_overlaps <- data.frame(
      chr1 = character(0), start1 = numeric(0), end1 = numeric(0),
      chr2 = character(0), start2 = numeric(0), end2 = numeric(0)
    )
  }
  
  final_bedpe <- rbind(merged_bedpe, non_overlaps)
  final_bedpe <- unique(final_bedpe)
  final_bedpe <- final_bedpe[order(final_bedpe$chr1, final_bedpe$start1), ]
  
  return(final_bedpe)
}

union_bedpe <- function(bedpe1, bedpe2, is_single_file = FALSE) {
  # Single file mode
  if (is_single_file) {
    
    current <- bedpe1[, 1:6]
    max_iterations <- max(1, nrow(current) - 1)
    
    # Iterate until no more merges occur
    for (iteration in 1:max_iterations) {
      previous_nrow <- nrow(current)
      current <- .merge_single_pass(current)
      if (nrow(current) == previous_nrow) break
    }
    
    return(current)
  }
  
  # Two file mode
  gi1 <- makeGInteractionsFromDataFrame(bedpe1)
  gi2 <- makeGInteractionsFromDataFrame(bedpe2)

  overlaps <- as.data.frame(suppressWarnings(findOverlaps(gi1, gi2, type = "any",
                                        ignore.strand = TRUE,
                                        minoverlap = 2L)))
  
  # If no overlaps, return combined unique rows
  if (nrow(overlaps) == 0) {
    final_bedpe <- unique(rbind(bedpe1[, 1:6], bedpe2[, 1:6]))
    colnames(final_bedpe) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
    final_bedpe <- final_bedpe[order(final_bedpe$chr1, final_bedpe$start1), ]
    return(final_bedpe)
  }
  
  ##################################################################
  ##                      Multi M-N Overlaps                      ##
  ##################################################################
  
  # identify multi m-n overlap intersection indices 
  dupQ <- overlaps$queryHits[duplicated(overlaps$queryHits)]
  dupS <- overlaps$subjectHits[duplicated(overlaps$subjectHits)]
  # vectorized tagging: TRUE=1, FALSE=0
  overlaps$x <- as.integer(overlaps$queryHits %in% dupQ)
  overlaps$y <- as.integer(overlaps$subjectHits  %in% dupS)
  # keep any row where x or y is 1
  multi_overlaps <- overlaps[overlaps$x != 0 | overlaps$y != 0, 1:2]

  if (nrow(multi_overlaps) > 0) {
    multi_overlaps$queryHits   <- paste0(multi_overlaps$queryHits,  "-1")
    multi_overlaps$subjectHits <- paste0(multi_overlaps$subjectHits, "-2")
    g               <- graph_from_data_frame(multi_overlaps, directed = T)
    clusters        <- igraph::components(g)
    cluster_members <- split(V(g)$name, clusters$membership)
    cluster_df      <- data.frame(
      cluster = rep(names(cluster_members), lengths(cluster_members)),
      index = unlist(cluster_members) )
    
    clusters     <- unique(cluster_df$cluster)
    multi_list <- vector("list", length(clusters))
    for (j in seq_along(clusters)) {
      clust_id <- clusters[j]
      df_clust <- cluster_df[cluster_df$cluster == clust_id, ]
      parts    <- do.call(rbind, strsplit(as.character(df_clust$index), "-", fixed = TRUE))
      df_clust$index <- as.numeric(parts[,1])
      df_clust$file  <- as.numeric(parts[,2])
      idx1  <- df_clust$index[df_clust$file == 1]
      idx2  <- df_clust$index[df_clust$file == 2]
      file1 <- bedpe1[idx1, 1:6]
      file2 <- bedpe2[idx2, 1:6]
      comb  <- rbind(file1, file2)
      multi_list[[j]] <- data.frame(
        chr1   = comb[1, 1],
        start1 = min(comb[,2]),
        end1   = max(comb[,3]),
        chr2   = comb[1, 4],
        start2 = min(comb[,5]),
        end2   = max(comb[,6])
      )
    }
    multi_unions <- do.call(rbind, multi_list)
  } else {
    multi_unions <- data.frame(
      chr1 = character(0), start1 = numeric(0), end1 = numeric(0),
      chr2 = character(0), start2 = numeric(0), end2 = numeric(0)
    )
  }

  #################################################################
  ##                     Simple 1-1 Overlaps                     ##
  #################################################################
  single_overlaps  <- overlaps[overlaps$x == 0 & overlaps$y == 0, 1:2]
  if (nrow(single_overlaps) > 0) {
    single_list     <- vector("list", nrow(single_overlaps))
    for (i in seq_len(nrow(single_overlaps))) {
      idx1 <- single_overlaps[i, 1]
      idx2 <- single_overlaps[i, 2]
      comb <- rbind(bedpe1[idx1, 1:6], bedpe2[idx2, 1:6])
      single_list[[i]] <- data.frame(
        chr1   = comb[1, 1],
        start1 = min(comb[,2]),
        end1   = max(comb[,3]),
        chr2   = comb[1, 4],
        start2 = min(comb[,5]),
        end2   = max(comb[,6])
      )
    }
    single_unions <- do.call(rbind, single_list)
  } else {
    single_unions <- data.frame(
      chr1 = character(0), start1 = numeric(0), end1 = numeric(0),
      chr2 = character(0), start2 = numeric(0), end2 = numeric(0)
    )
  }

  ##################################################################
  ##                         Non Overlaps                         ##
  ##################################################################
  
  non_overlaps_bedpe1 <- setdiff(1:nrow(bedpe1), unique(overlaps$queryHits))
  non_overlaps_bedpe2 <- setdiff(1:nrow(bedpe2), unique(overlaps$subjectHits))
  
  non_overlaps <- rbind(
    bedpe1[non_overlaps_bedpe1, 1:6, drop = FALSE],
    bedpe2[non_overlaps_bedpe2, 1:6, drop = FALSE])
  
  if (nrow(non_overlaps) > 0) {
    colnames(non_overlaps) <- c(
      "chr1", "start1", "end1",
      "chr2", "start2", "end2")
  } else {
    non_overlaps <- data.frame(
      chr1 = character(0), start1 = numeric(0), end1 = numeric(0),
      chr2 = character(0), start2 = numeric(0), end2 = numeric(0)
    )
  }
  
  final_bedpe <- rbind(
    multi_unions, 
    single_unions,
    non_overlaps)
  
  final_bedpe <- unique(final_bedpe)
  final_bedpe <- final_bedpe[order(final_bedpe$chr1, final_bedpe$start1), ]
  
  return(final_bedpe)
}

#################################################################
##                          Arguments                          ##
#################################################################

path_bedpe1 <- args[1]
path_bedpe2 <- args[2]

bedpe1 <- read.table(path_bedpe1, as.is = T)
bedpe2 <- read.table(path_bedpe2, as.is = T)

# Detect single-file mode
is_single_file <- (path_bedpe1 == path_bedpe2) || 
                  (nrow(bedpe1) == nrow(bedpe2) && all(bedpe1 == bedpe2))

union <- union_bedpe(bedpe1, bedpe2, is_single_file = is_single_file)

# For two-file mode, run iterative cleanup to catch new containments
if (!is_single_file) {
  union <- union_bedpe(union, union, is_single_file = TRUE)
}

for (i in 1:nrow(union)) { 
  cat(paste(union[i, ], collapse = "\t"), "\n", sep = "")
}
