suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(InteractionSet))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(S4Vectors))

options(stringsAsFactors = FALSE)
options(scipen = 999)
options(warn = -1 )


args           <- commandArgs( trailingOnly = TRUE )



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

# Internal function to split anchors
.splitAnchors <- function(df, a1Coords = seq(1, 3), a2Coords = seq(4, 6), starts.in.df.are.0based) {
  a1 <- df[a1Coords]
  colnames(a1) <- c("seqnames", "start", "end")
  a2 <- df[a2Coords]
  colnames(a2) <- c("seqnames", "start", "end")
  
  a1 <- makeGRangesFromDataFrame(df = a1, starts.in.df.are.0based = starts.in.df.are.0based)
  a2 <- makeGRangesFromDataFrame(df = a2, starts.in.df.are.0based = starts.in.df.are.0based)
  
  GInteractions(a1, a2)
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



union_bedpe <- function(bedpe1,bedpe2){
  
  # add offset of 1bp to not merge bedpes that are side-by-side
  bedpe2_offset     <- bedpe2
  offset            <- 1
  bedpe2_offset[,5] <- bedpe2_offset[,5] + offset 
  
  
  # covert to GI object
  gi1 <- makeGInteractionsFromDataFrame(bedpe1)
  gi2 <- makeGInteractionsFromDataFrame(bedpe2_offset)
  
  
  # perform intersection
  overlaps <- as.data.frame(findOverlaps(gi1, gi2, type = "any", ignore.strand = TRUE))
  if (nrow(overlaps) == 0) {
    cat("No overlaps found between BEDPE files.")
    stop()
  }
  
  
  # It's better now to handle the two kinds of intersections separately:
  #  a. simple 1-1 intersections
  #  b. multi  m-n intersections (from either side)
  
  
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

  # append the suffixes
  if (nrow(multi_overlaps) > 0) {
    multi_overlaps$queryHits   <- paste0(multi_overlaps$queryHits,  "-1")
    multi_overlaps$subjectHits <- paste0(multi_overlaps$subjectHits, "-2")
  }
  
  # assign unique ids to multi m-n overlaps
  g               <- graph_from_data_frame(multi_overlaps, directed = T)
  clusters        <- igraph::components(g)
  cluster_members <- split(V(g)$name, clusters$membership)
  cluster_df      <- data.frame(
    cluster = rep(names(cluster_members), lengths(cluster_members)),
    index = unlist(cluster_members) )
  
  # obtain new coordinates for multi m-n overlaps
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

    chr <- unique(c(comb[,1], comb[,4]))
    multi_list[[j]] <- data.frame(
      chr1   = chr,
      start1 = min(comb[,2]),
      end1   = max(comb[,3]),
      chr2   = chr,
      start2 = min(comb[,5]),
      end2   = max(comb[,6])
    )
  }

  multi_unions <- do.call(rbind, multi_list)
  
  #################################################################
  ##                     Simple 1-1 Overlaps                     ##
  #################################################################
  
  single_overlaps  <- overlaps[overlaps$x == 0 & overlaps$y == 0,1:2]
  single_list     <- vector("list", nrow(single_overlaps))

  for (i in seq_len(nrow(single_overlaps))) {
    idx1 <- single_overlaps[i,1]
    idx2 <- single_overlaps[i,2]
    comb <- rbind(bedpe1[idx1, ], bedpe2[idx2, ])

    chr <- unique(c(comb[,1], comb[,4]))
    single_list[[i]] <- data.frame(
      chr1   = chr,
      start1 = min(comb[,2]),
      end1   = max(comb[,3]),
      chr2   = chr,
      start2 = min(comb[,5]),
      end2   = max(comb[,6])
    )
  }

  single_unions <- do.call(rbind, single_list)

  ##################################################################
  ##                         Non Overlaps                         ##
  ##################################################################
  
  non_overlaps_bedpe1 <- setdiff(1:nrow(bedpe1), unique(overlaps$queryHits))
  non_overlaps_bedpe2 <- setdiff(1:nrow(bedpe2), unique(overlaps$subjectHits))
  
  non_overlaps <- rbind(
    bedpe1[non_overlaps_bedpe1,],
    bedpe2[non_overlaps_bedpe2,] )
  
  colnames(non_overlaps) <- c(
    "chr1", "start1", "end1",
    "chr2", "start2", "end2")
  
  final_bedpe <- rbind(
    multi_unions, 
    single_unions,
    non_overlaps)
  
  final_bedpe <- final_bedpe[order(
    final_bedpe$chr1,
    final_bedpe$start1),]
  
  return(final_bedpe)
  
}


#################################################################
##                          Arguments                          ##
#################################################################

path_bedpe1 <- args[1]
path_bedpe2 <- args[2]


bedpe1 <- read.table(path_bedpe1, as.is = T)
bedpe2 <- read.table(path_bedpe2, as.is = T)


union <- union_bedpe(
  bedpe1,
  bedpe2)


for( i in 1:nrow(union) ){ 
  try(cat( paste( union[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE) 
}
