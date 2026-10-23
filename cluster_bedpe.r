suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(igraph))

options(stringsAsFactors = FALSE)
options(scipen = 999)
options(warn = -1 )

args           <- commandArgs( trailingOnly = TRUE )



#################################################################
##                          Arguments                          ##
#################################################################


file_bedpe <- args[1]
flank      <- as.numeric(args[2])

bedpe <- read.table( file_bedpe, as.is = TRUE )

flag_meta <- FALSE

if (ncol(bedpe) < 6) {

  cat("Please use a minimum of 6 column bedpe, exiting....\n")
  q(save = "no")

} else if (ncol(bedpe) > 6) {

  flag_meta <- TRUE
  meta_cols <- bedpe[, 7:ncol(bedpe)]
  meta_cols <- data.frame(meta_cols, stringsAsFactors = FALSE)
  num_meta_cols <- ncol(meta_cols)
  bedpe <- bedpe[, 1:6]
}

colnames(bedpe) <- c(
  "chr1",
  "start1",
  "end1",
  "chr2",
  "start2",
  "end2")


# cis only; a trans pair would bridge two chromosomes into one cluster
n_trans <- sum(bedpe$chr1 != bedpe$chr2)

if (n_trans > 0) {
  cat(paste0("Found ", n_trans, " trans pairs. cluster_bedpe only handles cis pairs, exiting....\n"))
  q(save = "no")
}

# add flank
bedpe$start1 <- bedpe$start1 - flank
bedpe$start2 <- bedpe$start2 - flank

bedpe$end1 <- bedpe$end1 + flank
bedpe$end2 <- bedpe$end2 + flank

# sort rows; meta_cols is re-sorted below to match
ord   <- order(bedpe$chr1, bedpe$start1)
bedpe <- bedpe[ord, ]

if (flag_meta) {
  meta_cols <- meta_cols[ord, , drop = FALSE]
}

n <- nrow(bedpe)

# cluster loops that share an anchor
# one row per anchor; idx remembers which loop it came from
feet <- data.frame(
  chr   = c(bedpe$chr1,   bedpe$chr2),
  start = c(bedpe$start1, bedpe$start2),
  end   = c(bedpe$end1,   bedpe$end2),
  idx   = seq_len(2L * n),
  stringsAsFactors = FALSE)

# sort anchors by position
o  <- order(feet$chr, feet$start, feet$end)
fs <- feet[o, ]

# furthest end seen before this row, per chromosome
prev_max <- ave(
  as.numeric(fs$end), fs$chr,
  FUN = function(z) c(-Inf, cummax(z)[-length(z)]))

# number the runs; >= means anchors that only touch stay apart
run_of_foot         <- integer(2L * n)
run_of_foot[fs$idx] <- cumsum(fs$start >= prev_max)

# run number for each loop's left and right anchor
run1 <- run_of_foot[seq_len(n)]
run2 <- run_of_foot[seq_len(n) + n]

# each loop is an edge between its two runs
el <- cbind(run1, run2)
storage.mode(el) <- "integer"

cl <- components(graph_from_edgelist(el, directed = FALSE))$membership[run1]

# number clusters in order of appearance
bedpe$cluster <- paste0("cluster-", match(cl, unique(cl)))

# remove flank
bedpe$start1 <- bedpe$start1 + flank
bedpe$start2 <- bedpe$start2 + flank

bedpe$end1 <- bedpe$end1 - flank
bedpe$end2 <- bedpe$end2 - flank

if(flag_meta){
  bedpe <- cbind(
    bedpe[,1:6,],
    meta_cols,
    bedpe$cluster
  )
} else {
  bedpe <- cbind(
    bedpe[,1:6,],
    bedpe$cluster
  )
}


tryCatch(
  write.table(
    bedpe,
    "",
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE),
  error = function(e) {
    if (identical(conditionMessage(e), "ignoring SIGPIPE signal")) {
      quit(save = "no", status = 0)
    }
    stop(e)
  })
