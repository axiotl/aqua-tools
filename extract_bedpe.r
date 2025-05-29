suppressPackageStartupMessages(library(strawr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dbscan))
suppressPackageStartupMessages(library(dplyr))

options(stringsAsFactors = FALSE)
options(scipen = 999)
options(warn = -1 )



args           <- commandArgs( trailingOnly = TRUE )

analysis_type  <<- args[1]

zero_diag <- function(
    matrix,
    width ){
  
  if( nrow(matrix) == ncol(matrix) ){
    for ( j in 0:width ){
      for ( i in 1:(nrow(matrix)-j) ){
        matrix[i+j,i  ] <- 0
        matrix[i  ,i+j] <- 0
      }
    }
    return( matrix )
  } else {
    cat("zero_diag: matrix isn't square\n")
  }
}



if(analysis_type == "range"){
  
  #################################################################
  ##                          Arguments                          ##
  #################################################################
  
  path_hic   <- args[2]
  chr        <- as.character(args[3])
  start      <- as.numeric(args[4])
  end        <- as.numeric(args[5])
  bin_size   <- as.numeric(args[6])
  score      <- as.numeric(args[7])
  sample_dir <- args[8]
  radius     <- as.character(args[9])
  mode       <- as.character(args[10])
  min_dist   <- as.numeric(args[11])
  norm  <- "NONE"
  unit  <- "BP"
  
  ##################################################################
  ##                          Code-block                          ##
  ##################################################################
  
  if( any( grepl( "chr", chr, ignore.case = T ) ) == FALSE ){
    cat("\n  Please add 'chr' prefix for chromosome\n")
    quit(save="no")
  }
  
  ## 'chr' prefix handler:
  hic_A_chroms <- readHicChroms(path_hic)$name
  
  if( sum( grepl("chr", hic_A_chroms ) )  > 0 ){ flag_strip_chr <- FALSE }
  if( sum( grepl("chr", hic_A_chroms ) ) == 0 ){ flag_strip_chr <- TRUE  }
  
  if( flag_strip_chr ){
    interval_chr <- sub( "chr", "", chr )
  } else {
    interval_chr <- chr
  }
  
  interval_start <-   floor(start/bin_size)*bin_size
  interval_end   <- ceiling(  end/bin_size)*bin_size
  
  if((interval_end - interval_start)/bin_size < 10 ){
    cat("Range less than 10 bins? Expand your imagination! \n")
    cat("Please enter a range of at least 10 bins \n")
    q( save = "no" )
  }
  
  sparMat <- straw(
    norm,
    path_hic,
    paste( interval_chr, interval_start, interval_end, sep = ":"  ),
    paste( interval_chr, interval_start, interval_end, sep = ":"  ),
    unit,
    bin_size)
  
  if(nrow(sparMat)==0){
    q( save = "no" )
  }
  
  if (flag_strip_chr) {
    interval_chr <- paste0("chr", interval_chr)
  }
  
  tad_matrix <- matrix(
    data = 0,
    nrow = length(seq(interval_start,interval_end,bin_size)),
    ncol = length(seq(interval_start,interval_end,bin_size))
  )
  
  colnames(tad_matrix) <- seq(interval_start,interval_end,bin_size)
  rownames(tad_matrix) <- seq(interval_start,interval_end,bin_size)

  pos <- rownames(tad_matrix)

  # convert every sparMat$x and sparMat$y into integer row/col indices
  i1 <- match(as.character(sparMat[ , "x"]), pos)
  j1 <- match(as.character(sparMat[ , "y"]), pos)

  # bulk-assign counts into the upper triangle in one C-level call
  tad_matrix[cbind(i1, j1)] <- sparMat[ , "counts"]

  # mirror for lower triangle
  tad_matrix[cbind(j1, i1)] <- sparMat[ , "counts"]
  
  power_law_path <- list.files(
    sample_dir,
    pattern = paste0("inherentStats",".txt"),
    full.names = T)
  
  if(length(power_law_path) == 0){
    cat("Power laws unavailable for this sample! \n")
    q( save = "no" )
  }
  
  power_laws <- read.table(
    power_law_path,
    as.is = T,
    skip = 1 )
  
  colnames(power_laws) <- c("off","on","res") ; power_laws <- power_laws[! power_laws$off %in% "off",]
  
  power_laws$off  <- as.numeric(power_laws$off)
  power_laws$on   <- as.numeric(power_laws$on)
  power_laws$res  <- as.numeric(power_laws$res)
  power_laws      <- power_laws[power_laws$res == bin_size,]
  
  if(nrow(power_laws) == 0){
    cat("Power laws unavailable at this resolution! \n")
    q( save = "no" )
  }
  
  if(nrow(power_laws) < ncol(tad_matrix)){
    
    diff <- ncol(tad_matrix) - nrow(power_laws)
    
    power_laws <- rbind(
      power_laws,
      do.call(
        "rbind",
        replicate(diff, power_laws[nrow(power_laws),], simplify = FALSE)) )
  }
  
  background <- power_laws[,"off"]
  foreground <- power_laws[,"on" ]
  
  # Remove background
  tad_matrix_bg_removed <- tad_matrix

  # n = number of bins = nrow = ncol of tad_matrix
  n   <- ncol(tad_matrix)

  # compute offset for every [i,j] cell: 0 on diagonal
  off <- col(tad_matrix) - row(tad_matrix)

  # map offsets to background indices, clamp lower triangle to 1
  idx          <- off + 1
  idx[ idx < 1 ] <- 1

  # build a full matrix of background values
  bg_mat      <- matrix(background[idx], nrow = n, ncol = n)
  bg_mat[ off < 0 ] <- 0

  # subtract in one vectorized step
  tad_matrix_bg_removed <- tad_matrix - bg_mat
  
  # Standardize
  tad_matrix_sd <- tad_matrix_bg_removed

  # compute the per-offset denominator vector
  denom <- foreground - background

  # map each cell offset to the correct denom entry
  idx       <- off + 1
  idx[idx < 1] <- 1
  den_mat   <- matrix(denom[idx], nrow = n, ncol = n)
  den_mat[off < 0] <- 1        # so lower-triangle divisions are by 1

  # vectorized divide
  tad_matrix_sd <- tad_matrix_bg_removed / den_mat

  A <- tad_matrix_sd
  
  # zero out lower triangle
  A[lower.tri(A)] <- 0
  

  if(mode == "loop"){
    A <- zero_diag(A,1)
    A[A < score] <- 0
    
    if(sum(A)>0){
      
      A[A>0] <- 1
      
      
      df <- which(A > 0, arr.ind = T)
      
      for( i in 1:nrow(df)){
        
        loops <- data.frame(
          chr1   = as.character(interval_chr),
          start1 = as.numeric(rownames(A)[df[i,1]]),
          end1   = as.numeric(rownames(A)[df[i,1]]) + bin_size,
          chr2   = as.character(interval_chr),
          start2 = as.numeric(rownames(A)[df[i,2]]),
          end2   = as.numeric(rownames(A)[df[i,2]]) + bin_size,
          score  = round(A[df[i,1],df[i,2]],3),
          stringsAsFactors = FALSE)
        
        loops <- loops[loops$start2 - loops$start1 >= min_dist,]
        
        if(nrow(loops)>0){
          try(cat( paste( loops[,1:6], collapse = "\t"), "\n", sep = "" ), silent=TRUE)
        }
        
      }
      
    } else {
      trap <- 1 # do nothing
    }
    
  }
  
  if(mode == "glob"){
    A <- zero_diag(A,1)
    A[A < score] <- 0
    
    if(sum(A)>0){
      
      A[A>0] <- 1
      
      positions <- as.data.frame(which(A > 0, arr.ind = TRUE))
      positions <- positions[positions$row != positions$col,]
      
      if(nrow(positions > 0)){
        
        clust     <- dbscan(positions, eps = radius, minPts = 2)
        
        positions$cluster <- clust$cluster
        
        positions <- positions[positions$cluster != 0,]
        
        if(nrow(positions) > 0){
          
          num_clusters <- unique(positions$cluster)
          
          bedpe <- data.frame()
          
          for(cluster in num_clusters){
            
            indices <- positions[positions$cluster == cluster,]
            
            start1  <- as.numeric(rownames(A)[min(indices$row)])
            end1    <- as.numeric(rownames(A)[max(indices$row)])
            start2  <- as.numeric(colnames(A)[min(indices$col)])
            end2    <- as.numeric(colnames(A)[max(indices$col)])
            
            loops <- data.frame(
              chr1    = as.character(interval_chr),
              start1  = start1,
              end1    = end1 + bin_size,
              chr2    = as.character(interval_chr),
              start2  = start2,
              end2    = end2 + bin_size,
              #cluster = cluster,
              stringsAsFactors = FALSE)
            
            bedpe <- rbind(bedpe, loops)
          }
          
          midpoint_x <- bedpe$start1 + (bedpe$end1 - bedpe$start1)/2
          midpoint_y <- bedpe$start2 + (bedpe$end2 - bedpe$start2)/2
          
          
          # height of equilateral triangle = distance of centre of loop to diagonal
          bedpe$midpoint_dist_to_diag <- ceiling((sqrt(3)/2)*(midpoint_y - midpoint_x))
          
          bedpe <- bedpe[bedpe$midpoint_dist_to_diag >= min_dist,1:6]
          
          if(nrow(bedpe)>0){
            for(i in 1:nrow(bedpe)){try(cat( paste( bedpe[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE)}
          }
        }
        
      }
      
    } else {
      trap <- 1 # do nothing
    }
  }
  
  if(mode == "flare"){
    A <- zero_diag(A,1)
    
    A[A < score] <- 0
    
    if(sum(A)>0){
      
      A[A>0] <- 1
      
      
      df <- which(A > 0, arr.ind = T)
      real_df = data.frame(row = df[,"row"], col = df[,"col"], step = rownames(df))
      df_ordered <- real_df[order(real_df$row), ]
      
      
      
      # Combine rows where col differs by parameter radius
      result <- df_ordered %>%
        arrange(row, col) %>% # Ensure data is ordered by row and col
        group_by(row) %>%
        mutate(group = cumsum(c(1, diff(col) > radius))) %>% # Create groups where col differs by more than 2
        group_by(row, group) %>%
        summarise(
          start_col = min(col),      # Keep the smallest col in the group
          end_col = max(col),        # Keep the largest col in the group
          step = first(step),        # Keep the first step value
          combined_count = n(),       # Count how many rows were combined
          .groups = "drop" # Suppresses the warning
        ) %>%
        ungroup() %>%
        mutate(end_row = row) %>%    # Add an end_row column (equal to start row)
        select(row, start_col, end_col, step, combined_count, end_row) # Rearrange columns
      
      
      # Step 2: Combine rows where start_col and end_col are the same and row differs by less-than parameter radius
      final_result <- result %>%
        arrange(start_col, end_col, row) %>%
        group_by(start_col, end_col) %>%
        mutate(row_group = cumsum(c(1, diff(row) > radius))) %>%
        group_by(start_col, end_col, row_group) %>%
        summarise(
          start_row = min(row),
          end_row = max(row),
          start_col = first(start_col),
          end_col = first(end_col),
          step = first(step),
          combined_count = first(combined_count),          # Sum of combined rows
          combined_count_row = n(),                       # Number of rows combined in this step
          .groups = "drop" # Suppresses the warning
        ) %>%
        ungroup() %>%
        select(start_row, end_row, start_col, end_col, step, combined_count, combined_count_row)
      
      final_result = final_result %>% arrange(start_row)
      for( i in 1:nrow(final_result)){
        iter = final_result[i,]
        loops <- data.frame(
          chr1   = as.character(interval_chr),
          start1 = as.numeric(rownames(A)[iter$start_row]), #   as.numeric(rownames(A)[df[i,1]]),
          end1   = as.numeric(rownames(A)[iter$start_row]) + bin_size*iter$combined_count_row,
          chr2   = as.character(interval_chr),
          start2 = as.numeric(rownames(A)[iter$start_col]),
          end2   = as.numeric(rownames(A)[iter$start_col]) + bin_size*iter$combined_count,
          score  = round(mean(A[iter$start_row: (iter$start_row+ iter$combined_count_row-1),
                                iter$start_col: (iter$start_col+iter$combined_count-1)]),3),
          stringsAsFactors = FALSE)
        
        # loops <- loops[loops$start2 - loops$start1 >= min_dist,]
        
        if(nrow(loops)>0){
          try(cat( paste( loops[,1:6], collapse = "\t"), "\n", sep = "" ), silent=TRUE)
        }
        
      }
      
    } else {
      trap <- 1 # do nothing
    }
    
    
  }
  if(mode == "minimal"){
    
    A[A < score] <- 0
    
    if(sum(A) > 0){
      
      high_score = score + 3 * score / 7 # Default 1 for score = 0.7
      
      A[A > 0 & A < high_score ] <- 1
      A[A > high_score ] <- 2
      
      
      df <- which(A > 0, arr.ind = T)
      indices <- which(A > 0, arr.ind = TRUE)
      
      # Extract the values from A corresponding to these indices
      values <- A[indices]
      df <- data.frame(row = indices[, 1], col = indices[, 2], value = values)
      
      unique_elements <- unique(c(df$row, df$col))
      # Initialize a vector to store the sums
      # Count the number of 1s and 2s for each element
      # Initialize a data frame to store the results
      result_df <- data.frame(
        element = unique_elements,
        count_1 = 0,  # Column for counting 1s
        count_2 = 0   # Column for counting 2s
      )
      
      for (element in unique_elements) {
        # Count 1s where the element appears in the row
        count_1_row <- sum(df$value[df$row == element] == 1)
        # Count 1s where the element appears in the column
        count_1_col <- sum(df$value[df$col == element] == 1)
        # Total count of 1s for the element
        result_df$count_1[result_df$element == element] <- count_1_row + count_1_col
        # Count 2s where the element appears in the row
        count_2_row <- sum(df$value[df$row == element] == 2)
        # Count 2s where the element appears in the column
        count_2_col <- sum(df$value[df$col == element] == 2)
        # Total count of 2s for the element
        result_df$count_2[result_df$element == element] <- count_2_row + count_2_col
      }
      result_df$count_sum = result_df$count_1 + 2* result_df$count_2 #  %>% arrange()
      
      result_df = result_df %>% arrange(-count_sum)
      
      result_df$prop_count_sum = result_df$count_sum  / max(result_df$count_sum )
      
      # Assign diagonal values if the element is on the diagonal
      for (element in unique_elements) {
        if (element <= nrow(A) && element <= ncol(A)) {  # Ensure within matrix bounds
          result_df$diagonal_value[result_df$element == element] <- A[element, element]
        }
      }
      
      # Get the top-1 element
      top_1_element <- result_df$element[1]
      
      # Identify neighbors (2 before, 1 before, 1 after, 2 after)
      all_elements <- result_df$element
      top_1_index <- which(all_elements == top_1_element)
      neighbors <- c(top_1_element-2, top_1_element-1, top_1_element+1, top_1_element+2)  # all_elements[pmax(1, top_1_index - 2):pmin(length(all_elements), top_1_index + 2)]
      

      # remove edge cases end and start -----------------------------------------
      neighbors = neighbors[neighbors > 0 & neighbors <= nrow(A)]
      
      l_neigh = length(neighbors) +1
      
      # Initialize new columns
      result_df$direct_connection <- 0  # "1" or "2" if direct connection to top-1
      result_df$neighborhood_1 <- 0  # "1" if at least one connection with a neighbor for 1
      result_df$neighborhood_2 <- 0  # "1" if at least one connection with a neighbor for 2
      
      # Loop over all elements in result_df
      for (element in all_elements) {
        if (element == top_1_element) next  # Skip the top-1 itself
        
        # Check direct connection to top-1
        if (A[element, top_1_element] == 1 || A[top_1_element, element] == 1) {
          result_df$direct_connection[result_df$element == element] <- 1
        } else if (A[element, top_1_element] == 2 || A[top_1_element, element] == 2) {
          result_df$direct_connection[result_df$element == element] <- 2
        }
        
        # Check connections with neighbors
        for (neighbor in neighbors) {
          if (neighbor == element) next  # Skip itself
          
          if (A[element, neighbor] == 1 || A[neighbor, element] == 1) {
            result_df$neighborhood_1[result_df$element == element] <- result_df$neighborhood_1[result_df$element == element]+1
          }
          if (A[element, neighbor] == 2 || A[neighbor, element] == 2) {
            result_df$neighborhood_2[result_df$element == element] <- result_df$neighborhood_2[result_df$element == element]+1
          }
        }
      }
      
      
      
      # Omnicomprehensive score ----------------------------------------------------
      
      # A bit of magic numbers..
      result_df$score = 5*result_df$prop_count_sum + result_df$diagonal_value + result_df$direct_connection + result_df$neighborhood_1 / (2*l_neigh) + result_df$neighborhood_2 / l_neigh
      result_df$score[1] = result_df$score[1] + 4 ## comes out on top.  
      
      # Filter: filtering element with a neighbor with higher score, or with very low number of contacts
      # Identify elements to remove
      up_df <- result_df %>%
        rowwise() %>%
        mutate(has_higher_adjacent = any(result_df$score[result_df$element %in% c(element - 1, element - 2, element + 1, element + 2)] > score),
               has_double_higher_adjacent = any(result_df$score[result_df$element %in% c(element - 1, element - 2, element + 1, element + 2)] > 2*score),
               has_low_interactions = prop_count_sum < 0.2,
               has_very_low_interactions = prop_count_sum < 0.05) %>%
        ungroup() 
      
      
      # remove individual with either very low adj or very low contacts --------
      removed_df = up_df %>%
        filter(!has_double_higher_adjacent) %>%
        filter(!has_very_low_interactions)  %>%
        filter(! (has_higher_adjacent | has_low_interactions) ) ## This is the difference between High cut or normal
      
      # Extract interactions only between non-removed loci ----------------------------------------------------
      A <- zero_diag(A,1)
      smallerA = A[sort(removed_df$element),sort(removed_df$element)]
      
      if(sum(smallerA)>0){
        
        smallerA[smallerA>0] <- 1
        
        
        df <- which(smallerA > 0, arr.ind = T)
        
        for( i in 1:nrow(df)){
          
          loops <- data.frame(
            chr1   = as.character(interval_chr),
            start1 = as.numeric(rownames(smallerA)[df[i,1]]),
            end1   = as.numeric(rownames(smallerA)[df[i,1]]) + bin_size,
            chr2   = as.character(interval_chr),
            start2 = as.numeric(rownames(smallerA)[df[i,2]]),
            end2   = as.numeric(rownames(smallerA)[df[i,2]]) + bin_size,
            score  = round(smallerA[df[i,1],df[i,2]],3),
            stringsAsFactors = FALSE)
          
          loops <- loops[loops$start2 - loops$start1 >= min_dist,]
          
          if(nrow(loops)>0){
            try(cat( paste( loops[,1:6], collapse = "\t"), "\n", sep = "" ), silent=TRUE)
          }
          
        }
        
      } else {
        trap <- 1 # do nothing
      }
    }
  }
}

if(analysis_type == "TAD"){
  
  #################################################################
  ##                          Arguments                          ##
  #################################################################
  
  
  path_hic   <- args[2]
  path_tad   <- args[3]
  bin_size   <- as.numeric(args[4])
  score      <- as.numeric(args[5])
  sample_dir <- args[6]
  radius     <- as.numeric(args[7])
  mode       <- as.character(args[8])
  min_dist   <- as.numeric(args[9])
  
  if (!file.exists(path_tad)) {
    cat("TAD file not found \n")
    q(save = "no")
  }
  
  norm  <- "NONE"
  unit  <- "BP"
  
  ##################################################################
  ##                          Code-block                          ##
  ##################################################################
  
  tads    <- read.table(path_tad, as.is = T)[,1:3] ; colnames(tads) <- c("chr","start","end")
  tads$id <- paste(tads$chr,tads$start,tads$end,sep=":")
  
  
  if( any( grepl( "chr", tads$chr, ignore.case = T ) ) == FALSE ){
    cat("\n  Please add 'chr' prefix for chromosome in TAD file\n")
    quit(save="no")
  }
  
  ## 'chr' prefix handler:
  hic_A_chroms <- readHicChroms(path_hic)$name
  
  if( sum( grepl("chr", hic_A_chroms ) )  > 0 ){ flag_strip_chr <- FALSE }
  if( sum( grepl("chr", hic_A_chroms ) ) == 0 ){ flag_strip_chr <- TRUE  }
  
  
  power_law_path <- list.files(
    sample_dir,
    pattern = paste0("inherentStats",".txt"),
    full.names = T)
  
  if(length(power_law_path) == 0){
    cat("Power laws unavailable for this sample! \n")
    q( save = "no" )
  }
  
  power_laws <- read.table(
    power_law_path,
    as.is = T,
    skip = 1 )
  
  colnames(power_laws) <- c("off","on","res") ; power_laws <- power_laws[! power_laws$off %in% "off",]
  
  power_laws$off  <- as.numeric(power_laws$off)
  power_laws$on   <- as.numeric(power_laws$on)
  power_laws$res  <- as.numeric(power_laws$res)
  power_laws      <- power_laws[power_laws$res == bin_size,]
  
  
  
  for(tad in 1:nrow(tads)){
    
    if( flag_strip_chr ){
      interval_chr <- sub( "chr", "", tads[tad,"chr"] )
    } else {
      interval_chr <- tads[tad,"chr"]
    }
    
    interval_start <- tads[tad,"start"]
    interval_end   <- tads[tad,  "end"]
    interval_tag   <- tads[tad,   "id"]
    
    interval_start <-   floor(interval_start/bin_size)*bin_size
    interval_end   <- ceiling(  interval_end/bin_size)*bin_size
    
    sparMat <- straw(
      norm,
      path_hic,
      paste( interval_chr, interval_start, interval_end, sep = ":"  ),
      paste( interval_chr, interval_start, interval_end, sep = ":"  ),
      unit,
      bin_size)
    
    if(nrow(sparMat)==0){
      next
    }
    
    if (flag_strip_chr) {
      interval_chr <- paste0("chr", interval_chr)
    }
    
    tad_matrix <- matrix(
      data = 0,
      nrow = length(seq(interval_start,interval_end,bin_size)),
      ncol = length(seq(interval_start,interval_end,bin_size))
    )
    
    colnames(tad_matrix) <- seq(interval_start,interval_end,bin_size)
    rownames(tad_matrix) <- seq(interval_start,interval_end,bin_size)
    
    pos <- rownames(tad_matrix)

    # convert every sparMat$x and sparMat$y into integer row/col indices
    i1 <- match(as.character(sparMat[ , "x"]), pos)
    j1 <- match(as.character(sparMat[ , "y"]), pos)

    # bulk-assign counts into the upper triangle in one C-level call
    tad_matrix[cbind(i1, j1)] <- sparMat[ , "counts"]

    # mirror for lower triangle
    tad_matrix[cbind(j1, i1)] <- sparMat[ , "counts"]
    
    if(nrow(power_laws) < ncol(tad_matrix)){
      
      diff <- ncol(tad_matrix) - nrow(power_laws)
      
      power_laws <- rbind(
        power_laws,
        do.call(
          "rbind",
          replicate(diff, power_laws[nrow(power_laws),], simplify = FALSE)) )
    }
    
    background <- power_laws[,"off"]
    foreground <- power_laws[,"on" ]
    
    # Remove background
    tad_matrix_bg_removed <- tad_matrix

    # n = number of bins = nrow = ncol of tad_matrix
    n   <- ncol(tad_matrix)

    # compute offset for every [i,j] cell: 0 on diagonal
    off <- col(tad_matrix) - row(tad_matrix)

    # map offsets to background indices, clamp lower triangle to 1
    idx          <- off + 1
    idx[ idx < 1 ] <- 1

    # build a full matrix of background values
    bg_mat      <- matrix(background[idx], nrow = n, ncol = n)
    bg_mat[ off < 0 ] <- 0

    # subtract in one vectorized step
    tad_matrix_bg_removed <- tad_matrix - bg_mat
    
    # Standardize
    tad_matrix_sd <- tad_matrix_bg_removed
    
    # compute the per-offset denominator vector
    denom <- foreground - background

    # map each cell offset to the correct denom entry
    idx       <- off + 1
    idx[idx < 1] <- 1
    den_mat   <- matrix(denom[idx], nrow = n, ncol = n)
    den_mat[off < 0] <- 1        # so lower-triangle divisions are by 1

    # vectorized divide
    tad_matrix_sd <- tad_matrix_bg_removed / den_mat

    A <- tad_matrix_sd
    
    # zero out lower triangle
    A[lower.tri(A)] <- 0
    
    
    if(mode == "loop"){
      A <- zero_diag(A,1)
      
      A[A < score] <- 0
      
      if(sum(A)>0){
        
        df <- which(A > 0, arr.ind = T)
        
        for( i in 1:nrow(df)){
          
          loops <- data.frame(
            chr1   = as.character(interval_chr),
            start1 = as.numeric(rownames(A)[df[i,1]]),
            end1   = as.numeric(rownames(A)[df[i,1]]) + bin_size,
            chr2   = as.character(interval_chr),
            start2 = as.numeric(rownames(A)[df[i,2]]),
            end2   = as.numeric(rownames(A)[df[i,2]]) + bin_size,
            score  = round(A[df[i,1],df[i,2]],3),
            tag    = interval_tag,
            stringsAsFactors = FALSE )
          
          loops <- loops[loops$start2 - loops$start1 >= min_dist,]
          
          
          if(nrow(loops)>0){
            try(cat( paste( loops[,1:6], collapse = "\t"), "\n", sep = "" ), silent=TRUE)
          }
          
          
        }
        
      } else {
        trap <- 1 # do nothing
      }
    }
    
    if(mode == "glob"){
      A <- zero_diag(A,1)
      
      A[A < score] <- 0
      
      if(sum(A)>0){
        
        A[A>0] <- 1
        
        positions <- as.data.frame(which(A > 0, arr.ind = TRUE))
        positions <- positions[positions$row != positions$col,]
        
        if(nrow(positions>0)){
          
          clust     <- dbscan(positions, eps = radius, minPts = 2)
          
          positions$cluster <- clust$cluster
          
          positions <- positions[positions$cluster != 0,]
          
          if(nrow(positions) > 0){
            
            num_clusters <- unique(positions$cluster)
            
            bedpe <- data.frame()
            
            for(cluster in num_clusters){
              
              indices <- positions[positions$cluster == cluster,]
              
              start1  <- as.numeric(rownames(A)[min(indices$row)])
              end1    <- as.numeric(rownames(A)[max(indices$row)])
              start2  <- as.numeric(colnames(A)[min(indices$col)])
              end2    <- as.numeric(colnames(A)[max(indices$col)])
              
              loops <- data.frame(
                chr1    = as.character(interval_chr),
                start1  = start1,
                end1    = end1 + bin_size,
                chr2    = as.character(interval_chr),
                start2  = start2,
                end2    = end2 + bin_size,
                #cluster = cluster,
                stringsAsFactors = FALSE)
              
              
              bedpe <- rbind(bedpe, loops)
            }
            
            midpoint_x <- bedpe$start1 + (bedpe$end1 - bedpe$start1)/2
            midpoint_y <- bedpe$start2 + (bedpe$end2 - bedpe$start2)/2
            
            
            # height of equilateral triangle = distance of centre of loop to diagonal
            bedpe$midpoint_dist_to_diag <- ceiling((sqrt(3)/2)*(midpoint_y - midpoint_x))
            
            bedpe <- bedpe[bedpe$midpoint_dist_to_diag >= min_dist,1:6]
            
            if(nrow(bedpe)>0){
              for(i in 1:nrow(bedpe)){try(cat( paste( bedpe[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE)}
            }
          }
          
        }
        
        
      } else {
        trap <- 1 # do nothing
      }
    }
    
    if(mode == "flare"){
      A <- zero_diag(A,1)
      
      A[A < score] <- 0
      
      if(sum(A)>0){
        
        A[A>0] <- 1
        
        
        df <- which(A > 0, arr.ind = T)
        real_df = data.frame(row = df[,"row"], col = df[,"col"], step = rownames(df))
        df_ordered <- real_df[order(real_df$row), ]
        
        
        
        # Combine rows where col differs by parameter radius
        result <- df_ordered %>%
          arrange(row, col) %>% # Ensure data is ordered by row and col
          group_by(row) %>%
          mutate(group = cumsum(c(1, diff(col) > radius))) %>% # Create groups where col differs by more than 2
          group_by(row, group) %>%
          summarise(
            start_col = min(col),      # Keep the smallest col in the group
            end_col = max(col),        # Keep the largest col in the group
            step = first(step),        # Keep the first step value
            combined_count = n(),       # Count how many rows were combined
            .groups = "drop" # Suppresses the warning
          ) %>%
          ungroup() %>%
          mutate(end_row = row) %>%    # Add an end_row column (equal to start row)
          select(row, start_col, end_col, step, combined_count, end_row) # Rearrange columns
        
        
        # Step 2: Combine rows where start_col and end_col are the same and row differs by less-than parameter radius
        final_result <- result %>%
          arrange(start_col, end_col, row) %>%
          group_by(start_col, end_col) %>%
          mutate(row_group = cumsum(c(1, diff(row) > radius))) %>%
          group_by(start_col, end_col, row_group) %>%
          summarise(
            start_row = min(row),
            end_row = max(row),
            start_col = first(start_col),
            end_col = first(end_col),
            step = first(step),
            combined_count = first(combined_count),          # Sum of combined rows
            combined_count_row = n(),                       # Number of rows combined in this step
            .groups = "drop" # Suppresses the warning
          ) %>%
          ungroup() %>%
          select(start_row, end_row, start_col, end_col, step, combined_count, combined_count_row)
        
        final_result = final_result %>% arrange(start_row)
        for( i in 1:nrow(final_result)){
          iter = final_result[i,]
          loops <- data.frame(
            chr1   = as.character(interval_chr),
            start1 = as.numeric(rownames(A)[iter$start_row]), #   as.numeric(rownames(A)[df[i,1]]),
            end1   = as.numeric(rownames(A)[iter$start_row]) + bin_size*iter$combined_count_row,
            chr2   = as.character(interval_chr),
            start2 = as.numeric(rownames(A)[iter$start_col]),
            end2   = as.numeric(rownames(A)[iter$start_col]) + bin_size*iter$combined_count,
            score  = round(mean(A[iter$start_row: (iter$start_row+ iter$combined_count_row-1),
                                  iter$start_col: (iter$start_col+iter$combined_count-1)]),3),
            stringsAsFactors = FALSE)
          
          # loops <- loops[loops$start2 - loops$start1 >= min_dist,]
          
          if(nrow(loops)>0){
            try(cat( paste( loops[,1:6], collapse = "\t"), "\n", sep = "" ), silent=TRUE)
          }
          
        }
        
      } else {
        trap <- 1 # do nothing
      }
      
      
    }
    if(mode == "minimal"){
      
      A[A < score] <- 0
      
      if(sum(A) > 0){
        
        high_score = score + 3 * score / 7 # Default 1 for score = 0.7
        
        A[A > 0 & A < high_score ] <- 1
        A[A > high_score ] <- 2
        
        
        df <- which(A > 0, arr.ind = T)
        indices <- which(A > 0, arr.ind = TRUE)
        
        # Extract the values from A corresponding to these indices
        values <- A[indices]
        df <- data.frame(row = indices[, 1], col = indices[, 2], value = values)
        
        unique_elements <- unique(c(df$row, df$col))
        # Initialize a vector to store the sums
        # Count the number of 1s and 2s for each element
        # Initialize a data frame to store the results
        result_df <- data.frame(
          element = unique_elements,
          count_1 = 0,  # Column for counting 1s
          count_2 = 0   # Column for counting 2s
        )
        
        for (element in unique_elements) {
          # Count 1s where the element appears in the row
          count_1_row <- sum(df$value[df$row == element] == 1)
          # Count 1s where the element appears in the column
          count_1_col <- sum(df$value[df$col == element] == 1)
          # Total count of 1s for the element
          result_df$count_1[result_df$element == element] <- count_1_row + count_1_col
          # Count 2s where the element appears in the row
          count_2_row <- sum(df$value[df$row == element] == 2)
          # Count 2s where the element appears in the column
          count_2_col <- sum(df$value[df$col == element] == 2)
          # Total count of 2s for the element
          result_df$count_2[result_df$element == element] <- count_2_row + count_2_col
        }
        result_df$count_sum = result_df$count_1 + 2* result_df$count_2 #  %>% arrange()
        
        result_df = result_df %>% arrange(-count_sum)
        
        result_df$prop_count_sum = result_df$count_sum  / max(result_df$count_sum )
        
        # Assign diagonal values if the element is on the diagonal
        for (element in unique_elements) {
          if (element <= nrow(A) && element <= ncol(A)) {  # Ensure within matrix bounds
            result_df$diagonal_value[result_df$element == element] <- A[element, element]
          }
        }
        
        # Get the top-1 element
        top_1_element <- result_df$element[1]
        
        # Identify neighbors (2 before, 1 before, 1 after, 2 after)
        all_elements <- result_df$element
        top_1_index <- which(all_elements == top_1_element)
        neighbors <- c(top_1_element-2, top_1_element-1, top_1_element+1, top_1_element+2)  # all_elements[pmax(1, top_1_index - 2):pmin(length(all_elements), top_1_index + 2)]
        
        neighbors = neighbors[neighbors > 0 & neighbors <= nrow(A)]
        
        l_neigh = length(neighbors) +1
        
        # Initialize new columns
        result_df$direct_connection <- 0  # "1" or "2" if direct connection to top-1
        result_df$neighborhood_1 <- 0  # "1" if at least one connection with a neighbor for 1
        result_df$neighborhood_2 <- 0  # "1" if at least one connection with a neighbor for 2
        
        # Loop over all elements in result_df
        for (element in all_elements) {
          if (element == top_1_element) next  # Skip the top-1 itself
          
          # Check direct connection to top-1
          if (A[element, top_1_element] == 1 || A[top_1_element, element] == 1) {
            result_df$direct_connection[result_df$element == element] <- 1
          } else if (A[element, top_1_element] == 2 || A[top_1_element, element] == 2) {
            result_df$direct_connection[result_df$element == element] <- 2
          }
          
          # Check connections with neighbors
          for (neighbor in neighbors) {
            if (neighbor == element) next  # Skip itself
            
            if (A[element, neighbor] == 1 || A[neighbor, element] == 1) {
              result_df$neighborhood_1[result_df$element == element] <- result_df$neighborhood_1[result_df$element == element]+1
            }
            if (A[element, neighbor] == 2 || A[neighbor, element] == 2) {
              result_df$neighborhood_2[result_df$element == element] <- result_df$neighborhood_2[result_df$element == element]+1
            }
          }
        }
        
        
        
        # Omnicomprehensive score ----------------------------------------------------
        
        # A bit of magic numbers..
        result_df$score = 5*result_df$prop_count_sum + result_df$diagonal_value + result_df$direct_connection + result_df$neighborhood_1 / (2*l_neigh) + result_df$neighborhood_2 / l_neigh
        result_df$score[1] = result_df$score[1] + 4 ## comes out on top.  
        
        # Filter: filtering element with a neighbor with higher score, or with very low number of contacts
        # Identify elements to remove
        up_df <- result_df %>%
          rowwise() %>%
          mutate(has_higher_adjacent = any(result_df$score[result_df$element %in% c(element - 1, element - 2, element + 1, element + 2)] > score),
                 has_double_higher_adjacent = any(result_df$score[result_df$element %in% c(element - 1, element - 2, element + 1, element + 2)] > 2*score),
                 has_low_interactions = prop_count_sum < 0.2,
                 has_very_low_interactions = prop_count_sum < 0.05) %>%
          ungroup() 
        
        
        # remove individual with either very low adj or very low contacts --------
        removed_df = up_df %>%
          filter(!has_double_higher_adjacent) %>%
          filter(!has_very_low_interactions)  %>%
          filter(! (has_higher_adjacent | has_low_interactions) ) ## This is the difference between High cut or normal
        
        # Extract interactions only between non-removed loci ----------------------------------------------------
        A <- zero_diag(A,1)
        smallerA = A[sort(removed_df$element),sort(removed_df$element)]
        
        if(sum(smallerA)>0){
          
          smallerA[smallerA>0] <- 1
          
          
          df <- which(smallerA > 0, arr.ind = T)
          
          for( i in 1:nrow(df)){
            
            loops <- data.frame(
              chr1   = as.character(interval_chr),
              start1 = as.numeric(rownames(smallerA)[df[i,1]]),
              end1   = as.numeric(rownames(smallerA)[df[i,1]]) + bin_size,
              chr2   = as.character(interval_chr),
              start2 = as.numeric(rownames(smallerA)[df[i,2]]),
              end2   = as.numeric(rownames(smallerA)[df[i,2]]) + bin_size,
              score  = round(smallerA[df[i,1],df[i,2]],3),
              stringsAsFactors = FALSE)
            
            loops <- loops[loops$start2 - loops$start1 >= min_dist,]
            
            if(nrow(loops)>0){
              try(cat( paste( loops[,1:6], collapse = "\t"), "\n", sep = "" ), silent=TRUE)
            }
            
          }
          
        } else {
          trap <- 1 # do nothing
        }
      }
    }
  }
}
