suppressPackageStartupMessages(library(strawr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dbscan))

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

  if( end - start > 5000000 ){
    cat("\n  Please keep loop search space within 5Mb\n")
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
  
  for( i in 1:nrow(sparMat)){
    tad_matrix[
      as.character(sparMat[i,"x"]),
      as.character(sparMat[i,"y"])] <- sparMat[i,"counts"]
    tad_matrix[
      as.character(sparMat[i,"y"]),
      as.character(sparMat[i,"x"])] <- sparMat[i,"counts"]
  }
  
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
  
  for ( j in 0:(ncol(tad_matrix) - 1) ){
    for ( i in 1:(nrow(tad_matrix)-j) ){
      tad_matrix_bg_removed[i  ,i+j  ] <-
        tad_matrix[i,i+j  ] - background[j+1]
    }
  }
  
  # Standardize
  tad_matrix_sd <- tad_matrix_bg_removed
  
  for ( j in 0:(ncol(tad_matrix_bg_removed) - 1) ){
    for ( i in 1:(nrow(tad_matrix_bg_removed)-j) ){
      tad_matrix_sd[i  ,i+j  ] <-
        tad_matrix_bg_removed[i ,i+j  ] / (foreground[j+1] - background[j+1])
    }
  }
  
  A <- tad_matrix_sd
  for(i in 2:nrow(A)){
    for(j in 1:(i-1)){
      A[i,j] <- 0
    }
  }
  
  A <- zero_diag(A,1)
  
  if(mode == "loop"){
    
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
    cat("Under development, use --mode loop\n")
    q(save = "no")
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
  
  
  if(mode == "flare"){
    cat("Under development, use loop\n")
    q(save = "no")
  }
  
  
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
    
    for( i in 1:nrow(sparMat)){
      tad_matrix[
        as.character(sparMat[i,"x"]),
        as.character(sparMat[i,"y"])] <- sparMat[i,"counts"]
      tad_matrix[
        as.character(sparMat[i,"y"]),
        as.character(sparMat[i,"x"])] <- sparMat[i,"counts"]
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
    
    for ( j in 0:(ncol(tad_matrix) - 1) ){
      for ( i in 1:(nrow(tad_matrix)-j) ){
        tad_matrix_bg_removed[i  ,i+j  ] <-
          tad_matrix[i,i+j  ] - background[j+1]
      }
    }
    
    # Standardize
    tad_matrix_sd <- tad_matrix_bg_removed
    
    for ( j in 0:(ncol(tad_matrix_bg_removed) - 1) ){
      for ( i in 1:(nrow(tad_matrix_bg_removed)-j) ){
        tad_matrix_sd[i  ,i+j  ] <-
          tad_matrix_bg_removed[i ,i+j  ] / (foreground[j+1] - background[j+1])
      }
    }
    
    A <- tad_matrix_sd
    for(i in 2:nrow(A)){
      for(j in 1:(i-1)){
        A[i,j] <- 0
      }
    }
    
    A <- zero_diag(A,1)
    
    if(mode == "loop"){
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
    
    
  }
  
}
