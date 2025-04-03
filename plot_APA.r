suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(parallel))

options(scipen = 999)

Args <- commandArgs(trailingOnly=T)

###########################################################################
###########################################################################
###                                                                     ###
###                                FUNCTIONS                            ###
###                                                                     ###
###########################################################################
###########################################################################


# Function to compute distance between bedpe feet
  compute_distance <- function(i, pairs) {
    if (pairs[i, 1] == pairs[i, 4]) {
      mp1 <- as.numeric(pairs[i, 2]) + (as.numeric(pairs[i, 3]) - as.numeric(pairs[i, 2])) / 2
      mp2 <- as.numeric(pairs[i, 5]) + (as.numeric(pairs[i, 6]) - as.numeric(pairs[i, 5])) / 2
      return(abs(mp2 - mp1))
    }
    return(NA)  # return NA for inter-chromosomal pairs
  }


###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if( length(Args) == 14 ) {

  matrix_file   <- Args[1]
  sample1       <- Args[2]
  max_cap       <- Args[3]
  max_cap_delta <- Args[4]
  out_file      <- Args[5]
  win_size      <- Args[6]
  res           <- Args[7]
  pairs         <- Args[8]
  flag_norm     <- Args[9]
  path_mgs      <- Args[10]
  flag_loop_norm <- Args[11]
  num_loops      <- as.numeric(Args[12])
  out_dir        <- Args[13]
  min_cap        <- Args[14]

  if( ! flag_norm %in% c( "blank", "none", "cpm", "aqua" ) ){
    cat("Norm should strictly be none, cpm or aqua in lower case \n")
    q( save = "no" )
  }


  pairs <- as.data.frame( read.table( pairs, as.is = TRUE ) )
  pairs <- pairs[ , 1:6 ]


  cat( sprintf( "sample1        <- %s\n", sample1        ) )
  cat( sprintf( "out_file       <- %s\n", out_file       ) )
  cat( sprintf( "bin_size       <- %s\n", res            ) ) 

  matrix_raw <- read.table(matrix_file, header=FALSE, sep="\t")

  ## CPM or AQuA factors
  
  mergeStats_A <- read.table( path_mgs, as.is = T )
  spikeVar     <- ncol(mergeStats_A)
  hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
  mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
  total1       <- sum( hg_total1 , mm_total1 )
  norm_factor1 <- 1000000 / total1
  aqua_factor1 <- hg_total1 / mm_total1


  cat(sprintf("norm_factor    <- %f\n", norm_factor1))
  cat(sprintf("aqua_factor    <- %f\n", aqua_factor1))
  

# set normalization factor for spike-in and non-spike-in data.
  # non-spike-in data defaults to cpm.
  # spike in data defaults to aqua.
  
  if(spikeVar == 1){
    if(flag_norm == "blank"){
      norm_factor1 <- norm_factor1
      aqua_factor1 <- 1
      flag_norm <- "cpm"
    } else if(flag_norm == "none"){
      norm_factor1 <- 1
      aqua_factor1 <- 1
    } else if(flag_norm == "cpm"){
      norm_factor1 <- norm_factor1
      aqua_factor1 <- 1
    } else if(flag_norm == "aqua"){
      cat("\n\n--norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
      norm_factor1 <- norm_factor1
      aqua_factor1 <- 1
      flag_norm <- "cpm"
    }
  }else if(spikeVar == 2){
    if(flag_norm == "blank"){
      norm_factor1 <- norm_factor1
      aqua_factor1 <- aqua_factor1
      flag_norm <- "aqua"
    } else if(flag_norm == "none"){
      norm_factor1 <- 1
      aqua_factor1 <- 1
    } else if(flag_norm == "cpm"){
      norm_factor1 <- norm_factor1
      aqua_factor1 <- 1
    }else if(flag_norm == "aqua"){
      norm_factor1 <- norm_factor1
      aqua_factor1 <- aqua_factor1
    }
  }

  cat(sprintf("normalization  <- %s\n", flag_norm      ) )
  cat(sprintf("loop norm      <- %s\n", flag_loop_norm    )) 

  # Multiply matrix by the norm and aqua factors
  matrix_norm <- matrix_raw * norm_factor1 * aqua_factor1
 
  
  if (flag_loop_norm == "TRUE") {
    matrix_norm <- matrix_norm / num_loops
  }

  write.table(matrix_norm, file = paste0(out_dir, "/matrix.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)

  # Check if max_cap is "no_cap" and set the cap accordingly
  if( max_cap == "no_cap" ) { 
    cap <- max(matrix_norm) 
  } else { 
    cap <- as.numeric(max_cap) 
  }

  # Check if min_cap is "no_cap" and set the min_cap accordingly
  if( min_cap == "no_cap" ) { 
    min_cap_value <- 0 
  } else { 
    min_cap_value <- as.numeric(min_cap) 
  }

  # Print the max_cap and min_cap values
  cat(sprintf("max_cap        <- %.2f\n", round(cap, 2)))
  cat(sprintf("min_cap        <- %.2f\n", round(min_cap_value, 2)))

  # Define colors and breakList
  breakList <- seq(min_cap_value, cap, by = (cap - min_cap_value) / 100)
  color_ramp_red  <-  colorRampPalette( c( "white", "red") )( 100 )


  ## define window sizes and resolution
  win_size <- as.numeric( win_size )
  res      <- as.numeric( res )
  res      <- res/1000


  ## define row and column names of APA plot
  cols <- c( paste( "-", rev( seq( res, res*win_size, res ) ), "kb", sep="" ), 
             "", 
             paste(           seq( res, res*win_size, res )  , "kb", sep="" ) 
            )
  
  rownames( matrix_norm) <- cols
  colnames( matrix_norm ) <- cols

  
  ## plotting
  p1 <- pheatmap(  matrix_norm,
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   border_color = NA, 
                   color = color_ramp_red, 
                   breaks = breakList, 
                   show_colnames = TRUE,
                   show_rownames = FALSE    
                )

  g <- arrangeGrob( p1[[4]], 
                    nrow = 1,
                    top = paste( sample1, " (n.loop=", num_loops, ")", sep = "" )
                  )

    ggsave( file = out_file, 
          g, 
          width = 4.5, 
          height = 4.5, 
          device = "pdf"
        ) 
  
  ## PATCH: https://github.com/tidyverse/ggplot2/issues/2787
  file.exists("Rplots.pdf")
  file.remove("Rplots.pdf")

  distance_between_feet <- c()
  
  # Use mclapply to parallelize the computation
  distance_between_feet <- mclapply(1:nrow(pairs), compute_distance, pairs = pairs, mc.cores = detectCores())
  
  # Unlist the result to get a vector
  distance_between_feet <- unlist(distance_between_feet)
  
  distance_between_feet <- distance_between_feet[!is.na(distance_between_feet)]

  distribution <- data.frame( distance = round(distance_between_feet) ) 

  if( nrow( distribution ) < 50 ){
    bins <- 10
  } else if( nrow( distribution ) <= 500 ){
    bins <- 30
  } else {
    bins <- 100
  }

  out_file_2 <- paste( dirname(out_file), "bedpe_size-distribution.pdf", sep = "/" )

  pdf( out_file_2,
       width = 7,
       height = 6,
       useDingbats = FALSE 
      )

  p2 <- ggplot( distribution, aes( x = distance ) ) + 
        geom_histogram( bins = bins ) + 
        xlab( "Distance between bedpe feet" ) + 
        ylab( "Count" ) +
        xlim( 0 , ceiling( quantile(distribution$distance)[4] ) )  +
        theme_classic()

    suppressWarnings( print( p2 ) )

    dev.off()
 
    cat("\nOutput folder created at", out_dir, "\n\n") 

}


###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if(length(Args) == 17) {

  matrix_fileA    <- Args[1]
  matrix_fileB    <- Args[2]
  sampleA         <- Args[3]
  sampleB         <- Args[4]
  max_cap         <- Args[5]
  max_cap_delta   <- Args[6]
  out_file        <- Args[7]
  win_size        <- Args[8]
  res             <- Args[9]
  pairs           <- Args[10]
  flag_norm       <- Args[11]
  path_mgsA       <- Args[12]
  path_mgsB       <- Args[13]
  flag_loop_norm  <- Args[14]
  num_loops       <- as.numeric(Args[15])
  out_dir         <- Args[16]
  min_cap         <- Args[17]


  pairs <- as.data.frame( read.table( pairs, as.is = TRUE ) )
  pairs <- pairs[ , 1:6 ]


  cat("\n")
  cat(sprintf("sampleA         <- %s\n", sampleA        ))
  cat(sprintf("sampleB         <- %s\n", sampleB        ))
  cat(sprintf("out_file        <- %s\n", out_file       ))


  matrix_rawA <- read.table(matrix_fileA, header=FALSE, sep="\t")
  matrix_rawB <- read.table(matrix_fileB, header=FALSE, sep="\t")

  ## CPM or AQuA factors for Sample A
  mergeStats_A  <- read.table(path_mgsA, as.is = T)
  spikeVar_A    <- ncol(mergeStats_A)
  hg_total1     <- as.numeric(mergeStats_A["valid_interaction_rmdup", 1])
  mm_total1     <- as.numeric(mergeStats_A["valid_interaction_rmdup", 2])
  total1        <- sum(hg_total1, mm_total1)
  norm_factor1  <- 1000000 / total1
  aqua_factor1  <- hg_total1 / mm_total1

  ## CPM or AQuA factors for Sample B
  mergeStats_B  <- read.table(path_mgsB, as.is = T)
  spikeVar_B    <- ncol(mergeStats_B)
  hg_total2     <- as.numeric(mergeStats_B["valid_interaction_rmdup", 1])
  mm_total2     <- as.numeric(mergeStats_B["valid_interaction_rmdup", 2])
  total2        <- sum(hg_total2, mm_total2)
  norm_factor2  <- 1000000 / total2
  aqua_factor2  <- hg_total2 / mm_total2

  if( ! flag_norm %in% c("blank",  "none", "cpm", "aqua" ) ){
    cat("Norm should strictly be none, cpm or aqua in lower case \n")
    q( save = "no" )
  }

  # set normalization factor for spike-in and non-spike-in data.
  # non-spike-in data defaults to cpm
  # spike-in data defaults to aqua
  
  
  if( spikeVar_A == 1 || spikeVar_B == 1){
    if (flag_norm == "blank"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
      flag_norm <- "cpm"
    } else if (flag_norm == "none"){
      norm_factor1 <- 1 ; aqua_factor1 <- 1
      norm_factor2 <- 1 ; aqua_factor2 <- 1
    } else if (flag_norm == "cpm"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
    } else if (flag_norm == "aqua"){
      cat("\n\n--norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
      flag_norm <- "cpm"
    }
  } else if (spikeVar_A == 2 || spikeVar_B == 2){
    if (flag_norm == "blank"){
      norm_factor1 <- norm_factor1
      aqua_factor1 <- aqua_factor1
      norm_factor2 <- norm_factor2
      aqua_factor2 <- aqua_factor2
      flag_norm <- "aqua"
    } else if (flag_norm == "none"){
      norm_factor1 <- 1 ; aqua_factor1 <- 1
      norm_factor2 <- 1 ; aqua_factor2 <- 1
    } else if (flag_norm == "cpm"){
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
    } else if (flag_norm == "aqua"){
      norm_factor1 <- norm_factor1
      aqua_factor1 <- aqua_factor1
      norm_factor2 <- norm_factor2
      aqua_factor2 <- aqua_factor2
    }
  }

  cat(sprintf("norm_factor1    <- %f\n", norm_factor1 ))
  cat(sprintf("aqua_factor1    <- %f\n", aqua_factor1 ))
  cat(sprintf("norm_factor2    <- %f\n", norm_factor2 ))
  cat(sprintf("aqua_factor2    <- %f\n", aqua_factor2 ))
  cat(sprintf("normalization   <- %s\n", flag_norm      ))
  cat(sprintf("loop norm       <- %s\n", flag_loop_norm    ))

  # For Sample A
  matrix_normA <- matrix_rawA * norm_factor1 * aqua_factor1

  # For Sample B
  matrix_normB <- matrix_rawB * norm_factor2 * aqua_factor2

  if (flag_loop_norm == "TRUE") {
    matrix_normA <- matrix_normA / num_loops
    matrix_normB <- matrix_normB / num_loops
  }

  delta <- matrix_normB - matrix_normA

  write.table(delta, file = paste0(out_dir, "/matrix.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)

  if( max_cap       == "no_cap" ) { 
    cap <- max( matrix_normA, matrix_normB ) 
  } else { cap <- as.numeric(max_cap)  }

  if( min_cap == "no_cap" ) { 
    min_cap_value <- 0 
  } else { min_cap_value <- as.numeric(min_cap) }

  if( max_cap_delta == "no_cap" ) { 
    delta_cap <- max( abs(delta) ) 
  } else { delta_cap <- as.numeric(max_cap_delta) }


  ## define colors
  breakList       <- seq(min_cap_value, cap, by = (cap - min_cap_value) / 100)
  breakList_delta <- seq(-delta_cap, delta_cap, by = delta_cap/50 )

  cat(sprintf("max_cap          <- %.2f\n", cap))
  cat(sprintf("min_cap          <- %.2f\n", min_cap_value))
  cat(sprintf("max_cap_delta    <- %.2f\n", delta_cap))

  color_ramp_red <- colorRampPalette(c(     "white",                      "red"))(100)
  color_ramp_vio <- colorRampPalette(c("dodgerblue", "white", "mediumvioletred"))(100)


  ## define window sizes and resolution
  win_size <- as.numeric(win_size)
  res      <- as.numeric(res)
  res      <- res/1000


  ## define row and column names of APA plot
  cols <- c( paste( "-", rev( seq( res, res*win_size, res ) ), "kb", sep="" ), 
             "", 
             paste(           seq( res, res*win_size, res )  , "kb", sep="" ) 
            )
  
  rownames(matrix_normA) <- cols ; rownames(matrix_normB) <- cols ; rownames(delta) <- cols
  colnames(matrix_normA) <- cols ; colnames(matrix_normB) <- cols ; colnames(delta) <- cols


  ## plotting
  p1 <- pheatmap(  matrix_normA, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   border_color = NA, 
                   color = color_ramp_red, 
                   breaks = breakList, 
                   show_colnames = TRUE, 
                   show_rownames = FALSE     
                )

  p2 <- pheatmap(  matrix_normB, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   border_color = NA, 
                   color = color_ramp_red, 
                   breaks = breakList, 
                   show_colnames = TRUE, 
                   show_rownames = FALSE      
                )

  p3 <- pheatmap( delta, 
                  cluster_rows = FALSE, 
                  cluster_cols = FALSE, 
                  border_color = NA, 
                  color = color_ramp_vio, 
                  breaks = breakList_delta, 
                  show_colnames = TRUE, 
                  show_rownames = FALSE 
                )

  g <- arrangeGrob( p1[[4]], p2[[4]], p3[[4]],
                    nrow = 1,
                    top = paste(sampleA, ", ", sampleB, ", and delta", " (n.loop=",  num_loops, ")", sep = "" )
                  )

  ggsave( file = out_file, 
          g, 
          width = 14, 
          height = 4.5, 
          device = "pdf"
        )
  
  ## PATCH: https://github.com/tidyverse/ggplot2/issues/2787
  file.exists("Rplots.pdf")
  file.remove("Rplots.pdf")

  ## plotting distribution of distances between bedpe feet
  distance_between_feet <- c()

  # Use mclapply to parallelize the computation
  distance_between_feet <- mclapply(1:nrow(pairs), compute_distance, pairs = pairs, mc.cores = detectCores())

  # Unlist the result to get a vector
  distance_between_feet <- unlist(distance_between_feet)

  distance_between_feet <- distance_between_feet[!is.na(distance_between_feet)]

  distribution <- data.frame( distance = round(distance_between_feet) ) 
   
  if( nrow( distribution ) < 50 ){
    bins <- 10
  } else if( nrow( distribution ) <= 500 ){
    bins <- 30
  } else {
    bins <- 100
  }

  out_file_2 <- paste( dirname(out_file), "bedpe_size-distribution.pdf", sep = "/" )

  pdf( out_file_2,
       width = 10,
       height = 10,
       useDingbats = FALSE 
      )

  p2 <- ggplot( distribution, aes( x = distance ) ) + 
        geom_histogram( bins = bins, color = "#757575" ) + 
        xlab( "Distance between bedpe feet" ) + 
        ylab( "Count" ) +
        xlim( 0 , ceiling( quantile(distribution$distance)[4] ) )  +
        theme_classic()

  suppressWarnings( print( p2 ) )

  dev.off()

  cat("\nOutput folder created at", out_dir, "\n\n")  

}
