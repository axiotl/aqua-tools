suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

options(scipen = 999)

Args <- commandArgs(trailingOnly=T)

###########################################################################
###########################################################################
###                                                                     ###
###                                FUNCTIONS                            ###
###                                                                     ###
###########################################################################
###########################################################################


# compute distance between bedpe feet
  compute_distance <- function(i, pairs) {
    if (pairs[i, 1] == pairs[i, 4]) {
      mp1 <- as.numeric(pairs[i, 2]) + (as.numeric(pairs[i, 3]) - as.numeric(pairs[i, 2])) / 2
      mp2 <- as.numeric(pairs[i, 5]) + (as.numeric(pairs[i, 6]) - as.numeric(pairs[i, 5])) / 2
      return(abs(mp2 - mp1))
    }
    return(NA)  # return NA for inter-chromosomal pairs
  }

# calculate peak and peak to corner scores
compute_peak_scores <- function(apa, res) {
  apa_m <- as.matrix(apa)
  n     <- nrow(apa_m)
  mid   <- (n + 1) %/% 2
  
  if (res <= 5000) {
    q <- 3
  } else {
    q <- 6
  }
  
  # central pixel
  peak <- as.numeric(apa_m[mid, mid])
  
  # corner regions (q x q blocks)
  corners <- list(
    p2ll = apa_m[(n - q + 1):n,        1:q],
    p2ul = apa_m[1:q,                  1:q],
    p2ur = apa_m[1:q,         (n - q + 1):n],
    p2lr = apa_m[(n - q + 1):n, (n - q + 1):n]
  )
  
  # compute ratios
  scores <- lapply(corners, function(region) peak / mean(region, na.rm = TRUE))
  
  # add central peak value
  scores$peak <- peak
  
  return(scores)
}

# Create peak and corner annotations
build_box_annotations <- function(matrix_norm, res, peak_scores) {
  n   <- nrow(matrix_norm)
  mid <- (n + 1) %/% 2
  q   <- if (res <= 5000) 3 else 6

  region_names <- c("P2LL", "P2UL", "P2UR", "P2LR", "peak")
  scores <- c(
    peak_scores$p2ll,
    peak_scores$p2ul,
    peak_scores$p2ur,
    peak_scores$p2lr,
    peak_scores$peak
  )

  # Conditional label formatting
  label <- ifelse(
    region_names == "peak",
    sprintf("peak: %.2f", scores),
    sprintf("%s:\n%.2f", region_names, scores)
  )

  box_data <- data.frame(
    region = region_names,
    xmin   = c(1, 1, n - q + 1, n - q + 1, mid) - 0.5,
    xmax   = c(q, q, n, n, mid) + 0.5,
    ymin   = n + 1 - c(n, q, q, n, mid) - 0.5,
    ymax   = n + 1 - c(n - q + 1, 1, 1, n - q + 1, mid) + 0.5,
    label  = label
  )

  # Add label positions
  box_data$x_label <- (box_data$xmin + box_data$xmax) / 2
  box_data$y_label <- (box_data$ymin + box_data$ymax) / 2

  # Shift peak label
  box_data$x_label[box_data$region == "peak"] <- box_data$x_label[box_data$region == "peak"] + 3
  box_data$y_label[box_data$region == "peak"] <- box_data$y_label[box_data$region == "peak"] + 1.2

  return(box_data)
}


# Helper function for scale bar format
get_label_fmt <- function(cap) {
  if (cap < 1) {
    scales::label_number(accuracy = 0.01, big.mark = "")
  } else if (cap < 500) {
    scales::label_number(accuracy = 0.1, big.mark = "")
  } else {
    scales::label_number(accuracy = 1, big.mark = "")
  }
}



###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if( length(Args) == 15 ) {

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
  apa_scores     <- Args[15]

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
  breakList <- seq(min_cap_value, cap, by = (cap - min_cap_value) / 5)
  color_ramp_red  <-  colorRampPalette( c( "white", "red") )( 100 )

  # define window sizes and resolution
  win_size <- as.numeric( win_size )
  res      <- as.numeric( res )
  res_div      <- res/1000

  # define row and column names of APA plot
  cols <- c( paste( "-", rev( seq( res_div, res_div*win_size, res_div ) ), "kb", sep="" ), 
             "", 
             paste(           seq( res_div, res_div*win_size, res_div )  , "kb", sep="" ) 
            )
  cols[cols == ""] <- "0"
  
  rownames( matrix_norm) <- cols
  colnames( matrix_norm ) <- cols

  matrix_long <- as.data.frame(matrix_norm) %>%
  mutate(row = seq_len(nrow(matrix_norm))) %>% 
  pivot_longer(-row, names_to = "col", values_to = "value") %>%
  mutate(
    col = factor(col, levels = cols),
    row = factor(row, levels = seq_len(nrow(matrix_norm))) 
  )

  if (apa_scores){
    peak_scores <- compute_peak_scores(matrix_norm, res)
    box_data    <- build_box_annotations(matrix_norm, res, peak_scores)
    cat(sprintf("\nAPA scores for %s\n", sample1))
    cat(sprintf("peak           <- %.3f\n", peak_scores$peak))
    cat(sprintf("P2LL           <- %.3f\n", peak_scores$p2ll))
    cat(sprintf("P2UL           <- %.3f\n", peak_scores$p2ul))
    cat(sprintf("P2UR           <- %.3f\n", peak_scores$p2ur))
    cat(sprintf("P2LR           <- %.3f\n", peak_scores$p2lr))
  }

  p1_gg <- ggplot(matrix_long, aes(x = col, y = row, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = color_ramp_red, breaks = breakList, limits = c(min_cap_value, cap),
                        oob = scales::squish,  # <-- clamp values to limits without altering matrix
                        labels = get_label_fmt(cap),
                        name = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(matrix_long$row)), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.key.height = unit(1.2, "cm")
    )

  # Add APA score annotations if flag is TRUE
  if (apa_scores) {
    p1_gg <- p1_gg +
      geom_rect(data = box_data,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                color = "black", fill = NA, linewidth = 0.6) +
      geom_text(data = box_data,
                aes(x = x_label, y = y_label, label = label),
                inherit.aes = FALSE,
                color = "black", size = 3)
  }

  pdf(out_file, width = 4.5, height = 4.5, useDingbats = FALSE)
  grid.draw(arrangeGrob(
      p1_gg,
      nrow = 1,
      top = paste(sample1, " (n.loop=",  num_loops, ")", sep = "")
  ))
  dev.off()

  distance_between_feet <- c()
  distance_between_feet <- lapply(1:nrow(pairs), compute_distance, pairs = pairs)

  # Unlist the result to get a vector
  distance_between_feet <- unlist(distance_between_feet)
  
  distance_between_feet <- distance_between_feet[!is.na(distance_between_feet)]

  distribution <- data.frame( distance = round(distance_between_feet) ) 

  if( nrow( distribution ) < 50 ){
    bins <- 10
  } else if( nrow( distribution ) <= 1000 ){
    bins <- 30
  } else {
    bins <- 50
  }

  out_file_2 <- paste( dirname(out_file), "bedpe_size-distribution.pdf", sep = "/" )

  pdf( out_file_2,
       width = 7,
       height = 6,
       useDingbats = FALSE 
      )

  p2 <- ggplot( distribution, aes( x = distance ) ) + 
        geom_histogram(bins = bins, color = "#757575" ) + 
        xlab( "Distance between bedpe feet" ) + 
        ylab( "Count" ) +
        xlim( 0 , ceiling( quantile(distribution$distance)[4] ) )  +
        theme_classic()

  suppressWarnings(suppressMessages(print(p2)))

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

if(length(Args) == 18) {

  matrix_file_A    <- Args[1]
  matrix_file_B    <- Args[2]
  sampleA         <- Args[3]
  sampleB         <- Args[4]
  max_cap         <- Args[5]
  max_cap_delta   <- Args[6]
  out_file        <- Args[7]
  win_size        <- Args[8]
  res             <- Args[9]
  pairs           <- Args[10]
  flag_norm       <- Args[11]
  path_mgs_A      <- Args[12]
  path_mgs_B      <- Args[13]
  flag_loop_norm  <- Args[14]
  num_loops       <- as.numeric(Args[15])
  out_dir         <- Args[16]
  min_cap         <- Args[17]
  apa_scores      <- Args[18]


  pairs <- as.data.frame( read.table( pairs, as.is = TRUE ) )
  pairs <- pairs[ , 1:6 ]


  cat("\n")
  cat(sprintf("sample A        <- %s\n", sampleA        ))
  cat(sprintf("sample B        <- %s\n", sampleB        ))
  cat(sprintf("out_file        <- %s\n", out_file       ))
  cat(sprintf("bin_size        <- %s\n", res            )) 


  matrix_raw_A <- read.table(matrix_file_A, header=FALSE, sep="\t")
  matrix_raw_B <- read.table(matrix_file_B, header=FALSE, sep="\t")

  ## CPM or AQuA factors for Sample A
  mergeStats_A  <- read.table(path_mgs_A, as.is = T)
  spikeVar_A    <- ncol(mergeStats_A)
  hg_total1     <- as.numeric(mergeStats_A["valid_interaction_rmdup", 1])
  mm_total1     <- as.numeric(mergeStats_A["valid_interaction_rmdup", 2])
  total1        <- sum(hg_total1, mm_total1)
  norm_factor1  <- 1000000 / total1
  aqua_factor1  <- hg_total1 / mm_total1

  ## CPM or AQuA factors for Sample B
  mergeStats_B  <- read.table(path_mgs_B, as.is = T)
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
  matrix_norm_A <- matrix_raw_A * norm_factor1 * aqua_factor1

  # For Sample B
  matrix_norm_B <- matrix_raw_B * norm_factor2 * aqua_factor2

  if (flag_loop_norm == "TRUE") {
    matrix_norm_A <- matrix_norm_A / num_loops
    matrix_norm_B <- matrix_norm_B / num_loops
  }

  delta <- matrix_norm_B - matrix_norm_A

  write.table(delta, file = paste0(out_dir, "/matrix.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)

  # Set max cap
  if( max_cap == "no_cap" ) { 
    cap <- max(matrix_norm_A, matrix_norm_B)
  } else { 
    cap <- as.numeric(max_cap) 
  }

  # Determine min cap
  if( min_cap == "no_cap" ) { 
    min_cap_value <- 0 
  } else { min_cap_value <- as.numeric(min_cap) }

  if( max_cap_delta == "no_cap" ) { 
    delta_cap <- max( abs(delta) ) 
  } else { delta_cap <- as.numeric(max_cap_delta) }


  ## define colors
  breakList       <- seq(min_cap_value, cap, by = (cap - min_cap_value) / 5)
  breakList_delta <- seq(-delta_cap, delta_cap, by = delta_cap / 3 )

  cat(sprintf("max_cap         <- %.2f\n", cap))
  cat(sprintf("min_cap         <- %.2f\n", min_cap_value))
  cat(sprintf("max_cap_delta   <- %.2f\n", delta_cap))


  color_ramp_red <- colorRampPalette(c(     "white",                      "red"))(100)
  color_ramp_vio <- colorRampPalette(c("dodgerblue", "white", "mediumvioletred"))(100)


  # define window sizes and resolution
  win_size <- as.numeric(win_size)
  res      <- as.numeric(res)
  res_div      <- res/1000


  # define row and column names of APA plot
  cols <- c( paste( "-", rev( seq( res_div, res_div*win_size, res_div ) ), "kb", sep="" ), 
             "", 
             paste(           seq( res_div, res_div*win_size, res_div )  , "kb", sep="" ) 
            )
  cols[cols == ""] <- "0"
  
  rownames(matrix_norm_A) <- cols ; rownames(matrix_norm_B) <- cols ; rownames(delta) <- cols
  colnames(matrix_norm_A) <- cols ; colnames(matrix_norm_B) <- cols ; colnames(delta) <- cols

  matrix_long_A <- as.data.frame(matrix_norm_A) %>%
  mutate(row = seq_len(nrow(matrix_norm_A))) %>% 
  pivot_longer(-row, names_to = "col", values_to = "value") %>%
  mutate(
    col = factor(col, levels = cols),
    row = factor(row, levels = seq_len(nrow(matrix_norm_A))) 
  )

  matrix_long_B <- as.data.frame(matrix_norm_B) %>%
  mutate(row = seq_len(nrow(matrix_norm_B))) %>% 
  pivot_longer(-row, names_to = "col", values_to = "value") %>%
  mutate(
    col = factor(col, levels = cols),
    row = factor(row, levels = seq_len(nrow(matrix_norm_B))) 
  )

  matrix_long_D <- as.data.frame(delta) %>%
  mutate(row = seq_len(nrow(delta))) %>% 
  pivot_longer(-row, names_to = "col", values_to = "value") %>%
  mutate(
    col = factor(col, levels = cols),
    row = factor(row, levels = seq_len(nrow(delta))) 
  )

  if (apa_scores){
    peak_scores_A <- compute_peak_scores(matrix_norm_A, res)
    peak_scores_B <- compute_peak_scores(matrix_norm_B, res)

    box_data_A <- build_box_annotations(matrix_norm_A, res, peak_scores_A)
    box_data_B <- build_box_annotations(matrix_norm_B, res, peak_scores_B)

    cat(sprintf("\nAPA scores for %s\n", sampleA))
    cat(sprintf("peak           <- %.3f\n", peak_scores_A$peak))
    cat(sprintf("P2LL           <- %.3f\n", peak_scores_A$p2ll))
    cat(sprintf("P2UL           <- %.3f\n", peak_scores_A$p2ul))
    cat(sprintf("P2UR           <- %.3f\n", peak_scores_A$p2ur))
    cat(sprintf("P2LR           <- %.3f\n", peak_scores_A$p2lr))

    cat(sprintf("\nAPA scores for %s\n", sampleB))
    cat(sprintf("peak           <- %.3f\n", peak_scores_B$peak))
    cat(sprintf("P2LL           <- %.3f\n", peak_scores_B$p2ll))
    cat(sprintf("P2UL           <- %.3f\n", peak_scores_B$p2ul))
    cat(sprintf("P2UR           <- %.3f\n", peak_scores_B$p2ur))
    cat(sprintf("P2LR           <- %.3f\n", peak_scores_B$p2lr))
  }

  p1_gg <- ggplot(matrix_long_A, aes(x = col, y = row, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = color_ramp_red, breaks = breakList, limits = c(min_cap_value, cap),
                        oob = scales::squish,  # <-- clamp values to limits without altering matrix
                        labels = get_label_fmt(cap),
                        name = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(matrix_long_A$row)), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.key.height = unit(1.2, "cm")
    )

  # Add APA score annotations if flag is TRUE
  if (apa_scores) {
    p1_gg <- p1_gg +
      geom_rect(data = box_data_A,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                color = "black", fill = NA, linewidth = 0.6) +
      geom_text(data = box_data_A,
                aes(x = x_label, y = y_label, label = label),
                inherit.aes = FALSE,
                color = "black", size = 3)
  }

  p2_gg <- ggplot(matrix_long_B, aes(x = col, y = row, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = color_ramp_red, breaks = breakList, limits = c(min_cap_value, cap),
                        oob = scales::squish,  # <-- clamp values to limits without altering matrix
                        labels = get_label_fmt(cap),
                        name = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(matrix_long_B$row)), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.key.height = unit(1.2, "cm")
    )

  # Add APA score annotations if flag is TRUE
  if (apa_scores) {
    p2_gg <- p2_gg +
      geom_rect(data = box_data_B,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                color = "black", fill = NA, linewidth = 0.6) +
      geom_text(data = box_data_B,
                aes(x = x_label, y = y_label, label = label),
                inherit.aes = FALSE,
                color = "black", size = 3)
  }

  p3_gg <- ggplot(matrix_long_D, aes(x = col, y = row, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = color_ramp_vio, breaks = breakList_delta, limits = c(-delta_cap, delta_cap),
                        oob = scales::squish,  # <-- clamp values to limits without altering matrix
                        labels = get_label_fmt(cap),
                        name = NULL) +
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(limits = rev(levels(matrix_long_D$row)), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.key.height = unit(1.2, "cm")
  )

  # Add APA score annotations if flag is TRUE
  if (apa_scores) {

    # Get center box coordinates
    n_delta <- nrow(delta)
    mid     <- (n_delta + 1) %/% 2
    peak_value <- delta[mid, mid]

    # Make a one-row box data frame for the center only
    delta_peak_box <- data.frame(
      xmin   = mid - 0.5,
      xmax   = mid + 0.5,
      ymin   = mid - 0.5,
      ymax   = mid + 0.5,
      x_label = mid + 3,      # nudge right
      y_label = mid + 1.2,    # nudge up
      label  = sprintf("peak: %.2f", peak_value)
    )

    p3_gg <- p3_gg +
      geom_rect(data = delta_peak_box,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                color = "black", fill = NA, linewidth = 0.6) +
      geom_text(data = delta_peak_box,
                aes(x = x_label, y = y_label, label = label),
                inherit.aes = FALSE,
                color = "black", size = 3)
  }

  pdf(out_file, width = 14, height = 4.5, useDingbats = FALSE)
  grid.draw(arrangeGrob(
      p1_gg, p2_gg, p3_gg,
      nrow = 1,
      top = paste(sampleA, ", ", sampleB, ", and delta", " (n.loop=",  num_loops, ")", sep = "")
  ))
  dev.off()
  
  ## plotting distribution of distances between bedpe feet
  distance_between_feet <- c()
  distance_between_feet <- lapply(1:nrow(pairs), compute_distance, pairs = pairs)

  # Unlist the result to get a vector
  distance_between_feet <- unlist(distance_between_feet)

  distance_between_feet <- distance_between_feet[!is.na(distance_between_feet)]

  distribution <- data.frame( distance = round(distance_between_feet) ) 
   
  if( nrow( distribution ) < 50 ){
    bins <- 10
  } else if( nrow( distribution ) <= 1000 ){
    bins <- 30
  } else {
    bins <- 50
  }

  out_file_2 <- paste( dirname(out_file), "bedpe_size-distribution.pdf", sep = "/" )

  pdf( out_file_2,
       width = 10,
       height = 10,
       useDingbats = FALSE 
      )

  p2 <- ggplot( distribution, aes( x = distance ) ) + 
        geom_histogram(bins = bins, color = "#757575" ) + 
        xlab( "Distance between bedpe feet" ) + 
        ylab( "Count" ) +
        xlim( 0 , ceiling( quantile(distribution$distance)[4] ) )  +
        theme_classic()

  suppressWarnings(suppressMessages(print(p2)))

  dev.off()

  cat("\nOutput folder created at", out_dir, "\n\n")  

}
