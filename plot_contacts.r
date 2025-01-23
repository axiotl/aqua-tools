suppressPackageStartupMessages(library(HiCcompare))
suppressPackageStartupMessages(library(strawr))

options(scipen = 999)
options(warn = -1 )

data_dir <- "~/lab-data"

Args <- commandArgs(trailingOnly=T)

analysis_type <- Args[1]
flag_inter    <- Args[2]

#################################################################
##                          Functions                          ##
#################################################################

calculate_w <- function(
    plot_width, 
    plot_height, 
    genomic_range,
    resolution ){
  
  range_length <- genomic_range[2] - genomic_range[1]
  aspect_ratio <- plot_height / plot_width
  
  # calculate maximum number of bins that can fit on the page
  # +1 ensures 1 bin buffer on the right for y-axis profiles
  
  max_bins_h <- floor(range_length / resolution)+1 
  w <- aspect_ratio * plot_width / max_bins_h
  
  return(w)
}


calculate_cap <- function(
    denMat, 
    max_cap, 
    quant_cut, 
    flag_inter = "FALSE") {
  
  # Function to calculate skewness
  calculate_skewness <- function(values) {
    n <- length(values)
    mean_val <- mean(values)
    sd_val <- sd(values)
    return ((n / ((n - 1) * (n - 2))) * sum(((values - mean_val) / sd_val) ^ 3))
  }

  if (all(denMat == 0)) {
    cat(sprintf("\nNo contact values found for this interval at resolution %d\n", bin_size))
    quit()
  }
  
  # Check for 'max_cap' and 'quant_cut' conditions
  if (max_cap != "none") {
    numeric_max_cap <- as.numeric(max_cap)
    
    if (numeric_max_cap < 0) {
      cat("\nMax_cap values cannot be negative. Switching to its absolute value.\n\n")
      return(abs(numeric_max_cap))
    } else {
      return(numeric_max_cap)
    }
  } else if (quant_cut != 1.00) {
    non_zero_values <- denMat[denMat != 0]
    quantile_value <- quantile(non_zero_values, probs = c(quant_cut), na.rm = TRUE)
    calculated_cap <- max(quantile_value)
    
    if (calculated_cap < 0) {
      cat("\nThe calculated cap value based on quant_cut is negative. Switching to its absolute value.\n\n")
      return(abs(calculated_cap))
    } else {
      return(calculated_cap)
    }
  } else {
    # Check for interchromosomal plot for no zero diag
    if (flag_inter == "FALSE") { 
      modified_mat <- zero_diag(denMat, 4)
      non_zero_values <- modified_mat[modified_mat != 0]
    } else {
      non_zero_values <- denMat[denMat != 0]
    }
    
    # Calculate skewness
    skewness <- calculate_skewness(non_zero_values)
    
    # Define thresholds for skewness to apply scaling factors
    skewness_threshold1 <- 1  # For moderate skewness
    skewness_threshold2 <- 2  # For high skewness
    
    # Determine the scaling factor based on skewness
    scaling_factor <- ifelse(abs(skewness) > skewness_threshold2, 3,
                             ifelse(abs(skewness) > skewness_threshold1, 2, 1))
    
    # Calculate the 90th percentile of the absolute values
    high_percentile <- quantile(abs(non_zero_values), probs = 0.90, na.rm = TRUE)
    
    # Apply scaling factor to high percentile
    cap_value <- high_percentile * scaling_factor
    
    return(abs(cap_value))
  }
}

#################################################################
##                          Intra-Chr                          ##
#################################################################

get_diagonals <- function( 
    matrix, 
    num_off_diagonals ){
  
  I <- nrow( matrix )
  
  off_diagonal_matrix <- matrix( nrow = I , ncol = num_off_diagonals )
  
  for( i in 1:num_off_diagonals ){
    for( j in 1:I ){
      if( (j + i) > I ){ next() }
      off_diagonal_matrix[ j , i ] <- matrix[ j , j + i ]
    }
  }
  
  off_diagonal_matrix[ is.na( off_diagonal_matrix ) ] <- 0
  colnames( off_diagonal_matrix ) <- paste( "off_diagonal_", 1:num_off_diagonals, sep = "" )
  return( off_diagonal_matrix )
}

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

draw_scale <- function( ){
  
  if (analysis_type == "single_sample"){
    break_steps <- seq(1, length( breakList ), by = 4)
  } else {
    break_steps <- seq(1, length( breakList), by = 8)
  }
  
  wid <- 1
  i <- 1
  
  for (j in 1:length( break_steps )){
    color <- rownames(C[ order( abs( C[,"breaks"] - breakList[break_steps[j]] ), decreasing = FALSE ), ])[1]
    rect(
      (j-1)*wid, top-(i-1)*wid +2*wid,
      j   *wid, top- i   *wid +2*wid,
      col = color, border = NA )
  }
  
  text( 
    0, 
    top+3*wid, 
    labels = format(round(breakList[1],4)), 
    cex = 0.4, 
    col = "#888888", 
    pos = 4, 
    offset = 0 )
  
  text( 
    j*wid, 
    top+3*wid, 
    labels = format( round(breakList[length( breakList )],4), nsmall = 2 ), 
    cex = 0.4, 
    col = "#888888",
    pos = 2,
    offset = 0 )
  
}

draw_scale_inh <- function(
    steps, 
    breaks){
  
  wid <- 1
  i <- 1
  break_steps <- seq( 1, length( breaks ), by = steps )
  
  for (j in 1:length( break_steps )){
    color <- rownames(
      C[ order( abs( C[,"breaks"] - breaks[break_steps[j]] ),
                decreasing = FALSE ), ])[1]
    rect(
      (j-1)*wid, top-(i-1)*wid +2*wid,
      j   *wid, top- i   *wid +2*wid,
      col = color, border = NA )
  }
  
  text(
    0  ,
    top+3*wid,
    labels = format( round(breaks[1],4), nsmall = 2),
    cex = 0.4,
    col = "#888888",
    pos = 4,
    offset = 0 )
  
  text(
    j*wid,
    top+3*wid,
    labels = format( round(breaks[length( breaks )],4), nsmall = 2 ),
    cex = 0.4,
    col = "#888888",
    pos = 2,
    offset = 0 )
  
}

draw_feature_2 <- function( 
    chr, 
    up, 
    do, 
    txt, 
    col, 
    sh, 
    flag_text, 
    offset ){
  
  text_margin <- 1
  
  if (missing(offset)){
    offset_horizontal <- 8.0
  } else {
    offset_horizontal <- 8.0 * offset
  }
  
  offset_vertical   <- 0.6
  
  str <- ( up - interval_start ) / interval_len * w * (ncol(A)-1)
  end <- ( do - interval_start ) / interval_len * w * (ncol(A)-1)
  
  lines(
    c(    str-sh,     end-sh),
    c(top-str, top-end),
    col = col, lwd = 0.5 )
  
  if ( flag_text ){
    
    txt_coord <- sprintf( "%s:%d-%d", chr, up, do )
    
    text(
      str - sh - offset_horizontal, 
      top - str - text_margin,       
      labels = txt_coord,
      col = col, offset = 0.1, pos = 4, cex = 0.2
    )

    text(
      str - sh - offset_horizontal,                
      top - str - text_margin + offset_vertical, 
      labels = txt,
      col = col, offset = 0.1, pos = 4, cex = 0.28
    )

    # Dotted line connecting text to the start of the line
    segments(
      str - sh - (0.5 * offset_horizontal), # Text horizontal position
      top - str,                            # Same y-coordinate as the line
      str - sh,                             # Line horizontal position
      top - str,                            # Same y-coordinate as the line
      col = "#B2B2B2", lwd = 0.5, lty = "dotted"
    )
  }
}

draw_bed <- function(
    bed_file, 
    color, 
    depth, 
    flag_text) {
  
  offset <- bin_size/2	
  if (!file.exists(bed_file)) {
    cat("draw_bed error: file not found\n")
    stop()
  }
  
  bed <- read.table(bed_file, as.is = TRUE)
  colnames(bed) <- c("chr", "start", "end")
  
  if (nrow(bed) > 0) {
    
    # make sure there's a fourth column for text
    if (ncol(bed) == 3 && flag_text) {
      flag_text <- FALSE
      cat("draw_bed error: no additional columns to print as text\n")}
    
    # subset the bed file to the visualized interval
    bed <- bed[bed[,"chr"] == interval_chr, ]
    bed <- bed[bed[,"start"] < interval_end & bed[,"end"] > interval_start, ]
    
    # get index of bins for each element in bed
    bin_index <- floor(bed[,"start"] / bin_size)
    bed <- bed[order(bin_index),]
    
    # get counts inside each bin
    bin_counts <- table(bin_index)
    
    # find genes that fall within each bin
    for (i in seq_along(bin_counts)) {
      bin_names <- names(bin_counts)[i]
      gene_index <- which(bin_index == bin_names)
      
      # retrieve gene information from bed using index
      for (j in seq_along(gene_index)) {
        b <- bed[gene_index[j],]
        if(flag_text) { btxt <- b[4] } else { btxt <- "" }
        
        # print in place if one gene/bin
        if (length(gene_index) == 1) {	
          draw_feature_2(b$chr, b$start, b$end, btxt, color, depth, flag_text)
        }    
        
        # offset labels for two genes in the same bin
        if (length(gene_index) == 2) {
          if (j %% 2 == 0) {
            # even numbered rows
            b[2] <- b[2] + offset
            draw_feature_2(b$chr, b$start, b$end, btxt, color, depth, flag_text, length(gene_index))
          } else {									        # odd numbered rows
            b[2] <- b[2] - offset
            draw_feature_2(b$chr, b$start, b$end, btxt, color, depth, flag_text, length(gene_index))
          }
        }
      }
      if (length(gene_index) > 2) {
        btxt <- paste0(bed[gene_index,4], collapse = ", ")
        draw_feature_2(b$chr, b$start, b$end, btxt, color, depth, flag_text, length(gene_index))
      } 
    }
  } else {
    cat("draw_bed error: empty bed file\n")
  }
}

draw_contacts <- function( 
    A, 
    C, 
    left_side, 
    bedpe ){
  
  for ( i in 1:nrow(A) ){
    for (j in i:ncol(A)){
      color <- rownames(C[ order( abs( C[,"breaks"] - A[i,j] ), decreasing = FALSE ), ])[1]
      
      rect(
        left_side + (j-1)*w, top-(i-1)*w,
        left_side +  j   *w, top- i   *w,
        col = color, border = NA )
    }
    text(
      (j+1)*w, top-(i-0.5)*w,
      labels = i, cex = 0.03, col = "#dddddd" )
  }
  
  if( bedpe != "FALSE" ){
    
    bedpe_coordinates <- data.frame(
      x1 = c(), y1 = c(),
      x2 = c(), y2 = c() )
    
    for( i in 1:nrow(pairs)){
      
      bedpe_coordinates[i,"x1"] <- which( pairs[i,"start2_bin"] == colnames(A) )
      bedpe_coordinates[i,"y1"] <- which( pairs[i,"end1_bin"  ] == rownames(A) )
      
      bedpe_coordinates[i,"x2"] <- which( pairs[i,"end2_bin"  ] == colnames(A) )
      bedpe_coordinates[i,"y2"] <- which( pairs[i,"start1_bin"] == rownames(A) )
    }
    
    
    for( i in 1:nrow(bedpe_coordinates)){
      
      if( bedpe_coordinates[i,"x1"] != bedpe_coordinates[i,"y2"] &&
          bedpe_coordinates[i,"y1"] != bedpe_coordinates[i,"x2"]){
        
        x1 <- bedpe_coordinates[i,"x1"]
        y1 <- bedpe_coordinates[i,"y1"]
        x2 <- bedpe_coordinates[i,"x2"] - 1
        y2 <- bedpe_coordinates[i,"y2"] - 1
        
        rect(
          (x1-1)*w, top-(y1-1)*w,
          x2   *w, top- y2   *w,
          border = bedpe_color )
        
      } else if( bedpe_coordinates[i,"x1"] == bedpe_coordinates[i,"y2"] &&
                 bedpe_coordinates[i,"y1"] == bedpe_coordinates[i,"x2"] ){
        
        x1 <-       w*(bedpe_coordinates[i,"x1"] - 1)
        y1 <- top - w*(bedpe_coordinates[i,"y2"] - 1)
        x2 <-       w*(bedpe_coordinates[i,"x2"] - 1)
        y2 <- top - w*(bedpe_coordinates[i,"y2"] - 1)
        x3 <-       w*(bedpe_coordinates[i,"x2"] - 1)
        y3 <- top - w*(bedpe_coordinates[i,"y1"] - 1)
        
        polygon(
          c(x1,x2,x3),
          c(y1,y2,y3),
          border = bedpe_color )
      }
      
    }
  }
}

draw_title <- function(
    interval_chr, 
    interval_start, 
    interval_end) {
  
  # Construct the title based on analysis type
  if (analysis_type == "single_sample") {
    sample_name <- basename(sample_dir)
    title <- sprintf(
      "Sample: %s\n\nRange: %s:%d-%d | Genome Build: %s | Resolution: %s",
      sample_name,
      interval_chr,
      interval_start,
      interval_end,
      genome_build,
      bin_size
    )
  } else if (analysis_type == "two_sample") {
    sample_nameA <- basename(sample_dirA)
    sample_nameB <- basename(sample_dirB)
    title <- sprintf(
      "Delta = \n%s - \n%s\n\nRange: %s:%d-%d | Genome Build: %s | Resolution: %s",
      sample_nameB, sample_nameA,
      interval_chr,
      interval_start,
      interval_end,
      genome_build,
      bin_size
    )
  }
  text(-2, top + 2, labels = title, col = "#888888", pos = 4, cex = 0.8)
}

draw_profiles <- function( ){
  
  if(analysis_type == "single_sample"){
    
    pos_profile_gap <- 1.1	
    
    ## diagonal profiles
    
    off_diagonals <- get_diagonals( A, 3 )
    
    off_diagonals_positive <- off_diagonals ; off_diagonals_positive[ off_diagonals_positive < 0 ] <- 0
    
    means_positive  <- apply( off_diagonals_positive, 1, mean ) ; means_positive <- c( means_positive, 0 )
    
    means_positive  <- smooth.spline( means_positive )$y
    means_positive  <- means_positive /   max( means_positive )
    
    for ( k in 1:I ){
      for (j in k:J ){
        
        x_source <-       (k-0.5)*w - pos_profile_gap
        y_source <- top - (k-0.5)*w - pos_profile_gap
        
        x_target <- x_source - means_positive[k] #* profile_scale_factor
        y_target <- y_source - means_positive[k] #* profile_scale_factor				            
        lines( c( x_source , x_target ) , c( y_source , y_target ), col = contact_color, lwd = 0.8 )
        
        if( means_positive[k] == max( means_positive ) ){
          text( x_target - 2*w, y_target - 2*w, labels = format(means_positive[k], nsmall = 3), cex = 0.2, col = "#888888", pos = 4, offset = 0 )
        }
      }
    }
    
    pos_profile_gap       <- 0.25
    
    ## x-axis profiles
    
    x_profile <- c()
    
    for( j in 3:J ){
      
      observations <- A[, j]
      observations[observations < 0] <- 0
      x_profile[(j-2)]  <- mean( observations )
      
    }
    
    x_profile <- smooth.spline(x_profile)$y
    x_profile  <- x_profile /  max( x_profile ) 
    
    for( j in 3:J ){
      
      x_source <-       (j-0.5)*w
      y_source <- top + pos_profile_gap
      x_target <- x_source
      y_target <- y_source + x_profile[(j-2)] #* profile_scale_factor
      
      lines( c( x_source , x_target ) , c( y_source , y_target ), col = contact_color, lwd = 0.8 )
      
      if( x_profile[(j-2)] == max( x_profile ) ){
        text( x_target, y_target + 2*w, labels = format( x_profile[(j-2)], nsmall = 3 ), cex = 0.2, col = "#888888", pos = 4, offset = 0 )
      }
    }
    
    
    ## y-axis profiles
    
    y_profile <- c()
    
    for( i in 1:(I-3) ){
      observations  <- A[i ,]
      observations[observations < 0] <- 0
      y_profile[i]  <- mean( observations )
    }
    
    y_profile <- smooth.spline(y_profile)$y
    y_profile  <- y_profile /  max( y_profile ) 
    
    for( i in 1:(I-3) ){
      x_source <- I*w   +  pos_profile_gap
      y_source <- top   -  (i-0.5)*w
      
      x_target <- x_source + y_profile[i]
      y_target <- y_source  #* profile_scale_factor
      
      lines( c( x_source , x_target ) , c( y_source , y_target ), col = contact_color, lwd = 0.8 )
      
      if( y_profile[i] == max( y_profile ) ){
        
        text( x_target + 2*w, y_target , labels = format( y_profile[i], nsmall = 2 ), cex = 0.2, col = "#888888", pos = 4, offset = 0 )
      }
    }
  } else { 
    pos_profile_gap <- 1.1
    ## diagonal profiles
    
    off_diagonals <- get_diagonals( A, 3 )
    
    max <-  max( abs( min( off_diagonals )),
                 abs( max( off_diagonals )) )
    max <- max / 30
    
    off_diagonals_positive <- off_diagonals ; off_diagonals_positive[ off_diagonals_positive < 0 ] <- 0
    off_diagonals_negative <- off_diagonals ; off_diagonals_negative[ off_diagonals_negative > 0 ] <- 0
    
    means_positive  <- apply(  off_diagonals_positive, 1, mean ) #; means_positive <- c( means_positive, 0 )
    means_negative  <- apply( -off_diagonals_negative, 1, mean ) #; means_negative <- c( means_negative, 0 )
    
    means_positive  <- smooth.spline( means_positive, spar = 0.3 )$y
    check_if_uniq   <- suppressWarnings( unique(means_positive) )
    
    if( check_if_uniq != 0 ){means_positive <- means_positive / max}
    
    means_negative  <- smooth.spline( means_negative, spar = 0.3 )$y
    check_if_uniq   <- suppressWarnings( unique(means_negative) )
    
    if( check_if_uniq != 0 ){means_negative <- means_negative / max}
    
    
    profile_scale_factor <- 0.5
    
    
    for ( k in 1:I ){
      for (j in k:J ){
        
        if( length(unique(means_positive)) > 1 ){
          
          x_source <-       (k-0.2)*w - pos_profile_gap
          y_source <- top - (k-0.2)*w - pos_profile_gap
          
          x_target <- x_source - means_positive[k] * profile_scale_factor
          y_target <- y_source - means_positive[k] * profile_scale_factor
          
          lines( c( x_source , x_target ) , c( y_source , y_target ), col = color_pos, lwd = 0.8 )
          
          if( means_positive[k] == max( means_positive ) ){
            text( x_target - 2*w, y_target - 2*w, labels = format(means_positive[k], nsmall = 3), cex = 0.2, col = "#888888", pos = 4, offset = 0 )
          }
          
        }
        
        if( length(unique(means_negative)) > 1 ){
          
          x_source <-       (k-0.8)*w - pos_profile_gap
          y_source <- top - (k-0.8)*w - pos_profile_gap
          
          x_target <- x_source - means_negative[k] * profile_scale_factor
          y_target <- y_source - means_negative[k] * profile_scale_factor
          
          lines( c( x_source , x_target ) , c( y_source , y_target ), col = color_neg, lwd = 0.8 )
          
          if( means_negative[k] == max( means_negative ) ){
            text( x_target - 2*w, y_target - 2*w, labels = format(means_negative[k], nsmall = 3), cex = 0.2, col = "#888888", pos = 4, offset = 0 )
          }
          
        }
        
      }
    }
    
    pos_profile_gap       <- 0.25
    
    
    A_pos <- zero_diag( A , 2)
    A_pos[ A_pos < 0 ] <- 0
    
    A_neg <- zero_diag( A , 2)
    A_neg[ A_neg > 0 ] <- 0
    
    ## x-axis profiles
    
    x_profile_pos <- apply( A_pos, 2, mean )
    x_profile_pos <- smooth.spline(x_profile_pos, spar = 0.3)$y
    
    check_if_uniq <- suppressWarnings( unique(x_profile_pos) )
    
    if( check_if_uniq != 0 ){x_profile_pos <- x_profile_pos /  max( x_profile_pos ) }
    
    x_profile_neg <- apply( -A_neg, 2, mean )
    x_profile_neg <- smooth.spline(x_profile_neg, spar = 0.3)$y
    
    check_if_uniq <- suppressWarnings( unique(x_profile_neg) )
    
    if( check_if_uniq != 0 ){x_profile_neg <- x_profile_neg / max( x_profile_neg ) }
    
    for( j in 3:J ){
      
      if( length(unique(x_profile_pos)) > 1 ){
        x_source <-       (j-0.2)*w
        y_source <- top + pos_profile_gap
        
        x_target <- x_source
        y_target <- y_source + x_profile_pos[(j-2)] #* (profile_scale_factor*0.5)
        
        lines( c( x_source , x_target ) , c( y_source , y_target ), col = color_pos, lwd = 0.8 )
        
        if( x_profile_pos[(j-2)] == max( x_profile_pos ) ){
          text( x_target, y_target + 2*w, labels = format( x_profile_pos[(j-2)], nsmall = 3 ), cex = 0.2, col = "#888888", pos = 4, offset = 0 )
        }
      }
      
      if( length(unique(x_profile_neg)) > 1 ){
        x_source <-       (j-0.8)*w
        y_source <- top + pos_profile_gap
        
        x_target <- x_source
        y_target <- y_source + abs(x_profile_neg[(j-2)]) #* (profile_scale_factor*0.5)
        
        lines( c( x_source , x_target ) , c( y_source , y_target ), col = color_neg, lwd = 0.8 )
        
        if( x_profile_neg[(j-2)] == max( x_profile_neg ) ){
          text( x_target, y_target + 2*w, labels = format( x_profile_neg[(j-2)], nsmall = 3 ), cex = 0.2, col = "#888888", pos = 4, offset = 0 )
        }
      }
    }
    
    
    ## y-axis profiles
    
    y_profile_pos <- apply( A_pos, 1, mean )
    y_profile_pos <- smooth.spline(y_profile_pos, spar = 0.3)$y
    
    check_if_uniq <- suppressWarnings( unique(y_profile_pos) )
    
    if( check_if_uniq != 0 ){y_profile_pos <- y_profile_pos /   max( y_profile_pos ) }
    
    y_profile_neg <- apply( -A_neg, 1, mean )
    y_profile_neg <- smooth.spline(y_profile_neg, spar = 0.3)$y
    
    check_if_uniq <- suppressWarnings( unique(y_profile_neg) )
    
    if( check_if_uniq != 0 ){y_profile_neg <- y_profile_neg /  max( y_profile_neg ) }
    
    for( i in 1:(I-3) ){
      
      if( length(unique(y_profile_pos)) > 1 ){
        
        x_source <- I*w   +  pos_profile_gap
        y_source <- top   -  (i-0.2)*w
        
        x_target <- x_source + y_profile_pos[i] #* (profile_scale_factor*0.5)
        y_target <- y_source
        
        lines( c( x_source , x_target ) , c( y_source , y_target ), col = color_pos, lwd = 0.8 )
        
        if( y_profile_pos[i] == max( y_profile_pos ) ){
          text( x_target, y_target + 2*w, labels = format( y_profile_pos[i], nsmall = 3 ), cex = 0.2, col = "#888888", pos = 4, offset = 0 )
        }
      }
      
      if( length(unique(y_profile_neg)) > 1 ){
        
        x_source <- I*w   +  pos_profile_gap
        y_source <- top   -  (i-0.8)*w
        
        x_target <- x_source + abs(y_profile_neg[i]) #* (profile_scale_factor*0.5)
        y_target <- y_source
        
        lines( c( x_source , x_target ) , c( y_source , y_target ), col = color_neg, lwd = 0.8 )
        
        if( y_profile_neg[i] == max( y_profile_neg ) ){
          text( x_target, y_target + 2*w, labels = format( y_profile_neg[i], nsmall = 3 ), cex = 0.2, col = "#888888", pos = 4, offset = 0 )
        }
      }
    }
  }
}

#################################################################
##                          Inter-Chr                          ##
#################################################################

calculate_plot_dimensions <- function(interval_start1, interval_end1, 
                                      interval_start2, interval_end2, 
                                      bin_size, 
                                      max_dimension = 100) {
  # Calculate the number of bins for each interval
  bins1 <- ceiling((interval_end1 - interval_start1) / bin_size)
  bins2 <- ceiling((interval_end2 - interval_start2) / bin_size)
  
  # Calculate the aspect ratio
  aspect_ratio <- bins2 / bins1
  
  # Determine plot dimensions
  if (aspect_ratio > 1) {
    # Taller than wide
    plot_height <- max_dimension
    plot_width <- max_dimension / aspect_ratio
  } else {
    # Wider than tall or square
    plot_width <- max_dimension
    plot_height <- max_dimension * aspect_ratio
  }
  
  return(list(width = plot_width, height = plot_height))
}


draw_feature_inter <- function(
    chr, up, do,
    gene_name,
    col,
    flag_text,
    w, top, axis,
    interval_start, interval_len,
    plot_width, plot_height,
    depth) {
  
  # Clip the start and end positions to the interval boundaries
  clipped_start <- max(up, interval_start)
  clipped_end <- min(do, interval_start + interval_len)
  
  genomic_region_label <- sprintf("%s:%d-%d", chr, clipped_start, clipped_end)
  
  calculate_position <- function(pos) {
    relative_position <- (pos - interval_start) / interval_len
    return(relative_position * ifelse(axis == "x", plot_width, plot_height))
  }
  
  str <- calculate_position(clipped_start)
  end <- calculate_position(clipped_end)
  
  x_adjustment <- 2 * (plot_height / 100)
  y_adjustment <- 2 * (plot_height / 100)
  
  if (axis == "x") {
    y_position <- top - depth  # Adjust y-position based on depth for x-axis
    lines(c(str, end), c(y_position, y_position), col = col, lwd = 0.5)
    if (flag_text && gene_name != "") {
      label <- sprintf("%s (%s)", gene_name, genomic_region_label)
      label_x <- end + x_adjustment
      label_y <- y_position + y_adjustment
      text(label_x, label_y, label, col = col, cex = 0.2 * (plot_height / 100), srt = 45, pos = 3)
      
      # Add dotted line
      segments(
        end, y_position,
        end, label_y,
        col = "#B2B2B2", lwd = 0.5, lty = "dotted"
      )
    }
  } else {  # Axis == "y"
    x_position <- depth  # Use depth for x-position on y-axis
    y_start <- top - str
    y_end <- top - end
    lines(c(x_position, x_position), c(y_start, y_end), col = col, lwd = 0.5)
    if (flag_text && gene_name != "") {
      label <- sprintf("%s (%s)", gene_name, genomic_region_label)
      
      # Adjusts position of text label
      text_offset <- x_adjustment * 0.5  # Reduced offset for text
      label_x <- x_position - text_offset 
      label_y <- y_end
      text(label_x, label_y, label, col = col, cex = 0.2 * (plot_height / 100), pos = 2, adj = 1)
      
      # Add dotted line
      dotted_line_end <- x_position - x_adjustment  # Keep the original length for dotted line
      segments(
        x_position, y_end,
        dotted_line_end, y_end,
        col = "#B2B2B2", lwd = 0.5, lty = "dotted"
      )
    }
  }
}

draw_bed_inter <- function(
    bed_file, 
    color, 
    depth,  
    flag_text, 
    w, top, axis, 
    interval_start1, interval_len1, interval_start2, interval_len2,
    plot_width, plot_height) {
  
  if (!file.exists(bed_file)) {
    cat("bed file not found\n")
    stop()
  }

  bed <- read.table(bed_file, as.is = TRUE, header = FALSE)
  
  # Assign column names based on the number of columns
  if (ncol(bed) == 3) {
    colnames(bed) <- c("chr", "start", "end")
  } else if (ncol(bed) >= 4) {
    colnames(bed) <- c("chr", "start", "end", "name")
  } else {
    cat("bed file must have at least 3 columns\n")
    stop()
  }
  
  if (nrow(bed) > 0) {
    # Check if there's a fourth column for text
    if (ncol(bed) == 3 && flag_text) {
      flag_text <- FALSE
      cat("bed file warning: no additional columns to print as text\n")
    }
    
    # Separate beds for each axis
    bed_x <- bed[bed$chr == interval_chr1, ]
    bed_y <- bed[bed$chr == interval_chr2, ]
    
    process_axis <- function(bed_data, axis, interval_start, interval_len) {
      
      # Filter genes within the interval
      bed_data <- bed_data[bed_data$start < interval_start + interval_len & bed_data$end > interval_start, ]
      
      # Assign bins
      bed_data$bin <- floor(bed_data$start / bin_size)
      
      # Count genes per bin
      bin_counts <- table(bed_data$bin)
      
      # Process each bin
      for (bin in names(bin_counts)) {
        bin_genes <- bed_data[bed_data$bin == as.numeric(bin), ]
        
        if (nrow(bin_genes) == 1) {
          b <- bin_genes[1, ]
          gene_name <- if(flag_text && ncol(bed) >= 4) b$name else ""
          draw_feature_inter(b$chr, b$start, b$end, gene_name, color, flag_text, w, top, axis, interval_start, interval_len, plot_width, plot_height, depth)
        } else {
          # For bins with 2 or more genes, apply offsetting
          offset <- bin_size / 2
          for (i in 1:nrow(bin_genes)) {
            b <- bin_genes[i, ]
            gene_name <- if(flag_text && ncol(bed) >= 4) b$name else ""
            if (nrow(bin_genes) == 2) {
              # For exactly 2 genes, offset in opposite directions
              offset_direction <- if(i == 1) -1 else 1
            } else {
              # For 3 or more genes, alternate offsetting
              offset_direction <- if(i %% 2 == 1) -1 else 1
            }
            offset_start <- b$start + offset_direction * offset
            offset_end <- b$end + offset_direction * offset
            draw_feature_inter(b$chr, offset_start, offset_end, gene_name, color, flag_text, w, top, axis, interval_start, interval_len, plot_width, plot_height, depth)
          }
        }
      }
    }
    
    # Process X-axis
    process_axis(bed_x, "x", interval_start1, interval_len1)
    
    # Process Y-axis
    process_axis(bed_y, "y", interval_start2, interval_len2)
    
  } else {
    cat("draw_bed_inter error: empty bed file\n")
  }
}

draw_contacts_inter <- function(
    A, 
    C, 
    left_side, 
    bedpe,
    plot_width,
    plot_height,
    top) {  
  
  # Draw the contact matrix
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A))  {
      color <- rownames(C[order(abs(C[,"breaks"] - A[i,j]), decreasing = FALSE),])[1]
      rect(
        left_side + (j-1)*plot_width/ncol(A), top - (i-1)*plot_height/nrow(A),
        left_side + j*plot_width/ncol(A), top - i*plot_height/nrow(A),
        col = color, border = NA)
    }
  }
  
  if (bedpe != "FALSE") {
    bedpe_coordinates <- data.frame(
      x1 = numeric(nrow(pairs)), y1 = numeric(nrow(pairs)),
      x2 = numeric(nrow(pairs)), y2 = numeric(nrow(pairs))
    )
    
    for (i in 1:nrow(pairs)) {
      bedpe_coordinates$x1[i] <- which(as.character(pairs$start1_bin[i]) == colnames(A))
      bedpe_coordinates$y1[i] <- which(as.character(pairs$start2_bin[i]) == rownames(A))
      bedpe_coordinates$x2[i] <- which(as.character(pairs$end1_bin[i]) == colnames(A))
      bedpe_coordinates$y2[i] <- which(as.character(pairs$end2_bin[i]) == rownames(A))
    }
    
    # Draw bedpe highlights
    for (i in 1:nrow(bedpe_coordinates)) {

      x1 <- bedpe_coordinates$x1[i] - 1 # - 1 to match contact plot bin indices
      y1 <- bedpe_coordinates$y1[i] - 1
      x2 <- bedpe_coordinates$x2[i] - 1
      y2 <- bedpe_coordinates$y2[i] - 1
      
      if (length(x1) > 0 && length(y1) > 0 && length(x2) > 0 && length(y2) > 0) {
        rect(
          left_side + (x1) * plot_width/ncol(A), 
          top - (y1) * plot_height/nrow(A),
          left_side + x2 * plot_width/ncol(A), 
          top - y2 * plot_height/nrow(A),
          border = bedpe_color, col = NA, lwd = 1)
      }
    }
  }
}


draw_title_inter <- function(
    interval_chr1, interval_start1, interval_end1, 
    interval_chr2, interval_start2, interval_end2) {
  if (analysis_type == "single_sample") {
    sample_name <- basename(sample_dir)
    title <- sprintf(
      "Sample: %s\n\nRegions: %s:%d:%d - %s:%d:%d\n\nGenome Build: %s | Resolution: %s",
      sample_name,
      interval_chr1,
      interval_start1,
      interval_end1,
      interval_chr2,
      interval_start2,
      interval_end2,
      genome_build,
      bin_size
    )
  } else if (analysis_type == "two_sample") {
    sample_nameA <- basename(sample_dirA)
    sample_nameB <- basename(sample_dirB)
    title <- sprintf(
      "Samples: %s and %s\n\nRegions: %s:%d:%d - %s:%d:%d\n\nGenome Build: %s | Resolution: %s",
      sample_nameA, sample_nameB,
      interval_chr1,
      interval_start1,
      interval_end1,
      interval_chr2,
      interval_start2,
      interval_end2,
      genome_build,
      bin_size
    )
  }
  
  # Draw the title
  text(-2, top + 2, labels = title, col = "#888888", pos = 4, cex = 0.8)
}

draw_profiles_inter <- function( ) {
  
  # Profiles along the x-axis
  x_profile <- apply(A, 2, mean) # Average contacts for each column (x-axis)
  x_profile <- x_profile / max(x_profile) # Normalization
  
  for (j in 1:ncol(A)) {
    x_source <- (j-0.5)*w
    y_source <- top + 0.5 # Gap above the top of the plot
    x_target <- x_source
    y_target <- y_source + x_profile[j]
    
    lines(c(x_source, x_target), c(y_source, y_target), col = contact_color, lwd = 0.8)
  }
  
  # Profiles along the y-axis
  y_profile <- apply(A, 1, mean) # Average contacts for each row (y-axis)
  y_profile <- y_profile / max(y_profile) # Normalization
  
  for (i in 1:nrow(A)) {
    x_source <- -0.5 # Gap to the left of the plot
    y_source <- top - (i-0.5)*w
    x_target <- x_source - y_profile[i]
    y_target <- y_source
    
    lines(c(x_source, x_target), c(y_source, y_target), col = contact_color, lwd = 0.8)
  }
}

draw_profiles_inter_two_sample <- function() {
  # Normalize and smooth the profiles for both positive and negative parts
  normalize_and_smooth <- function(profile) {
    profile_pos <- pmax(profile, 0)
    profile_neg <- pmax(-profile, 0)
    
    if (max(profile_pos) == 0) {
      profile_pos_smoothed <- rep(0, length(profile_pos))
    } else {
      profile_pos_smoothed <- smooth.spline(profile_pos / max(profile_pos), spar = 0.3)$y
    }
    
    if (max(profile_neg) == 0) {
      profile_neg_smoothed <- rep(0, length(profile_neg))
    } else {
      profile_neg_smoothed <- smooth.spline(profile_neg / max(profile_neg), spar = 0.3)$y
    }
    
    list(positive = profile_pos_smoothed, negative = profile_neg_smoothed)
  }
  
  bar_width <- 0.1 * w  # Width of each bar
  bar_gap <- 0.05 * w   # Additional gap between bars
  
  # Profiles along the x-axis
  x_profile <- apply(A, 2, mean)
  x_profile_smoothed <- normalize_and_smooth(x_profile)
  
  for (j in 1:ncol(A)) {
    x_source_center <- (j-0.5)*w  # Center of the bin
    
    # Calculate positions for positive and negative bars
    x_source_pos <- x_source_center - bar_width / 2 - bar_gap / 2
    x_source_neg <- x_source_center + bar_width / 2 + bar_gap / 2
    
    y_source_pos <- top + 0.5  # Position above the top of the plot for positive
    y_source_neg <- top + 0.5  # Same position for negative, but different x_source
    
    # Positive profile
    x_target_pos <- x_source_pos
    y_target_pos <- y_source_pos + x_profile_smoothed$positive[j]
    lines(c(x_source_pos, x_target_pos), c(y_source_pos, y_target_pos), col = color_pos, lwd = 0.8)
    
    # Negative profile
    x_target_neg <- x_source_neg
    y_target_neg <- y_source_neg + x_profile_smoothed$negative[j]
    lines(c(x_source_neg, x_target_neg), c(y_source_neg, y_target_neg), col = color_neg, lwd = 0.8)
  }
  
  # Profiles along the y-axis
  y_profile <- apply(A, 1, mean)
  y_profile_smoothed <- normalize_and_smooth(y_profile)
  
  for (i in 1:nrow(A)) {
    y_source_center <- top - (i-0.5)*w  # Center of the bin
    
    # Base position to the left of the plot, same for both positive and negative
    x_source_base <- -0.5
    
    # Calculate positions for positive and negative bars
    x_target_pos <- x_source_base - y_profile_smoothed$positive[i]  # Position for positive bar
    x_target_neg <- x_source_base - y_profile_smoothed$negative[i]  # Position for negative bar
    
    # Positive profile
    lines(c(x_source_base, x_target_pos), c(y_source_center, y_source_center), col = color_pos, lwd = 0.8)
    
    # Negative profile
    lines(c(x_source_base, x_target_neg), c(y_source_center + 3*bar_gap, y_source_center + 3*bar_gap), col = color_neg, lwd = 0.8)
  }
}



underplot_message <- "\n\nParameters provided may generate poor visualization.\nConsider increasing interval range or decreasing resolution\nfor improved visualization.\n\n"

overplot_message_1 <-"\n\nParameters provided may generate poor visualization.\nConsider decreasing interval range or increasing resolution\nfor improved visualization.\n\n"

overplot_message_2 <-"\n\nParameters provided may generate poor visualization.\nConsider decreasing interval range for improved visualization.\n\n"

###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if(flag_inter == "FALSE"){
  if( analysis_type == "single_sample" ){
    
    ## Parameters
    
    path_hic           <- Args[3]
    path_mgs           <- Args[4]
    norm               <- Args[5]
    unit               <- Args[6]
    bin_size           <- as.numeric( Args[7] )
    prefix             <- Args[8]
    range_text         <- Args[9]
    out_dir            <- Args[10]
    interval_chr       <- Args[11]
    interval_start     <- as.numeric( Args[12] )
    interval_end       <- as.numeric( Args[13] )
    genome_build       <- Args[14]
    flag_profiles      <- as.logical( Args[15] )
    contact_color      <- sprintf( "#%s", Args[16] )
    ann_default        <- Args[17]
    ann_custom         <- Args[18]
    ann_custom_colors  <- Args[19]
    max_cap            <- Args[20]
    flag_norm          <- Args[21]
    quant_cut          <- as.numeric( Args[22] )
    output_name        <- Args[23]
    bedpe              <- Args[24]
    bedpe_color        <- sprintf( "#%s", Args[25] )
    inherent           <- as.logical(Args[26])
    sample_dir         <- Args[27]
    flag_w             <- Args[28]
    inh_col_floor      <- sprintf( "#%s", Args[29] )
    inh_col_off        <- sprintf( "#%s", Args[30] )
    inh_col_on         <- sprintf( "#%s", Args[31] )
    inh_col_ceil       <- sprintf( "#%s", Args[32] )
    flag_matrix        <- Args[33]

    interval_len <- interval_end - interval_start
    
    if( any( grepl( "chr", interval_chr, ignore.case = T ) ) == FALSE ){
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save="no")
    }
    
    if (max_cap != "none" && quant_cut != 1) {
      cat("\nPlease provide either max_cap or quant_cut, not both!\n")
      quit(save="no")
    }
    
    if (flag_matrix == "FALSE") {
      
      if (bedpe != "FALSE") {
        pairs <- read.table( bedpe, as.is = TRUE)[,1:6]
        colnames(pairs) <- c("chr1","start1","end1","chr2","start2","end2")
        
        if( sum( grepl( "chr", pairs$chr1 ) ) == 0 ) { pairs$chr1 <- paste0( "chr", pairs$chr1 ) }
        if( sum( grepl( "chr", pairs$chr2 ) ) == 0 ) { pairs$chr2 <- paste0( "chr", pairs$chr2 ) }
        
        # subset pairs to interval length
        pairs <- pairs[ pairs[, "chr1"] == interval_chr, ]
        pairs <- pairs[ 
          pairs[,"end2"] < interval_end & 
            pairs[,"start1"] > interval_start & 
            pairs[,"end1"] < interval_end & 
            pairs[,"start2"] > interval_start, ]
        
        if (nrow(pairs) == 0) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          pairs <- pairs[ pairs$chr1 %in% c( paste0( rep("chr",22), 1:22 ), "chrX", "chrY" ) , ]
          
          # convert coordinates to bins
          pairs$start1_bin <- as.integer( floor(   pairs$start1 / bin_size ) * bin_size )
          pairs$end1_bin   <- as.integer( ceiling( pairs$end1   / bin_size ) * bin_size )
          pairs$start2_bin <- as.integer( floor(   pairs$start2 / bin_size ) * bin_size )
          pairs$end2_bin   <- as.integer( ceiling( pairs$end2   / bin_size ) * bin_size )
          
          
          for( i in 1:nrow(pairs) ){
            
            if( pairs[i,"start2_bin"] < pairs[i,"start1_bin"] ){
              a <- pairs[i,"start1_bin"] ; b <- pairs[i,"end1_bin"]
              c <- pairs[i,"start2_bin"] ; d <- pairs[i,"end2_bin"]
              
              pairs[i,"start1_bin"] <- c
              pairs[i,"end1_bin"  ] <- d
              pairs[i,"start2_bin"] <- a
              pairs[i,"end2_bin"  ] <- b
            }
            
            if( pairs[i,"end1_bin"] == pairs[i,"start1_bin"] ){
              pairs[i,"end1_bin"] <- pairs[i,"end1_bin"] + bin_size
            }
            if( pairs[i,"end2_bin"] == pairs[i,"start2_bin"] ){
              pairs[i,"end2_bin"] <- pairs[i,"end2_bin"] + bin_size
            }
            
            if( pairs[i,"chr1"] != pairs[i,"chr2"] ){
              cat("Please provide intra-chromosomal bedpes\n") ; quit(save="no")
            }
          }
        }
      }
    }
    
    cat(      "\n  interval\n")
    cat(sprintf("  interval_chr:   %s\n", interval_chr   ))
    cat(sprintf("  interval_start: %d\n", interval_start ))
    cat(sprintf("  interval_end:   %d\n", interval_end   ))
    
    ## 'chr' prefix handler:
    
    hic_A_chroms <- readHicChroms(path_hic)$name
    
    if( sum( grepl("chr", hic_A_chroms ) )  > 0 ){ flag_strip_chr <- FALSE }
    if( sum( grepl("chr", hic_A_chroms ) ) == 0 ){ flag_strip_chr <- TRUE  }
    
    if( flag_strip_chr ){
      new_interval_chr <- sub( "chr", "", interval_chr )
    } else {
      new_interval_chr <- interval_chr
    }
    
    
    ## CPM or AQuA factors
    
    mergeStats_A <- read.table( path_mgs, as.is = T )
    spikeVar <- ncol(mergeStats_A)
    hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
    mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
    total1       <- sum( hg_total1 , mm_total1 )
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1
    
    if( ! flag_norm %in% c( "blank", "none", "cpm", "aqua" ) ){
      cat("Norm should strictly be none, cpm or aqua in lower case \n")
      q( save = "no" )
    }
    
    # set normalization factor for spike-in and non-spike-in data.
    # non-spike-in data defaults to cpm.
    # spike in data defaults to aqua.
    
    if(spikeVar == 1){
      if(flag_norm == "blank"){
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      } else if(flag_norm == "none"){
        norm_factor1 <- 1
        aqua_factor1 <- 1
      } else if(flag_norm == "cpm"){
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      } else if(flag_norm == "aqua"){
        cat("\n\nError: --norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      }
    }else if(spikeVar == 2){
      if(flag_norm == "blank"){
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
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
    
    
    cat(      "\n  factors\n")
    cat(sprintf("  norm_factor: %f\n", norm_factor1 ))
    cat(sprintf("  aqua_factor: %f\n", aqua_factor1 ))
    
    ## Process matrix:
    
    cat(      "\n  straw parameters\n")
    cat(sprintf("       norm: %s\n", norm     ))
    cat(sprintf("        hic: %s\n", path_hic ))
    cat(sprintf("   interval: %s\n", paste( new_interval_chr, interval_start, interval_end, sep = ":"  ) ))
    cat(sprintf("       unit: %s\n", unit     ))
    cat(sprintf("   bin_size: %s\n", bin_size ))
    
    
    sparMat1 <- straw(
      norm,
      path_hic,
      paste( new_interval_chr, interval_start, interval_end, sep = ":"  ),
      paste( new_interval_chr, interval_start, interval_end, sep = ":"  ),
      unit,
      bin_size )

    if( nrow(sparMat1) == 0 ){
      cat("\n No contact values found for this interval\n")
      quit()
    }

    
    if(isFALSE(inherent)){

      sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1

      denMat1 <- matrix(
          data = 0,
          nrow = length(seq(interval_start, interval_end, bin_size)),
          ncol = length(seq(interval_start, interval_end, bin_size))
      )

      rownames(denMat1) <- seq(interval_start, interval_end, bin_size)
      colnames(denMat1) <- seq(interval_start, interval_end, bin_size)

      # Populate denMat1
      for (i in 1:nrow(sparMat1)) {
          x_index <- as.character(sparMat1[i, "x"])
          y_index <- as.character(sparMat1[i, "y"])
            
          if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
              denMat1[x_index, y_index] <- sparMat1[i, "counts"]
          }
      }
      

      if( ncol(denMat1) >= 5){
        
        off_diagonal_max <- max( zero_diag( denMat1 , 4 ) )
        denMat1[denMat1 > off_diagonal_max] <- off_diagonal_max
        
      } else{
        
        denMat1 <- denMat1
      }


      if(flag_matrix){ 
        sample_name <- basename(sample_dir)
        out_file <- paste0(sample_name, "_matrix.txt")
        write.table(denMat1, file = out_file, row.names = FALSE, col.names = FALSE)
        cat("\n", out_file, " output here:", paste0(out_dir, "/", out_file), "\n")
        quit(save = "no")
      }
      
      # Get cap value
      cap <- calculate_cap(denMat1, max_cap, quant_cut)
      cap <- round(cap, 3)
      
      # Print the cap
      cat(sprintf("   cap: %s\n", cap ))
      
      # Cap the values in the original matrix
      denMat1[denMat1 > cap] <- cap
      denMat1max <- denMat1
      
      color_ramp <- colorRampPalette(c("white", contact_color))(101)
      
      breakList = seq( 0, cap, by = cap/100 )
      
      C <- matrix( ncol = 2, nrow = length(breakList), data = 0)
      rownames( C ) <- color_ramp
      colnames( C ) <- c("rank","breaks")
      C[,1] <- 1:nrow(C)
      C[,2] <- breakList
      
      
      A <- denMat1max
      I <- nrow( A )
      J <- ncol( A )
      
      ## Plot:
      
      pdf(
        output_name,
        width = 8.5, height = 11 )
      
      par(family = "Helvetica")
      par(omi = rep(0.5,4))
      par(mai = rep(0.5,4))
      par(bg  = "#eeeeee")
      
      
      plot( NULL,
            xlim = c(-10,100), ylim = c(1,140),
            xlab = NA,       ylab = NA,
            xaxt = "n",      yaxt = "n",
            bty = "n",       asp = 1 )
      
      
      if(flag_w == "blank"){
        w  <- calculate_w(100, 100, c(interval_start, interval_end), 
                          bin_size)
      }
      
      if (flag_w != "blank"){
        w <- as.numeric(flag_w)
        max_w <- calculate_w(100, 100, c(interval_start, interval_end), 
                             bin_size)
        if (w > max_w){
          cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
        }
      }
      
      if (w > 2 && bin_size >= 5000){cat(underplot_message)}
      if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}
      if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}
      
      cat(sprintf("plotting bin width: %s\n", round(w,3)))
      top        <- 138
      
      draw_title(interval_chr, interval_start, interval_end)
      
      top <- top - 10
      
      draw_scale( )
      
      top <- top - 20
      
      draw_contacts( A , C , left_side = 0, bedpe )
      
      no_def_ann <- FALSE
      no_cus_ann <- FALSE
      
      if(isFALSE(flag_profiles)){
        if( ann_default ){
          
          if( genome_build == "hg19" ){
            
            cat( "\nDrawing default annotations\n" )
            
            bed_tss <- paste( data_dir,"/",genome_build,"/reference/TSS.3000.bed", sep = "" )
            
            draw_bed( bed_tss , "#ef0000", 0.5, TRUE  )
            legend(-2, 125, legend="transcription start sites", col = "#ef0000", lwd = 1, 
                   box.lty = 0, cex =0.8, text.col ="#888888")
          }
          
          if( genome_build == "mm10" ){
            
            cat( "\nDrawing default annotations\n" )
            
            bed_tss <- paste( data_dir,"/",genome_build,"/reference/GENCODE_TSSs_200bp_mm10.bed", sep = "" )
            
            print(bed_tss)
            
            draw_bed( bed_tss , "#ef0000", 0.5, TRUE  )
            legend(-2, 125, legend="transcription start sites", col = "#ef0000", lwd = 1, 
                   box.lty = 0, cex =0.8, text.col ="#888888")
          }
          
          if( genome_build == "hg38" ){
            
            cat( "\nDrawing default annotations\n" )
            
            bed_tss <- paste( data_dir,"/",genome_build,"/reference/GENCODE_TSSs_hg38.bed", sep = "" )
            bed_enh <- paste( data_dir,"/",genome_build,"/reference/ENCODE3_cCRE-enhancers_hg38.bed", sep = "" )
            bed_cpg <- paste( data_dir,"/",genome_build,"/reference/UCSC_CpG-islands_hg38.bed", sep = "" )
            
            draw_bed( bed_tss , "#ef0000", 0.5, TRUE  )
            draw_bed( bed_enh , "#ffd700", 0.7, FALSE )
            draw_bed( bed_cpg , "#008949",   1, FALSE )
            
            legend(-2, 125, legend=c("transcription start sites", "ENCODE 3 cCRE enhancers", "CpG sites"), 
                   col = c("#ef0000", "#ffd700", "#008949"), lwd = 1, box.lty = 0, cex =0.8, text.col ="#888888")
          }
          
        } else {
          no_def_ann <- TRUE
        }
        
        if(ann_custom != "NONE") {
          # Split the annotation files string into vector
          annotation_files <- unlist(strsplit(ann_custom, " "))
          
          # Define different colors and values for each file
          default_colors <- c("#00829d", "#8622b7")  # First and second file colors
          track_depths <- c(1.5, 2.0)  # Different depths for first and second file
          
          # Process each annotation file separately
          for(i in seq_along(annotation_files)) {
              current_file <- annotation_files[i]
              
              if(ann_custom_colors == "NONE") {
                  cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                  draw_bed(current_file, default_colors[i], track_depths[i], FALSE)
                  
              } else if(ann_custom_colors != "NONE") {
                  # Split the colors string into vector if colors were provided
                  colors <- unlist(strsplit(ann_custom_colors, " "))
                  if(length(colors) >= i) {
                      cat(sprintf("\nDrawing custom annotations for file %s\n", current_file))
                      draw_bed(current_file, sprintf("#%s", colors[i]), track_depths[i], FALSE)
                  } else {
                      # Fall back to default color if not enough colors provided
                      cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                      draw_bed(current_file, default_colors[i], track_depths[i], FALSE)
                  }
              }
          }
          
        } else {
            no_cus_ann <- TRUE
        }
      }
      
      if( flag_profiles ){
        
        draw_profiles( )
        
      }
      
      device_off <- dev.off()
      
    }
    
    if(isTRUE(inherent)){
      
      cat(sprintf("   inh_col_floor: %s\n", inh_col_floor ))
      cat(sprintf("   inh_col_off:   %s\n", inh_col_off   ))
      cat(sprintf("   inh_col_on:    %s\n", inh_col_on    ))
      cat(sprintf("   inh_col_ceil:  %s\n", inh_col_ceil  ))
      
      denMat1 <- matrix(
          data = 0,
          nrow = length(seq(interval_start, interval_end, bin_size)),
          ncol = length(seq(interval_start, interval_end, bin_size))
      )

      rownames(denMat1) <- seq(interval_start, interval_end, bin_size)
      colnames(denMat1) <- seq(interval_start, interval_end, bin_size)

      # Populate denMat1
      for (i in 1:nrow(sparMat1)) {
          x_index <- as.character(sparMat1[i, "x"])
          y_index <- as.character(sparMat1[i, "y"])
            
          if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
              denMat1[x_index, y_index] <- sparMat1[i, "counts"]
          }
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
      
      if(nrow(power_laws) < ncol(denMat1)){
        
        diff <- ncol(denMat1) - nrow(power_laws)
        
        power_laws <- rbind(
          power_laws,
          do.call(
            "rbind",
            replicate(diff, power_laws[nrow(power_laws),], simplify = FALSE)) )
      }
      
      background <- power_laws[,"off"]
      foreground <- power_laws[,"on" ]
      
      # Remove background
      tad_matrix_bg_removed <- denMat1
      
      for ( j in 0:(ncol(denMat1) - 1) ){
        for ( i in 1:(nrow(denMat1)-j) ){
          tad_matrix_bg_removed[i  ,i+j  ] <-
            denMat1[i,i+j  ] - background[j+1]
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
      
      if(flag_matrix){ 
        sample_name <- basename(sample_dir)
        out_file <- paste0(sample_name, "_matrix.txt")
        write.table(A, file = out_file, row.names = FALSE, col.names = FALSE)
        cat("\n", out_file, " output here:", paste0(out_dir, "/", out_file), "\n")
        quit(save = "no")
      }
      
      I <- nrow( A )
      J <- ncol( A )
      
      # Max_cap
      if( max_cap != "none" ){
        cat("\n\nParameter --max_cap is not applicable when inherent = TRUE.\nContinuing without --max_cap...\n\n")
      }
      if( quant_cut != 1 ){
        cat("\n\nParameter --quant_cut is not applicable when inherent = TRUE.\nContinuing without --quant_cut...\n\n")
      }
      
      ## Plot:
      
      pdf(
        output_name,
        width = 8.5, height = 11 )
      
      par(family = "Helvetica")
      par(omi = rep(0.5,4))
      par(mai = rep(0.5,4))
      par(bg  = "#eeeeee")
      
      plot( NULL,
            xlim = c(-10,100), ylim = c(1,140),
            xlab = NA,       ylab = NA,
            xaxt = "n",      yaxt = "n",
            bty = "n",       asp = 1 )
      
      
      if(flag_w == "blank"){
        w  <- calculate_w(100, 100, c(interval_start, interval_end), 
                          bin_size)
      }
      
      if (flag_w != "blank"){
        w <- as.numeric(flag_w)
        max_w <- calculate_w(100, 100, c(interval_start, interval_end), 
                             bin_size)
        if (w > max_w){
          cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
        }
      }
      
      if (w > 2 && bin_size >= 5000){cat(underplot_message)}
      if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}
      if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}
      
      cat(sprintf("plotting bin width: %s\n", round(w,3)))
      
      top        <- 138
      
      colors_n       <- 10
      colors_0_1     <- colorRampPalette(c(inh_col_off , inh_col_on  ))( colors_n )
      colors_1_2     <- colorRampPalette(c(inh_col_on  , inh_col_ceil))( colors_n )
      #colors_2_3     <- colorRampPalette(c(inh_col_ceil, inh_col_ceil))( colors_n )
      colors         <- c( inh_col_floor, colors_0_1, colors_1_2 )
      
      breaks <- c(
        (1:colors_n)/colors_n,
        (1:colors_n)/colors_n+1 )
      #(1:colors_n)/colors_n+2)
      breaks <- c( 0, breaks )
      
      
      C <- matrix( ncol = 2, nrow = length(breaks), data = 0)
      rownames( C ) <- colors
      colnames( C ) <- c("rank","breaks")
      C[,1] <- 1:nrow(C)
      C[,2] <- breaks
      
      
      draw_title(interval_chr, interval_start, interval_end)
      top <- top - 10
      draw_scale_inh(steps = 1, breaks)
      top <- top - 20
      
      draw_contacts( A , C , left_side = 0, bedpe )
      
      no_def_ann <- FALSE
      no_cus_ann <- FALSE
      
      if(isFALSE(flag_profiles)){
        if( ann_default ){
          
          if( genome_build == "hg19" ){
            
            cat( "\nDrawing default annotations\n" )
            
            bed_tss <- paste( data_dir,"/",genome_build,"/reference/TSS.3000.bed", sep = "" )
            
            draw_bed( bed_tss , "#ef0000", 0.5, TRUE  )
            legend(-2, 125, legend="transcription start sites", col = "#ef0000", lwd = 1, 
                   box.lty = 0, cex =0.8, text.col ="#888888")
          }
          
          if( genome_build == "mm10" ){
            
            cat( "\nDrawing default annotations\n" )
            
            bed_tss <- paste( data_dir,"/",genome_build,"/reference/GENCODE_TSSs_200bp_mm10.bed", sep = "" )
            
            draw_bed( bed_tss , "#ef0000", 0.5, TRUE  )
            legend(-2, 125, legend="transcription start sites", col = "#ef0000", lwd = 1, 
                   box.lty = 0, cex =0.8, text.col ="#888888")
          }
          
          if( genome_build == "hg38" ){
            
            cat( "\nDrawing default annotations\n" )
            
            bed_tss <- paste( data_dir,"/",genome_build,"/reference/GENCODE_TSSs_hg38.bed", sep = "" )
            bed_enh <- paste( data_dir,"/",genome_build,"/reference/ENCODE3_cCRE-enhancers_hg38.bed", sep = "" )
            bed_cpg <- paste( data_dir,"/",genome_build,"/reference/UCSC_CpG-islands_hg38.bed", sep = "" )
            
            draw_bed( bed_tss , "#ef0000", 0.5, TRUE  )
            draw_bed( bed_enh , "#ffd700", 0.7, FALSE )
            draw_bed( bed_cpg , "#008949",   1, FALSE )
            
            legend(-2, 125, legend=c("transcription start sites", "ENCODE 3 cCRE enhancers", "CpG sites"), 
                   col = c("#ef0000", "#ffd700", "#008949"), lwd = 1, box.lty = 0, cex =0.8, text.col ="#888888")
          }
          
        } else {
          no_def_ann <- TRUE
        }
        
        if(ann_custom != "NONE") {
          # Split the annotation files string into vector
          annotation_files <- unlist(strsplit(ann_custom, " "))
          
          # Define different colors and values for each file
          default_colors <- c("#00829d", "#8622b7")  # First and second file colors
          track_depths <- c(1.5, 2.0)  # Different depths for first and second file
          
          # Process each annotation file separately
          for(i in seq_along(annotation_files)) {
              current_file <- annotation_files[i]
              
              if(ann_custom_colors == "NONE") {
                  cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                  draw_bed(current_file, default_colors[i], track_depths[i], FALSE)
                  
              } else if(ann_custom_colors != "NONE") {
                  # Split the colors string into vector if colors were provided
                  colors <- unlist(strsplit(ann_custom_colors, " "))
                  if(length(colors) >= i) {
                      cat(sprintf("\nDrawing custom annotations for file %s\n", current_file))
                      draw_bed(current_file, sprintf("#%s", colors[i]), track_depths[i], FALSE)
                  } else {
                      # Fall back to default color if not enough colors provided
                      cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                      draw_bed(current_file, default_colors[i], track_depths[i], FALSE)
                  }
              }
          }
          
        } else {
            no_cus_ann <- TRUE
        }
      }
      
      if( flag_profiles ){
        
        draw_profiles( )
        
      }
      
      device_off <- dev.off()
      
    }
  }
  
}

###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if(flag_inter == "FALSE"){
  if( analysis_type == "two_sample" ){
    
    ## Parameters
    
    path_hic_A         <- Args[3]
    path_hic_B         <- Args[4]
    path_mgs_A         <- Args[5]
    path_mgs_B         <- Args[6]
    norm               <- Args[7]
    unit               <- Args[8]
    bin_size           <- as.numeric( Args[9] )
    prefix             <- Args[10]
    range_text         <- Args[11]
    out_dir            <- Args[12]
    interval_chr       <- Args[13]
    interval_start     <- as.numeric( Args[14] )
    interval_end       <- as.numeric( Args[15] )
    genome_build       <- Args[16]
    flag_profiles      <- as.logical( Args[17] )
    contact_color      <- Args[18]
    ann_default        <- Args[19]
    ann_custom         <- Args[20]
    ann_custom_colors  <- Args[21]
    max_cap            <- Args[22]
    flag_norm          <- Args[23]
    quant_cut          <- as.numeric( Args[24] )
    output_name        <- Args[25]
    bedpe              <- Args[26]
    bedpe_color        <- sprintf( "#%s", Args[27] )
    flag_w             <- Args[28]
    sample_dirA        <- Args[29]
    sample_dirB        <- Args[30]
    flag_matrix        <- Args[31]
    inh_col_floor      <- sprintf( "#%s", Args[32] )
    inh_col_off        <- sprintf( "#%s", Args[33] )
    inh_col_on         <- sprintf( "#%s", Args[34] )
    inh_col_ceil       <- sprintf( "#%s", Args[35] )
    inherent           <- as.logical(Args[36])

    interval_len <- interval_end - interval_start
    
    if( any( grepl( "chr", interval_chr, ignore.case = T ) ) == FALSE ){
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit()
    }
    
    if (flag_matrix == "FALSE") {
      
      if (bedpe != "FALSE") {
        pairs <- read.table( bedpe, as.is = TRUE)[,1:6]
        colnames(pairs) <- c("chr1","start1","end1","chr2","start2","end2")
        
        if( sum( grepl( "chr", pairs$chr1 ) ) == 0 ) { pairs$chr1 <- paste0( "chr", pairs$chr1 ) }
        if( sum( grepl( "chr", pairs$chr2 ) ) == 0 ) { pairs$chr2 <- paste0( "chr", pairs$chr2 ) }
        
        pairs <- pairs[ pairs[,"chr1"] == interval_chr , ]
        pairs <- pairs[ pairs[,"end2"]           < interval_end
                        & pairs[,"start1"] > interval_start
                        & pairs[, "end1"]  < interval_end
                        & pairs[,"start2"] > interval_start, ]
        
        if (nrow(pairs) == 0) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          pairs <- pairs[ pairs$chr1 %in% c( paste0( rep("chr",22), 1:22 ), "chrX", "chrY" ) , ]
          
          # convert coordinates to bins
          pairs$start1_bin <- as.integer( floor(   pairs$start1 / bin_size ) * bin_size )
          pairs$end1_bin   <- as.integer( ceiling( pairs$end1   / bin_size ) * bin_size )
          pairs$start2_bin <- as.integer( floor(   pairs$start2 / bin_size ) * bin_size )
          pairs$end2_bin   <- as.integer( ceiling( pairs$end2   / bin_size ) * bin_size )
          
          
          for( i in 1:nrow(pairs) ){
            
            if( pairs[i,"start2_bin"] < pairs[i,"start1_bin"] ){
              a <- pairs[i,"start1_bin"] ; b <- pairs[i,"end1_bin"]
              c <- pairs[i,"start2_bin"] ; d <- pairs[i,"end2_bin"]
              
              pairs[i,"start1_bin"] <- c
              pairs[i,"end1_bin"  ] <- d
              pairs[i,"start2_bin"] <- a
              pairs[i,"end2_bin"  ] <- b
            }
            
            if( pairs[i,"end1_bin"] == pairs[i,"start1_bin"] ){
              pairs[i,"end1_bin"] <- pairs[i,"end1_bin"] + bin_size
            }
            if( pairs[i,"end2_bin"] == pairs[i,"start2_bin"] ){
              pairs[i,"end2_bin"] <- pairs[i,"end2_bin"] + bin_size
            }
            
            if( pairs[i,"chr1"] != pairs[i,"chr2"] ){
              cat("Please provide intra-chromosomal bedpes\n") ; quit(save="no")
            }
          }
        }
      }
    }
    
    cat(      "\n  interval\n")
    cat(sprintf("  interval_chr:   %s\n", interval_chr   ))
    cat(sprintf("  interval_start: %d\n", interval_start ))
    cat(sprintf("  interval_end:   %d\n", interval_end   ))
    
    ## 'chr' prefix handler:
    
    hic_A_chroms <- readHicChroms(path_hic_A)$name
    hic_B_chroms <- readHicChroms(path_hic_B)$name
    
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
    
    
    if( flag_strip_chr == "yes" ){
      interval_chr_A <- sub( "chr", "", interval_chr )
      interval_chr_B <- sub( "chr", "", interval_chr )
    } else if( flag_strip_chr == "B" ){
      interval_chr_A <- interval_chr
      interval_chr_B <- sub( "chr", "", interval_chr )
    } else if( flag_strip_chr == "A" ){
      interval_chr_A <- sub( "chr", "", interval_chr )
      interval_chr_B <- interval_chr
    } else if( flag_strip_chr == "no" ){
      interval_chr_A <- interval_chr
      interval_chr_B <- interval_chr
    }
    
    
    ## CPM or AQuA factors
    
    mergeStats_A <- read.table( path_mgs_A, as.is = T )
    spikeVar_A    <- ncol(mergeStats_A)
    hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
    mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
    total1       <- sum( hg_total1 , mm_total1 )
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1
    
    mergeStats_B <- read.table( path_mgs_B, as.is = T )
    spikeVar_B    <- ncol(mergeStats_B)
    hg_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 1 ])
    mm_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 2 ])
    total2       <- sum( hg_total2 , mm_total2 )
    norm_factor2 <- 1000000 / total2
    aqua_factor2 <- hg_total2 / mm_total2
    
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
      } else if (flag_norm == "none"){
        norm_factor1 <- 1 ; aqua_factor1 <- 1
        norm_factor2 <- 1 ; aqua_factor2 <- 1
      } else if (flag_norm == "cpm"){
        norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
        norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
      } else if (flag_norm == "aqua"){
        cat("\n\nError: --norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
        norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
        norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
      }
    } else if (spikeVar_A == 2 || spikeVar_B == 2){
      if (flag_norm == "blank"){
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- aqua_factor2
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
    
    cat(      "\n  factors\n")
    cat(sprintf("  norm_factor1: %f\n", norm_factor1 ))
    cat(sprintf("  norm_factor2: %f\n", norm_factor2 ))
    cat(sprintf("  aqua_factor1: %f\n", aqua_factor1 ))
    cat(sprintf("  aqua_factor2: %f\n", aqua_factor2 ))
    
    ## Process matrix:
    
    cat(      "\n  straw parameters\n")
    cat(sprintf("       norm: %s\n", norm       ))
    cat(sprintf("      hic_A: %s\n", path_hic_A ))
    cat(sprintf("      hic_B: %s\n", path_hic_B ))
    cat(sprintf("   interval: %s\n", paste( interval_chr, interval_start, interval_end, sep = ":"  ) ))
    cat(sprintf("       unit: %s\n", unit       ))
    cat(sprintf("   bin_size: %s\n", bin_size   ))
    
    
    sparMat1 <- straw( norm,
                       path_hic_A,
                       paste( interval_chr_A, interval_start, interval_end, sep = ":"  ),
                       paste( interval_chr_A, interval_start, interval_end, sep = ":"  ),
                       unit,
                       bin_size)
    
    if( nrow(sparMat1) == 0 ){
      cat("\n No contact values found for this interval in sample A\n")
      quit()
    }
    
    sparMat2 <- straw( norm,
                       path_hic_B,
                       paste( interval_chr_B, interval_start, interval_end, sep = ":"  ),
                       paste( interval_chr_B, interval_start, interval_end, sep = ":"  ),
                       unit,
                       bin_size)
    
    if( nrow(sparMat2) == 0 ){
      cat("\n No contact values found for this interval in sample B\n")
      quit()
    }

    if(isFALSE(inherent)){
      
      sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1
      sparMat2$counts <- sparMat2$counts * norm_factor2 * aqua_factor2

      denMat1 <- matrix(
          data = 0,
          nrow = length(seq(interval_start, interval_end, bin_size)),
          ncol = length(seq(interval_start, interval_end, bin_size))
      )

      rownames(denMat1) <- seq(interval_start, interval_end, bin_size)
      colnames(denMat1) <- seq(interval_start, interval_end, bin_size)

      denMat2 <- matrix(
          data = 0,
          nrow = length(seq(interval_start, interval_end, bin_size)),
          ncol = length(seq(interval_start, interval_end, bin_size))
      )

      rownames(denMat2) <- seq(interval_start, interval_end, bin_size)
      colnames(denMat2) <- seq(interval_start, interval_end, bin_size)

      # Populate denMat1
      for (i in 1:nrow(sparMat1)) {
          x_index <- as.character(sparMat1[i, "x"])
          y_index <- as.character(sparMat1[i, "y"])
            
          if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
              denMat1[x_index, y_index] <- sparMat1[i, "counts"]
          }
      }

      # Populate denMat2
      for (i in 1:nrow(sparMat2)) {
          x_index <- as.character(sparMat2[i, "x"])
          y_index <- as.character(sparMat2[i, "y"])
            
          if (x_index %in% rownames(denMat2) && y_index %in% colnames(denMat2)) {
              denMat2[x_index, y_index] <- sparMat2[i, "counts"]
          }
      }

      
      if (nrow(denMat1) != nrow(denMat2)) {
        intersect <- intersect(rownames(denMat1),rownames(denMat2))
        denMat1 <- denMat1[intersect,intersect]
        denMat2 <- denMat2[intersect,intersect]
      }
      
      denDelta  <- denMat2 - denMat1
      
      if(flag_matrix){ 
        sample_nameA <- basename(sample_dirA)
        sample_nameB <- basename(sample_dirB)
        out_file <- paste0(sample_nameA,"_",sample_nameB,"_matrix.txt")
        write.table(denDelta, file = out_file, row.names = FALSE, col.names = FALSE)
        cat("\n", out_file, " output here:", paste0(out_dir, "/", out_file), "\n")
        quit(save = "no")
      }
      
      if( max_cap != "none" && quant_cut != 1 ){
        cat("\n Please provide either max_cap or quant_cut, not both!\n")
        quit()
      }
      
      # Max_cap and quant_cut handling
      cap <- calculate_cap(denDelta, max_cap, quant_cut)
      cap <- round(cap, 3)
      
      denDelta[denDelta > cap] <- cap
      
      cat(sprintf("        cap: %s\n", cap ))
      
      if( isFALSE( grepl( "-", contact_color ) ) ){
        cat("\nPlease separate the two delta colors with '-' only! \n")
        quit()
      }
      
      color_neg <- sprintf( "#%s", unlist( strsplit(contact_color,"-") )[1] )
      color_pos <- sprintf( "#%s", unlist( strsplit(contact_color,"-") )[2] )
      
      breakList = seq(
        -cap,
        cap,
        by = cap/100)
      
      color_ramp <- colorRampPalette(c(color_neg, "white", color_pos))(length(breakList))
      
      C <- matrix( ncol = 2, nrow = length(breakList), data = 0)
      rownames( C ) <- color_ramp
      colnames( C ) <- c("rank","breaks")
      C[,1] <- 1:nrow(C)
      C[,2] <- breakList
      
      
      A <- denDelta
      I <- nrow( A )
      J <- ncol( A )
      
      ## Plot:
      
      pdf(
        output_name,
        width = 8.5, height = 11,
        onefile = T)
      
      par(family = "Helvetica")
      par(omi = rep(0.5,4))
      par(mai = rep(0.5,4))
      par(bg  = "#eeeeee")
      
      plot( NULL,
            xlim = c(-10,100), ylim = c(1,140),
            xlab = NA,       ylab = NA,
            xaxt = "n",      yaxt = "n",
            bty = "n",       asp = 1 )
      
      
      if(flag_w == "blank"){
        w  <- calculate_w(100, 100, c(interval_start, interval_end), 
                          bin_size)
      }
      
      if (flag_w != "blank"){
        w <- as.numeric(flag_w)
        max_w <- calculate_w(100, 100, c(interval_start, interval_end), 
                             bin_size)
        if (w > max_w){
          cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
        }
      }
      
      if (w > 2 && bin_size >= 5000){cat(underplot_message)}
      if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}        
      if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}
      
      
      cat(sprintf("plotting bin width: %s\n", round(w,3)))
      
      top        <- 138
      
      draw_title(interval_chr, interval_start, interval_end)
      
      top <- top - 10
      
      draw_scale( )
      
      top <- top - 20
      
      draw_contacts( A , C , left_side = 0, bedpe )
      
      no_def_ann <- FALSE
      no_cus_ann <- FALSE
      
      if(isFALSE(flag_profiles)){
        if( ann_default ){
          
          if( genome_build == "hg19" ){
            cat( "\nDrawing default annotations\n" )
            
            bed_tss <- paste( data_dir,"/",genome_build,"/reference/TSS.3000.bed", sep = "" )
            
            draw_bed( bed_tss , "#ef0000", 0.5, TRUE  )
          }
          
          if( genome_build == "mm10" ){
            
            cat( "\nDrawing default annotations\n" )
            
            bed_tss <- paste( data_dir,"/",genome_build,"/reference/GENCODE_TSSs_200bp_mm10.bed", sep = "" )
            
            draw_bed( bed_tss , "#ef0000", 0.5, TRUE  )
            legend(-2, 125, legend="transcription start sites", col = "#ef0000", lwd = 1, 
                   box.lty = 0, cex =0.8, text.col ="#888888")
          }
          
          if( genome_build == "hg38" ){
            
            cat( "\nDrawing default annotations\n" )
            
            bed_tss <- paste( data_dir,"/",genome_build,"/reference/GENCODE_TSSs_hg38.bed", sep = "" )
            bed_enh <- paste( data_dir,"/",genome_build,"/reference/ENCODE3_cCRE-enhancers_hg38.bed", sep = "" )
            bed_cpg <- paste( data_dir,"/",genome_build,"/reference/UCSC_CpG-islands_hg38.bed", sep = "" )
            
            draw_bed( bed_tss , "#ef0000", 0.5, TRUE  )
            draw_bed( bed_enh , "#ffd700", 0.7, FALSE )
            draw_bed( bed_cpg , "#008949",   1, FALSE )
            
            legend(-2, 125, legend=c("transcription start sites", "ENCODE 3 cCRE enhancers", "CpG sites"), 
                   col = c("#ef0000", "#ffd700", "#008949"), lwd = 1, box.lty = 0, cex =0.8, text.col ="#888888")
          }
          
        } else {
          no_def_ann <- TRUE
        }
        
        if(ann_custom != "NONE") {
          # Split the annotation files string into vector
          annotation_files <- unlist(strsplit(ann_custom, " "))
          
          # Define different colors and values for each file
          default_colors <- c("#00829d", "#8622b7")  # First and second file colors
          track_depths <- c(1.5, 2.0)  # Different depths for first and second file
          
          # Process each annotation file separately
          for(i in seq_along(annotation_files)) {
              current_file <- annotation_files[i]
              
              if(ann_custom_colors == "NONE") {
                  cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                  draw_bed(current_file, default_colors[i], track_depths[i], FALSE)
                  
              } else if(ann_custom_colors != "NONE") {
                  # Split the colors string into vector if colors were provided
                  colors <- unlist(strsplit(ann_custom_colors, " "))
                  if(length(colors) >= i) {
                      cat(sprintf("\nDrawing custom annotations for file %s\n", current_file))
                      draw_bed(current_file, sprintf("#%s", colors[i]), track_depths[i], FALSE)
                  } else {
                      # Fall back to default color if not enough colors provided
                      cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                      draw_bed(current_file, default_colors[i], track_depths[i], FALSE)
                  }
              }
          }
          
        } else {
            no_cus_ann <- TRUE
        }
      }
      
      if( flag_profiles ){
        
        draw_profiles( )
      }
      device_off <- dev.off()
    }
    
    if(isTRUE(inherent)){
      
      cat("\n--inherent TRUE is currently only available for one sample tests only.\n")
      
    }
  }
}


###########################################################################
###########################################################################
###                                                                     ###
###                   INTER-CHROMOSOMAL ONE-SAMPLE                      ###
###                                                                     ###
###########################################################################
###########################################################################

if(flag_inter == "TRUE"){
  if( analysis_type == "single_sample" ){
    
    path_hic           <- Args[3]
    path_mgs           <- Args[4]
    norm               <- Args[5]
    unit               <- Args[6]
    bin_size           <- as.numeric( Args[7] )
    prefix             <- Args[8]
    range_text         <- Args[9]
    out_dir            <- Args[10]
    interval_chr       <- Args[11]
    interval_start     <- as.numeric( Args[12] )
    interval_end       <- as.numeric( Args[13] )
    genome_build       <- Args[14]
    flag_profiles      <- as.logical( Args[15] )
    contact_color      <- sprintf( "#%s", Args[16] )
    ann_default        <- Args[17]
    ann_custom         <- Args[18]
    ann_custom_colors  <- Args[19]
    max_cap            <- Args[20]
    flag_norm          <- Args[21]
    quant_cut          <- as.numeric( Args[22] )
    output_name        <- Args[23]
    bedpe              <- Args[24]
    bedpe_color        <- sprintf( "#%s", Args[25] )
    inherent           <- as.logical(Args[26])
    sample_dir         <- Args[27]
    flag_w             <- Args[28]
    inh_col_floor      <- sprintf( "#%s", Args[29] )
    inh_col_off        <- sprintf( "#%s", Args[30] )
    inh_col_on         <- sprintf( "#%s", Args[31] )
    inh_col_ceil       <- sprintf( "#%s", Args[32] )
    flag_matrix        <- Args[33]
    
    # Split the range_text at underscores
    range_components <- unlist(strsplit(range_text, "_"))
    
    # Extract the individual components
    interval_chr1 <- range_components[1]
    interval_start1 <- as.numeric(range_components[2])
    interval_end1 <- as.numeric(range_components[3])
    interval_chr2 <- range_components[4]
    interval_start2 <- as.numeric(range_components[5])
    interval_end2 <- as.numeric(range_components[6])

    interval_len1 <- interval_end1 - interval_start1
    interval_len2 <- interval_end2 - interval_start2

    dimensions <- calculate_plot_dimensions(interval_start1, interval_end1, 
                                        interval_start2, interval_end2, 
                                        bin_size)

    plot_width <- dimensions$width
    plot_height <- dimensions$height
    
    if( any( grepl( "chr", interval_chr1, ignore.case = T ) ) == FALSE ){
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save="no")
    }
    if( any( grepl( "chr", interval_chr2, ignore.case = T ) ) == FALSE ){
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save="no")
    }
    
    if (max_cap != "none" && quant_cut != 1) {
      cat("\nPlease provide either max_cap or quant_cut, not both!\n")
      quit(save="no")
    }
    
    
    ## 'chr' prefix handler:
    hic_A_chroms <- readHicChroms(path_hic)$name
    
    if( sum( grepl("chr", hic_A_chroms ) )  > 0 ){ flag_strip_chr <- FALSE }
    if( sum( grepl("chr", hic_A_chroms ) ) == 0 ){ flag_strip_chr <- TRUE  }
    
    if(flag_strip_chr){
      new_interval_chr1 <- sub("chr", "", interval_chr1)
      new_interval_chr2 <- sub("chr", "", interval_chr2)
    } else {
      new_interval_chr1 <- interval_chr1
      new_interval_chr2 <- interval_chr2
    }
    
    if (flag_matrix == "FALSE") {
      if (bedpe != "FALSE") {
        pairs <- read.table(bedpe, as.is = TRUE, header = FALSE)[, 1:6]
        colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
        
        # Ensure chromosome names have 'chr' prefix
        pairs$chr1 <- ifelse(grepl("chr", pairs$chr1), pairs$chr1, paste0("chr", pairs$chr1))
        pairs$chr2 <- ifelse(grepl("chr", pairs$chr2), pairs$chr2, paste0("chr", pairs$chr2))
        
        # Filter the BEDPE pairs
        pairs <- pairs[
          ((pairs$chr1 == new_interval_chr1 & pairs$start1 >= interval_start1 & pairs$end1 <= interval_end1 &
              pairs$chr2 == new_interval_chr2 & pairs$start2 >= interval_start2 & pairs$end2 <= interval_end2) |
             (pairs$chr1 == new_interval_chr2 & pairs$start1 >= interval_start2 & pairs$end1 <= interval_end2 &
                pairs$chr2 == new_interval_chr1& pairs$start2 >= interval_start1 & pairs$end2 <= interval_end1)
          ), ]
        
        if (nrow(pairs) == 0) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          # Reorder columns if feet are in swapped position
          for (i in 1:nrow(pairs)) {
            if (pairs[i, "chr1"] != new_interval_chr1) {
              pairs[i, c("chr1", "start1", "end1", "chr2", "start2", "end2")] <- 
                pairs[i, c("chr2", "start2", "end2", "chr1", "start1", "end1")]
            }
          }
          
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
          
          # Convert coordinates to bins
          pairs$start1_bin <- as.integer( floor(pairs$start1 / bin_size) * bin_size )
          pairs$end1_bin   <- as.integer( ceiling(pairs$end1 / bin_size) * bin_size )
          pairs$start2_bin <- as.integer( floor(pairs$start2 / bin_size) * bin_size )
          pairs$end2_bin   <- as.integer( ceiling(pairs$end2 / bin_size) * bin_size )
          
          # Apply adjustments to bin ranges
          for (i in 1:nrow(pairs)) {
            if (pairs[i, "end1_bin"] == pairs[i, "start1_bin"]) {
              pairs[i, "end1_bin"] <- pairs[i, "end1_bin"] + bin_size
            }
            if (pairs[i, "end2_bin"] == pairs[i, "start2_bin"]) {
              pairs[i, "end2_bin"] <- pairs[i, "end2_bin"] + bin_size
            }
          }
        }
      }
    }
    
    
    
    
    ## CPM or AQuA factors
    
    mergeStats_A <- read.table( path_mgs, as.is = T )
    spikeVar <- ncol(mergeStats_A)
    hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
    mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
    total1       <- sum( hg_total1 , mm_total1 )
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1
    
    if( ! flag_norm %in% c( "blank", "none", "cpm", "aqua" ) ){
      cat("Norm should strictly be none, cpm or aqua in lower case \n")
      q( save = "no" )
    }
    
    # set normalization factor for spike-in and non-spike-in data.
    # non-spike-in data defaults to cpm.
    # spike in data defaults to aqua.
    
    if(spikeVar == 1){
      if(flag_norm == "blank"){
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      } else if(flag_norm == "none"){
        norm_factor1 <- 1
        aqua_factor1 <- 1
      } else if(flag_norm == "cpm"){
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      } else if(flag_norm == "aqua"){
        cat("\n\nError: --norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
      }
    }else if(spikeVar == 2){
      if(flag_norm == "blank"){
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
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
    
    cat(      "\n  factors\n")
    cat(sprintf("  norm_factor1: %f\n", norm_factor1 ))
    cat(sprintf("  aqua_factor1: %f\n", aqua_factor1 ))
    
    ## Process matrix:
    
    cat(      "\n  straw parameters\n")
    cat(sprintf("       norm: %s\n", norm     ))
    cat(sprintf("        hic: %s\n", path_hic ))
    cat(sprintf("   interval: %s\n", paste( new_interval_chr1, interval_start1, interval_end1, sep = ":"  ) ))
    cat(sprintf("   interval: %s\n", paste( new_interval_chr2, interval_start2, interval_end2, sep = ":"  ) ))
    cat(sprintf("       unit: %s\n", unit     ))
    cat(sprintf("   bin_size: %s\n", bin_size ))
    
    sparMat1 <- straw(
      norm,
      path_hic,
      paste( new_interval_chr1, interval_start1, interval_end1, sep = ":"  ),
      paste( new_interval_chr2, interval_start2, interval_end2, sep = ":"  ),
      unit,
      bin_size )
    
    
    if( nrow(sparMat1) == 0 ){
      cat("\n No contact values found for this interval\n")
      quit()
    }
    
    sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1
    
    denMat1 <- matrix(
      data = 0,
      nrow = length(seq(interval_start1, interval_end1, bin_size)),
      ncol = length(seq(interval_start2, interval_end2, bin_size))
    )

    rownames(denMat1) <- seq(interval_start1, interval_end1, bin_size)
    colnames(denMat1) <- seq(interval_start2, interval_end2, bin_size)

    # Populate denMat1
    for (i in 1:nrow(sparMat1)) {
      x_index <- as.character(sparMat1[i, "x"])
      y_index <- as.character(sparMat1[i, "y"])
      
      if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
        denMat1[x_index, y_index] <- sparMat1[i, "counts"]
      }
    }

    denMat1    <- t(denMat1)
    
    if(flag_matrix){ 
      sample_name <- basename(sample_dir)
      out_file <- paste0(sample_name, "_matrix.txt")
      write.table(denMat1, file = out_file, row.names = FALSE, col.names = FALSE)
      cat("\n", out_file, " output here:", paste0(out_dir, "/", out_file), "\n")
      quit(save = "no")
    }
    
    cap <- calculate_cap(denMat1, max_cap, quant_cut, flag_inter)
    cap <- round(cap, 3)
    
    # Print the cap
    cat(sprintf("        cap: %s\n", cap ))
    
    # Cap the values in the original matrix
    denMat1[denMat1 > cap] <- cap
    denMat1max <- denMat1
    
    color_ramp <- colorRampPalette(c("white", contact_color))(101)
    
    breakList = seq( 0, cap, by = cap/100 )
    
    C <- matrix( ncol = 2, nrow = length(breakList), data = 0)
    rownames( C ) <- color_ramp
    colnames( C ) <- c("rank","breaks")
    C[,1] <- 1:nrow(C)
    C[,2] <- breakList
    
    A <- denMat1max
    I <- nrow( A )
    J <- ncol( A )
    
    ## Plot:
    
    pdf(
      output_name,
      width = 8.5, height = 11 )
    
    par(family = "Helvetica")
    par(omi = rep(0.5,4))
    par(mai = rep(0.5,4))
    par(bg  = "#eeeeee")
    
    
    plot( NULL,
          xlim = c(-10,100), ylim = c(1,140),
          xlab = NA,       ylab = NA,
          xaxt = "n",      yaxt = "n",
          bty = "n",       asp = 1 )
    
    
    if(flag_w == "blank"){
      
      # Call calculate_w for each interval
      w_interval1 <- calculate_w(plot_width, plot_height, c(interval_start1, interval_end1), bin_size)
      w_interval2 <- calculate_w(plot_width, plot_height, c(interval_start2, interval_end2), bin_size)

      
      # Choose the smaller of the two
      w <- min(w_interval1, w_interval2)
    }
    
    if (flag_w != "blank"){
      w <- as.numeric(flag_w)
      w_interval1 <- calculate_w(plot_width, plot_height, c(interval_start1, interval_end1), bin_size)
      w_interval2 <- calculate_w(plot_width, plot_height, c(interval_start2, interval_end2), bin_size)
      max_w <- min(w_interval1, w_interval2)
      if (w > max_w){
        cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
      }
    }
    
    if (w > 2 && bin_size >= 5000){cat(underplot_message)}
    if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}
    if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}
    
    cat(sprintf("plotting bin width: %s\n", round(w,3)))
    top        <- 135
    
    draw_title_inter(new_interval_chr1, interval_start1, interval_end1, new_interval_chr2, interval_start2, interval_end2)
    
    top <- top - 12
    
    draw_scale( )
    
    top <- top - 13

    draw_contacts_inter(A, C, left_side = 0, bedpe, plot_width, plot_height, top)
    
    if(isFALSE(flag_profiles)){
      if(ann_default) {
        if(genome_build == "hg19") {
          cat("\nDrawing default annotations for hg19\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/TSS.3000.bed", sep = "")
          
          # Draw annotations for chromosome on the x-axis
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
          
          # Draw annotations for chromosome on the y-axis
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "y", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
        }
        
        if(genome_build == "hg38") {
          cat("\nDrawing default annotations for hg38\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/GENCODE_TSSs_hg38.bed", sep = "")
          bed_enh <- paste(data_dir, "/", genome_build, "/reference/ENCODE3_cCRE-enhancers_hg38.bed", sep = "")
          bed_cpg <- paste(data_dir, "/", genome_build, "/reference/UCSC_CpG-islands_hg38.bed", sep = "")

        # Draw annotations for chromosome on the x-axis
            draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
            draw_bed_inter(bed_enh, "#ffd700", -0.2, FALSE, w, top, "x", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
            draw_bed_inter(bed_cpg, "#008949", -0.3, FALSE, w, top, "x", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)

            # Draw annotations for chromosome on the y-axis
            draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "y", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
            draw_bed_inter(bed_enh, "#ffd700", -0.2, FALSE, w, top, "y", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
            draw_bed_inter(bed_cpg, "#008949", -0.3, FALSE, w, top, "y", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
        }
      }
      
      if(ann_custom != "NONE") {
            # Split the annotation files string into vector
            annotation_files <- unlist(strsplit(ann_custom, " "))
            
            # Define different colors and values for each file
            default_colors <- c("#00829d", "#8622b7")
            track_depths <- c(-0.5, -0.8)  # Different depths for first and second file
            
            # Process each annotation file separately
            for(i in seq_along(annotation_files)) {
                current_file <- annotation_files[i]
                
                if(ann_custom_colors == "NONE") {
                    cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                    color <- default_colors[i]
                } else {
                    # Split the colors string into vector if colors were provided
                    colors <- unlist(strsplit(ann_custom_colors, " "))
                    if(length(colors) >= i) {
                        cat(sprintf("\nDrawing custom annotations for file %s\n", current_file))
                        color <- sprintf("#%s", colors[i])
                    } else {
                        cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                        color <- default_colors[i]
                    }
                }
                
                # Draw annotations for both x and y axes
                draw_bed_inter(current_file, color, track_depths[i], FALSE, w, top, "x", 
                            interval_start1, interval_len1, interval_start2, interval_len2, 
                            plot_width, plot_height)
                draw_bed_inter(current_file, color, track_depths[i], FALSE, w, top, "y",
                            interval_start1, interval_len1, interval_start2, interval_len2, 
                            plot_width, plot_height)
            }
      } else {
            no_cus_ann <- TRUE
      }
    }
    
    if( flag_profiles ){
      draw_profiles_inter()
    }
    
    device_off <- dev.off()
  }
}


###########################################################################
###########################################################################
###                                                                     ###
###                   INTER-CHROMOSOMAL TWO-SAMPLE                      ###
###                                                                     ###
###########################################################################
###########################################################################


if(flag_inter == "TRUE"){
  if( analysis_type == "two_sample" ){
    
    path_hic_A         <- Args[3]
    path_hic_B         <- Args[4]
    path_mgs_A         <- Args[5]
    path_mgs_B         <- Args[6]
    norm               <- Args[7]
    unit               <- Args[8]
    bin_size           <- as.numeric( Args[9] )
    prefix             <- Args[10]
    range_text         <- Args[11]
    out_dir            <- Args[12]
    interval_chr       <- Args[13]
    interval_start     <- as.numeric( Args[14] )
    interval_end       <- as.numeric( Args[15] )
    genome_build       <- Args[16]
    flag_profiles      <- as.logical( Args[17] )
    contact_color      <- Args[18]
    ann_default        <- Args[19]
    ann_custom         <- Args[20]
    ann_custom_colors  <- Args[21]
    max_cap            <- Args[22]
    flag_norm          <- Args[23]
    quant_cut          <- as.numeric( Args[24] )
    output_name        <- Args[25]
    bedpe              <- Args[26]
    bedpe_color        <- sprintf( "#%s", Args[27] )
    flag_w             <- Args[28]
    sample_dirA        <- Args[29]
    sample_dirB        <- Args[30]
    flag_matrix        <- Args[31]
    
    # Split the range_text at underscores
    range_components <- unlist(strsplit(range_text, "_"))
    
    # Extract the individual components
    interval_chr1 <- range_components[1]
    interval_start1 <- as.numeric(range_components[2])
    interval_end1 <- as.numeric(range_components[3])
    interval_chr2 <- range_components[4]
    interval_start2 <- as.numeric(range_components[5])
    interval_end2 <- as.numeric(range_components[6])
    
    
    if( any( grepl( "chr", interval_chr1, ignore.case = T ) ) == FALSE ){
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save="no")
    }
    if( any( grepl( "chr", interval_chr2, ignore.case = T ) ) == FALSE ){
      cat("\n  Please add 'chr' prefix for chromosome\n")
      quit(save="no")
    }
    
    if (max_cap != "none" && quant_cut != 1) {
      cat("\nPlease provide either max_cap or quant_cut, not both!\n")
      quit(save="no")
    }
    
    interval_len1 <- interval_end1 - interval_start1
    interval_len2 <- interval_end2 - interval_start2

    dimensions <- calculate_plot_dimensions(interval_start1, interval_end1, 
                                    interval_start2, interval_end2, 
                                    bin_size)

    plot_width <- dimensions$width
    plot_height <- dimensions$height
    
    ## 'chr' prefix handler:
    
    hic_A_chroms <- readHicChroms(path_hic_A)$name
    hic_B_chroms <- readHicChroms(path_hic_B)$name
    
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
    
    
    if( flag_strip_chr == "yes" ){
      interval_chr_A <- sub( "chr", "", interval_chr1 )
      interval_chr_B <- sub( "chr", "", interval_chr2 )
    } else if( flag_strip_chr == "B" ){
      interval_chr_A <- interval_chr1
      interval_chr_B <- sub( "chr", "", interval_chr2 )
    } else if( flag_strip_chr == "A" ){
      interval_chr_A <- sub( "chr", "", interval_chr1 )
      interval_chr_B <- interval_chr2
    } else if( flag_strip_chr == "no" ){
      interval_chr_A <- interval_chr1
      interval_chr_B <- interval_chr2
    }
    
    if (flag_matrix == "FALSE") {
      if (bedpe != "FALSE") {
        pairs <- read.table(bedpe, as.is = TRUE, header = FALSE)[, 1:6]
        colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
        
        # Ensure chromosome names have 'chr' prefix
        pairs$chr1 <- ifelse(grepl("chr", pairs$chr1), pairs$chr1, paste0("chr", pairs$chr1))
        pairs$chr2 <- ifelse(grepl("chr", pairs$chr2), pairs$chr2, paste0("chr", pairs$chr2))
        
        # Filter the BEDPE pairs
        pairs <- pairs[
          ((pairs$chr1 == interval_chr_A & pairs$start1 >= interval_start1 & pairs$end1 <= interval_end1 &
              pairs$chr2 == interval_chr_B & pairs$start2 >= interval_start2 & pairs$end2 <= interval_end2) |
             (pairs$chr1 == interval_chr_B & pairs$start1 >= interval_start2 & pairs$end1 <= interval_end2 &
                pairs$chr2 == interval_chr_A & pairs$start2 >= interval_start1 & pairs$end2 <= interval_end1)
          ), ]
        
        if (nrow(pairs) == 0) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          
          # Reorder columns if feet are in swapped position
          for (i in 1:nrow(pairs)) {
            if (pairs[i, "chr1"] != interval_chr_A) {
              pairs[i, c("chr1", "start1", "end1", "chr2", "start2", "end2")] <- 
                pairs[i, c("chr2", "start2", "end2", "chr1", "start1", "end1")]
            }
          }
          
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
          
          # Convert coordinates to bins
          pairs$start1_bin <- as.integer( floor(pairs$start1 / bin_size) * bin_size )
          pairs$end1_bin   <- as.integer( ceiling(pairs$end1 / bin_size) * bin_size )
          pairs$start2_bin <- as.integer( floor(pairs$start2 / bin_size) * bin_size )
          pairs$end2_bin   <- as.integer( ceiling(pairs$end2 / bin_size) * bin_size )
          
          # Apply adjustments to bin ranges
          for (i in 1:nrow(pairs)) {
            if (pairs[i, "end1_bin"] == pairs[i, "start1_bin"]) {
              pairs[i, "end1_bin"] <- pairs[i, "end1_bin"] + bin_size
            }
            if (pairs[i, "end2_bin"] == pairs[i, "start2_bin"]) {
              pairs[i, "end2_bin"] <- pairs[i, "end2_bin"] + bin_size
            }
          }
        }
      }
    }
    
    ## CPM or AQuA factors
    
    mergeStats_A <- read.table( path_mgs_A, as.is = T )
    spikeVar_A    <- ncol(mergeStats_A)
    hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
    mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
    total1       <- sum( hg_total1 , mm_total1 )
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1
    
    mergeStats_B <- read.table( path_mgs_B, as.is = T )
    spikeVar_B    <- ncol(mergeStats_B)
    hg_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 1 ])
    mm_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 2 ])
    total2       <- sum( hg_total2 , mm_total2 )
    norm_factor2 <- 1000000 / total2
    aqua_factor2 <- hg_total2 / mm_total2
    
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
      } else if (flag_norm == "none"){
        norm_factor1 <- 1 ; aqua_factor1 <- 1
        norm_factor2 <- 1 ; aqua_factor2 <- 1
      } else if (flag_norm == "cpm"){
        norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
        norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
      } else if (flag_norm == "aqua"){
        cat("\n\nError: --norm cannot be aqua for non-spike-in samples.\nPlease use cpm or none. Continuing with cpm...\n\n")
        norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
        norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
      }
    } else if (spikeVar_A == 2 || spikeVar_B == 2){
      if (flag_norm == "blank"){
        norm_factor1 <- norm_factor1
        aqua_factor1 <- aqua_factor1
        norm_factor2 <- norm_factor2
        aqua_factor2 <- aqua_factor2
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
    
    cat(      "\n  factors\n")
    cat(sprintf("  norm_factor1: %f\n", norm_factor1 ))
    cat(sprintf("  norm_factor2: %f\n", norm_factor2 ))
    cat(sprintf("  aqua_factor1: %f\n", aqua_factor1 ))
    cat(sprintf("  aqua_factor2: %f\n", aqua_factor2 ))
    
    ## Process matrix:
    
    cat(      "\n  straw parameters\n")
    cat(sprintf("       norm: %s\n", norm       ))
    cat(sprintf("      hic_A: %s\n", path_hic_A ))
    cat(sprintf("      hic_B: %s\n", path_hic_B ))
    cat(sprintf("   interval: %s\n", paste( interval_chr_A, interval_start1, interval_end1, sep = ":"  ) ))
    cat(sprintf("   interval: %s\n", paste( interval_chr_B, interval_start2, interval_end2, sep = ":"  ) ))
    cat(sprintf("       unit: %s\n", unit       ))
    cat(sprintf("   bin_size: %s\n", bin_size   ))
    
    
    sparMat1 <- straw( norm,
                       path_hic_A,
                       paste( interval_chr_A, interval_start1, interval_end1, sep = ":"  ),
                       paste( interval_chr_B, interval_start2, interval_end2, sep = ":"  ),
                       unit,
                       bin_size)
    
    if( nrow(sparMat1) == 0 ){
      cat("\n No contact values found for this interval in sample A\n")
      quit()
    }
    
    sparMat2 <- straw( norm,
                       path_hic_B,
                       paste( interval_chr_A, interval_start1, interval_end1, sep = ":"  ),
                       paste( interval_chr_B, interval_start2, interval_end2, sep = ":"  ),
                       unit,
                       bin_size)
    
    if( nrow(sparMat2) == 0 ){
      cat("\n No contact values found for this interval in sample B\n")
      quit()
    }
    
    sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1
    sparMat2$counts <- sparMat2$counts * norm_factor2 * aqua_factor2
    
    denMat1 <- matrix(
      data = 0,
      nrow = length(seq(interval_start1, interval_end1, bin_size)),
      ncol = length(seq(interval_start2, interval_end2, bin_size))
    )
    rownames(denMat1) <- seq(interval_start1, interval_end1, bin_size)
    colnames(denMat1) <- seq(interval_start2, interval_end2, bin_size)
    
    denMat2 <- matrix(
      data = 0,
      nrow = length(seq(interval_start1, interval_end1, bin_size)),
      ncol = length(seq(interval_start2, interval_end2, bin_size))
    )
    rownames(denMat2) <- seq(interval_start1, interval_end1, bin_size)
    colnames(denMat2) <- seq(interval_start2, interval_end2, bin_size)
    
    # Populate denMat1
    for (i in 1:nrow(sparMat1)) {
      x_index <- as.character(sparMat1[i, "x"])
      y_index <- as.character(sparMat1[i, "y"])
      
      if (x_index %in% rownames(denMat1) && y_index %in% colnames(denMat1)) {
        denMat1[x_index, y_index] <- sparMat1[i, "counts"]
      }
    }
    
    # Populate denMat2
    for (i in 1:nrow(sparMat2)) {
      x_index <- as.character(sparMat2[i, "x"])
      y_index <- as.character(sparMat2[i, "y"])
      
      if (x_index %in% rownames(denMat2) && y_index %in% colnames(denMat2)) {
        denMat2[x_index, y_index] <- sparMat2[i, "counts"]
      }
    }
    
    denMat1    <- t(denMat1)
    denMat2    <- t(denMat2)
    
    if (nrow(denMat1) != nrow(denMat2)) {
      intersect <- intersect(rownames(denMat1),rownames(denMat2))
      denMat1 <- denMat1[intersect,intersect]
      denMat2 <- denMat2[intersect,intersect]
    }
    
    denDelta  <- denMat2 - denMat1
    
    if(flag_matrix){ 
      sample_nameA <- basename(sample_dirA)
      sample_nameB <- basename(sample_dirB)
      out_file <- paste0(sample_nameA,"_",sample_nameB,"_matrix.txt")
      write.table(denDelta, file = out_file, row.names = FALSE, col.names = FALSE)
      cat("\n", out_file, " output here:", paste0(out_dir, "/", out_file), "\n")
      quit(save = "no")
    }
    
    
    # Max_cap and quant_cut handling
    cap <- calculate_cap(denDelta, max_cap, quant_cut, flag_inter)
    cap <- round(cap, 3)
    
    cat(sprintf("       cap: %s\n", cap))
    
    if( isFALSE( grepl( "-", contact_color ) ) ){
      cat("\nPlease separate the two delta colors with '-' only! \n")
      quit()
    }
    
    color_neg <- sprintf( "#%s", unlist( strsplit(contact_color,"-") )[1] )
    color_pos <- sprintf( "#%s", unlist( strsplit(contact_color,"-") )[2] )
    
    
    # Create the break list
    breakList <- seq(
      from = -cap,
      to = cap,
      by = cap/100
    )
    
    color_ramp <- colorRampPalette(c(color_neg, "white", color_pos))(length(breakList))
    
    C <- matrix( ncol = 2, nrow = length(breakList), data = 0)
    rownames( C ) <- color_ramp
    colnames( C ) <- c("rank","breaks")
    C[,1] <- 1:nrow(C)
    C[,2] <- breakList
    
    
    A <- denDelta
    I <- nrow( A )
    J <- ncol( A )
    
    ## Plot:
    
    pdf(
      output_name,
      width = 8.5, height = 11 )
    
    par(family = "Helvetica")
    par(omi = rep(0.5,4))
    par(mai = rep(0.5,4))
    par(bg  = "#eeeeee")
    
    
    plot( NULL,
          xlim = c(-10,100), ylim = c(1,140),
          xlab = NA,       ylab = NA,
          xaxt = "n",      yaxt = "n",
          bty = "n",       asp = 1 )
    
    
    if(flag_w == "blank"){
      
      # Call calculate_w for each interval
      w_interval1 <- calculate_w(100, 100, c(interval_start1, interval_end1), bin_size)
      w_interval2 <- calculate_w(100, 100, c(interval_start2, interval_end2), bin_size)
      
      # Choose the smaller of the two
      w <- min(w_interval1, w_interval2)
    }
    
    if (flag_w != "blank"){
      w <- as.numeric(flag_w)
      w_interval1 <- calculate_w(100, 100, c(interval_start1, interval_end1), bin_size)
      w_interval2 <- calculate_w(100, 100, c(interval_start2, interval_end2), bin_size)
      max_w <- min(w_interval1, w_interval2)
      if (w > max_w){
        cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
      }
    }
    
    if (w > 2 && bin_size >= 5000){cat(underplot_message)}
    if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}
    if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}
    
    cat(sprintf("plotting bin width: %s\n", round(w,3)))
    top        <- 135
    
    draw_title_inter(interval_chr_A, interval_start1, interval_end1, interval_chr_B, interval_start2, interval_end2)
    
    top <- top - 12
    
    draw_scale( )
    
    top <- top - 13

    draw_contacts_inter(A, C, left_side = 0, bedpe, plot_width, plot_height, top)
    
    if(isFALSE(flag_profiles)){
      if(ann_default) {
        if(genome_build == "hg19") {
          cat("\nDrawing default annotations for hg19\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/TSS.3000.bed", sep = "")
          
          # Draw annotations for chromosome on the x-axis
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
          
          # Draw annotations for chromosome on the y-axis
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "y", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
        }
        
        if(genome_build == "hg38") {
          cat("\nDrawing default annotations for hg38\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/GENCODE_TSSs_hg38.bed", sep = "")
          bed_enh <- paste(data_dir, "/", genome_build, "/reference/ENCODE3_cCRE-enhancers_hg38.bed", sep = "")
          bed_cpg <- paste(data_dir, "/", genome_build, "/reference/UCSC_CpG-islands_hg38.bed", sep = "")

          # Draw annotations for chromosome on the x-axis
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
          draw_bed_inter(bed_enh, "#ffd700", -0.2, FALSE, w, top, "x", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
          draw_bed_inter(bed_cpg, "#008949", -0.3, FALSE, w, top, "x", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
          
          # Draw annotations for chromosome on the y-axis
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "y", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
          draw_bed_inter(bed_enh, "#ffd700", -0.2, FALSE, w, top, "y", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
          draw_bed_inter(bed_cpg, "#008949", -0.3, FALSE, w, top, "y", interval_start1, interval_len1, interval_start2, interval_len2, plot_width, plot_height)
        }
      }
      
      if(ann_custom != "NONE") {
            # Split the annotation files string into vector
            annotation_files <- unlist(strsplit(ann_custom, " "))
            
            # Define different colors and values for each file
            default_colors <- c("#00829d", "#8622b7")
            track_depths <- c(-0.5, -0.8)  # Different depths for first and second file
            
            # Process each annotation file separately
            for(i in seq_along(annotation_files)) {
                current_file <- annotation_files[i]
                
                if(ann_custom_colors == "NONE") {
                    cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                    color <- default_colors[i]
                } else {
                    # Split the colors string into vector if colors were provided
                    colors <- unlist(strsplit(ann_custom_colors, " "))
                    if(length(colors) >= i) {
                        cat(sprintf("\nDrawing custom annotations for file %s\n", current_file))
                        color <- sprintf("#%s", colors[i])
                    } else {
                        cat(sprintf("\nUsing default color %s for custom annotation file %s\n", default_colors[i], current_file))
                        color <- default_colors[i]
                    }
                }
                
                # Draw annotations for both x and y axes
                draw_bed_inter(current_file, color, track_depths[i], FALSE, w, top, "x", 
                            interval_start1, interval_len1, interval_start2, interval_len2, 
                            plot_width, plot_height)
                draw_bed_inter(current_file, color, track_depths[i], FALSE, w, top, "y",
                            interval_start1, interval_len1, interval_start2, interval_len2, 
                            plot_width, plot_height)
            }
      } else {
            no_cus_ann <- TRUE
      }
    }
    
    if( flag_profiles ){
      draw_profiles_inter_two_sample()
    }
    
    device <- dev.off()
  }
}
