suppressPackageStartupMessages(library(strawr))

options(scipen = 999)
options(help_type = "text")
options(warn = -1 )

source(file.path(Sys.getenv("AQUA_TOOLS_DIR", file.path(Sys.getenv("HOME"), "aqua_tools")), "aqua_functions.r"))

data_dir <- Sys.getenv("LAB_DATA_DIR", file.path(Sys.getenv("HOME"), "lab-data"))

Args <- commandArgs(trailingOnly=T)

analysis_type <- Args[1]


#################################################################
##                          Functions                          ##
#################################################################

calculate_w <- function(
    plot_width,
    plot_height,
    n_bins ){

  aspect_ratio <- plot_height / plot_width
  w <- aspect_ratio * plot_width / n_bins

  return(w)
}

# Dense upper-triangle matrix on the full bin grid, labeled by bin start
sparse_to_grid <- function( sparMat, grid ){

  denMat <- matrix( 0, nrow = length(grid), ncol = length(grid),
                    dimnames = list( as.character(grid), as.character(grid) ) )

  i    <- match( sparMat$x, grid )
  j    <- match( sparMat$y, grid )
  keep <- !is.na(i) & !is.na(j)
  denMat[ cbind( i[keep], j[keep] ) ] <- sparMat$counts[keep]

  return(denMat)
}

calculate_cap <- function(
  denMat,
  max_cap,
  quant_cut,
  flag_inter = FALSE
  ) {
  calculate_skewness <- function(values) {
    n <- length(values)
    if (n < 2) {
      return(0)
    }
    mean_val <- mean(values)
    sd_val <- sd(values)
    if (sd_val == 0) {
      return(0)
    }
    return((n / ((n - 1) * (n - 2))) * sum(((values - mean_val) / sd_val)^3))
  }

  if (all(denMat == 0)) {
    cat(sprintf("\nNo contact values found for this interval at resolution %d\n", bin_size))
    quit(save = "no", status = 1)
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
    if (isFALSE(flag_inter)) {
      modified_mat <- zero_diag(denMat, 4)
      non_zero_values <- modified_mat[modified_mat != 0]
    } else {
      non_zero_values <- denMat[denMat != 0]
    }

    if (length(non_zero_values) == 0) {
      cat("\nNot enough data to plot. Please decrease the resolution or increase the range\n")
      quit(save = "no", status = 1)
    }

    # Calculate skewness
    skewness <- calculate_skewness(non_zero_values)

    # Define thresholds for skewness to apply scaling factors
    skewness_threshold1 <- 1 # For moderate skewness
    skewness_threshold2 <- 2 # For high skewness

    # Determine the scaling factor based on skewness
    scaling_factor <- ifelse(abs(skewness) > skewness_threshold2, 3,
      ifelse(abs(skewness) > skewness_threshold1, 2, 1)
    )

    # Calculate the 90th percentile of the absolute values
    high_percentile <- quantile(abs(non_zero_values), probs = 0.90, na.rm = TRUE)

    # Apply scaling factor to high percentile
    cap_value <- high_percentile * scaling_factor

    return(abs(cap_value))
  }
}

draw_scale <- function(
    inherent = FALSE,
    steps = NULL,
    breaks = NULL ) {
  wid <- 1
  i <- 1

  if (!inherent) {
    # Determine steps and breaks based on analysis_type
    steps <- if (analysis_type == "single_sample") 4 else 8
    breaks <- breakList
  }

  break_steps <- seq(1, length(breaks), by = steps)

  for (j in 1:length(break_steps)) {
    color <- rownames(
      C[order(abs(C[,"breaks"] - breaks[break_steps[j]]), decreasing = FALSE), ])[1]
    rect(
      (j-1)*wid, top-(i-1)*wid +2*wid,
      j   *wid, top- i   *wid +2*wid,
      col = color, border = NA
    )
  }

  text(
    0,
    top + 3*wid,
    labels = format(round(breaks[1], 4), nsmall = 2),
    cex = 0.4,
    col = "#888888",
    pos = 4,
    offset = 0
  )

  text(
    j * wid,
    top + 3*wid,
    labels = format(round(breaks[length(breaks)], 4), nsmall = 2),
    cex = 0.4,
    col = "#888888",
    pos = 2,
    offset = 0
  )
}

draw_title <- function(
    interval_chr,
    interval_start,
    interval_end) {

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

# Arrow under the viewpoint's diagonal bin, tip centered on the bin, plus its
# coordinate label
draw_viewpoint_marker <- function( bin_idx, left_side, label ){

  cex_lab <- 0.45
  text_h  <- strheight( "0", cex = cex_lab )
  arrow_h <- text_h * 1.2
  half_w  <- arrow_h * 0.6
  gap     <- text_h * 1.5   # space between the heatmap and the marker
  gap_lab <- text_h * 0.6   # space between the arrow and its label

  # Highest y a box can reach and still be clear of the drawn cells
  clear_y <- function( x_right ){
    col <- ceiling( (x_right - left_side) / w )
    if( col <= 0 ) return( Inf )
    top - col * w - gap
  }

  # Arrow: centered on the bin, just below it, dropped further if needed
  x_tip <- left_side + (bin_idx - 0.5) * w
  y_tip <- min( top - bin_idx * w - gap, clear_y( x_tip + half_w ) )

  # thin leader from the arrow tip up to the viewpoint bin, so the exact bin is
  # clear even when bins are much narrower than the arrow
  segments( x_tip, y_tip, x_tip, top - bin_idx * w, col = "#000000", lwd = 0.3 )

  polygon( c(x_tip, x_tip - half_w, x_tip + half_w),
           c(y_tip, y_tip - arrow_h, y_tip - arrow_h),
           col = "#aaaaaa", border = "#000000", lwd = 0.3 )

  # Label left of the arrow, level with it
  # xpd = NA lets it extend into the left margin
  text( x_tip - half_w - gap_lab, y_tip - arrow_h / 2, labels = label,
        adj = c(1, 0.5), cex = cex_lab, col = "#888888", xpd = NA )
}

# Min/max labels and a thin axis line to the left of the profile
draw_profile_axis <- function( x, y_lo, y_hi, lab_lo, lab_hi ){

  cex_ax <- 0.2
  min_sep <- strheight( "0", cex = cex_ax ) * 1.6

  if( y_hi > y_lo ){
    segments( x, y_lo, x, y_hi, col = "#aaaaaa", lwd = 0.1 )
  }

  if( y_hi - y_lo < min_sep ){
    mid  <- (y_lo + y_hi) / 2
    y_lo <- mid - min_sep / 2
    y_hi <- mid + min_sep / 2
  }

  x_lab <- x - strwidth( "0", cex = cex_ax )
  text( x_lab, y_lo, labels = lab_lo, adj = c(1, 0.5), col = "#aaaaaa", cex = cex_ax, xpd = NA )
  text( x_lab, y_hi, labels = lab_hi, adj = c(1, 0.5), col = "#aaaaaa", cex = cex_ax, xpd = NA )
}

draw_viewpoints <- function( A, C, left_side, viewpoint_start, analysis_type, height = "blank"){

  viewpoint_bin     <- as.integer( interval_to_bins( viewpoint_start, viewpoint_start + 1, bin_size )$first * bin_size )
  viewpoint_band    <- as.character(c( viewpoint_bin - width*bin_size , viewpoint_bin + width*bin_size ))
  viewpoint_label <- sprintf("%s:%d", viewpoint_chr, viewpoint_start)

  if (!all(viewpoint_band %in% colnames(A))) {
    cat("\n Viewpoint probably too close to range border, try increasing range \n")
    q(status = 1, save = "no")
  }

  viewpoint_bin_idx <- which(viewpoint_bin == colnames(A))
  v4C_slice        <- as.data.frame(A[ , which(colnames(A) == viewpoint_band[1]):(which(colnames(A) == viewpoint_band[2])) ])
  v4C_slice        <- as.data.frame(apply( v4C_slice, 1, mean )) ; colnames(v4C_slice) <- "contact_freq"
  v4C_slice$start  <- as.numeric(rownames(v4C_slice))

  v4C_slice_spline <- as.data.frame(cbind(v4C_slice$start, v4C_slice$contact_freq))

  y_pos <- -8

  # Define max vertical space for histogram
  max_plot_height <- ifelse(analysis_type == "single_sample", 3, 4)


  # Determine height_scale, either user-supplied or calculated
  if (height == "blank") {
    # Scale bars against the heatmap's color scale
    scale_ref    <- if (isTRUE(inherent)) 2 else cap
    height_scale <- max_plot_height / scale_ref
    cat(sprintf("  height: %.2f\n", height_scale))
  } else {
    height_scale <- as.numeric(height)
    cat(sprintf("  height: %.2f\n", height_scale))
  }

  if( analysis_type == "single_sample" ){
    for ( i in 1:nrow(A) ){
      for (j in i:ncol(A)){
        color <- rownames(C[ order( abs( C[,"breaks"] - A[i,j] ), decreasing = FALSE ), ])[1]

        rect(
          left_side + (j-1)*w, top-(i-1)*w,
          left_side +  j   *w, top- i   *w,
          col = color, border = NA
        )

        if( i == viewpoint_bin_idx && j == viewpoint_bin_idx ){
          draw_viewpoint_marker( viewpoint_bin_idx, left_side, viewpoint_label )
        }

        if( i == viewpoint_bin_idx && j == ncol(A) ){
          side <- left_side
          for( k in 1:length(v4C_slice_spline$V2)){
            x0 <- side + (k-1)*w
            x1 <- (side + k*w)-0.1
            y0 <- top - y_pos

            bar_h <- v4C_slice_spline$V2[k] * height_scale
            y1 <- y0 + min( max(bar_h, 0.1), max_plot_height )

            if (inherent) {
              value <- v4C_slice_spline$V2[k]
              color <- rownames(C[order(abs(C[, "breaks"] - value), decreasing = FALSE), ])[1]
            } else {
              color <- "red"
            }

            rect(x0, y0, x1, y1, border = NA, col = color)

            if( k == 1 ){
              draw_profile_axis(
                x0 - 1,
                y0 + min( max( min(v4C_slice_spline$V2) * height_scale, 0.1 ), max_plot_height ),
                y0 + min( max( max(v4C_slice_spline$V2) * height_scale, 0.1 ), max_plot_height ),
                round(min(v4C_slice_spline$V2), 2),
                round(max(v4C_slice_spline$V2), 2)
              )
            }

            if( k == viewpoint_bin_idx){
              points( x0+(x1-x0)/2, y1+1, cex = 0.3, pch = 6, bg = "#aaaaaa", col = "#000000" )
              points( x0+(x1-x0)/2, y0-1, cex = 0.3, pch = 2, bg = "#aaaaaa", col = "#000000" )
            }
          }
        }
      }
    }
  }

  if( analysis_type == "two_sample" ){
    for ( i in 1:nrow(A) ){
      for (j in i:ncol(A)){
        color <- rownames(C[ order( abs( C[,"breaks"] - A[i,j] ), decreasing = FALSE ), ])[1]

        rect(
          left_side + (j-1)*w, top-(i-1)*w,
          left_side +  j   *w, top- i   *w,
          col = color, border = NA
        )

        if( i == viewpoint_bin_idx && j == viewpoint_bin_idx ){
          draw_viewpoint_marker( viewpoint_bin_idx, left_side, viewpoint_label )
        }

        if( i == viewpoint_bin_idx && j == ncol(A) ){
          side <- left_side
          coordinates <- matrix(nrow=nrow(A), ncol=4) ; colnames(coordinates) <- c("x0","y0","x1","y1")

          for (k in 1:length(v4C_slice_spline$V2)) {
            x0 <- side + (k-1) * w
            x1 <- (side + k * w)-0.1
            y <- top - y_pos

            adjusted_value <- v4C_slice_spline$V2[k]
            bar_h <- max( min( adjusted_value * height_scale, max_plot_height ), -max_plot_height )

            if (adjusted_value > 0) {
              y0 <- y
              y1 <- y + bar_h
              color <- "mediumvioletred"
            } else {
              y0 <- y + bar_h
              y1 <- y
              color <- "dodgerblue"
            }

            rect(x0, y0, x1, y1, border = NA, col = color)

            coordinates[k,"x0"] <- x0
            coordinates[k,"y0"] <- y0
            coordinates[k,"x1"] <- x1
            coordinates[k,"y1"] <- y1

            if( k == viewpoint_bin_idx){
              points( x0+(x1-x0)/2, y1+1, cex = 0.3, pch = 6, bg = "#aaaaaa", col = "#000000" )
              points( x0+(x1-x0)/2, y0-1, cex = 0.3, pch = 2, bg = "#aaaaaa", col = "#000000" )
            }
          }
draw_profile_axis(
            min(coordinates[,"x0"]) - 1,
            min(coordinates[,"y0"]),
            max(coordinates[,"y1"]),
            round(min(v4C_slice_spline$V2), 2),
            round(max(v4C_slice_spline$V2), 2)
          )
        }
      }
    }
  }
}


underplot_message <- "\n\nParameters provided may generate poor visualization.\nConsider increasing interval range or decreasing resolution.\n\n"
overplot_message_1 <-"\n\nParameters provided may generate poor visualization.\nConsider decreasing interval range or increasing resolution.\n\n"
overplot_message_2 <-"\n\nParameters provided may generate poor visualization.\nConsider decreasing interval range.\n\n"


###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if( analysis_type == "single_sample" ){

  path_hic        <- Args[2]
  path_mgs        <- Args[3]
  norm_method     <- Args[4]
  unit            <- Args[5]
  bin_size        <- as.numeric(Args[6])
  output_name     <- Args[7]
  genome_build    <- Args[8]
  flag_norm       <- Args[9]
  interval_chr    <- Args[10]
  interval_start  <- as.numeric(Args[11])
  interval_end    <- as.numeric(Args[12])
  quant_cut       <- as.numeric(Args[13])
  max_cap         <- Args[14]
  width           <- as.numeric(Args[15])
  height          <- Args[16]
  viewpoint_chr   <- Args[17]
  viewpoint_start <- as.numeric(Args[18])
  sample_dir      <- Args[19]
  inherent        <- as.logical(Args[20])
  inh_col_floor   <- sprintf( "#%s", Args[21] )
  inh_col_off     <- sprintf( "#%s", Args[22] )
  inh_col_on      <- sprintf( "#%s", Args[23] )
  inh_col_ceil    <- sprintf( "#%s", Args[24] )

  if( viewpoint_chr != interval_chr ){
    cat("Viewpoint locus cannot be out of range\n")
    q(status = 1)
  }

  cat("Parameters\n")


  cat("\n")
  cat(sprintf("hic file            %s\n", path_hic          ))
  cat(sprintf("bin size            %s\n", bin_size          ))
  cat(sprintf("v4C_viewpoint_chr   %s\n", viewpoint_chr   ))
  cat(sprintf("v4C_viewpoint_start %d\n", viewpoint_start ))


  hic_A_chroms <- readHicChroms(path_hic)$name

  if( sum( grepl("chr", hic_A_chroms ) )  > 0 ){ flag_strip_chr <- FALSE }
  if( sum( grepl("chr", hic_A_chroms ) ) == 0 ){ flag_strip_chr <- TRUE  }

  if( flag_strip_chr ){
    interval_chr <- sub( "chr", "", interval_chr )
  } else {
    interval_chr <- interval_chr
  }

  mergeStats_A <- read.table( path_mgs, as.is = T )
  spikeVar     <- ncol(mergeStats_A)
  hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
  mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
  total1       <- sum( hg_total1 , mm_total1 )
  norm_factor1 <- 1000000 / total1
  aqua_factor1 <- hg_total1 / mm_total1

  if( ! flag_norm %in% c("blank", "none", "cpm", "aqua" ) ){
    cat("\nNorm should strictly be none, cpm or aqua in lower case \n")
    q(status = 1)
  }

  # set normalization factor for spike and nonspike data
  # non-spike-in samples default to cpm
  # spike-in samples default to aqua

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
    } else if(flag_norm == "aqua"){
      norm_factor1 <- norm_factor1
      aqua_factor1 <- aqua_factor1
    }
  }

  if(isTRUE(inherent)){
    flag_norm <- "inherent"
  }

  cat(sprintf("normalization       %s\n", flag_norm         ))
  cat(       "\nfactors\n")
  cat(sprintf("  norm_factor:  %f\n", norm_factor1 ))
  cat(sprintf("  aqua_factor:  %f\n", aqua_factor1 ))


  b_range     <- interval_to_bins(interval_start, interval_end, bin_size)
  range_grid  <- bins_to_grid(b_range, bin_size)
  range_query <- bins_to_query(interval_chr, b_range, bin_size)

  sparMat1 <- straw( norm_method,
                     path_hic,
                     range_query,
                     range_query,
                     unit,
                     bin_size)

  if( nrow(sparMat1) == 0 ){
    cat("\nNo contact values found for this interval\n")
    q(status = 1)
  }

  if(isFALSE(inherent)){
    sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1

    denMat1 <- sparse_to_grid(sparMat1, range_grid)

    if( max_cap != "none" && quant_cut != 1 ){
      cat("\nPlease provide either max_cap or quant_cut, not both!\n")
      q(status = 1)
    }

    cap <- calculate_cap(denMat1, max_cap, quant_cut)
    cap <- round(cap, 3)

    cat(sprintf("  cap: %s\n", cap ))

    # Cap the values in the original matrix
    denMat1[denMat1 > cap] <- cap
    denMat1max <- denMat1

    # make symmetrical
    denMat1max[lower.tri(denMat1max)] <- t(denMat1max)[lower.tri(denMat1max)]

    # plot
    color_ramp <- colorRampPalette(c("white", "#FF0000"))(101)

    breakList = seq( 0, cap, by = cap/100 )

    C <- matrix( ncol = 2, nrow = length(breakList), data = 0)
    rownames( C ) <- color_ramp
    colnames( C ) <- c("rank","breaks")
    C[,1] <- 1:nrow(C)
    C[,2] <- breakList


    test_con <- tryCatch(file(output_name, open = "ab"), error = function(e) NULL)
    if (is.null(test_con)) {
      stop(sprintf(
        "Cannot write to '%s' -- the file may be open in another application. Please close it and try again.",
        output_name
      ))
    }
    close(test_con)

    pdf(
      output_name,
      width = 8.5, height = 11
    )

    par(bg  = "#eeeeee")

    plot( NULL,
          xlim = c(1,100), ylim = c(1,140),
          xlab = NA,       ylab = NA,
          xaxt = "n",      yaxt = "n",
          bty = "n",       asp = 1
    )

    par(omi = rep(0.5,4))
    par(mai = rep(0.5,4))
    par(bg  = "#eeeeee")

    w  <- calculate_w(100, 100, b_range$n)

    if (w > 2 && bin_size >= 5000){cat(underplot_message)}
    if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}
    if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}

    top   <- 134
    scale <- 1

    draw_title(interval_chr, interval_start, interval_end)

    top <- top - 10

    draw_scale(inherent=FALSE)

    top <- top - 20

    draw_viewpoints( denMat1max , C , left_side = 0, viewpoint_start, analysis_type, height)

    invisible(dev.off())

  }
  if(isTRUE(inherent)){
    cat(sprintf("  inh_col_floor: %s\n", inh_col_floor ))
    cat(sprintf("  inh_col_off:   %s\n", inh_col_off   ))
    cat(sprintf("  inh_col_on:    %s\n", inh_col_on    ))
    cat(sprintf("  inh_col_ceil:  %s\n", inh_col_ceil  ))

    denMat1 <- sparse_to_grid(sparMat1, range_grid)

    power_law_path <- list.files(
      sample_dir,
      pattern = paste0("inherentStats",".txt"),
      full.names = T)

    if(length(power_law_path) == 0){
      cat("\nPower laws unavailable for this sample! \n")
      q(status = 1)
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

    if (nrow(power_laws) == 0) {
      cat("\nPower laws unavailable at this resolution! \n")
      q(status = 1)
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

    # make symmetrical
    A[lower.tri(A)] <- t(A)[lower.tri(A)]

    if( max_cap != "none" ){
      cat("\n\nParameter --max_cap is not applicable when inherent = TRUE.\nContinuing without --max_cap...\n\n")
    }
    if( quant_cut != 1 ){
      cat("\n\nParameter --quant_cut is not applicable when inherent = TRUE.\nContinuing without --quant_cut...\n\n")
    }


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

    test_con <- tryCatch(file(output_name, open = "ab"), error = function(e) NULL)
    if (is.null(test_con)) {
      stop(sprintf(
        "Cannot write to '%s' -- the file may be open in another application. Please close it and try again.",
        output_name
      ))
    }
    close(test_con)

    pdf(
      output_name,
      width = 8.5, height = 11
    )

    par(bg  = "#eeeeee")

    plot( NULL,
          xlim = c(1,100), ylim = c(1,140),
          xlab = NA,       ylab = NA,
          xaxt = "n",      yaxt = "n",
          bty = "n",       asp = 1
    )

    par(omi = rep(0.5,4))
    par(mai = rep(0.5,4))
    par(bg  = "#eeeeee")

    w  <- calculate_w(100, 100, b_range$n)

    if (w > 2 && bin_size >= 5000){cat(underplot_message)}
    if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}
    if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}

    top   <- 134
    scale <- 1

    draw_title(interval_chr, interval_start, interval_end)

    top <- top - 10

    draw_scale(inherent=TRUE, steps=1, breaks=breaks)

    top <- top - 20

    draw_viewpoints( A , C , left_side = 0, viewpoint_start, analysis_type, height)

    invisible(dev.off())

  }
}




###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if( analysis_type == "two_sample" ){

  path_hic_A        <- Args[2]
  path_mgs_A        <- Args[3]
  path_hic_B        <- Args[4]
  path_mgs_B        <- Args[5]
  norm_method       <- Args[6]
  unit              <- Args[7]
  bin_size          <- as.numeric(Args[8])
  output_name       <- Args[9]
  genome_build      <- Args[10]
  flag_norm         <- Args[11]
  interval_chr         <- Args[12]
  interval_start       <- as.numeric(Args[13])
  interval_end         <- as.numeric(Args[14])
  quant_cut         <- as.numeric(Args[15])
  max_cap           <- Args[16]
  width             <- as.numeric(Args[17])
  height            <- Args[18]
  viewpoint_chr     <- Args[19]
  viewpoint_start   <- as.numeric(Args[20])
  sample_dirA       <- Args[21]
  sample_dirB       <- Args[22]
  inherent          <- Args[23]
  inh_col_floor     <- sprintf( "#%s", Args[24] )
  inh_col_off       <- sprintf( "#%s", Args[25] )
  inh_col_on        <- sprintf( "#%s", Args[26] )
  inh_col_ceil      <- sprintf( "#%s", Args[27] )


  if( viewpoint_chr != interval_chr ){
    cat("Viewpoint locus cannot be out of range\n")
    q(status = 1)
  }

  cat("Parameters\n")


  cat("\n")
  cat(sprintf("hic file A          %s\n", path_hic_A        ))
  cat(sprintf("hic file B          %s\n", path_hic_B        ))
  cat(sprintf("bin size            %s\n", bin_size          ))
  cat(sprintf("v4C_viewpoint_chr   %s\n", viewpoint_chr   ))
  cat(sprintf("v4C_viewpoint_start %d\n", viewpoint_start ))

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
  } else if( flag_strip_chr == "A" ){
    interval_chr_A <- interval_chr
    interval_chr_B <- sub( "chr", "", interval_chr )
  } else if( flag_strip_chr == "B" ){
    interval_chr_A <- sub( "chr", "", interval_chr )
    interval_chr_B <- interval_chr
  } else if( flag_strip_chr == "no" ){
    interval_chr_A <- interval_chr
    interval_chr_B <- interval_chr
  }

  mergeStats_A <- read.table( path_mgs_A, as.is = T )
  spikeVar_A   <- ncol(mergeStats_A)
  hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
  mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
  total1       <- sum( hg_total1 , mm_total1 )
  norm_factor1 <- 1000000 / total1
  aqua_factor1 <- hg_total1 / mm_total1

  mergeStats_B <- read.table( path_mgs_B, as.is = T )
  spikeVar_B   <- ncol(mergeStats_B)
  hg_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 1 ])
  mm_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 2 ])
  total2       <- sum( hg_total2 , mm_total2 )
  norm_factor2 <- 1000000 / total2
  aqua_factor2 <- hg_total2 / mm_total2

  if( ! flag_norm %in% c("blank", "none", "cpm", "aqua" ) ){
    cat("\nNorm should strictly be none, cpm or aqua in lower case \n")
    q(status = 1)
  }

  # set normalization factor for spike and non-spike data.
  # non-spike-in samples default to cpm.
  # spike in samples default to aqua.

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
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- aqua_factor1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- aqua_factor2
      flag_norm <- "aqua"
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
  cat(sprintf("normalization       %s\n", flag_norm         ))
  cat(      "\nfactors\n")
  cat(sprintf("  norm_factor1: %f\n", norm_factor1 ))
  cat(sprintf("  norm_factor2: %f\n", norm_factor2 ))
  cat(sprintf("  aqua_factor1: %f\n", aqua_factor1 ))
  cat(sprintf("  aqua_factor2: %f\n", aqua_factor2 ))


  b_range      <- interval_to_bins(interval_start, interval_end, bin_size)
  range_grid   <- bins_to_grid(b_range, bin_size)
  range_query_A <- bins_to_query(interval_chr_A, b_range, bin_size)
  range_query_B <- bins_to_query(interval_chr_B, b_range, bin_size)

  sparMat1 <- straw( norm_method,
                     path_hic_A,
                     range_query_A,
                     range_query_A,
                     unit,
                     bin_size)

  if( nrow(sparMat1) == 0 ){
    cat("\nNo contact values found for this interval in sample A\n")
    q(status = 1)
  }

  sparMat2 <- straw( norm_method,
                     path_hic_B,
                     range_query_B,
                     range_query_B,
                     unit,
                     bin_size)

  if( nrow(sparMat2) == 0 ){
    cat("\nNo contact values found for this interval in sample B\n")
    q(status = 1)
  }

  sparMat1$counts <- sparMat1$counts * norm_factor1 * aqua_factor1
  sparMat2$counts <- sparMat2$counts * norm_factor2 * aqua_factor2

  # convert sparse to dense on the shared bin grid (upper triangle)
  denMat1 <- sparse_to_grid(sparMat1, range_grid)
  denMat2 <- sparse_to_grid(sparMat2, range_grid)

  denDelta  <- denMat2 - denMat1

  if( max_cap != "none" && quant_cut != 1 ){
    cat("\nPlease provide either max_cap or quant_cut, not both!\n")
    q(status = 1)
  }

  # Max_cap and quant_cut handling
  cap <- calculate_cap(denDelta, max_cap, quant_cut)
  cap <- round(cap, 3)

  denDelta[denDelta >  cap] <-  cap
  denDelta[denDelta < -cap] <- -cap

  cat(sprintf("  cap: %s\n", cap ))

  # make symmetrical
  denDelta[lower.tri(denDelta)] <- t(denDelta)[lower.tri(denDelta)]

  # plot
  color_neg <- "#1E90FF"
  color_pos <- "#C71585"

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

  test_con <- tryCatch(file(output_name, open = "ab"), error = function(e) NULL)
  if (is.null(test_con)) {
    stop(sprintf(
      "Cannot write to '%s' -- the file may be open in another application. Please close it and try again.",
      output_name
    ))
  }
  close(test_con)

  pdf(
    output_name,
    width = 8.5, height = 11
  )

  par(omi = rep(0.5,4))
  par(mai = rep(0.5,4))
  par(bg  = "#eeeeee")

  plot( NULL,
        xlim = c(1,100), ylim = c(1,140),
        xlab = NA,       ylab = NA,
        xaxt = "n",      yaxt = "n",
        bty = "n",       asp = 1
  )

  w  <- calculate_w(100, 100, b_range$n)

  if (w > 2 && bin_size >= 5000){cat(underplot_message)}
  if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}
  if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}

  top   <- 134
  scale <- 1

  draw_title(interval_chr, interval_start, interval_end)

  top <- top - 10

  draw_scale( )

  top <- top - 20

  draw_viewpoints( denDelta , C , left_side = 0, viewpoint_start, analysis_type, height)

  invisible(dev.off())
}
