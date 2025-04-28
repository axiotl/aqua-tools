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
    A){
  
  aspect_ratio <- plot_height / plot_width
  
  diagonal_bins <- ceiling(sqrt(nrow(A)^2 + nrow(A)^2))
  w <- (aspect_ratio * plot_width / diagonal_bins)
  
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
    if (n < 2) {
      return(0)
    }
    mean_val <- mean(values)
    sd_val <- sd(values)
    if (sd_val == 0) {
      return(0)
    }
    return((n / ((n - 1) * (n - 2))) * sum(((values - mean_val) / sd_val) ^ 3))
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
    
    if (length(non_zero_values) == 0) {
      cat("\nNot enough data to plot. Please decrease the resolution or increase the range\n")
      quit(save = "no", status = 1)
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


is_gene_in_bedpe <- function(gene_start, gene_end, pairs) {
  if (!is.data.frame(pairs) || nrow(pairs) == 0) {
    return(FALSE)
  }
  
  for (p in seq_len(nrow(pairs))) {
    # Use the original coordinates
    region1_start <- pairs$start1[p]
    region1_end   <- pairs$end1[p]
    region2_start <- pairs$start2[p]
    region2_end   <- pairs$end2[p]
    
    # Check if the gene overlaps region 1
    if (gene_end >= region1_start && gene_start <= region1_end) {
      return(TRUE)
    }
    # Check if the gene overlaps region 2
    if (gene_end >= region2_start && gene_start <= region2_end) {
      return(TRUE)
    }
  }
  return(FALSE)
}


draw_feature_2 <- function(chr, start_pos, end_pos, label_txt, color, y_offset, flag_text, gene_cex=gene_cex, coord_cex=coord_cex) {

  scale_factor <- w * (ncol(A)-1) / interval_len # ncol(A)-1 removes annotation drift, related to strawr indexing 

  # Compute the x positions (in plot coordinates) corresponding to genomic start and end
  x_start <- ((start_pos - interval_start) * scale_factor) - (w*0.5) # shift half w to beginning of bin
  x_end   <- ((end_pos - interval_start) * scale_factor) - (w*0.5)

  y_pos <- 50 - y_offset

  # Draw line for gene span
  lines(c(x_start, x_end), rep(y_pos, 2), col = color, lwd = 0.5)

  if (flag_text) {
    composite_label <- sprintf("%s:%d-%d", chr, start_pos, end_pos)

    # Add gene coordinates
    text(x_start, y_pos, labels = composite_label, col = color,
         offset = -0.07, pos = 2, cex = coord_cex, srt = 45)
    # Add gene name in larger text size
    text(x_start, y_pos, labels = label_txt, col = color,
         offset = 0.07, pos = 2, cex = gene_cex, srt = 45)

    # Draw a dotted vertical line connecting the annotation to the plot
    segments(x_start, y_pos, x_start, 45, col = "#B2B2B2", lwd = 0.2, lty = "dotted")
  } else {
    # Dotted lines for custom annotations (when flag text = FALSE)
    segments(x_start, y_pos, x_start, 45, col = "#B2B2B2", lwd = 0.2, lty = "dotted")
    segments(x_end, y_pos, x_end, 45, col = "#B2B2B2", lwd = 0.2, lty = "dotted")
  }
}


draw_bed <- function(bed_file, color, depth, flag_text) {
  if (!file.exists(bed_file)) {
    stop("draw_bed error: file not found")
  }
  
  bed <- read.table(bed_file, as.is = TRUE, header = FALSE)
  colnames(bed) <- c("chr", "start", "end")
  
  # Subset the BED features to visualized genomic interval
  bed <- bed[bed$chr == interval_chr & bed$start < interval_end & bed$end > interval_start, ]
  if (nrow(bed) == 0) {
    warning("draw_bed warning: no features in the specified interval")
    return()
  }
  
  # Clip the feature boundaries to visualized interval
  bed$start <- pmax(bed$start, interval_start)
  bed$end   <- pmin(bed$end, interval_end)
  
  if (!flag_text) {
    for (i in seq_len(nrow(bed))) {
      b <- bed[i, ]
      draw_feature_2(b$chr, b$start, b$end, label_txt = "",
                     color = color, y_offset = depth, flag_text = FALSE)
    }
    return()
  }
  
  # Gene annotations
  bin_index <- floor(bed$start / bin_size)
  bed <- bed[order(bin_index), ]
  bin_counts <- table(bin_index)
  
  for (bin in names(bin_counts)) {
    feature_indices <- which(bin_index == bin)
    for (j in seq_along(feature_indices)) {
      b <- bed[feature_indices[j], ]
      # Adjust vertical offset for multiple genes per bin
      current_offset <- depth - (j - 1) * (-4)
      
      label_txt <- if (ncol(b) >= 4) b[[4]] else ""
      
      # Set gene name and coordinate sizes
      if (exists("pairs") && nrow(pairs) > 0 &&
          is_gene_in_bedpe(b$start, b$end, pairs)) {
        gene_cex  <- 0.25   
        coord_cex <- 0.08
        colour <- "#ef0000"
      } else {
        gene_cex  <- 0.10
        coord_cex <- 0.05
        colour <- "#aaaaaa"
      }
      
      draw_feature_2(b$chr, b$start, b$end, label_txt, color=colour, current_offset,
                     flag_text = TRUE, gene_cex = gene_cex, coord_cex = coord_cex)
    }
  }
}


draw_contacts <- function(A, C, left_side, bedpe) {
  start_x <- 0
  start_y <- 90

  # Compute diamond center coordinates given matrix indices
  get_center <- function(i, j, w) {
    cx <- start_x + (j - 1) * w + (i - 1) * w
    cy <- start_y + (i - 1) * (-w) + j * w
    c(cx, cy)
  }

  # Compute diamond vertices from the center
  get_vertices <- function(center, w) {
    list(
      left   = c((center[1] - w) / 2, center[2] / 2),
      top    = c(center[1] / 2, (center[2] + w) / 2),
      right  = c((center[1] + w) / 2, center[2] / 2),
      bottom = c(center[1] / 2, (center[2] - w) / 2)
    )
  }

  # Draw a diamond given its center
  draw_diamond <- function(center, w, color, highlight = FALSE) {
    verts <- get_vertices(center, w)
    x_coords <- c(verts$left[1], verts$top[1], verts$right[1], verts$bottom[1])
    y_coords <- c(verts$left[2], verts$top[2], verts$right[2], verts$bottom[2])
    polygon(x_coords, y_coords, col = color, border = NA)
    if (highlight) {
      polygon(x_coords, y_coords, col = color, border = bedpe_color, lwd = 0.5)
    }
  }

  # Draw all diamonds
  bin_coordinates <- colnames(A)
  n_rows <- nrow(A)
  n_cols <- ncol(A)

  for (i in seq_len(n_rows)) {
    for (j in i:n_cols) {
      color <- rownames(C[order(abs(C[,"breaks"] - A[i, j]), decreasing = FALSE), ])[1]
      center <- get_center(i, j, w)
      draw_diamond(center, w, color, highlight = FALSE)
    }

    # Set annotation interval based on bin size
    annotation_interval <- ifelse(bin_size == 1000, 20, ifelse(bin_size %in% c(5000, 10000, 25000, 50000), 10, ifelse(bin_size >= 100000, 2, 10)))

    for (j in seq(1, n_cols, by = annotation_interval)) {
      bin_label <- bin_coordinates[j]
      label_x <- (start_x + (j - 1) * w)-(0.5*w)
      label_y <- 44.25
      text(label_x, label_y, labels = bin_label, cex = 0.05, col = "#8B8B8B")
      segments(x0 = label_x, y0 = 45, x1 = label_x, y1 = label_y + 0.25,
               col = "#B2B2B2", lty = "dotted", lwd = 0.2)
    }
  }
  # Draw BEDPE highlight polygons
  if (bedpe != "FALSE") {
    for (p in seq_len(nrow(pairs))) {

      i_start <- pairs$i[p]
      i_end   <- pairs$i_end[p]
      j_start <- pairs$j[p]
      j_end   <- pairs$j_end[p]

      # Check if this BEDPE feature is within a single bin
      if (i_start == i_end && j_start == j_end) {
        # Single bin highlight: draw the diamond's border for that bin
        center <- get_center(i_start, j_start, w)
        draw_diamond(center, w, color = NA, highlight = TRUE)
      } else {
        # Multi-bin highlight: get the vertices for the four outer diamonds
        left_center   <- get_center(i_start, j_start, w)
        left_v        <- get_vertices(left_center, w)$left

        top_center    <- get_center(i_start, j_end, w)
        top_v         <- get_vertices(top_center, w)$top

        right_center  <- get_center(i_end, j_end, w)
        right_v       <- get_vertices(right_center, w)$right

        bottom_center <- get_center(i_end, j_start, w)
        bottom_v      <- get_vertices(bottom_center, w)$bottom

        polygon(x = c(left_v[1], top_v[1], right_v[1], bottom_v[1]),
                y = c(left_v[2], top_v[2], right_v[2], bottom_v[2]),
                border = bedpe_color, col = NA, lwd = 0.5)
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

  # +bin_size fixes annotation drift 
  calculate_position <- function(pos) {
    relative_position <- (pos - interval_start) / (interval_len + bin_size)
    return(relative_position * ifelse(axis == "x", plot_width, plot_height))
  }

  str <- calculate_position(clipped_start)
  end <- calculate_position(clipped_end)
  
  x_adjustment <- 2 * (plot_height / 100)
  y_adjustment <- 4 * (plot_height / 100)
  
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
      dotted_line_end <- x_position - x_adjustment
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


draw_contacts_inter <- function(A, C, left_side, bedpe, plot_width, plot_height, top) {
  
  n_rows <- nrow(A)
  n_cols <- ncol(A)
  cell_width  <- plot_width / n_cols
  cell_height <- plot_height / n_rows
  
  # Draw contact matrix
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      color <- rownames(C[order(abs(C[,"breaks"] - A[i, j]), decreasing = FALSE), ])[1]
      rect(
        left_side + (j - 1) * cell_width, 
        top - (i - 1) * cell_height,
        left_side + j * cell_width, 
        top - i * cell_height,
        col = color, 
        border = NA
      )
    }
  }
  
  # Annotate the x-axis (columns) with labels aligned at the start (left edge) of the bin
  annotation_interval <- ifelse(bin_size == 1000, 20, ifelse(bin_size %in% c(5000, 10000, 25000, 50000), 10, ifelse(bin_size >= 100000, 2, 10)))
  
  bottom_y <- top - n_rows * cell_height
  for (j in seq(1, n_cols, by = annotation_interval)) {
    bin_label <- colnames(A)[j]
    # Place the label at the left edge of the cell
    label_x <- left_side + (j - 1) * cell_width
    # Position the label just below the contact matrix
    label_y <- bottom_y - 0.5 * cell_height
    text(label_x, label_y, labels = bin_label, cex = 0.05, col = "#8B8B8B", adj = 0)
    segments(x0 = label_x, y0 = bottom_y, x1 = label_x, y1 = label_y + 0.25,
             col = "#B2B2B2", lty = "dotted", lwd = 0.2)
  }
  
  # Annotate the y-axis (rows) with labels aligned at the top edge of the bin
  for (i in seq(1, n_rows, by = annotation_interval)) {
    bin_label <- rownames(A)[i]
    # Place the label at the top edge of the row
    label_y <- top - (i - 1) * cell_height
    # Position the label just to the left of the contact matrix
    label_x <- left_side - 0.5 * cell_width
    text(label_x, label_y, labels = bin_label, cex = 0.05, col = "#8B8B8B", adj = 1)
    segments(x0 = left_side, y0 = label_y, x1 = label_x + 0.25, y1 = label_y,
             col = "#B2B2B2", lty = "dotted", lwd = 0.2)
  }
  
  # Draw BEDPE highlights
  if (bedpe != "FALSE") {
    for (p in seq_len(nrow(pairs))) {
      x1 <- left_side + (pairs$i[p] - 1) * cell_width
      y1 <- top - (pairs$j[p] - 1) * cell_height
      x2 <- left_side + (pairs$i_end[p]) * cell_width
      y2 <- top - (pairs$j_end[p]) * cell_height
      
      # Draw bedpe highlight
      rect(
        x1, y1, x2, y2,
        border = bedpe_color, 
        col = NA,
        lwd = 1
      )
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
    
    if (flag_matrix == "FALSE") {
      
      if (bedpe != "FALSE") {
        pairs <- read.table( bedpe, as.is = TRUE)[,1:6]
        colnames(pairs) <- c("chr1","start1","end1","chr2","start2","end2")
        
        if( sum( grepl( "chr", pairs$chr1 ) ) == 0 ) { pairs$chr1 <- paste0( "chr", pairs$chr1 ) }
        if( sum( grepl( "chr", pairs$chr2 ) ) == 0 ) { pairs$chr2 <- paste0( "chr", pairs$chr2 ) }
        
        # Adjust end coordinates for half open binning
        pairs$end1 <- ifelse(pairs$end1 > pairs$start1, pairs$end1 - 1, pairs$end1)
        pairs$end2 <- ifelse(pairs$end2 > pairs$start2, pairs$end2 - 1, pairs$end2)
        
        # subset pairs to interval length
        pairs <- pairs[ pairs[, "chr1"] == interval_chr, ]
        pairs <- pairs[ 
          pairs[,"end2"] <= interval_end & 
            pairs[,"start1"] >= interval_start & 
            pairs[,"end1"] <= interval_end & 
            pairs[,"start2"] >= interval_start, ]
        
        if (nrow(pairs) == 0) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          pairs <- pairs[ pairs$chr1 %in% c( paste0( rep("chr",22), 1:22 ), "chrX", "chrY" ) , ]
          
          # convert coordinates to bins
          pairs$start1_bin <- as.integer( floor(   pairs$start1 / bin_size ) * bin_size )
          pairs$end1_bin   <- as.integer( floor( pairs$end1   / bin_size ) * bin_size )
          pairs$start2_bin <- as.integer( floor(   pairs$start2 / bin_size ) * bin_size )
          pairs$end2_bin   <- as.integer( floor( pairs$end2   / bin_size ) * bin_size )
          
          
          for( i in 1:nrow(pairs) ){
            
            if( pairs[i,"start2_bin"] < pairs[i,"start1_bin"] ){
              a <- pairs[i,"start1_bin"] ; b <- pairs[i,"end1_bin"]
              c <- pairs[i,"start2_bin"] ; d <- pairs[i,"end2_bin"]
              
              pairs[i,"start1_bin"] <- c
              pairs[i,"end1_bin"  ] <- d
              pairs[i,"start2_bin"] <- a
              pairs[i,"end2_bin"  ] <- b
            }
            
            if( pairs[i,"chr1"] != pairs[i,"chr2"] ){
              cat("Please provide intra-chromosomal bedpes\n") ; quit(save="no")
            }
          }
          
          # Calculate plotting coordinates and add 1 for 1-based indexing
          pairs$i <- (pairs$start1_bin - interval_start) / bin_size + 1 
          pairs$j <- (pairs$start2_bin - interval_start) / bin_size + 1
          pairs$i_end <- (pairs$end1_bin - interval_start) / bin_size + 1
          pairs$j_end <- (pairs$end2_bin - interval_start) / bin_size + 1
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

    if(isTRUE(inherent)){
        flag_norm <- "inherent"
    }
    
    cat(      "\n  factors\n")
    cat(sprintf("  norm_factor: %f\n", norm_factor1 ))
    cat(sprintf("  aqua_factor: %f\n", aqua_factor1 ))
    
    ## Process matrix:
    
    cat(      "\n  straw parameters\n")
    cat(sprintf("       norm: %s\n", flag_norm  ))
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
        w  <- calculate_w(100, 140, A)
      }
      
      if (flag_w != "blank"){
        w <- as.numeric(flag_w)
        max_w <- calculate_w(100, 140, A)
        if (w > max_w){
          cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
          w <- max_w
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

      base_depth <- 6.00
      step_depth <- 0.25
      
      if(isFALSE(flag_profiles)){

        legend_labels <- c()
        legend_colors <- c()

        # Default Annotation (TSS)
        if (ann_default) {
          cat("\nDrawing default annotations\n")
          # figure out how many custom annotation files there are
          n_custom <- if (ann_custom != "NONE") {
          length(unlist(strsplit(ann_custom, " ")))
        } else {
          0
        }
          
        # compute dynamic TSS depth
        if (n_custom > 0) {
          tss_depth <- base_depth + step_depth * n_custom + step_depth
          cat("\n")
        } else {
          # when no customs, bring closer to the diagonal
          tss_depth <- base_depth + (2 * step_depth)
        }
        
        if (genome_build == "hg19") {
          bed_tss <- file.path(data_dir, genome_build, "reference", "GENCODE_TSSs_hg19.bed")
        } else if (genome_build == "mm10") {
          bed_tss <- file.path(data_dir, genome_build, "reference", "GENCODE_TSSs_200bp_mm10.bed")
        } else if (genome_build == "hg38") {
          bed_tss <- file.path(data_dir, genome_build, "reference", "GENCODE_TSSs_hg38.bed")
        }
          
        draw_bed(bed_tss, "#ef0000", tss_depth, TRUE)
        legend_labels <- c(legend_labels, "transcription start sites")
        legend_colors <- c(legend_colors, "#ef0000")

        } else {
          no_def_ann <- TRUE
        }
        
        # Custom Annotations
        if (ann_custom != "NONE") {
          # Split the annotation files string into a vector
         annotation_files <- unlist(strsplit(ann_custom, " "))
         n_files <- length(annotation_files)

          # default base colors
          base_colors <- c("#E69F00",  # orange 
                  "#377EB8",  # blue  
                  "#4DAF4A",  # green
                  "#C77CFF",  # light purple    
                  "#984EA3")  # purple

          # build dynamic depths
          track_depths <- seq(from = base_depth,
                              by   = step_depth,
                              length.out = n_files)

          # build a color vector of exactly n_files
          if (ann_custom_colors == "NONE") {
            # ramp base colors out to n_files
            final_colors <- colorRampPalette(base_colors)(n_files)
          } else {
            # split & sanitize any user-supplied ones
            custom_colors <- strsplit(ann_custom_colors, " ")[[1]]
            custom_colors <- sapply(custom_colors,
                                      function(x) if (substr(x,1,1)=="#") x else paste0("#",x))
            if (length(custom_colors) >= n_files) {
              final_colors <- custom_colors[1:n_files]
            } else {
              # fill the rest by ramping the base colors
              extra <- n_files - length(custom_colors)
              fill  <- colorRampPalette(base_colors)(extra)
              final_colors <- c(custom_colors, fill)
            }
          }

          # Loop and pick final_colors[i] every time
          for (i in seq_along(annotation_files)) {
            file  <- annotation_files[i]
            color <- final_colors[i]
            depth <- track_depths[i]

            cat(sprintf("Drawing custom annotation for file %s with color %s \n", basename(file), color))
            draw_bed(file, color, depth, FALSE)
            legend_labels <- c(legend_labels, basename(file))
            legend_colors <- c(legend_colors, color)
          }
        } else {
          no_cus_ann <- TRUE
        }
          
        # Draw legend
        if (length(legend_labels) > 0) {
          legend(-2, 125, legend = legend_labels, col = legend_colors, lwd = 1,
                  box.lty = 0, cex = 0.8, text.col = "#888888")
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
      
      if( max_cap != "none" ){
        cat("\nParameter --max_cap is not applicable when inherent = TRUE\nContinuing without --max_cap...\n\n")
      }
      if( quant_cut != 1 ){
        cat("\nParameter --quant_cut is not applicable when inherent = TRUE\nContinuing without --quant_cut...\n\n")
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
            xlim = c(0,100), ylim = c(0,140),
            xlab = NA,       ylab = NA,
            xaxt = "n",      yaxt = "n",
            bty = "n",       asp = 1 )
      
      
      if(flag_w == "blank"){
        w  <- calculate_w(100, 140, A)
      }
      
      if (flag_w != "blank"){
        w <- as.numeric(flag_w)
        max_w  <- calculate_w(100, 140, A)
        
        if (w > max_w){
          cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
          w <- max_w
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
      colors         <- c( inh_col_floor, colors_0_1, colors_1_2 )
      
      breaks <- c(
        (1:colors_n)/colors_n,
        (1:colors_n)/colors_n+1 )
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

      base_depth <- 6.00
      step_depth <- 0.25
      
      if(isFALSE(flag_profiles)){

        legend_labels <- c()
        legend_colors <- c()

        # Default Annotation (TSS)
        if (ann_default) {
          cat("\nDrawing default annotations\n")
          # figure out how many custom annotation files there are
          n_custom <- if (ann_custom != "NONE") {
          length(unlist(strsplit(ann_custom, " ")))
        } else {
          0
        }
          
        # compute dynamic TSS depth
        if (n_custom > 0) {
          tss_depth <- base_depth + step_depth * n_custom + step_depth
          cat("\n")
        } else {
          # when no customs, bring closer to the diagonal
          tss_depth <- base_depth + (2 * step_depth)
        }
        
        if (genome_build == "hg19") {
          bed_tss <- file.path(data_dir, genome_build, "reference", "GENCODE_TSSs_hg19.bed")
        } else if (genome_build == "mm10") {
          bed_tss <- file.path(data_dir, genome_build, "reference", "GENCODE_TSSs_200bp_mm10.bed")
        } else if (genome_build == "hg38") {
          bed_tss <- file.path(data_dir, genome_build, "reference", "GENCODE_TSSs_hg38.bed")
        }
          
        draw_bed(bed_tss, "#ef0000", tss_depth, TRUE)
        legend_labels <- c(legend_labels, "transcription start sites")
        legend_colors <- c(legend_colors, "#ef0000")

        } else {
          no_def_ann <- TRUE
        }
        
        # Custom Annotations
        if (ann_custom != "NONE") {
          # Split the annotation files string into a vector
         annotation_files <- unlist(strsplit(ann_custom, " "))
         n_files <- length(annotation_files)

          # default base colors
          base_colors <- c("#E69F00",  # orange 
                  "#377EB8",  # blue  
                  "#4DAF4A",  # green
                  "#C77CFF",  # light purple    
                  "#984EA3")  # purple

          # build dynamic depths
          track_depths <- seq(from = base_depth,
                              by   = step_depth,
                              length.out = n_files)

          # build a color vector of exactly n_files
          if (ann_custom_colors == "NONE") {
            # ramp base colors out to n_files
            final_colors <- colorRampPalette(base_colors)(n_files)
          } else {
            # split & sanitize any user-supplied ones
            custom_colors <- strsplit(ann_custom_colors, " ")[[1]]
            custom_colors <- sapply(custom_colors,
                                      function(x) if (substr(x,1,1)=="#") x else paste0("#",x))
            if (length(custom_colors) >= n_files) {
              final_colors <- custom_colors[1:n_files]
            } else {
              # fill the rest by ramping the base colors
              extra <- n_files - length(custom_colors)
              fill  <- colorRampPalette(base_colors)(extra)
              final_colors <- c(custom_colors, fill)
            }
          }

          # Loop and pick final_colors[i] every time
          for (i in seq_along(annotation_files)) {
            file  <- annotation_files[i]
            color <- final_colors[i]
            depth <- track_depths[i]

            cat(sprintf("Drawing custom annotation for file %s with color %s \n", basename(file), color))
            draw_bed(file, color, depth, FALSE)
            legend_labels <- c(legend_labels, basename(file))
            legend_colors <- c(legend_colors, color)
          }
        } else {
          no_cus_ann <- TRUE
        }
          
        # Draw legend
        if (length(legend_labels) > 0) {
          legend(-2, 125, legend = legend_labels, col = legend_colors, lwd = 1,
                  box.lty = 0, cex = 0.8, text.col = "#888888")
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
        
        # Adjust end coordinates for half open binning
        pairs$end1 <- ifelse(pairs$end1 > pairs$start1, pairs$end1 - 1, pairs$end1)
        pairs$end2 <- ifelse(pairs$end2 > pairs$start2, pairs$end2 - 1, pairs$end2)
        
        # subset pairs to interval length
        pairs <- pairs[ pairs[, "chr1"] == interval_chr, ]
        pairs <- pairs[ 
          pairs[,"end2"] <= interval_end & 
            pairs[,"start1"] >= interval_start & 
            pairs[,"end1"] <= interval_end & 
            pairs[,"start2"] >= interval_start, ]
        
        if (nrow(pairs) == 0) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          pairs <- pairs[ pairs$chr1 %in% c( paste0( rep("chr",22), 1:22 ), "chrX", "chrY" ) , ]
          
          # convert coordinates to bins
          pairs$start1_bin <- as.integer( floor(   pairs$start1 / bin_size ) * bin_size )
          pairs$end1_bin   <- as.integer( floor( pairs$end1   / bin_size ) * bin_size )
          pairs$start2_bin <- as.integer( floor(   pairs$start2 / bin_size ) * bin_size )
          pairs$end2_bin   <- as.integer( floor( pairs$end2   / bin_size ) * bin_size )
          
          
          for( i in 1:nrow(pairs) ){
            
            if( pairs[i,"start2_bin"] < pairs[i,"start1_bin"] ){
              a <- pairs[i,"start1_bin"] ; b <- pairs[i,"end1_bin"]
              c <- pairs[i,"start2_bin"] ; d <- pairs[i,"end2_bin"]
              
              pairs[i,"start1_bin"] <- c
              pairs[i,"end1_bin"  ] <- d
              pairs[i,"start2_bin"] <- a
              pairs[i,"end2_bin"  ] <- b
            }
            
            if( pairs[i,"chr1"] != pairs[i,"chr2"] ){
              cat("Please provide intra-chromosomal bedpes\n") ; quit(save="no")
            }
          }
          
          # Calculate plotting coordinates and add 1 for 1-based indexing
          pairs$i <- (pairs$start1_bin - interval_start) / bin_size + 1 
          pairs$j <- (pairs$start2_bin - interval_start) / bin_size + 1
          pairs$i_end <- (pairs$end1_bin - interval_start) / bin_size + 1
          pairs$j_end <- (pairs$end2_bin - interval_start) / bin_size + 1

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
    
    cat(      "\n  factors\n")
    cat(sprintf("  norm_factor1: %f\n", norm_factor1 ))
    cat(sprintf("  norm_factor2: %f\n", norm_factor2 ))
    cat(sprintf("  aqua_factor1: %f\n", aqua_factor1 ))
    cat(sprintf("  aqua_factor2: %f\n", aqua_factor2 ))
    
    ## Process matrix:
    
    cat(      "\n  straw parameters\n")
    cat(sprintf("       norm: %s\n", flag_norm  ))
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
        w  <- calculate_w(100, 140, A)
      }
      
      if (flag_w != "blank"){
        w <- as.numeric(flag_w)
        max_w <- calculate_w(100, 140, A)
        
        if (w > max_w){
          cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
          w <- max_w
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
      
      base_depth <- 6.00
      step_depth <- 0.25

      if(isFALSE(flag_profiles)){

        legend_labels <- c()
        legend_colors <- c()

        # Default Annotation (TSS)
        if (ann_default) {
          cat("\nDrawing default annotations\n")
          # figure out how many custom annotation files there are
          n_custom <- if (ann_custom != "NONE") {
          length(unlist(strsplit(ann_custom, " ")))
        } else {
          0
        }
            
        # compute dynamic TSS depth
        if (n_custom > 0) {
          tss_depth <- base_depth + step_depth * n_custom + step_depth
          cat("\n")
        } else {
          # when no customs, bring closer to the diagonal
          tss_depth <- base_depth + (2 * step_depth)
        }
          
        if (genome_build == "hg19") {
          bed_tss <- file.path(data_dir, genome_build, "reference", "GENCODE_TSSs_hg19.bed")
        } else if (genome_build == "mm10") {
          bed_tss <- file.path(data_dir, genome_build, "reference", "GENCODE_TSSs_200bp_mm10.bed")
        } else if (genome_build == "hg38") {
          bed_tss <- file.path(data_dir, genome_build, "reference", "GENCODE_TSSs_hg38.bed")
        }
            
        draw_bed(bed_tss, "#ef0000", tss_depth, TRUE)
        legend_labels <- c(legend_labels, "transcription start sites")
        legend_colors <- c(legend_colors, "#ef0000")

        } else {
            no_def_ann <- TRUE
        }
          
        # Custom Annotations
        if (ann_custom != "NONE") {
          # Split the annotation files string into a vector
          annotation_files <- unlist(strsplit(ann_custom, " "))
          n_files <- length(annotation_files)

          # default base colors
          base_colors <- c("#E69F00",  # orange 
                "#377EB8",  # blue  
                "#4DAF4A",  # green
                "#C77CFF",  # light purple    
                "#984EA3")  # purple

          # build dynamic depths
          track_depths <- seq(from = base_depth,
                              by   = step_depth,
                              length.out = n_files)

          # build a color vector of exactly n_files
          if (ann_custom_colors == "NONE") {
            # ramp base colors out to n_files
            final_colors <- colorRampPalette(base_colors)(n_files)
          } else {
            # split & sanitize any user-supplied ones
            custom_colors <- strsplit(ann_custom_colors, " ")[[1]]
            custom_colors <- sapply(custom_colors,
                                      function(x) if (substr(x,1,1)=="#") x else paste0("#",x))
            if (length(custom_colors) >= n_files) {
              final_colors <- custom_colors[1:n_files]
            } else {
              # fill the rest by ramping the base colors
              extra <- n_files - length(custom_colors)
              fill  <- colorRampPalette(base_colors)(extra)
              final_colors <- c(custom_colors, fill)
            }
          }

          # Loop and pick final_colors[i] every time
          for (i in seq_along(annotation_files)) {
            file  <- annotation_files[i]
            color <- final_colors[i]
            depth <- track_depths[i]

            cat(sprintf("Drawing custom annotation for file %s with color %s \n", basename(file), color))
            draw_bed(file, color, depth, FALSE)
            legend_labels <- c(legend_labels, basename(file))
            legend_colors <- c(legend_colors, color)
          }
        } else {
          no_cus_ann <- TRUE
        }
          
        # Draw legend
        if (length(legend_labels) > 0) {
          legend(-2, 125, legend = legend_labels, col = legend_colors, lwd = 1,
                box.lty = 0, cex = 0.8, text.col = "#888888")
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
        
        # Make sure chromosome names have 'chr' prefix
        pairs$chr1 <- ifelse(grepl("chr", pairs$chr1), pairs$chr1, paste0("chr", pairs$chr1))
        pairs$chr2 <- ifelse(grepl("chr", pairs$chr2), pairs$chr2, paste0("chr", pairs$chr2))
        
        # Filter the BEDPE pairs using the appropriate interval variables
        pairs <- pairs[
          (
            (pairs$chr1 == new_interval_chr1 & pairs$start1 >= interval_start1 & pairs$end1 <= interval_end1 &
               pairs$chr2 == new_interval_chr2 & pairs$start2 >= interval_start2 & pairs$end2 <= interval_end2) |
              (pairs$chr1 == new_interval_chr2 & pairs$start1 >= interval_start2 & pairs$end1 <= interval_end2 &
                 pairs$chr2 == new_interval_chr1 & pairs$start2 >= interval_start1 & pairs$end2 <= interval_end1)
          ), 
        ]
        
        if (nrow(pairs) == 0) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          # Reorder columns if intervals are swapped
          for (idx in 1:nrow(pairs)) {
            if (pairs[idx, "chr1"] != new_interval_chr1) {
              pairs[idx, c("chr1", "start1", "end1", "chr2", "start2", "end2")] <- 
                pairs[idx, c("chr2", "start2", "end2", "chr1", "start1", "end1")]
            }
          }
          
          # Reset column names (if needed)
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
          
          # Adjust end coordinates for half open binning
          pairs$end1 <- ifelse(pairs$end1 > pairs$start1, pairs$end1 - 1, pairs$end1)
          pairs$end2 <- ifelse(pairs$end2 > pairs$start2, pairs$end2 - 1, pairs$end2)
          
          # Convert coordinates to bins
          pairs$start1_bin <- as.integer(floor(pairs$start1 / bin_size) * bin_size)
          pairs$end1_bin   <- as.integer(floor(pairs$end1 / bin_size) * bin_size)
          pairs$start2_bin <- as.integer(floor(pairs$start2 / bin_size) * bin_size)
          pairs$end2_bin   <- as.integer(floor(pairs$end2 / bin_size) * bin_size)
          
          # Calculate plotting coordinates
          # Use interval_start1 for the first set and interval_start2 for the second set
          pairs$i <- (pairs$start1_bin - interval_start1) / bin_size + 1 
          pairs$j <- (pairs$start2_bin - interval_start2) / bin_size + 1
          pairs$i_end <- (pairs$end1_bin - interval_start1) / bin_size + 1
          pairs$j_end <- (pairs$end2_bin - interval_start2) / bin_size + 1

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
    
    # set normalization factor for spike-in and non-spike-in data.
    # non-spike-in data defaults to cpm.
    # spike in data defaults to aqua.
    
    if(spikeVar == 1){
      if(flag_norm == "blank"){
        norm_factor1 <- norm_factor1
        aqua_factor1 <- 1
        flag_nor <- "cpm"
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
    
    cat(      "\n  factors\n")
    cat(sprintf("  norm_factor1: %f\n", norm_factor1 ))
    cat(sprintf("  aqua_factor1: %f\n", aqua_factor1 ))
    
    ## Process matrix:
    
    cat(      "\n  straw parameters\n")
    cat(sprintf("       norm: %s\n", flag_norm ))
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
      w <- calculate_w(plot_width, plot_height, A)
    }
    
    if (flag_w != "blank"){
      w <- as.numeric(flag_w)
      max_w <- calculate_w(plot_width, plot_height, A)
      if (w > max_w){
        cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
        w <- max_w
      }
    }
    
    if (w > 2 && bin_size >= 5000){cat(underplot_message)}
    if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}
    if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}
    
    cat(sprintf("plotting bin width: %s\n", round(w,3)))
    top        <- 140
    
    draw_title_inter(new_interval_chr1, interval_start1, interval_end1, new_interval_chr2, interval_start2, interval_end2)
    
    top <- top - 12
    
    draw_scale( )
    
    top <- top - 30
    
    draw_contacts_inter(A, C, left_side = 0, bedpe, plot_width, plot_height, top)
    
    if (isFALSE(flag_profiles)) {
      # Initialize legend elements
      legend_labels <- c()
      legend_colors <- c()
      
      # Default Annotations
      if (ann_default) {
        if (genome_build == "hg19") {
          cat("\nDrawing default annotations for hg19\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/GENCODE_TSSs_hg19.bed", sep = "")
          
          # Draw annotations for x- and y-axes
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          draw_bed_inter(bed_tss, "#ef0000", -0.2, TRUE, w, top, "y", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          
          legend_labels <- c(legend_labels, "transcription start sites")
          legend_colors <- c(legend_colors, "#ef0000")
        }
        if (genome_build == "hg38") {
          cat("\nDrawing default annotations for hg38\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/GENCODE_TSSs_hg38.bed", sep = "")
          
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "y", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          
          legend_labels <- c(legend_labels, "transcription start sites")
          legend_colors <- c(legend_colors, "#ef0000")
        }
        if (genome_build == "mm10") {
          cat("\nDrawing default annotations for mm10\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/GENCODE_TSSs_200bp_mm10.bed", sep = "")
          
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "y", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          
          legend_labels <- c(legend_labels, "transcription start sites")
          legend_colors <- c(legend_colors, "#ef0000")
        }
      } else {
        no_def_ann <- TRUE
      }
      
      # Custom Annotations
      if (ann_custom != "NONE") {
        # Split the annotation files string into a vector
        annotation_files <- unlist(strsplit(ann_custom, " "))
        
        # Define default colors and track depths
        default_colors <- c("#ef0000", "#00829d", "#8622b7", "#4C8219")
        track_depths <- c(-0.5, -0.8, -1.1, -1.4)
        
        for (i in seq_along(annotation_files)) {
          current_file <- annotation_files[i]
          
          # Decide which color to use
          if (ann_custom_colors == "NONE") {
            # Use default color
            this_color <- default_colors[(i - 1) %% length(default_colors) + 1]
            cat(sprintf("\nUsing default color %s for custom annotation file %s", this_color, current_file))
          } else {
            # Split the provided colors and make sure a '#' is included
            colors <- unlist(strsplit(ann_custom_colors, " "))
            if (length(colors) >= i) {
              # Prepend '#' if necessary
              this_color <- if (substr(colors[i], 1, 1) == "#") colors[i] else paste0("#", colors[i])
            } else {
              this_color <- default_colors[(i - 1) %% length(default_colors) + 1]
            }
            cat(sprintf("\nDrawing custom annotations for file %s with color %s", current_file, this_color))
          }
          
          # Draw annotations for both x and y axes
          draw_bed_inter(current_file, this_color, track_depths[(i - 1) %% length(track_depths) + 1], FALSE, w, top, "x",
                         interval_start1, interval_len1, interval_start2, interval_len2,
                         plot_width, plot_height)
          draw_bed_inter(current_file, this_color, track_depths[(i - 1) %% length(track_depths) + 1], FALSE, w, top, "y",
                         interval_start1, interval_len1, interval_start2, interval_len2,
                         plot_width, plot_height)
          
          # Update legend elements
          legend_labels <- c(legend_labels, basename(current_file))
          legend_colors <- c(legend_colors, this_color)
        }
      } else {
        no_cus_ann <- TRUE
      }
      
      # Draw the legend
      if (length(legend_labels) > 0) {
        legend(-2, 127, legend = legend_labels, col = legend_colors, lwd = 1,
               box.lty = 0, cex = 0.8, text.col = "#888888")
      }
    }
    
    if (flag_profiles) {
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
      new_interval_chr1 <- sub( "chr", "", interval_chr1 )
      new_interval_chr2 <- sub( "chr", "", interval_chr2 )
    } else if( flag_strip_chr == "B" ){
      new_interval_chr1 <- interval_chr1
      new_interval_chr2 <- sub( "chr", "", interval_chr2 )
    } else if( flag_strip_chr == "A" ){
      new_interval_chr1 <- sub( "chr", "", interval_chr1 )
      new_interval_chr2 <- interval_chr2
    } else if( flag_strip_chr == "no" ){
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
        
        # Filter the BEDPE pairs using the appropriate interval variables
        pairs <- pairs[
          (
            (pairs$chr1 == new_interval_chr1 & pairs$start1 >= interval_start1 & pairs$end1 <= interval_end1 &
               pairs$chr2 == new_interval_chr2 & pairs$start2 >= interval_start2 & pairs$end2 <= interval_end2) |
              (pairs$chr1 == new_interval_chr2 & pairs$start1 >= interval_start2 & pairs$end1 <= interval_end2 &
                 pairs$chr2 == new_interval_chr1 & pairs$start2 >= interval_start1 & pairs$end2 <= interval_end1)
          ), 
        ]
        
        if (nrow(pairs) == 0) {
          cat("\nThere are no bedpe pairs in the specified range. Continuing...\n")
          bedpe <- "FALSE"
        } else {
          # Reorder columns if intervals are swapped
          for (idx in 1:nrow(pairs)) {
            if (pairs[idx, "chr1"] != new_interval_chr1) {
              pairs[idx, c("chr1", "start1", "end1", "chr2", "start2", "end2")] <- 
                pairs[idx, c("chr2", "start2", "end2", "chr1", "start1", "end1")]
            }
          }
          
          # Reset column names
          colnames(pairs) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
          
          # Adjust end coordinates for half open binning
          pairs$end1 <- ifelse(pairs$end1 > pairs$start1, pairs$end1 - 1, pairs$end1)
          pairs$end2 <- ifelse(pairs$end2 > pairs$start2, pairs$end2 - 1, pairs$end2)
          
          # Convert coordinates to bins
          pairs$start1_bin <- as.integer(floor(pairs$start1 / bin_size) * bin_size)
          pairs$end1_bin   <- as.integer(floor(pairs$end1 / bin_size) * bin_size)
          pairs$start2_bin <- as.integer(floor(pairs$start2 / bin_size) * bin_size)
          pairs$end2_bin   <- as.integer(floor(pairs$end2 / bin_size) * bin_size)
          
          # Calculate plotting coordinates
          pairs$i <- (pairs$start1_bin - interval_start1) / bin_size + 1 
          pairs$j <- (pairs$start2_bin - interval_start2) / bin_size + 1
          pairs$i_end <- (pairs$end1_bin - interval_start1) / bin_size + 1
          pairs$j_end <- (pairs$end2_bin - interval_start2) / bin_size + 1

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
    
    cat(      "\n  factors\n")
    cat(sprintf("  norm_factor1: %f\n", norm_factor1 ))
    cat(sprintf("  norm_factor2: %f\n", norm_factor2 ))
    cat(sprintf("  aqua_factor1: %f\n", aqua_factor1 ))
    cat(sprintf("  aqua_factor2: %f\n", aqua_factor2 ))
    
    ## Process matrix:
    
    cat(      "\n  straw parameters\n")
    cat(sprintf("       norm: %s\n", flag_norm  ))
    cat(sprintf("      hic_A: %s\n", path_hic_A ))
    cat(sprintf("      hic_B: %s\n", path_hic_B ))
    cat(sprintf("   interval: %s\n", paste( interval_chr1, interval_start1, interval_end1, sep = ":"  ) ))
    cat(sprintf("   interval: %s\n", paste( interval_chr2, interval_start2, interval_end2, sep = ":"  ) ))
    cat(sprintf("       unit: %s\n", unit       ))
    cat(sprintf("   bin_size: %s\n", bin_size   ))
    
    
    sparMat1 <- straw( norm,
                       path_hic_A,
                       paste( interval_chr1, interval_start1, interval_end1, sep = ":"  ),
                       paste( interval_chr2, interval_start2, interval_end2, sep = ":"  ),
                       unit,
                       bin_size)
    
    if( nrow(sparMat1) == 0 ){
      cat("\n No contact values found for this interval in sample A\n")
      quit()
    }
    
    sparMat2 <- straw( norm,
                       path_hic_B,
                       paste( interval_chr1, interval_start1, interval_end1, sep = ":"  ),
                       paste( interval_chr2, interval_start2, interval_end2, sep = ":"  ),
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
      w <- calculate_w(plot_width, plot_height, A)
    }
    
    if (flag_w != "blank"){
      w <- as.numeric(flag_w)
      max_w <- calculate_w(plot_width, plot_height, A)
      if (w > max_w){
        cat(sprintf("\n\nUser-supplied -w value is too large for the given range and/or resolution.\nThe maximum -w acceptable for these parameters is %f\n\n", max_w))
        w <- max_w
      }
    }
    
    if (w > 2 && bin_size >= 5000){cat(underplot_message)}
    if (w < 0.26 && bin_size <= 5000){cat(overplot_message_1)}
    if (w < 0.26 && bin_size > 5000){cat(overplot_message_2)}
    
    cat(sprintf("plotting bin width: %s\n", round(w,3)))
    top        <- 140
    
    draw_title_inter(interval_chr1, interval_start1, interval_end1, interval_chr2, interval_start2, interval_end2)
    
    top <- top - 12
    
    draw_scale( )
    
    top <- top - 30
    
    draw_contacts_inter(A, C, left_side = 0, bedpe, plot_width, plot_height, top)
    
    if (isFALSE(flag_profiles)) {
      # Initialize legend elements
      legend_labels <- c()
      legend_colors <- c()
      
      # Default Annotations
      if (ann_default) {
        if (genome_build == "hg19") {
          cat("\nDrawing default annotations for hg19\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/GENCODE_TSSs_hg19.bed", sep = "")
          
          # Draw annotations for x- and y-axes
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          draw_bed_inter(bed_tss, "#ef0000", -0.2, TRUE, w, top, "y", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          
          legend_labels <- c(legend_labels, "transcription start sites")
          legend_colors <- c(legend_colors, "#ef0000")
        }
        if (genome_build == "hg38") {
          cat("\nDrawing default annotations for hg38\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/GENCODE_TSSs_hg38.bed", sep = "")
          
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "y", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          
          legend_labels <- c(legend_labels, "transcription start sites")
          legend_colors <- c(legend_colors, "#ef0000")
        }
        if (genome_build == "mm10") {
          cat("\nDrawing default annotations for mm10\n")
          bed_tss <- paste(data_dir, "/", genome_build, "/reference/GENCODE_TSSs_200bp_mm10.bed", sep = "")
          
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "x", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          draw_bed_inter(bed_tss, "#ef0000", -0.1, TRUE, w, top, "y", 
                         interval_start1, interval_len1, interval_start2, interval_len2, 
                         plot_width, plot_height)
          
          legend_labels <- c(legend_labels, "transcription start sites")
          legend_colors <- c(legend_colors, "#ef0000")
        }
      } else {
        no_def_ann <- TRUE
      }
      
      # Custom Annotations
      if (ann_custom != "NONE") {
        # Split the annotation files string into a vector
        annotation_files <- unlist(strsplit(ann_custom, " "))
        
        # Define default colors and track depths
        default_colors <- c("#ef0000", "#00829d", "#8622b7", "#4C8219")
        track_depths <- c(-0.5, -0.8, -1.1, -1.4)
        
        for (i in seq_along(annotation_files)) {
          current_file <- annotation_files[i]
          
          # Decide which color to use
          if (ann_custom_colors == "NONE") {
            this_color <- default_colors[(i - 1) %% length(default_colors) + 1]
            cat(sprintf("\nUsing default color %s for custom annotation file %s", this_color, current_file))
          } else {
            # Split the provided colors and ensure a '#' is included
            colors <- unlist(strsplit(ann_custom_colors, " "))
            if (length(colors) >= i) {
              this_color <- if (substr(colors[i], 1, 1) == "#") colors[i] else paste0("#", colors[i])
            } else {
              this_color <- default_colors[(i - 1) %% length(default_colors) + 1]
            }
            cat(sprintf("\nDrawing custom annotations for file %s with color %s", current_file, this_color))
          }
          
          # Draw annotations for both x and y axes
          draw_bed_inter(current_file, this_color, track_depths[(i - 1) %% length(track_depths) + 1], FALSE, w, top, "x",
                         interval_start1, interval_len1, interval_start2, interval_len2,
                         plot_width, plot_height)
          draw_bed_inter(current_file, this_color, track_depths[(i - 1) %% length(track_depths) + 1], FALSE, w, top, "y",
                         interval_start1, interval_len1, interval_start2, interval_len2,
                         plot_width, plot_height)
          
          # Update legend elements using the base name of the file
          legend_labels <- c(legend_labels, basename(current_file))
          legend_colors <- c(legend_colors, this_color)
        }
      } else {
        no_cus_ann <- TRUE
      }
      
      # Draw the legend if there are any annotation labels
      if (length(legend_labels) > 0) {
        legend(-2, 127, legend = legend_labels, col = legend_colors, lwd = 1,
               box.lty = 0, cex = 0.8, text.col = "#888888")
      }
    }
    
    if (flag_profiles) {
      draw_profiles_inter()
    }
    
    device_off <- dev.off()
  }
}
