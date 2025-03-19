options(stringsAsFactors = FALSE)

################################################################################
#                                 Functions                                    #
################################################################################


read_file_or_return_empty <- function(file_path) {
  if (file.exists(file_path) && file.info(file_path)$size > 0) {
    lines <- readLines(file_path)
    max_num_cols <- max(sapply(strsplit(lines, "\t"), length))
    data <- lapply(lines, function(line) {
      fields <- unlist(strsplit(line, "\t"))
      length(fields) <- max_num_cols
      fields[fields == ""] <- "."
      fields[is.na(fields)] <- "."
      return(fields)
    })
    data <- do.call(rbind, data)
    return(data.frame(data, stringsAsFactors = FALSE))
  } else {
    return(data.frame())
  }
}


sort_chromosomes <- function(df) {
  return(order(df[,1], as.numeric(df[,2])))
}


check_chr_prefix <- function(bedpe_data) {
  valid_chr_1 <- all(grepl("^chr", bedpe_data[[1]]))
  valid_chr_4 <- all(grepl("^chr", bedpe_data[[4]]))
  
  if (!valid_chr_1 || !valid_chr_4) {
    cat("\nThe 1st and/or 4th columns do not start with 'chr'\n")
    quit(save="no")
  }
}


process_display_mode <- function(data, original_bedpe_list, flag_mode, print_env) {
  
  process_row <- function(i) {
    if (flag_mode == "coordinate") {
      id <- as.numeric(data$ID[i])
      subset_row_str <- original_bedpe_list[id]
      rest_of_row <- data[i, (start_col_index:ncol(data)),]
      
      # make new column with IDs (A-B)
      bed_ID_x_last_char <- substr(data$bed_ID.x[i], nchar(data$bed_ID.x[i]), nchar(data$bed_ID.x[i]))
      bed_ID_y_last_char <- substr(data$bed_ID.y[i], nchar(data$bed_ID.y[i]), nchar(data$bed_ID.y[i]))
      new_column_value <- paste0(bed_ID_x_last_char, "-", bed_ID_y_last_char)
      rest_of_row <- cbind(rest_of_row, new_column_value)
      
      cols_to_remove <- c("ID", "bed_ID.x", "bed_ID.y", "chr.y", "start.y", "end.y")
      row <- lapply(c(subset_row_str, rest_of_row), as.character)
      row <- row[!names(row) %in% cols_to_remove]
      row_key <- do.call(paste, c(row, sep = "\t"))
      if (!row_key %in% print_env$printed_rows) {
        write(row_key, stdout(), append = TRUE)
        print_env$printed_rows <- c(print_env$printed_rows, row_key)
      }
    } else if (flag_mode == "bed_ID") {
      id <- as.numeric(data$ID[i])
      bed_ID_x <- as.character(data$bed_ID.x[i])
      bed_ID_y <- as.character(data$bed_ID.y[i])
      subset_row_str <- original_bedpe_list[id]
      row_key <- paste(subset_row_str, bed_ID_x, bed_ID_y, sep = "\t")
      if (!row_key %in% print_env$printed_rows) {
        write(row_key, stdout(), append = TRUE)
        print_env$printed_rows <- c(print_env$printed_rows, row_key)
      }
    }
  }
  
  invisible(lapply(1:nrow(data), process_row))
}


################################################################################
#                               Arguments                                      #
################################################################################

args <- commandArgs( trailingOnly = TRUE )

path_foot_1     <- args[1]
path_foot_2     <- args[2]
path_bedpe      <- args[3]
print_bed       <- as.logical(args[4])
print_bool      <- as.logical(args[5])
absence         <- as.logical(args[6])
flag_mode       <- args[7]
two_bed         <- as.logical(args[8])

if (flag_mode == "blank"){
  flag_mode <- "coordinate"
}

start_col_index <- 6

################################################################################
#                               Initializing                                   #
################################################################################

foot1 <- read_file_or_return_empty(path_foot_1)
foot2 <- read_file_or_return_empty(path_foot_2)


if (nrow(foot1) == 0 && nrow(foot2) == 0) {
  if (!absence){
    cat("No intersections found on either foot of the pair.\n")
  } else if (absence){
    cat("No absences found on either foot of the pair.\n")
  }
  quit(save = "no")
}

original_bedpe <-  read_file_or_return_empty(path_bedpe)

check_chr_prefix(original_bedpe)

if (absence && print_bed) {
  column_names <- c("chr", "start", "end")
} else if (absence || print_bool) {
  column_names <- c("chr", "start", "end", "ID")
} else {
  if(!two_bed){
    if (nrow(foot1) > 0) {
      foot1$bed_ID <- "bed_A" 
      current_cols <- names(foot1)
      new_order <- c(current_cols[1:4], "bed_ID", current_cols[start_col_index:ncol(foot1)-1])
      foot1 <- foot1[, new_order]
      column_names <- c("chr", "start", "end", "ID", "bed_ID")
    }
    if (nrow(foot2) > 0) {
      foot2$bed_ID <- "bed_A" 
      current_cols <- names(foot2)
      new_order <- c(current_cols[1:4], "bed_ID", current_cols[start_col_index:ncol(foot2)-1])
      foot2 <- foot2[, new_order]
      column_names <- c("chr", "start", "end", "ID", "bed_ID")
    }
  } else {
    column_names <- c("chr", "start", "end", "ID", "bed_ID")
  }
}


if (nrow(foot1) > 0) {
  colnames(foot1) <- column_names
}

if (nrow(foot2) > 0) {
  colnames(foot2) <- column_names
}

if (!absence && !print_bool) {
  if (ncol(foot1) > length(column_names) || ncol(foot2) > length(column_names)){
    
    max_cols <- max(ncol(foot1), ncol(foot2))
    
    # Pad foot1 and foot2 with NA columns if necessary
    if (nrow(foot1) > 0) {
      if (ncol(foot1) < max_cols) {
        foot1[, (ncol(foot1) + 1):max_cols] <- "."
      }
    }
    if (nrow(foot2) > 0) {
      if (ncol(foot2) < max_cols) {
        foot2[, (ncol(foot2) + 1):max_cols] <- "."
      }
    }
    
    total_cols <- max_cols - length(column_names)
    bed_x_col_names <- paste("bed", seq_len(total_cols), "x", sep = "_")
    bed_y_col_names <- paste("bed", seq_len(total_cols), "y", sep = "_")
    foot1_names <- c(column_names, bed_x_col_names)
    foot2_names <- c(column_names, bed_y_col_names)
    
    if (nrow(foot1) > 0) {
      colnames(foot1) <- foot1_names
    }
    if (nrow(foot2) > 0) {
      colnames(foot2) <- foot2_names
    }
  }
}


################################################################################
#                               Code Block                                     #
################################################################################

# absence block
if (isTRUE(absence)) {
  if (nrow(foot1) > 0 && nrow(foot2) > 0) {
    if (isFALSE(print_bed)) {
      non_intersect_ids <- intersect(foot1$ID, foot2$ID)
      if (length(non_intersect_ids) > 0) {
        invisible(apply(original_bedpe[non_intersect_ids, ], 1, function(x) cat(paste(x, collapse = "\t"), "\n", sep = "")))
      } else {
        cat("\nThere are no non-intersecting rows in the BEDPE file.\n")
        quit(save="no")
      }
    } else if (isTRUE(print_bed)) {
      
      row_to_string <- function(row) {
        paste(row, collapse = "\t")
      }
      
      # Convert each dataframe to a character vector
      foot1_strings <- apply(foot1, 1, row_to_string)
      foot2_strings <- apply(foot2, 1, row_to_string)
      
      # Find matching rows
      matching_indices <- match(foot1_strings, foot2_strings, nomatch = 0)
      
      # Get the matching rows from foot1
      common_rows <- foot1[matching_indices > 0, , drop = FALSE]
      
      if (nrow(common_rows) > 0) {
        invisible(apply(common_rows, 1, function(x) {
        write(paste(x, collapse = "\t"), stdout(), append = TRUE)}))
      } else {
        write("There are no non-intersecting BED regions.", stdout())
        quit(save="no")
      }
    }
  } else {
    cat("\nThere are no non-intersecting rows.\n")
    quit(save="no")
  }
}

# intersection block
if (isFALSE(absence)) {
  if (isTRUE(print_bed)) {
    if (nrow(foot1) > 0) {
      bed_columns_foot1 <- foot1[, (start_col_index:ncol(foot1))]
      unique_bed_columns_foot1 <- unique(bed_columns_foot1)
      invisible(apply(unique_bed_columns_foot1, 1, function(x) {
        write(paste(x, collapse = "\t"), stdout(), append = TRUE)
      }))
    }
    if (nrow(foot2) > 0) {
      bed_columns_foot2 <- foot2[, (start_col_index:ncol(foot2))]
      unique_bed_columns_foot2 <- unique(bed_columns_foot2)
      invisible(apply(unique_bed_columns_foot2, 1, function(x) {
        write(paste(x, collapse = "\t"), stdout(), append = TRUE)
      }))
    }
    
  } else if (isTRUE(print_bool)) {
    
    original_bedpe$foot1_X <- FALSE 
    original_bedpe$foot2_X <- FALSE 
    original_bedpe[foot1$ID, "foot1_X"] <- TRUE  
    original_bedpe[foot2$ID, "foot2_X"] <- TRUE 
    write.table(original_bedpe, file = "", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  } else {
    
    # display_mode block
    
    if (nrow(foot1) > 0) {
      foot1 <- foot1[sort_chromosomes(foot1), ]
    }
    if (nrow(foot2) > 0) {
      foot2 <- foot2[sort_chromosomes(foot2), ]
    }
    
    if (nrow(foot1) > 0 && nrow(foot2) > 0) {
      intersections <- suppressWarnings(merge(foot1, foot2, by = "ID", all = TRUE))
    } else if (nrow(foot1) > 0) {
      intersections <- foot1
    } else if (nrow(foot2) > 0) {
      intersections <- foot2
    }
    
    intersections  <- replace(intersections, is.na(intersections), ".")
    intersections <- intersections[order(as.numeric(intersections$ID)), ]
    
    # list for faster subsetting
    original_bedpe_list <- apply(original_bedpe, 1, function(row) paste(row, collapse = "\t"))
    
    if (nrow(foot1) > 0 && nrow(foot2) > 0) {
      chunk_size <- 10000  
      num_chunks <- ceiling(nrow(intersections) / chunk_size)  
      
      for (j in 1:num_chunks) {
        
        start_row <- (j - 1) * chunk_size + 1
        end_row <- min(j * chunk_size, nrow(intersections))
        
        int_subset <- intersections[start_row:end_row, ]
        
        # env for storing row keys to unique on-the-fly
        print_env <- new.env()
        print_env$printed_rows <- character(0)
        
        process_display_mode(int_subset, original_bedpe_list, flag_mode, print_env)
      }
    } else {
      # cases where only one foot has intersections
      if (flag_mode == "bed_ID") {
        for (i in 1:nrow(intersections)) {
          id <- as.numeric(intersections$ID[i])
          bed_ID <- as.character(intersections$bed_ID[i])
          subset_row_str <- original_bedpe_list[id]
          
          # Check if foot1 or foot2 has data
          if (nrow(foot1) > 0 && nrow(foot2) == 0) {
            # foot1 has data, dot goes on the right
            row_key <- paste(subset_row_str, bed_ID, ".", sep = "\t")
          } else if (nrow(foot1) == 0 && nrow(foot2) > 0) {
            # foot2 has data, dot goes on the left
            row_key <- paste(subset_row_str, ".", bed_ID, sep = "\t")
          }
          write(row_key, stdout(), append = TRUE)
        }
      } else {  # coordinate mode
        for (i in 1:nrow(intersections)) {
          id <- as.numeric(intersections$ID[i])
          subset_row_str <- original_bedpe_list[id]
          bed_ID <- as.character(intersections$bed_ID[i])
          bed_ID_last_char <- substr(bed_ID, nchar(bed_ID), nchar(bed_ID))
          
          rest_of_row <- intersections[i, 4:ncol(intersections), drop = FALSE]
          if (nrow(foot1) > 0 && nrow(foot2) == 0) {
            # foot1 has data, append dot to the right of the letter
            suffix <- paste0(bed_ID_last_char, "-.")
          } else if (nrow(foot1) == 0 && nrow(foot2) > 0) {
            # foot2 has data, append dot to the left of the letter
            suffix <- paste0(".-", bed_ID_last_char)
          } 
          
          rest_of_row <- cbind(rest_of_row, suffix)
          
          cols_to_remove <- c("ID", "bed_ID")
          row <- lapply(c(subset_row_str, rest_of_row), as.character)
          row <- row[!names(row) %in% cols_to_remove]
          row_key <- do.call(paste, c(row, sep = "\t"))
          write(row_key, stdout(), append = TRUE)
        }
      }
    }
  }
}
