suppressPackageStartupMessages( library(strawr  ) )
suppressPackageStartupMessages( library(parallel) )

cores <- detectCores() - 1

options(scipen = 999)

#################################################################
##                          Functions                          ##
#################################################################

prefix_filter_bedpe <- function( bedpe ){
  
  chr_A <- bedpe[1]
  chr_B <- bedpe[4]
  
  if( chr_A != chr_B ){ cat("At least one entry is inter-chromosomal. Only intra-interactions supported") ; quit(save="no") }
  
  start1_bin <- as.numeric(bedpe[7])
  end1_bin   <- as.numeric(bedpe[8])
  
  if( end1_bin > start1_bin ){ end1_bin <- end1_bin - bin_size }
  
  start2_bin <- as.numeric(bedpe[9 ])
  end2_bin   <- as.numeric(bedpe[10])
  
  if( end2_bin > start2_bin ){ end2_bin <- end2_bin - bin_size }
  
  midpoint_1_bin <- as.numeric(bedpe[11])
  midpoint_2_bin <- as.numeric(bedpe[12])
  
  if( end1_bin == start1_bin ){ midpoint_1_bin <- start1_bin }
  if( end2_bin == start2_bin ){ midpoint_2_bin <- start2_bin }
  
  if(flag_strip_chr == "yes"){
    
    chr_L <- sub( "chr", "", chr_A )
    chr_R <- sub( "chr", "", chr_B )
    
    txt_L <- sprintf( "%s:%d:%d", chr_L, start1_bin, end1_bin )
    txt_R <- sprintf( "%s:%d:%d", chr_R, start2_bin, end2_bin )
    
    obj <- list( txt_L = txt_L,
                 txt_R = txt_R,
                 midpoint_1_bin = midpoint_1_bin,
                 midpoint_2_bin = midpoint_2_bin )
  }
  
  if(flag_strip_chr == "no"){
    
    chr_L <- chr_A
    chr_R <- chr_B
    
    txt_L <- sprintf( "%s:%d:%d", chr_L, start1_bin, end1_bin )
    txt_R <- sprintf( "%s:%d:%d", chr_R, start2_bin, end2_bin )
    
    obj <- list( txt_L = txt_L,
                 txt_R = txt_R,
                 midpoint_1_bin = midpoint_1_bin,
                 midpoint_2_bin = midpoint_2_bin ) 
  }  
  
  if(flag_strip_chr == "A"){
    
    sample_A_chr_L <- chr_A
    sample_A_chr_R <- chr_B
    
    sample_A_txt_L <- sprintf( "%s:%d:%d", sample_A_chr_L, start1_bin, end1_bin )
    sample_A_txt_R <- sprintf( "%s:%d:%d", sample_A_chr_R, start2_bin, end2_bin )
    
    sample_B_chr_L <- sub("chr", "", chr_A )
    sample_B_chr_R <- sub("chr", "", chr_B )
    
    sample_B_txt_L <- sprintf( "%s:%d:%d", sample_B_chr_L, start1_bin, end1_bin )
    sample_B_txt_R <- sprintf( "%s:%d:%d", sample_B_chr_R, start2_bin, end2_bin )
    
    obj <- list( sample_A_txt_L = sample_A_txt_L,
                 sample_A_txt_R = sample_A_txt_R,
                 sample_B_txt_L = sample_B_txt_L,
                 sample_B_txt_R = sample_B_txt_R,
                 midpoint_1_bin = midpoint_1_bin,
                 midpoint_2_bin = midpoint_2_bin )
    
  }
  
  if(flag_strip_chr == "B"){
    
    sample_B_chr_L <- chr_A
    sample_B_chr_R <- chr_B
    
    sample_B_txt_L <- sprintf( "%s:%d:%d", sample_B_chr_L, start1_bin, end1_bin )
    sample_B_txt_R <- sprintf( "%s:%d:%d", sample_B_chr_R, start2_bin, end2_bin )
    
    sample_A_chr_L <- sub("chr","", chr_A)
    sample_A_chr_R <- sub("chr","", chr_B)
    
    sample_A_txt_L <- sprintf( "%s:%d:%d", sample_A_chr_L, start1_bin, end1_bin )
    sample_A_txt_R <- sprintf( "%s:%d:%d", sample_A_chr_R, start2_bin, end2_bin )
    
    obj <- list( sample_A_txt_L = sample_A_txt_L,
                 sample_A_txt_R = sample_A_txt_R,
                 sample_B_txt_L = sample_B_txt_L,
                 sample_B_txt_R = sample_B_txt_R,
                 midpoint_1_bin = midpoint_1_bin,
                 midpoint_2_bin = midpoint_2_bin )
    
  }
  
  return( obj )
  
}

obtain_one_sample_counts <- function( mat_1, txt_L, txt_R, midpoint_1, midpoint_2 ){

  if( isTRUE(flag_fix)  && shrink_wrap == "FALSE" && split == "FALSE" ){
    
    if( nrow( mat_1) == 1 ){
      
      count_1 <- mat_1[1,"counts"] 
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
    } else if( nrow( mat_1 ) > 1 ){ 
      
      dim_1 <- seq( 
        as.numeric(unlist(strsplit(txt_L,":"))[2]),
        as.numeric(unlist(strsplit(txt_L,":"))[3]),
        bin_size )
      
      dim_2 <- seq( 
        as.numeric(unlist(strsplit(txt_R,":"))[2]),
        as.numeric(unlist(strsplit(txt_R,":"))[3]),
        bin_size )
      
      sparse_mat_1 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_1 ) <- dim_1
      rownames( sparse_mat_1 ) <- dim_2
      
      for( i in 1:nrow(mat_1) ){
        
        row   <- mat_1[ i , "x"      ]
        col   <- mat_1[ i , "y"      ]
        count <- mat_1[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_1 ) && 
            as.character(col) %in% colnames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_1 ) && 
                   as.character(col) %in% rownames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_1[ is.na(sparse_mat_1) ] <- 0
      }
      
      if( formula == "center"){
        
        count_1  <- sparse_mat_1[ as.character(midpoint_2) , as.character(midpoint_1) ]
        
      } else if( formula == "max"){
        
        count_1  <- max( sparse_mat_1, na.rm = T )
        
      } else if( formula == "sum"){
        
        count_1 <- sum( sparse_mat_1, na.rm = T )
        
      } else if( formula == "average"){
        
        count_1 <- mean( sparse_mat_1, na.rm = T )
      }
      
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
    } else if ( nrow( mat_1 ) == 0 ){ 
      
      count_1 <- 0
      
    } else { 
      
      count_1 <- "*"
      
    }
    
    return(count_1)
    
  }
  
  if( isFALSE(flag_fix) && shrink_wrap == "FALSE" && split == "FALSE" ){
    
    if( nrow( mat_1) == 1 ){
      
      count_1 <- mat_1[1,"counts"] 
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
      new_start1 <- as.numeric(mat_1[ 1 , "x" ])
      new_end1   <- new_start1 + bin_size
      
      new_start2 <- as.numeric(mat_1[ 1 , "y" ])
      new_end2   <- new_start2 + bin_size
      
      return_obj <- list( 
        feet1  = c(new_start1, new_end1),
        feet2  = c(new_start2, new_end2),
        count1 = count_1 )
      
    } else if( nrow( mat_1 ) > 1 ){ 
      
      dim_1 <- seq(          
        as.numeric(unlist(strsplit(txt_L,":"))[2]),         
        as.numeric(unlist(strsplit(txt_L,":"))[3]),         
        bin_size )
      dim_2 <- seq(          
        as.numeric(unlist(strsplit(txt_R,":"))[2]),         
        as.numeric(unlist(strsplit(txt_R,":"))[3]),         
        bin_size )
      
      sparse_mat_1 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_1 ) <- dim_1
      rownames( sparse_mat_1 ) <- dim_2
      
      for( i in 1:nrow(mat_1) ){
        
        row   <- mat_1[ i , "x"      ]
        col   <- mat_1[ i , "y"      ]
        count <- mat_1[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_1 ) && 
            as.character(col) %in% colnames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_1 ) && 
                   as.character(col) %in% rownames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_1[ is.na(sparse_mat_1) ] <- 0
      }
      
      if( formula == "center"){
        
        count_1    <- sparse_mat_1[ as.character(midpoint_2) , as.character(midpoint_1) ]
        count_1    <- count_1 * norm_factor1 * aqua_factor1
        count_1    <- round( count_1 , 3 )
        
        new_start1 <- as.numeric(midpoint_2)
        new_end1   <- new_start1 + bin_size
        
        new_start2 <- as.numeric(midpoint_1)
        new_end2   <- new_start2 + bin_size
        
        return_obj <- list( 
          feet1  = c(new_start1, new_end1),
          feet2  = c(new_start2, new_end2),
          count1 = count_1 )
        
      } else if( formula == "max"){
        
        count_1    <- max( sparse_mat_1, na.rm = T )
        count_1    <- count_1 * norm_factor1 * aqua_factor1
        count_1    <- round( count_1 , 3 )
        
        new_start1 <- rownames(sparse_mat_1)[ as.integer(which( max(sparse_mat_1) == sparse_mat_1, arr.ind = T )[1,1]) ]
        new_start1 <- as.numeric(new_start1)
        new_end1   <- new_start1 + bin_size
        
        new_start2 <- colnames(sparse_mat_1)[ as.integer(which( max(sparse_mat_1) == sparse_mat_1, arr.ind = T )[1,2]) ]
        new_start2 <- as.numeric(new_start2)
        new_end2   <- new_start2 + bin_size
        
        return_obj <- list( 
          feet1  = c(new_start1, new_end1),
          feet2  = c(new_start2, new_end2),
          count1 = count_1 )
        
      } else if( formula == "sum"){
        
        count_1    <- sum( sparse_mat_1, na.rm = T )
        count_1    <- count_1 * norm_factor1 * aqua_factor1
        count_1    <- round( count_1 , 3 )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size 
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size 
          
        return_obj <- list( 
          feet1  = c(new_start1, new_end1),
          feet2  = c(new_start2, new_end2),
          count1 = count_1 )
        
      } else if( formula == "average"){
        
        count_1    <- mean( sparse_mat_1, na.rm = T )
        count_1    <- count_1 * norm_factor1 * aqua_factor1
        count_1    <- round( count_1 , 3 )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size 
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
          
        return_obj <- list( 
          feet1  = c(new_start1, new_end1),
          feet2  = c(new_start2, new_end2),
          count1 = count_1 )
        
      }
      
    } else if ( nrow( mat_1 ) == 0 ){ 
      
      count_1 <- 0
      
      new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
      new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
      
      new_end1 <- new_end1 + bin_size
      
      new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
      new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
      
      new_end2 <- new_end2 + bin_size
      
      return_obj <- list( 
        feet1  = c(new_start1, new_end1),
        feet2  = c(new_start2, new_end2),
        count1 = count_1 )
      
    } else { 
      
      count_1 <- "*"
      
      new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
      new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
      
      new_end1 <- new_end1 + bin_size
      
      new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
      new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
      
      new_end2 <- new_end2 + bin_size
      
      return_obj <- list( 
        feet1  = c(new_start1, new_end1),
        feet2  = c(new_start2, new_end2),
        count1 = count_1 )
      
    }
    
    return(return_obj)
    
  }
  
  if( shrink_wrap != "FALSE" && split == "FALSE" ){
    
    shrink_wrap <- as.numeric(shrink_wrap)
    
    if( nrow( mat_1) == 1 ){
      
      count_1 <- mat_1[1,"counts"] 
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
      new_start1 <- as.numeric(mat_1[ 1 , "x" ])
      new_end1   <- new_start1 + bin_size
      
      new_start2 <- as.numeric(mat_1[ 1 , "y" ])
      new_end2   <- new_start2 + bin_size
      
      return_obj <- list( 
        feet1  = c(new_start1, new_end1),
        feet2  = c(new_start2, new_end2),
        count1 = count_1 )
      
    } else if( nrow( mat_1) > 1 ){
      
      dim_1 <- seq(          
        as.numeric(unlist(strsplit(txt_L,":"))[2]),         
        as.numeric(unlist(strsplit(txt_L,":"))[3]),         
        bin_size )
      dim_2 <- seq(          
        as.numeric(unlist(strsplit(txt_R,":"))[2]),         
        as.numeric(unlist(strsplit(txt_R,":"))[3]),         
        bin_size )
      
      sparse_mat_1 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_1 ) <- dim_1
      rownames( sparse_mat_1 ) <- dim_2
      
      for( i in 1:nrow(mat_1) ){
        
        row   <- mat_1[ i , "x"      ]
        col   <- mat_1[ i , "y"      ]
        count <- mat_1[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_1 ) && 
            as.character(col) %in% colnames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_1 ) && 
                   as.character(col) %in% rownames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(col) , as.character(row) ] <- count
          
        }
      }
      
      sparse_mat_1[ is.na(sparse_mat_1) ] <- 0
      
      if( rownames(sparse_mat_1)[1] > colnames(sparse_mat_1)[1] ){
        sparse_mat_1 <- t(sparse_mat_1)
      }
      
       if( txt_L == txt_R ){
        sparse_mat_1 <- zero_diag(sparse_mat_1, 2)
      }
      
      tuple_col <- apply( sparse_mat_1, 2, function(x){if(max(x) >= shrink_wrap){return(TRUE)}else{return(FALSE)}})
      tuple_row <- apply( sparse_mat_1, 1, function(x){if(max(x) >= shrink_wrap){return(TRUE)}else{return(FALSE)}})
      
      idx_row <- min(which(tuple_row)):max(which(tuple_row))
      idx_col <- min(which(tuple_col)):max(which(tuple_col))
      
      if( txt_L == txt_R ){
        idx_row <- min(c(idx_row,idx_col)):max(c(idx_row,idx_col))
        idx_col <- min(c(idx_row,idx_col)):max(c(idx_row,idx_col))
      }
      
      sparse_mat_1 <- sparse_mat_1[ idx_row , idx_col, drop = FALSE ] 
      
      new_start1 <- as.numeric(rownames(sparse_mat_1)[1]) ; new_end1 <- as.numeric(rownames(sparse_mat_1)[length(rownames(sparse_mat_1))]) + bin_size
      new_start2 <- as.numeric(colnames(sparse_mat_1)[1]) ; new_end2 <- as.numeric(colnames(sparse_mat_1)[length(colnames(sparse_mat_1))]) + bin_size
      
      if( formula == "center"){
        
        midpoint_2 <- (as.numeric(new_start1) + as.numeric(new_end1))/2
        midpoint_2 <- floor( midpoint_2/bin_size )*bin_size
        
        midpoint_1 <- (as.numeric(new_start2) + as.numeric(new_end2))/2
        midpoint_1 <- floor( midpoint_1/bin_size )*bin_size
        
        count_1    <- sparse_mat_1[ as.character(midpoint_2) , as.character(midpoint_1) ]
        count_1    <- count_1 * norm_factor1 * aqua_factor1
        count_1    <- round( count_1 , 3 )

      } else if( formula == "max"){
        
        count_1    <- max( sparse_mat_1, na.rm = T )
        count_1    <- count_1 * norm_factor1 * aqua_factor1
        count_1    <- round( count_1 , 3 )
      
      } else if( formula == "sum"){
        
        count_1    <- sum( sparse_mat_1, na.rm = T )
        count_1    <- count_1 * norm_factor1 * aqua_factor1
        count_1    <- round( count_1 , 3 )
        
      } else if( formula == "average"){
        
        count_1    <- mean( sparse_mat_1, na.rm = T )
        count_1    <- count_1 * norm_factor1 * aqua_factor1
        count_1    <- round( count_1 , 3 )
      }
      
      return_obj <- list( 
        feet1  = c(new_start1, new_end1),
        feet2  = c(new_start2, new_end2),
        count1 = count_1 )
      
      
    }
    
    return(return_obj)
  }
  
  if( shrink_wrap == "FALSE" && split != "FALSE" ){
    
    split <- as.numeric(split)
    
    if( nrow( mat_1) == 1 ){
      
      if( flag_strip_chr == "yes"){ chr <- paste( "chr", unlist(strsplit(txt_L,":"))[1], sep="" ) }
      if( flag_strip_chr == "no" ){ chr <-               unlist(strsplit(txt_L,":"))[1]           }
      
      count_1 <- mat_1[1,"counts"] 
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
      new_start1 <- as.numeric(mat_1[ 1 , "x" ])
      new_end1   <- new_start1 + bin_size
      
      new_start2 <- as.numeric(mat_1[ 1 , "y" ])
      new_end2   <- new_start2 + bin_size
      
      return_obj              <- list()
      return_obj[[1]]         <- list()
      return_obj[[1]]$chr     <- chr
      return_obj[[1]]$feet1   <- c(new_start1, new_end1)
      return_obj[[1]]$feet2   <- c(new_start2, new_end2)
      return_obj[[1]]$count_1 <- count_1
      
    } else if( nrow( mat_1) > 1 ){
      
      if( flag_strip_chr == "yes"){ chr <- paste( "chr", unlist(strsplit(txt_L,":"))[1], sep="" ) }
      if( flag_strip_chr == "no" ){ chr <-               unlist(strsplit(txt_L,":"))[1]           }
      
      dim_1 <- seq(          
        as.numeric(unlist(strsplit(txt_L,":"))[2]),         
        as.numeric(unlist(strsplit(txt_L,":"))[3]),         
        bin_size )
      dim_2 <- seq(          
        as.numeric(unlist(strsplit(txt_R,":"))[2]),         
        as.numeric(unlist(strsplit(txt_R,":"))[3]),         
        bin_size )
      
      sparse_mat_1 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_1 ) <- dim_1
      rownames( sparse_mat_1 ) <- dim_2
      
      for( i in 1:nrow(mat_1) ){
        
        row   <- mat_1[ i , "x"      ]
        col   <- mat_1[ i , "y"      ]
        count <- mat_1[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_1 ) && 
            as.character(col) %in% colnames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_1 ) && 
            as.character(col) %in% rownames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(col) , as.character(row) ] <- count
          
        }
      }
      
      sparse_mat_1[ is.na(sparse_mat_1) ] <- 0
      
      if( rownames(sparse_mat_1)[1] > colnames(sparse_mat_1)[1] ){
        sparse_mat_1 <- t(sparse_mat_1)
      }
      
      if( txt_L == txt_R ){
        sparse_mat_1 <- zero_diag( sparse_mat_1, 2 )
      }
      
      return_obj <- get_clusters( sparse_mat_1, split, padding, chr, bin_size, formula, norm_factor1, aqua_factor1 )
      
    }
    
    return(return_obj)
  }
  
}

call_one_sample_straw_bedpe <- function( bedpe_list, flag_fix, formula, norm_factor1, aqua_factor1, shrink_wrap, split ){
  
  txt_L <- as.character( bedpe_list[1]$txt_L )
  txt_R <- as.character( bedpe_list[2]$txt_R )
  
  midpoint_1 <- as.character( bedpe_list[3]$midpoint_1_bin )
  midpoint_2 <- as.character( bedpe_list[4]$midpoint_2_bin )
  
  mat_1 <- straw( norm, hic_A, txt_L, txt_R, "BP", bin_size )
  
  obtain_one_sample_counts( mat_1, txt_L, txt_R, midpoint_1, midpoint_2 )
}

obtain_two_sample_counts <- function( mat_1, mat_2, txt_L, txt_R, midpoint_1, midpoint_2 ){

  if( isTRUE(flag_fix) ){
    
    if( nrow( mat_1 ) == 1 && nrow( mat_2 ) == 1 ){ 
      
      count_1 <- mat_1[1,"counts"] 
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
      count_2 <- mat_2[1,"counts"] 
      count_2 <- count_2 * norm_factor2 * aqua_factor2
      count_2 <- round( count_2 , 3 )
      
    } else if ( nrow( mat_1 ) >= 1 && nrow( mat_2 ) == 0 ){ 
      
      dim_1 <- seq(         
        as.numeric(unlist(strsplit(txt_L,":"))[2]),         
        as.numeric(unlist(strsplit(txt_L,":"))[3]),         
        bin_size )
      dim_2 <- seq(          
        as.numeric(unlist(strsplit(txt_R,":"))[2]),         
        as.numeric(unlist(strsplit(txt_R,":"))[3]),         
        bin_size )
      
      sparse_mat_1 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_1 ) <- dim_1
      rownames( sparse_mat_1 ) <- dim_2
      
      for( i in 1:nrow(mat_1) ){
        
        row   <- mat_1[ i , "x"      ]
        col   <- mat_1[ i , "y"      ]
        count <- mat_1[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_1 ) && 
            as.character(col) %in% colnames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_1 ) && 
                   as.character(col) %in% rownames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_1[ is.na(sparse_mat_1) ] <- 0
      }
      
      if( formula == "center"){
        
        count_1    <- sparse_mat_1[ as.character(midpoint_2) , as.character(midpoint_1) ]
        
      } else if( formula == "max"){
        
        count_1    <- max( sparse_mat_1, na.rm = T )
        
        
      } else if( formula == "sum"){
        
        count_1 <- sum( sparse_mat_1, na.rm = T )
        
        
      } else if( formula == "average"){
        
        count_1 <- mean( sparse_mat_1, na.rm = T )
        
      }
      
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
      count_2 <- 0
      
    } else if ( nrow( mat_1 ) == 0 && nrow( mat_2 ) >= 1 ){ 
      
      count_1 <- 0
      
      dim_1 <- seq(                     
        as.numeric(unlist(strsplit(txt_L,":"))[2]),                    
        as.numeric(unlist(strsplit(txt_L,":"))[3]),                    
        bin_size )
      dim_2 <- seq(                     
        as.numeric(unlist(strsplit(txt_R,":"))[2]),                    
        as.numeric(unlist(strsplit(txt_R,":"))[3]),                    
        bin_size )
      
      sparse_mat_2 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_2 ) <- dim_1
      rownames( sparse_mat_2 ) <- dim_2
      
      for( i in 1:nrow(mat_2) ){
        
        row   <- mat_2[ i , "x"      ]
        col   <- mat_2[ i , "y"      ]
        count <- mat_2[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_2 ) && 
            as.character(col) %in% colnames( sparse_mat_2 )) {
          
          sparse_mat_2[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_2 ) && 
                   as.character(col) %in% rownames( sparse_mat_2 )) {
          
          sparse_mat_2[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_2[ is.na(sparse_mat_2) ] <- 0
      }
      
      if( formula == "center"){
        
        count_2    <- sparse_mat_2[ as.character(midpoint_2) , as.character(midpoint_1) ]
        
      } else if( formula == "max"){
        
        count_2    <- max( sparse_mat_2, na.rm = T )
        
        
      } else if( formula == "sum"){
        
        count_2 <- sum( sparse_mat_2, na.rm = T )
        
        
      } else if( formula == "average"){
        
        count_2 <- mean( sparse_mat_2, na.rm = T )
        
      }
      
      count_2 <- count_2 * norm_factor2 * aqua_factor2
      count_2 <- round( count_2 , 3 )
      
    } else if ( nrow( mat_1 ) > 1  && nrow( mat_2 ) >  1 ){
      
      dim_1 <- seq(          
        as.numeric(unlist(strsplit(txt_L,":"))[2]),         
        as.numeric(unlist(strsplit(txt_L,":"))[3]),         
        bin_size )
      dim_2 <- seq(          
        as.numeric(unlist(strsplit(txt_R,":"))[2]),         
        as.numeric(unlist(strsplit(txt_R,":"))[3]),         
        bin_size )
      
      sparse_mat_1 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_1 ) <- dim_1
      rownames( sparse_mat_1 ) <- dim_2
      
      for( i in 1:nrow(mat_1) ){
        
        row   <- mat_1[ i , "x"      ]
        col   <- mat_1[ i , "y"      ]
        count <- mat_1[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_1 ) && 
            as.character(col) %in% colnames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_1 ) && 
                   as.character(col) %in% rownames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_1[ is.na(sparse_mat_1) ] <- 0
      }
      if( formula == "center"){
        
        count_1    <- sparse_mat_1[ as.character(midpoint_2) , as.character(midpoint_1) ]
        
      } else if( formula == "max"){
        
        count_1    <- max( sparse_mat_1, na.rm = T )
        
        
      } else if( formula == "sum"){
        
        count_1 <- sum( sparse_mat_1, na.rm = T )
        
        
      } else if( formula == "average"){
        
        count_1 <- mean( sparse_mat_1, na.rm = T )
        
      }
      
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
      dim_1 <- seq(                     
        as.numeric(unlist(strsplit(txt_L,":"))[2]),                    
        as.numeric(unlist(strsplit(txt_L,":"))[3]),                    
        bin_size )
      dim_2 <- seq(                     
        as.numeric(unlist(strsplit(txt_R,":"))[2]),                    
        as.numeric(unlist(strsplit(txt_R,":"))[3]),                    
        bin_size )
      
      sparse_mat_2 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_2 ) <- dim_1
      rownames( sparse_mat_2 ) <- dim_2
      
      for( i in 1:nrow(mat_2) ){
        
        row   <- mat_2[ i , "x"      ]
        col   <- mat_2[ i , "y"      ]
        count <- mat_2[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_2 ) && 
            as.character(col) %in% colnames( sparse_mat_2 )) {
          
          sparse_mat_2[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_2 ) && 
                   as.character(col) %in% rownames( sparse_mat_2 )) {
          
          sparse_mat_2[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_2[ is.na(sparse_mat_2) ] <- 0
      }
      
      if( formula == "center"){
        
        count_2    <- sparse_mat_2[ as.character(midpoint_2) , as.character(midpoint_1) ]
        
      } else if( formula == "max"){
        
        count_2    <- max( sparse_mat_2, na.rm = T )
        
        
      } else if( formula == "sum"){
        
        count_2 <- sum( sparse_mat_2, na.rm = T )
        
        
      } else if( formula == "average"){
        
        count_2 <- mean( sparse_mat_2, na.rm = T )
        
      }
      
      count_2 <- count_2 * norm_factor2 * aqua_factor2
      count_2 <- round( count_2 , 3 )
      
    } else if ( nrow( mat_1 ) == 0 && nrow( mat_2 ) == 0 ){ 
      
      count_1 <- 0
      count_2 <- 0
      
    } else { 
      
      count_1 <- "*"
      count_2 <- "*"
      
    }
    
    obj <- list( 
      count_1 = count_1,
      count_2 = count_2 )
    
  }
  
  if( isFALSE(flag_fix) ){
    
    if( nrow( mat_1 ) == 1 && nrow( mat_2 ) == 1 ){ 
      
      count_1 <- mat_1[1,"counts"] 
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
      count_2 <- mat_2[1,"counts"] 
      count_2 <- count_2 * norm_factor2 * aqua_factor2
      count_2 <- round( count_2 , 3 )
      
      new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
      new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
      
      new_end1 <- new_end1 + bin_size 
      
      new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
      new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
      
      new_end2 <- new_end2 + bin_size 
      
    } 
    else if ( nrow( mat_1 ) >= 1 && nrow( mat_2 ) == 0 ){ 
      
      dim_1 <- seq(          
        as.numeric(unlist(strsplit(txt_L,":"))[2]),         
        as.numeric(unlist(strsplit(txt_L,":"))[3]),         
        bin_size )
      dim_2 <- seq(          
        as.numeric(unlist(strsplit(txt_R,":"))[2]),         
        as.numeric(unlist(strsplit(txt_R,":"))[3]),         
        bin_size )
      
      sparse_mat_1 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_1 ) <- dim_1
      rownames( sparse_mat_1 ) <- dim_2
      
      for( i in 1:nrow(mat_1) ){
        
        row   <- mat_1[ i , "x"      ]
        col   <- mat_1[ i , "y"      ]
        count <- mat_1[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_1 ) && 
            as.character(col) %in% colnames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_1 ) && 
                   as.character(col) %in% rownames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_1[ is.na(sparse_mat_1) ] <- 0
      }
      
      if( formula == "center"){
        
        count_1    <- sparse_mat_1[ as.character(midpoint_2) , as.character(midpoint_1) ]
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size 
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size 
        
      } else if( formula == "max"){
        
        count_1    <- max( sparse_mat_1, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size 
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size 
        
      } else if( formula == "sum"){
        
        count_1    <- sum( sparse_mat_1, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size 
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size 
        
      } else if( formula == "average"){
        
        count_1    <- mean( sparse_mat_1, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size 
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size 
        
      }
      
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
      count_2 <- 0
      
    } 
    else if ( nrow( mat_1 ) == 0 && nrow( mat_2 ) >= 1 ){ 
      
      count_1 <- 0
      
      dim_1 <- seq(                     
        as.numeric(unlist(strsplit(txt_L,":"))[2]),                    
        as.numeric(unlist(strsplit(txt_L,":"))[3]),                    
        bin_size )
      dim_2 <- seq(                     
        as.numeric(unlist(strsplit(txt_R,":"))[2]),                    
        as.numeric(unlist(strsplit(txt_R,":"))[3]),                    
        bin_size )
      
      sparse_mat_2 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_2 ) <- dim_1
      rownames( sparse_mat_2 ) <- dim_2
      
      for( i in 1:nrow(mat_2) ){
        
        row   <- mat_2[ i , "x"      ]
        col   <- mat_2[ i , "y"      ]
        count <- mat_2[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_2 ) && 
            as.character(col) %in% colnames( sparse_mat_2 )) {
          
          sparse_mat_2[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_2 ) && 
                   as.character(col) %in% rownames( sparse_mat_2 )) {
          
          sparse_mat_2[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_2[ is.na(sparse_mat_2) ] <- 0
      }
      
      if( formula == "center"){
        
        count_2    <- sparse_mat_2[ midpoint_2 , midpoint_1]
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      } else if( formula == "max"){
        
        count_2    <- max( sparse_mat_2, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size 
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      } else if( formula == "sum"){
        
        count_2    <- sum( sparse_mat_2, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size 
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      } else if( formula == "average"){
        
        count_2    <- mean( sparse_mat_2, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      }
      
      count_2 <- count_2 * norm_factor2 * aqua_factor2
      count_2 <- round( count_2 , 3 )
      
    } 
    else if ( nrow( mat_1 ) > 1  && nrow( mat_2 ) >  1 ){
      
      dim_1 <- seq(          
        as.numeric(unlist(strsplit(txt_L,":"))[2]),         
        as.numeric(unlist(strsplit(txt_L,":"))[3]),         
        bin_size )
      dim_2 <- seq(          
        as.numeric(unlist(strsplit(txt_R,":"))[2]),         
        as.numeric(unlist(strsplit(txt_R,":"))[3]),         
        bin_size )
      
      sparse_mat_1 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_1 ) <- dim_1
      rownames( sparse_mat_1 ) <- dim_2
      
      for( i in 1:nrow(mat_1) ){
        
        row   <- mat_1[ i , "x"      ]
        col   <- mat_1[ i , "y"      ]
        count <- mat_1[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_1 ) && 
            as.character(col) %in% colnames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_1 ) && 
                   as.character(col) %in% rownames( sparse_mat_1 )) {
          
          sparse_mat_1[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_1[ is.na(sparse_mat_1) ] <- 0
      }
      
      if( formula == "center"){
        
        count_1    <- sparse_mat_1[ as.character(midpoint_2) , as.character(midpoint_1) ]
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      } else if( formula == "max"){
        
        count_1    <- max( sparse_mat_1, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      } else if( formula == "sum"){
        
        count_1    <- sum( sparse_mat_1, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      } else if( formula == "average"){
        
        count_1    <- mean( sparse_mat_1, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      }
      
      count_1 <- count_1 * norm_factor1 * aqua_factor1
      count_1 <- round( count_1 , 3 )
      
      dim_1 <- seq(                     
        as.numeric(unlist(strsplit(txt_L,":"))[2]),                    
        as.numeric(unlist(strsplit(txt_L,":"))[3]),                    
        bin_size )
      dim_2 <- seq(                     
        as.numeric(unlist(strsplit(txt_R,":"))[2]),                    
        as.numeric(unlist(strsplit(txt_R,":"))[3]),                    
        bin_size )
      
      sparse_mat_2 <- matrix( ncol = length(dim_1), nrow = length(dim_2) )
      
      colnames( sparse_mat_2 ) <- dim_1
      rownames( sparse_mat_2 ) <- dim_2
      
      for( i in 1:nrow(mat_2) ){
        
        row   <- mat_2[ i , "x"      ]
        col   <- mat_2[ i , "y"      ]
        count <- mat_2[ i , "counts" ]
        
        if( as.character(row) %in% rownames( sparse_mat_2 ) && 
            as.character(col) %in% colnames( sparse_mat_2 )) {
          
          sparse_mat_2[ as.character(row) , as.character(col) ] <- count
          
        } else if( as.character(row) %in% colnames( sparse_mat_2 ) && 
                   as.character(col) %in% rownames( sparse_mat_2 )) {
          
          sparse_mat_2[ as.character(col) , as.character(row) ] <- count
          
        }
        
        sparse_mat_2[ is.na(sparse_mat_2) ] <- 0
      }
      
      if( formula == "center"){
        
        count_2    <- sparse_mat_2[ midpoint_2 , midpoint_1]
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      } else if( formula == "max"){
        
        count_2    <- max( sparse_mat_2, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      } else if( formula == "sum"){
        
        count_2    <- sum( sparse_mat_2, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      } else if( formula == "average"){
        
        count_2    <- mean( sparse_mat_2, na.rm = T )
        
        new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
        new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
        
        new_end1 <- new_end1 + bin_size
        
        new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
        new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
        
        new_end2 <- new_end2 + bin_size
        
      }
      
      count_2 <- count_2 * norm_factor2 * aqua_factor2
      count_2 <- round( count_2 , 3 )
      
    } 
    else if ( nrow( mat_1 ) == 0 && nrow( mat_2 ) == 0 ){ 
      
      count_1 <- 0
      count_2 <- 0
      
      new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
      new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
      
      new_end1 <- new_end1 + bin_size
      
      new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
      new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
      
      new_end2 <- new_end2 + bin_size
      
    } 
    else { 
      
      count_1 <- "*"
      count_2 <- "*"
      
      new_start1 <- as.numeric(unlist(strsplit(txt_L,":"))[2])
      new_end1   <- as.numeric(unlist(strsplit(txt_L,":"))[3])
      
      new_end1 <- new_end1 + bin_size 
      
      new_start2 <- as.numeric(unlist(strsplit(txt_R,":"))[2])
      new_end2   <- as.numeric(unlist(strsplit(txt_R,":"))[3])
      
      new_end2 <- new_end2 + bin_size
      
    }
    
    obj <- list( 
      feet1  = c(new_start1, new_end1),
      feet2  = c(new_start2, new_end2),
      count1 = count_1,
      count2 = count_2 )
  }
  
  return(obj)
}

call_two_sample_straw_bedpe <- function( bedpe_list, flag_fix, formula, norm_factor1, aqua_factor1, norm_factor2, aqua_factor2 ){
  
  if( length( unlist( bedpe_list ) ) == 4 ){
    
    txt_L <- as.character( bedpe_list[1]$txt_L )
    txt_R <- as.character( bedpe_list[2]$txt_R )
    
    midpoint_1 <- as.character( bedpe_list[3]$midpoint_1_bin )
    midpoint_2 <- as.character( bedpe_list[4]$midpoint_2_bin )
    
    mat_1 <- straw( norm, hic_A, txt_L, txt_R, "BP", bin_size )
    mat_2 <- straw( norm, hic_B, txt_L, txt_R, "BP", bin_size )
    
    ret <- obtain_two_sample_counts( mat_1, mat_2, txt_L, txt_R, midpoint_1, midpoint_2 )
  }
  
  if( length( unlist( bedpe_list ) ) == 6 ){
    
    sample_A_txt_L <- as.character( bedpe_list[1]$sample_A_txt_L )
    sample_A_txt_R <- as.character( bedpe_list[2]$sample_A_txt_R )
    sample_B_txt_L <- as.character( bedpe_list[3]$sample_B_txt_L )
    sample_B_txt_R <- as.character( bedpe_list[4]$sample_B_txt_R )

    midpoint_1 <- as.character( bedpe_list[3]$midpoint_1_bin )
    midpoint_2 <- as.character( bedpe_list[4]$midpoint_2_bin )
    
    mat_1 <- straw( norm, hic_A, sample_A_txt_L, sample_A_txt_R, "BP", bin_size )
    mat_2 <- straw( norm, hic_B, sample_B_txt_L, sample_B_txt_R, "BP", bin_size )
    
    ret <- obtain_two_sample_counts( mat_1, mat_2, txt_L, txt_R, midpoint_1, midpoint_2 )
    
  }
  
  return(ret)
}



# auxiliary functions for --split :- 

get_neighbours <- function( cell_row, cell_col, radius ){
  
  neighbours <- matrix( data = NA, nrow = radius*8, ncol = 2 )
  
  vec_1 <- paste( seq( -radius, radius, 1 ), -radius, sep = "," )
  vec_2 <- paste( seq( -radius, radius, 1 ),  radius, sep = "," )
  
  vec_3 <- paste( -radius, seq( -radius, radius, 1 ), sep = "," )
  vec_4 <- paste(  radius, seq( -radius, radius, 1 ), sep = "," )
    
  vec <- unique( c(vec_1, vec_2, vec_3, vec_4) )
  
  for( i in 1:length(vec) ){
    
    a <- cell_row + as.numeric(unlist(strsplit(vec[i],","))[1])
    b <- cell_col + as.numeric(unlist(strsplit(vec[i],","))[2])
    
    neighbours[i,1] <- a
    neighbours[i,2] <- b
  }
  
  return(neighbours)
}

get_clusters <- function( sparse_mat_1, split, padding, chr, bin_size, formula, norm_factor1, aqua_factor1 ){
  
  radius <- padding
  
  split_indices     <- which( sparse_mat_1 >= split, arr.ind = T )
  split_indices     <- split_indices[order(split_indices[,2], split_indices[,1]),]
  split_indices_ids <- paste( split_indices[,1], split_indices[,2], sep = "-")
  
  boolean_matrix <- matrix( data = NA, nrow = nrow(sparse_mat_1), ncol = ncol(sparse_mat_1) )
  boolean_matrix[split_indices] <- 1
  
  rownames(boolean_matrix) <- 1:nrow(boolean_matrix)
  colnames(boolean_matrix) <- 1:ncol(boolean_matrix)
  
  # 1. Get the neighboring indices of all cells greater than split threshold 
  
  cell_clusters <- list()
  
  for(i in 1:nrow(split_indices)){
    
    
    cell_row <- split_indices[i,"row"]
    cell_col <- split_indices[i,"col"]
    
    
    cell_clusters[[i]] <- paste( cell_row, cell_col, sep = "-")
    
    for( j in 1:radius){
      
      neighbours <- get_neighbours( cell_row, cell_col, j )
      
      for( k in 1:nrow(neighbours)){
        
        neighbour_row <- neighbours[ k , 1 ]
        neighbour_col <- neighbours[ k , 2 ]
        
        if( paste( neighbour_row, neighbour_col, sep = "-") %in% split_indices_ids ){
          
          cell_clusters[[i]] <- c( cell_clusters[[i]], paste( neighbour_row, neighbour_col, sep = "-") )
        }    
      }
    }
  }
  
  # 2. Construct shared neighbor matrix 
  
  names(cell_clusters) <- split_indices_ids
  
  shared_neighbours <- matrix( data = NA, nrow = length(split_indices_ids), ncol = length(split_indices_ids) )
  
  colnames(shared_neighbours) <- split_indices_ids
  rownames(shared_neighbours) <- split_indices_ids  
  
  for( i in 1:length(cell_clusters) ){
    
    for( j in 1:length(cell_clusters[[i]])){
      
      row <- names(cell_clusters[i])
      col <- cell_clusters[[i]][j]
      
      shared_neighbours[row,col] <- 1
    }
  }
  

  # 3. Obtain cluster indices
  
  raw_clusters <- list()
  
  for( i in 1:ncol(shared_neighbours) ){
    
    vec <- !(is.na(shared_neighbours[,i]))
    vec <- which(vec == TRUE)
    
    raw_clusters[[i]] <- as.numeric(vec)
  }
  
  
  for( i in 1:length(raw_clusters) ){
    
    members <- c()
    
    vec_1 <- raw_clusters[[i]]
    
    for( j in (i+1):length(raw_clusters) ){
      
      if(j > length(split_indices_ids)){break}
      
      vec_2 <- raw_clusters[[j]]
      
      if( any( vec_2 %in% vec_1 )){
        
        vec_1 <- unique( c(vec_1, vec_2) )
        
        members <- c( members, j)
      }
    }
    
    raw_clusters[[i]] <- unique( c(raw_clusters[[i]] , as.numeric(unlist(raw_clusters[members] ))) )
    
    raw_clusters[members] <- NA
    
  }
  
  final_clusters <- Filter(Negate(anyNA),raw_clusters)
  
  for( i in 1:length(final_clusters) ){
    
    cluster <- split_indices_ids[ final_clusters[[i]] ]
    cluster <- as.numeric(unlist(strsplit(cluster,"-")))
    
    row_min <- min(cluster[seq(1,length(cluster),2)])
    row_max <- max(cluster[seq(1,length(cluster),2)])
    
    col_min <- min(cluster[seq(2,length(cluster),2)])
    col_max <- max(cluster[seq(2,length(cluster),2)])
    
    mat <- as.matrix(sparse_mat_1[ row_min:row_max , col_min:col_max ])
    
    if( formula == "center"){
      
      count_1    <- mat[ ceiling(nrow(mat)/2) , ceiling(ncol(mat)/2) ]
      count_1    <- count_1 * norm_factor1 * aqua_factor1
      count_1    <- round( count_1 , 3 )
      
    } else if( formula == "max"){
      
      count_1    <- max( mat, na.rm = T )
      count_1    <- count_1 * norm_factor1 * aqua_factor1
      count_1    <- round( count_1 , 3 )
      
    } else if( formula == "sum"){
      
      count_1    <- sum( mat, na.rm = T )
      count_1    <- count_1 * norm_factor1 * aqua_factor1
      count_1    <- round( count_1 , 3 )
      
    } else if( formula == "average"){
      
      count_1    <- mean( mat, na.rm = T )
      count_1    <- count_1 * norm_factor1 * aqua_factor1
      count_1    <- round( count_1 , 3 )
    }
    
    final_clusters[[i]]         <- list()
    final_clusters[[i]]$chr     <- chr
    final_clusters[[i]]$feet1   <- c( as.numeric(rownames(sparse_mat_1)[row_min]) - (radius*bin_size) , as.numeric(rownames(sparse_mat_1)[row_max]) + (radius*bin_size) + bin_size )
    final_clusters[[i]]$feet2   <- c( as.numeric(colnames(sparse_mat_1)[col_min]) - (radius*bin_size) , as.numeric(colnames(sparse_mat_1)[col_max]) + (radius*bin_size) + bin_size )
    #final_clusters[[i]]$feet1   <- c( as.numeric(rownames(sparse_mat_1)[row_min]) - 0 , as.numeric(rownames(sparse_mat_1)[row_max]) + bin_size )
    #final_clusters[[i]]$feet2   <- c( as.numeric(colnames(sparse_mat_1)[col_min]) - 0 , as.numeric(colnames(sparse_mat_1)[col_max]) + bin_size )
    final_clusters[[i]]$count_1 <- count_1
  }
  
  return(final_clusters)
  
}

zero_diag <- function( matrix, width ){

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

#################################################################
##                          Arguments                          ##
#################################################################

args        <- commandArgs( trailingOnly = TRUE )
path_pairs  <- args[1]
bin_size    <- as.numeric(args[2])
norm        <- "NONE"
flag_debug  <- FALSE

if( length(args) == 12 ){
  
  one_sample_analysis <-  TRUE
  two_sample_analysis <- FALSE
  
  hic_A              <-  args[3]
  path_mergeStats_A  <-  args[4]
  flag_norm          <-  args[5]
  num_loops          <-  as.numeric(args[6])
  formula            <-  args[7]
  flag_fix           <-  as.logical(args[8])
  split              <-  args[9]
  shrink_wrap        <-  args[10]
  padding            <-  as.numeric(args[11])
  expand             <-  as.numeric(args[12])
  
  hic_A_chroms      <-  readHicChroms(hic_A)$name
  
  if(formula == "mean"){ formula <- "average" }
}

if( length(args) == 14 ){
  
  one_sample_analysis <- FALSE
  two_sample_analysis <-  TRUE
  
  hic_A             <-  args[3]
  path_mergeStats_A <-  args[4]
  hic_B             <-  args[5]
  path_mergeStats_B <-  args[6]
  flag_norm         <-  args[7]
  num_loops         <-  as.numeric(args[8])
  formula           <-  args[9]
  flag_fix          <-  as.logical(args[10])
  split             <-  args[11]
  shrink_wrap       <-  args[12]
  padding           <-  as.numeric(args[13])
  expand            <-  as.numeric(args[14])
  
  hic_A_chroms      <- readHicChroms(hic_A)$name
  hic_B_chroms      <- readHicChroms(hic_B)$name
  
  if(formula == "mean"){ formula <- "average" }

  if(shrink_wrap != "FALSE"){ cat("--shrink_wrap can only be applied to one sample analysis \n") ; q(save="no") }
  if(split       != "FALSE"){ cat("--split can only be applied to one sample analysis \n")       ; q(save="no") }
}



if( num_loops <= 10000){ 
  chunk_size <- num_loops 
} else if( num_loops > 10000){
  chunk_size <- 10000
}

if( shrink_wrap != "FALSE" || split != "FALSE" ){
  flag_fix <- FALSE
}

if( expand < 0 ){
  cat("Invalid use of --expand \n") ; q(save="no") 
}

if( shrink_wrap != "FALSE" && split != "FALSE" ){ 
  cat("Please use --split or --shirnk_wrap one at a time \n") ; q(save="no") 
}


num_chunks <- ceiling( num_loops / chunk_size )

con <- file( description = path_pairs, open = "r" )
pairs <- read.table(con, nrows=chunk_size, as.is = T )

index <- 0
repeat {
  
  index <- index + 1
  
  if( ncol(pairs) < 6 || ncol(pairs) > 6 ){
    cat("Please use a 6 column bedpe only, exiting....\n")
    q( save = "no" )
  }
  
  colnames(pairs) <- c("chr1","start1","end1","chr2","start2","end2")
  
  if( sum( grepl( "chr", pairs$chr1 ) ) == 0 ) { pairs$chr1 <- paste0( "chr", pairs$chr1 ) }
  if( sum( grepl( "chr", pairs$chr2 ) ) == 0 ) { pairs$chr2 <- paste0( "chr", pairs$chr2 ) }
  
  pairs <- pairs[ pairs$chr1 %in% c( paste0( rep("chr",22), 1:22 ), "chrX", "chrY" ) , ]
  
  # convert coordinates to bins
  pairs$start1_bin <- as.integer( floor( pairs$start1 / bin_size ) * bin_size )
  pairs$end1_bin   <- as.integer( floor( pairs$end1   / bin_size ) * bin_size )
  pairs$start2_bin <- as.integer( floor( pairs$start2 / bin_size ) * bin_size )
  pairs$end2_bin   <- as.integer( floor( pairs$end2   / bin_size ) * bin_size )
  
  if( expand > 0 ){
    
    expand_size <- expand*bin_size
    
    for( i in 1:nrow(pairs) ){ 
      
      if( pairs[i,"start2_bin"] < pairs[i,"start1_bin"] ){
        a <- pairs[i,"start1_bin"] ; b <- pairs[i,"end1_bin"]
        c <- pairs[i,"start2_bin"] ; d <- pairs[i,"end2_bin"]
        
        pairs[i,"start1_bin"] <- c
        pairs[i,"end1_bin"  ] <- d
        pairs[i,"start2_bin"] <- a
        pairs[i,"end2_bin"  ] <- b
      }
      
      pairs[i,"start1_bin"] <- pairs[i,"start1_bin"] - expand_size
      pairs[i,"end1_bin"  ] <- pairs[i,"end1_bin"  ] + expand_size
      pairs[i,"start2_bin"] <- pairs[i,"start2_bin"] - expand_size
      pairs[i,"end2_bin"  ] <- pairs[i,"end2_bin"  ] + expand_size
      
      if( pairs[i,"end1_bin"] > pairs[i,"start2_bin"] ){
        
        pairs[i,"end1_bin"  ] <- pairs[i,"end1_bin"  ] - expand_size 
        pairs[i,"start2_bin"] <- pairs[i,"start2_bin"] + expand_size 
        
        new_expand <- (pairs[i,"start2_bin"] - pairs[i,"end1_bin"])/bin_size
        new_expand <- (new_expand/2) - 1
        
        pairs[i,"end1_bin"  ] <- pairs[i,"end1_bin"  ] + new_expand*bin_size
        pairs[i,"start2_bin"] <- pairs[i,"start2_bin"] - new_expand*bin_size
      }
    }
    
    pairs[i,"start1"] <- pairs[i,"start1_bin"]
    pairs[i,"end1"  ] <- pairs[i,"end1_bin"  ]
    pairs[i,"start2"] <- pairs[i,"start2_bin"]
    pairs[i,"end2"  ] <- pairs[i,"end2_bin"  ]
  }
  
  pairs$midpoint_1_bin <- as.integer( floor( (pairs$start1 + pairs$end1)/2   / bin_size ) * bin_size )
  pairs$midpoint_2_bin <- as.integer( floor( (pairs$start2 + pairs$end2)/2   / bin_size ) * bin_size )
  
  ##################################################################
  ##                          Code block                          ##
  ##################################################################
  
  if( one_sample_analysis ){
    
    mergeStats_A <- read.table( path_mergeStats_A, as.is = TRUE)
    hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
    mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
    total1       <- sum( hg_total1 , mm_total1 )
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1
    
    # add or remove chr prefix from bedpe  
    # if the .hic has 'chr' prefix
    if( sum( grepl("chr", hic_A_chroms ) ) > 0 ){
      flag_strip_chr <- "no"
    }
    # if the .hic has 'chr' prefix
    if( sum( grepl("chr", hic_A_chroms ) ) == 0 ){
      flag_strip_chr <- "yes"
    }
    
    if( ! flag_norm %in% c( "none", "cpm", "aqua" ) ){
      cat("Norm should strictly be none, cpm or aqua in lower case \n")
      q( save = "no" )
    }
    
    if( flag_norm == "none" ){
      norm_factor1 <- 1 ; aqua_factor1 <- 1
    } else if ( flag_norm == "cpm" ) {
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
    } else if ( flag_norm == "aqua" ) {
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- aqua_factor1
    }
    
    A  <-  apply( pairs, 1, prefix_filter_bedpe )
    B  <-  mclapply( A, call_one_sample_straw_bedpe, mc.cores = cores, flag_fix, formula, norm_factor1, aqua_factor1, shrink_wrap, split )
    
    if( isTRUE(flag_fix) ){
      
      C  <-  cbind( pairs[,1:6], do.call( "rbind", B ) )
      
    } else if( isFALSE(flag_fix) ){
      
      if( split == "FALSE" ){
        
        C  <-  as.data.frame(matrix(ncol=7)) ; colnames(C) <- c("chr1","start1","end1","chr2","start2","end2","count")
        
        for( i in 1:length(B)){
          
          C[ i , "chr1" ] <- pairs[ i , "chr1" ] ; C[ i , "start1" ] <- B[[i]]$feet1[1] ; C[ i , "end1" ] <- B[[i]]$feet1[2]
          C[ i , "chr2" ] <- pairs[ i , "chr2" ] ; C[ i , "start2" ] <- B[[i]]$feet2[1] ; C[ i , "end2" ] <- B[[i]]$feet2[2]
          C[ i , "count"] <- B[[i]]$count1
        }
        
      } else {
        
        C  <-  data.frame()
        
        for( i in 1:length(B)){
          
          length <- as.numeric(lengths(B[i]))
          
          data <- as.data.frame(matrix(ncol=7)) ; colnames(data) <- c("chr1","start1","end1","chr2","start2","end2","count")
          
          for( j in 1:length){

            data[ j , "chr1"   ] <- B[[as.character(i)]][[j]]$chr  
            data[ j , "start1" ] <- B[[as.character(i)]][[j]]$feet1[1] 
            data[ j , "end1"   ] <- B[[as.character(i)]][[j]]$feet1[2]
            
            data[ j , "chr2"   ] <- B[[as.character(i)]][[j]]$chr  
            data[ j , "start2" ] <- B[[as.character(i)]][[j]]$feet2[1] 
            data[ j , "end2"   ] <- B[[as.character(i)]][[j]]$feet2[2]
            
            data[ j , "count"  ] <- B[[as.character(i)]][[j]]$count_1
            
          }
          
          C <- rbind( C, data )
        }
        
      }
    }
    
    C[ , "chr1" ] <- as.character( C[ , "chr1" ] )
    C[ , "chr2" ] <- as.character( C[ , "chr2" ] )
    
    for( i in 1:nrow(C) ){ 
      
      if( C[ i , "start2" ] < C[ i , "start1" ] ){
        C[ i , ] <- C[ i , c(4:6,1:3,7) ]
      }
      
      try(cat( paste( C[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE) }
  }
  
  
  if( two_sample_analysis ){
    
    mergeStats_A <- read.table( path_mergeStats_A, as.is = TRUE)
    hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
    mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
    total1       <- sum( hg_total1 , mm_total1 )
    norm_factor1 <- 1000000 / total1
    aqua_factor1 <- hg_total1 / mm_total1
    
    mergeStats_B <- read.table( path_mergeStats_B, as.is = TRUE)
    hg_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 1 ])
    mm_total2    <- as.numeric(mergeStats_B[ "valid_interaction_rmdup" , 2 ])
    total2       <- sum( hg_total2 , mm_total2 )
    norm_factor2 <- 1000000 / total2
    aqua_factor2 <- hg_total2 / mm_total2
    
    ## add or remove chr prefix from bedpe  
    ## if both .hics have 'chr' prefix
    if( sum( grepl("chr", hic_A_chroms ) ) > 0 && sum( grepl("chr", hic_B_chroms ) ) > 0 ){
      flag_strip_chr <- "no"
    }
    ## if sample A .hic has 'chr' prefix but sample B doesn't
    if( sum( grepl("chr", hic_A_chroms ) ) > 0 && sum( grepl("chr", hic_B_chroms ) ) == 0 ){
      flag_strip_chr <- "A"
    }
    ## if sample B .hic has 'chr' prefix but sample A doesn't
    if( sum( grepl("chr", hic_B_chroms ) ) > 0 && sum( grepl("chr", hic_A_chroms ) ) == 0 ){
      flag_strip_chr <- "B"
    }
    ## if both .hics don't have 'chr' prefix
    if( sum( grepl("chr", hic_A_chroms ) ) == 0 && sum( grepl("chr", hic_B_chroms ) ) == 0 ){
      flag_strip_chr <- "yes"
    }
    
    if( ! flag_norm %in% c( "none", "cpm", "aqua" ) ){
      cat("Norm should strictly be none, cpm or aqua in lower case \n")
      q( save = "no" )
    }
    
    if( flag_norm == "none" ){
      norm_factor1 <- 1 ; aqua_factor1 <- 1
      norm_factor2 <- 1 ; aqua_factor2 <- 1
    } else if ( flag_norm == "cpm" ) {
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- 1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- 1
    } else if ( flag_norm == "aqua" ) {
      norm_factor1 <- norm_factor1 ; aqua_factor1 <- aqua_factor1
      norm_factor2 <- norm_factor2 ; aqua_factor2 <- aqua_factor2
    }
    
    A <- apply( pairs, 1, prefix_filter_bedpe )
    B <- mclapply( A, call_two_sample_straw_bedpe, mc.cores = cores, flag_fix, formula, norm_factor1, aqua_factor1, norm_factor2, aqua_factor2 )
    
    if( isTRUE(flag_fix) ){
      
      C  <-  cbind( pairs[,1:6], as.data.frame( matrix( unlist(B), nrow=length(B), byrow=TRUE ) ) )
      
    } else if( isFALSE(flag_fix) ){
      
      C  <-  as.data.frame(matrix(ncol=8)) ; colnames(C) <- c("chr1","start1","end1","chr2","start2","end2","count1","count2")
      
      for( i in 1:length(B)){
        C[ i , "chr1"  ] <- pairs[ i , "chr1" ] ; C[ i , "start1" ] <- B[[i]]$feet1[1] ; C[ i , "end1" ] <- B[[i]]$feet1[2]
        C[ i , "chr2"  ] <- pairs[ i , "chr2" ] ; C[ i , "start2" ] <- B[[i]]$feet2[1] ; C[ i , "end2" ] <- B[[i]]$feet2[2]
        C[ i , "count1"] <- B[[i]]$count1
        C[ i , "count2"] <- B[[i]]$count2
      }
    }
    
    C[ , "chr1" ] <- as.character( C[ , "chr1" ] )
    C[ , "chr2" ] <- as.character( C[ , "chr2" ] )
    
    for( i in 1:nrow(C) ){ 
      
      if( C[ i , "start2" ] < C[ i , "start1" ] ){
        C[ i , ] <- C[ i , c(4:6,1:3,7:8) ]
      }
      C$delta <- C[,8] - C[,7]
      
      try(cat( paste( C[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE) 
      
    }
    
  }
  
  if(index == num_chunks) break
  pairs <- read.table(con, nrows=chunk_size, as.is = T)
  
}

close(con)
