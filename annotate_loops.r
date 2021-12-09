library(strawr)

args        <- commandArgs( trailingOnly = TRUE )
path_pairs  <- args[1]
bin_size    <- as.numeric(args[2])
norm        <- "NONE"
flag_debug  <- FALSE

if( length(args) == 5 ){
  
  one_sample_analysis <-  TRUE
  two_sample_analysis <- FALSE
  
  hic_A             <- args[3]
  path_mergeStats_A <- args[4]
  flag_aqua         <- args[5]
}



if( length(args) == 7 ){
  
  one_sample_analysis <- FALSE
  two_sample_analysis <-  TRUE

  hic_A             <- args[3]
  path_mergeStats_A <- args[4]
  hic_B             <- args[5]
  path_mergeStats_B <- args[6]
  flag_aqua         <- args[7]
}

if( grepl("false", flag_aqua, ignore.case = TRUE) ){
  #cat("# no, we won't use aqua factors\n")
  flag_aqua <- FALSE
} else {
  #cat("# yes, we'll use aqua factors\n")
  flag_aqua <- TRUE
}



pairs <- read.table( path_pairs, as.is = TRUE )
pairs <- pairs[,1:6]

colnames(pairs) <- c("chr1","start1","end1","chr2","start2","end2")


# convert coordinates to bins
pairs$a1_bin <- (pairs$start1 + pairs$end1)/2
pairs$a2_bin <- (pairs$start2 + pairs$end2)/2
pairs$a1_bin <- as.integer( floor( pairs$a1_bin / bin_size ) * bin_size )
pairs$a2_bin <- as.integer( floor( pairs$a2_bin / bin_size ) * bin_size )


# colouring loops
if(two_sample_analysis){
  
  delta_max_value <- 2.5349 #magic number for colours
  break_list = seq(
    -delta_max_value, 
    delta_max_value, 
    by = delta_max_value/100)
  
  color_ramp_vio <- colorRampPalette(
    c("dodgerblue", "white", "mediumvioletred"))(length(break_list))
  
  C <- matrix( ncol = 2, nrow = length(break_list), data = 0)
  rownames( C ) <- color_ramp_vio
  colnames( C ) <- c("rank","breaks")
  C[,1] <- 1:nrow(C)
  C[,2] <- break_list
  
}

# colouring loops
if(one_sample_analysis){
  
  sample_A_max_value <- 25 #magic number for colours
  break_list = seq(
    0, 
    sample_A_max_value, 
    by = sample_A_max_value/100)
  
  color_ramp_vio <- colorRampPalette(
    c("dodgerblue", "white", "mediumvioletred"))(length(break_list))
  
  C <- matrix( ncol = 2, nrow = length(break_list), data = 0)
  rownames( C ) <- color_ramp_vio
  colnames( C ) <- c("rank","breaks")
  C[,1] <- 1:nrow(C)
  C[,2] <- break_list
  
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
  
 
  if( flag_aqua ){
    cat("# AQuA-CPM\n")
  } else {
    cat("# CPM\n")
    aqua_factor1 <- 1
    aqua_factor2 <- 1
  }
 
  for(i in 1:nrow(pairs)){
    
    # left and right chromosomes and bins
    chr_L <- sub("chr","", pairs[i,"chr1"])
    chr_R <- sub("chr","", pairs[i,"chr2"])
    bin_L <- pairs[i,"a1_bin"]
    bin_R <- pairs[i,"a2_bin"]
    txt_L <- sprintf( "%s:%d:%d", chr_L, bin_L, bin_L )
    txt_R <- sprintf( "%s:%d:%d", chr_R, bin_R, bin_R )


    # use straw to query the .hic file
    mat_1 <- straw( norm, hic_A, txt_L, txt_R, "BP", bin_size )
    mat_2 <- straw( norm, hic_B, txt_L, txt_R, "BP", bin_size )

    if ( nrow( mat_1 ) == 1 ){
      count_1 <- mat_1[1,"counts"]
    } else if (nrow( mat_1 ) == 0) {
      count_1 <- 0
    } else {
      cat("# straw anomaly\n")
    }

    if ( nrow( mat_2 ) == 1 ){
      count_2 <- mat_2[1,"counts"]
    } else if (nrow( mat_2 ) == 0) {
      count_2 <- 0
    } else {
      cat("# straw anomaly\n")
    }

    if( flag_debug ){
      cat(sprintf( 
        "# straw out count_1: %d count_2: %d\n", 
        count_1, count_2 ))
    }

    # normalize counts
    count_1 <- count_1 * norm_factor1 * aqua_factor1
    count_2 <- count_2 * norm_factor2 * aqua_factor2

    aqua_delta <- count_2 - count_1

    aqua_color <- rownames(
      C[ order( abs( C[,"breaks"] - aqua_delta ), decreasing = FALSE ), ])[1]
    
    try(cat(sprintf(
      "%s\t%d\t%d\t%s\t%d\t%d\t%6.4f\t%6.4f\t%6.4f\t%s\n",
      pairs[i,"chr1"], pairs[i,"start1"], pairs[i,"end1"],
      pairs[i,"chr2"], pairs[i,"start2"], pairs[i,"end2"],
      count_1, count_2, aqua_delta, aqua_color
    )), silent=TRUE)
  }


}





if(one_sample_analysis){

  mergeStats_A <- read.table( path_mergeStats_A, as.is = TRUE)
  hg_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 1 ])
  mm_total1    <- as.numeric(mergeStats_A[ "valid_interaction_rmdup" , 2 ])
  total1       <- sum( hg_total1 , mm_total1 )
  norm_factor1 <- 1000000 / total1
  aqua_factor1 <- hg_total1 / mm_total1


  if( flag_aqua ){
    cat("# AQuA-CPM\n")
  } else {
    cat("# CPM\n")
    aqua_factor1 <- 1
  }
 
  for(i in 1:nrow(pairs)){
    
    # left and right chromosomes and bins
    chr_L <- sub("chr","", pairs[i,"chr1"])
    chr_R <- sub("chr","", pairs[i,"chr2"])
    bin_L <- pairs[i,"a1_bin"]
    bin_R <- pairs[i,"a2_bin"]
    txt_L <- sprintf( "%s:%d:%d", chr_L, bin_L, bin_L )
    txt_R <- sprintf( "%s:%d:%d", chr_R, bin_R, bin_R )


    # use straw to query the .hic file
    mat_1 <- straw( norm, hic_A, txt_L, txt_R, "BP", bin_size )

    if ( nrow( mat_1 ) == 1 ){
      count_1 <- mat_1[1,"counts"]
    } else if (nrow( mat_1 ) == 0) {
      count_1 <- 0
    } else {
      cat("# straw anomaly\n")
    }

    if( flag_debug ){
      cat(sprintf( 
        "# straw out count_1: %d\n", 
        count_1 ))
    }

    # normalize counts
    count_1 <- count_1 * norm_factor1 * aqua_factor1

    aqua_color <- rownames(
      C[ order( abs( C[,"breaks"] - count_1 ), decreasing = FALSE ), ])[1]
    
    try(cat(sprintf(
      "%s\t%d\t%d\t%s\t%d\t%d\t%6.4f\t%s\n",
      pairs[i,"chr1"], pairs[i,"start1"], pairs[i,"end1"],
      pairs[i,"chr2"], pairs[i,"start2"], pairs[i,"end2"],
      count_1, aqua_color
    )), silent=TRUE)
  }

  
}
