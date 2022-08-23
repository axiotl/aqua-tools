suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

options(scipen = 999)

Args <- commandArgs(trailingOnly=T)



###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

## define arguments
if( length(Args) == 8 ) {


  matrix_file1  <- Args[1]
  sample1       <- Args[2]
  hard_cap      <- Args[3]
  prefix        <- Args[4]
  out_file      <- Args[5]
  win_size      <- Args[6]
  res           <- Args[7]
  pairs         <- Args[8]


  pairs <- as.data.frame( read.table( pairs, as.is = TRUE ) )
  pairs <- pairs[ , 1:6 ]


  cat( "\n" )
  cat( sprintf( "matrix_file1   <- %s\n", matrix_file1   ) )
  cat( sprintf( "sample1        <- %s\n", sample1        ) )
  cat( sprintf( "hard_cap       <- %s\n", hard_cap       ) )
  cat( sprintf( "prefix         <- %s\n", prefix         ) )
  cat( sprintf( "out_file       <- %s\n", out_file       ) )
  cat( sprintf( "win_size       <- %s\n", win_size       ) )
  cat( sprintf( "res            <- %s\n", res            ) )


  ## define input matix ( cpm / AQuA / AQuA --cpml )
  mat1  <- as.matrix(read.table(matrix_file1, head=F, sep="\t"))


  if( hard_cap == "no_cap" ) { max_value <- max( mat1 ) } else { max_value <- as.numeric( hard_cap ) }
  

  ## define colors
  breakList       <-  seq( 0, max_value, by = max_value/100 )
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
  
  rownames( mat1 ) <- cols
  colnames( mat1 ) <- cols

  
  ## actual plotting
  
  p1 <- pheatmap(  mat1,
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   border_color = NA, 
                   color = color_ramp_red, 
                   breaks = breakList, 
                   show_colnames = TRUE,
                   show_rownames = FALSE    
                )

  #grid.arrange(arrangeGrob(p1[[4]], nrow=1,top = paste(prefix, ": ", sample1,sep="")))
  g <- arrangeGrob( p1[[4]], 
                    nrow = 1,
                    top = paste( prefix, ": ", sample1, " (n.loop=", nrow(pairs), ")", sep = "" ) 
                  )
  #ggsave(file=out_file, g, width = 4.5, height = 4.5) #saves g

    ggsave( file = out_file, 
          g, 
          width = 4.5, 
          height = 4.5, 
          device = "pdf"
        ) #saves g
  
  ## PATCH: https://github.com/tidyverse/ggplot2/issues/2787
  file.exists("Rplots.pdf")
  file.remove("Rplots.pdf")


  ## plotting distribution of distances between bedpe feet

  ## NOTE: this for now is done using a for loop which for 
  ##       large .bedpe files will be slow. Need to create
  ##       an apply function and parallelize this in the 
  ##       near future

  distance_between_feet <- c()

  for( i in 1:nrow( pairs ) ){
    
    ## only consider those feet that are intra-chromosomal
    if( pairs[ i , 1 ] == pairs[ i , 4 ] ){

      mp1 <- as.numeric( pairs[ i , 2 ] ) + ( as.numeric( pairs[ i , 3 ] ) - as.numeric( pairs[ i , 2 ] ) ) / 2
      mp2 <- as.numeric( pairs[ i , 5 ] ) + ( as.numeric( pairs[ i , 6 ] ) - as.numeric( pairs[ i , 5 ] ) ) / 2

      distance_between_feet[ i ] <- abs( mp2 - mp1 )

    }

  }

  distance_between_feet <- distance_between_feet[ complete.cases( distance_between_feet ) ]

  distribution <- data.frame( distance = round(distance_between_feet) ) 

  if( nrow( distribution ) < 500 ){ bins <- 30  }
  if( nrow( distribution ) > 500 ){ bins <- 100 }

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

}



###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

## define arguments

if( length(Args) == 11 ) {



  matrix_file1   <-  Args[1]
  matrix_file2   <-  Args[2]
  sample1        <-  Args[3]
  sample2        <-  Args[4]
  hard_cap       <-  Args[5]
  hard_cap_delta <-  Args[6]
  prefix         <-  Args[7]
  out_file       <-  Args[8]
  win_size       <-  Args[9]
  res            <- Args[10]
  pairs          <- Args[11]

  pairs <- as.data.frame( read.table( pairs, as.is = TRUE ) )
  pairs <- pairs[ , 1:6 ]

  # Hey Berkley,
  # These variables are usually sent from the shell script.
  # I'm printing them out here so that you can load them into an R session.

  cat("\n")
  cat(sprintf("matrix_file1   <- %s\n", matrix_file1   ))
  cat(sprintf("matrix_file2   <- %s\n", matrix_file2   ))
  cat(sprintf("sample1        <- %s\n", sample1        ))
  cat(sprintf("sample2        <- %s\n", sample2        ))
  cat(sprintf("hard_cap       <- %s\n", hard_cap       ))
  cat(sprintf("hard_cap_delta <- %s\n", hard_cap_delta ))
  cat(sprintf("prefix         <- %s\n", prefix         ))
  cat(sprintf("out_file       <- %s\n", out_file       ))


  ## define input matix ( cpm / AQuA / AQuA --cpml )
  mat1  <- as.matrix(read.table(matrix_file1, head=F, sep="\t"))
  mat2  <- as.matrix(read.table(matrix_file2, head=F, sep="\t"))
  delta <- mat2 - mat1

  if( hard_cap       == "no_cap" ) { max_value <- max( mat1, mat2 ) } else { max_value <- as.numeric(hard_cap)  }
  if( hard_cap_delta == "no_cap" ) { delta_max <- max( abs(delta) ) } else { delta_max <- as.numeric(hard_cap_delta) }


  ## define colors
  breakList       <- seq(         0, max_value, by = max_value/100)
  breakList_delta <- seq(-delta_max, delta_max, by = delta_max/50 )

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
  
  rownames(mat1) <- cols ; rownames(mat2) <- cols ; rownames(delta) <- cols
  colnames(mat1) <- cols ; colnames(mat2) <- cols ; colnames(delta) <- cols


  ## actual plotting
  p1 <- pheatmap(  mat1, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   border_color = NA, 
                   color = color_ramp_red, 
                   breaks = breakList, 
                   show_colnames = TRUE, 
                   show_rownames = FALSE     
                )

  p2 <- pheatmap(  mat2, 
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

  #grid.arrange(arrangeGrob(p1[[4]],p2[[4]],p3[[4]],nrow=1,top = paste(prefix, ": ",sample1,", ",sample2,", and delta",sep="")))
  g <- arrangeGrob( p1[[4]], p2[[4]], p3[[4]],
                    nrow = 1,
                    top = paste(prefix, ": ", sample1, ", ", sample2, ", and delta", " (n.loop=", nrow(pairs), ")", sep = "" )
                  )

  ggsave( file = out_file, 
          g, 
          width = 14, 
          height = 4.5, 
          device = "pdf"
        ) #saves g
  
  ## PATCH: https://github.com/tidyverse/ggplot2/issues/2787
  file.exists("Rplots.pdf")
  file.remove("Rplots.pdf")

  ## plotting distribution of distances between bedpe feet

  ## NOTE: this for now is done using a for loop which for 
  ##       large .bedpe files will be slow. Need to create
  ##       an apply function and parallelize this in the 
  ##       near future

  distance_between_feet <- c()

  for( i in 1:nrow( pairs ) ){
    
    ## only consider those feet that are intra-chromosomal
    if( pairs[ i , 1 ] == pairs[ i , 4 ] ){

      mp1 <- as.numeric( pairs[ i , 2 ] ) + ( as.numeric( pairs[ i , 3 ] ) - as.numeric( pairs[ i , 2 ] ) ) / 2
      mp2 <- as.numeric( pairs[ i , 5 ] ) + ( as.numeric( pairs[ i , 6 ] ) - as.numeric( pairs[ i , 5 ] ) ) / 2

      distance_between_feet[ i ] <- abs( mp2 - mp1 )

    }

  }

  distance_between_feet <- distance_between_feet[ complete.cases( distance_between_feet ) ]

  distribution <- data.frame( distance = round(distance_between_feet) ) 

  if( nrow( distribution ) < 500 ){ bins <- 30  }
  if( nrow( distribution ) > 500 ){ bins <- 100 }

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

}

