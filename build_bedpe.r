
#################################################################
##                          Functions                          ##
#################################################################

## this function checks if the two .beds are contained
## within a given TAD (if any)
add_tad               <- function( x ){

  x <- cbind( x, rep( 0, nrow(x) ))
  colnames(x)[ncol(x)] <- "tad"

  for( i in 1:nrow(t)){
    sub <- which(
      x[,"end"  ]  > t[i,"start"] &
      x[,"start"]  < t[i,"end"  ] &
      x[,"chr"  ] == t[i,"chr"  ] )
    if( length(sub) > 0 ){
      x[sub,"tad"] <- t[i,"tad"]
    }
  }
  x <- x[x[,"tad"] != 0,]
  return(x)
}



#################################################################
##                          Arguments                          ##
#################################################################

args     <- commandArgs( trailingOnly = TRUE )
file_a   <- args[1]
file_b   <- args[2]
file_t   <- args[3]
min_dist <- args[4]
max_dist <- args[5]


a <- read.table( file_a, as.is = TRUE )[,1:3]
b <- read.table( file_b, as.is = TRUE )[,1:3]
colnames(a) <- c("chr","start","end")
colnames(b) <- c("chr","start","end")



##################################################################
##                          Code-block                          ##
##################################################################

if (file_t != "NULL" ){
  t <- read.table( file_t, as.is = TRUE )[,1:3]  
  colnames(t) <- c("chr","start","end")
}

min_dist <- as.numeric( min_dist )
max_dist <- as.numeric( max_dist )

## add TAD information to the two .beds
if ( file_t == "NULL" ){
  a <- cbind( a, rep( 0, nrow(a) )) 
  b <- cbind( b, rep( 0, nrow(b) )) 
  colnames(a)[ncol(a)] <- "tad"
  colnames(b)[ncol(b)] <- "tad"
} else {
  t <- cbind( t,  1:nrow(t) )
  colnames(t)[ncol(t)] <- "tad"  
  a <- add_tad(a)
  b <- add_tad(b)
}


if( nrow(b) > nrow(a)  ){
  bed_A <- b ; bed_B <- a ; flag_flipped_cols <- TRUE  ; rm(a,b)
} else {
  bed_A <- a ; bed_B <- b ; flag_flipped_cols <- FALSE ; rm(a,b)
}


for( i in 1:nrow(bed_A) ){
  
  vec <- which( abs( bed_B[ , "start" ]  -  bed_A[ i , "start"] ) >= min_dist & 
                abs( bed_B[ , "start" ]  -  bed_A[ i , "start"] ) <= max_dist &
                bed_B[ , "chr"   ] ==  bed_A[ i , "chr"  ] &
                bed_B[ , "tad"   ] ==  bed_A[ i , "tad"  ]
                )
  
  
  if( length(vec) > 1 ){
    
    df <- cbind( do.call("rbind", replicate( length( vec ), bed_A[ i , ], simplify = FALSE)),
                 bed_B[ vec , ] ) 
    
  } else if( length(vec) == 1 ){
    
    df <- cbind( bed_A[ i , ] , bed_B[ vec , ] )
    
  } else {
    
    next
    
  }
  
  if( flag_flipped_cols ){
    
    df[ , 1 ] <- as.character( df[ , 1 ] )
    df[ , 5 ] <- as.character( df[ , 5 ] )
    
    df <- df[ , c(5:7,1:3) ]
    
    for( i in 1:nrow(df) ){ try(cat( paste( df[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE) }
    
  } else {
    
    df[ , 1 ] <- as.character( df[ , 1 ] )
    df[ , 5 ] <- as.character( df[ , 5 ] )
    
    df <- df[ , c(1:3,5:7) ]
    
    for( i in 1:nrow(df) ){ try(cat( paste( df[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE) }
    
  }


}
