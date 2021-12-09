
args <- commandArgs( trailingOnly = TRUE )
file_a <- args[1]
file_b <- args[2]
file_t <- args[3]
drop_d <- args[4]


a <- read.table( file_a, as.is = TRUE )[,1:3]
b <- read.table( file_b, as.is = TRUE )[,1:3]
colnames(a) <- c("chr","start","end")
colnames(b) <- c("chr","start","end")


if (file_t != "NULL" ){
  t <- read.table( file_t, as.is = TRUE )[,1:3]  
  colnames(t) <- c("chr","start","end")
}

d <- as.numeric( drop_d )





add_tad <- function( x ){

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


if ( file_t == "NULL" ){
  a <- cbind( a, rep( 0, nrow(a) )) 
  b <- cbind( b, rep( 0, nrow(b) )) 
  colnames(a)[ncol(a)] <- "tad"
  colnames(b)[ncol(b)] <- "tad"
} else {
  # t <- cbind( t, sprintf("%TAD04d", 1:nrow(t)) )
  t <- cbind( t,  1:nrow(t) )
  colnames(t)[ncol(t)] <- "tad"  
  a <- add_tad(a)
  b <- add_tad(b)
}

A <- merge( a, b, suffixes = c("_a","_b"), by = "tad" )
A <- A[,-1]


for ( i in 1:nrow( A ) ){

  if ( A[i,"chr_a"] == A[i,"chr_b"] ){
    center_a <- ( A[i,"end_a"] + A[i,"start_a"] ) / 2
    center_b <- ( A[i,"end_b"] + A[i,"start_b"] ) / 2
    distance <- abs( center_a - center_b )
    if ( distance >= d ){
      try(cat( paste( A[i,], collapse = "\t"), "\n", sep = "" ),silent=TRUE)
    }
  } else {
      try(cat( paste( A[i,], collapse = "\t"), "\n", sep = "" ),silent=TRUE)
  }

}

