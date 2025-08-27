suppressPackageStartupMessages(library(dbscan))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(igraph))

options(stringsAsFactors = FALSE)
options(scipen = 999)
options(warn = -1 )

args           <- commandArgs( trailingOnly = TRUE )



#################################################################
##                          Arguments                          ##
#################################################################


file_bedpe <- args[1]
flank      <- as.numeric(args[2])
temp_path  <- args[3]


bedpe <- read.table( file_bedpe, as.is = TRUE )


flag_meta <- FALSE


if (ncol(bedpe) < 6) {
  
  cat("Please use a minimum of 6 column bedpe, exiting....\n")
  q(save = "no")
  
} else if (ncol(bedpe) > 6) {
  
  flag_meta <- TRUE
  
  meta_cols <- bedpe[, 7:ncol(bedpe)]
  
  meta_cols <- data.frame(meta_cols, stringsAsFactors = FALSE)
  
  num_meta_cols <- ncol(meta_cols)
  
  bedpe <- bedpe[, 1:6]
}


colnames(bedpe) <- c(
  "chr1",
  "start1",
  "end1",
  "chr2",
  "start2",
  "end2")


# add flank
bedpe$start1 <- bedpe$start1 - flank
bedpe$start2 <- bedpe$start2 - flank

bedpe$end1 <- bedpe$end1 + flank
bedpe$end2 <- bedpe$end2 + flank

# prepare for intersections
bedpe <- bedpe[order(bedpe$chr1,bedpe$start1),]
bedpe$id <- 1:nrow(bedpe)

feet1 <- bedpe[,c(1:3,7)] 
names(feet1) <- c("chr","start","end","id")

feet2 <- bedpe[,c(4:6,7)] 
names(feet2) <- c("chr","start","end","id")

bedpe_bed <- rbind(
  feet1,
  feet2 )

bedpe_bed <- bedpe_bed[order(
  bedpe_bed$chr,
  bedpe_bed$start ),]


write.table(
  bedpe_bed, 
  paste0(temp_path,"/","cluster_bedpe-tmp.bed"), 
  row.names = F, 
  col.names = F, 
  quote     = F, 
  sep       = "\t")


system(
  paste(
    "bedtools intersect -a ",
    paste0(temp_path,"/","cluster_bedpe-tmp.bed"),
    " -b ",
    paste0(temp_path,"/","cluster_bedpe-tmp.bed"),
    " -wa -wb | awk '{print $4,$8}' | awk '$1 != $2{print $0}' > ",
    paste0(temp_path,"/","pair_id.txt"),
    sep = ""
  )
)


pair_file <- paste0(temp_path,"/","pair_id.txt")

if(!file.exists(pair_file) || file.info(pair_file)$size == 0){
  
  bedpe$cluster <- paste0("cluster-", 1:nrow(bedpe))
  
} else {

  pair_id <- read.table(
    paste0(temp_path,"/","pair_id.txt"), 
    sep = " ")
  
  pair_id <- data.frame(
    x = pmin(pair_id$V1, pair_id$V2),
    y = pmax(pair_id$V1, pair_id$V2)
  )
  
  pair_id <- pair_id[order(pair_id$x, pair_id$y),]
  pair_id <- unique(pair_id)
  
  g               <- graph_from_edgelist(as.matrix(pair_id), directed = FALSE)
  components      <- components(g)
  pair_id$cluster <- components$membership[pair_id$x]
  pair_id$cluster <- as.integer(factor(pair_id$cluster))
  
  cluster_mapping <- rbind(
    pair_id %>% select(id = x, cluster),
    pair_id %>% select(id = y, cluster)
  ) %>% distinct()
  
  bedpe <- bedpe %>%
    full_join(cluster_mapping, by = "id")
  
  
  
  vec <- bedpe$cluster
  
  unique_counter <- 1
  index_map <- list()
  
  for (i in 1:length(vec)) {
    
    if (is.na(vec[i])) {
      
      vec[i]         <- unique_counter
      unique_counter <- unique_counter + 1
      
    } else {
      
      if (!is.null(index_map[[as.character(vec[i])]])) {
        vec[i] <- index_map[[as.character(vec[i])]]
        
      } else {
        
        index_map[[as.character(vec[i])]] <- unique_counter
        vec[i]         <- unique_counter
        unique_counter <- unique_counter + 1
        
      }
    }
  }
  
  bedpe$cluster <- paste0("cluster-",vec)
  
}

# remove flank
bedpe$start1 <- bedpe$start1 + flank
bedpe$start2 <- bedpe$start2 + flank

bedpe$end1 <- bedpe$end1 - flank
bedpe$end2 <- bedpe$end2 - flank

if(flag_meta){
  bedpe <- cbind(
    bedpe[,1:6,],
    meta_cols,
    bedpe$cluster
  )
} else {
  bedpe <- cbind(
    bedpe[,1:6,],
    bedpe$cluster
  )
}

for( i in 1:nrow(bedpe) ){ 
  try(cat( paste( bedpe[i,], collapse = "\t"), "\n", sep = "" ), silent=TRUE) 
}
