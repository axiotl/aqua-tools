#!/bin/bash

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

sample_sheet="$HOME/setup/sample_sheet.txt"

ctrlc_cont=0

function no_sigint {

    let ctrlc_count++
    if [[ $ctrlc_count == 1 ]]; then
        exit 1
    fi

}

trap no_sigint EXIT

function usage {
    echo -e "usage: "
    echo -e "  query_bedpe.sh \\"
    echo -e "    -G GENOME_BUILD \\"
    echo -e "    -P PATH_TO_BEDPE_YOU_WANT_TO_ANNOTATE \\"
    echo -e "    -A NAME_OF_FIRST_SAMPLE \\"
    echo -e "   [-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Annotate a bedpe file with AQuA normalized contact values"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -P|--bedpe           : Full path to the bedpe file you want to annotate, without headers"
    echo "    -G|--genome          : The genome build the sample(s) has been processed using: hg19, hg38, or mm10"
    echo "    -A|--sample1         : Name of the sample as it appears on the Tinkerbox"
    echo "  [    --hic1          ] : If not using the tinkerbox, directly specify the full path to the .hic"
    echo "  [ -B|--sample2       ] : For two sample delta contacts, name of the second sample"
    echo "  [    --hic2          ] : If not using the tinkerbox, directly specify the full path to the second .hic"
    echo "  [ -Q|--norm          ] : Which normalization to use: none, cpm, aqua, or abc in lower case"
    echo "  [ -R|--resolution    ] : Resolution in bp for contact value calculation. Default 5000"
    echo "  [ -f|--formula       ] : Arithmetic for contact values: center, max, average, or sum. Default center"
    echo "  [ -F|--fix           ] : If FALSE, reports new coordinates using center or max arithmetic. Default TRUE"
    echo "  [    --expand        ] : Expands 1D bedpe feet in both directions by given bin value. Default 0"
    echo "  [ -i|--inherent      ] : If TRUE, hic values transformed to inherent units. Default FALSE"
    echo "  [ -m|--preserve_meta ] : If TRUE, bedpe metadata columns will be preserved. Default TRUE"
    echo "  [ -h|--help          ] : Help message"
    exit;
}

if [ $# -lt 1 ]
    then
    usage
    exit
fi

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--bedpe")         set -- "$@" "-P" ;;
      "--sample1")       set -- "$@" "-A" ;;
      "--hic1")          set -- "$@" "-H" ;;
      "--hic2")          set -- "$@" "-I" ;;
      "--genome")        set -- "$@" "-G" ;;
      "--norm")          set -- "$@" "-Q" ;;
      "--sample2")       set -- "$@" "-B" ;;
      "--resolution")    set -- "$@" "-R" ;;
      "--formula")       set -- "$@" "-f" ;;
      "--fix")           set -- "$@" "-F" ;;
      "--expand")        set -- "$@" "-e" ;; 
      "--inherent")      set -- "$@" "-i" ;;
      "--preserve_meta") set -- "$@" "-m" ;;
      "--help")          set -- "$@" "-h" ;;
       *)                set -- "$@" "$arg"
  esac
done


Q="blank"
R=5000
f=center
F=TRUE
S=FALSE
s=FALSE
p=2
e=0
i=FALSE
m=TRUE
H=blank
I=blank

while getopts ":P:A:H:I:G:Q:B:R:f:F:e:i:m:h" OPT
do
    case $OPT in
  P) P=$OPTARG;;
  A) A=$OPTARG;;
  H) H=$OPTARG;;
  I) I=$OPTARG;;
  G) G=$OPTARG;;
  Q) Q=$OPTARG;;
  B) B=$OPTARG;;
  R) R=$OPTARG;;
  f) f=$OPTARG;;
  F) F=$OPTARG;;
  e) e=$OPTARG;;
  i) i=$OPTARG;;
  m) m=$OPTARG;;
  h) help ;;
  \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
  :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      exit 1
      ;;
    esac
done


# Do all necessary parameter checks
#----------------------------------

# First validate sample1 inputs (A and H)
if [[ -z "$A" && "$H" == "blank" ]]; then
    echo "No sample1 or hic1 provided. Exiting."
    usage
    exit 1
fi

if [[ ! -z "$A" && "$H" != "blank" ]]; then
    echo "Need either sample1 name or hic1 path, not both. Exiting."
    usage
    exit 1
fi

# Initialize as single sample analysis
two_sample=FALSE

# Now check sample2 inputs (B and I)
if [[ ! -z "$B" || "$I" != "blank" ]]; then
    # At least one of B or I is provided - validate combination
    if [[ ! -z "$B" && "$I" != "blank" ]]; then
        echo "Need either sample2 name or hic2 path, not both. Exiting."
        usage
        exit 1
    fi
    # Valid two sample case
    two_sample=TRUE
fi



if [[ -z $P || ! -f $P ]]; then
    echo "\nInvalid file specified in '$P'\n"
    usage
    exit 1
fi

if [[ "$i" == "TRUE" && ( "$Q" != "blank" && "$Q" != "none" ) ]];
then
    echo -e "\n\n--inherent is not compatible with aqua, cpm, or abc normalization. Choose either --norm $Q or --inherent TRUE\n"
    exit 1
fi

if [ -z "$G" ]; then
    echo "\nNo genome build provided. Exiting.\n"
    usage
    exit 1
fi


#----------------------------------


cat > /tmp/calculate_stats.r << 'EOF'
suppressPackageStartupMessages(library(strawr))
args <- commandArgs(trailingOnly = T)
hic_file <- args[1] 
sample   <- args[2]
genome   <- args[3]
chromosomes <- readHicChroms(hic_file)
chromosomes <- chromosomes[chromosomes[,"length"] > 46000000,"name"]
if("All" %in% chromosomes){
 chromosomes <- chromosomes[-(which(chromosomes == "All"))]}
total_counts <- 0
for(i in 1:length(chromosomes)){
 for(j in i:length(chromosomes)){
   try({
     mat <- straw("NONE", hic_file, chromosomes[i], chromosomes[j], "BP", 2500000)
     count <- sum(mat[,"counts"])
     total_counts <- total_counts + count
   }, silent = TRUE)  
 }
}
interaction_types <- c("valid_interaction", "valid_interaction_rmdup", "trans_interaction", "cis_interaction", "cis_shortRange", "cis_longRange")
values            <- c(total_counts, total_counts, "-", "-", "-", "-")
mergeStats <- data.frame(interaction_type = interaction_types, values = as.character(values))
header <- paste(sample, genome, sep = ".")
cat("\t",header, "\n")
for (row in 1:nrow(mergeStats)) {
 cat(paste(mergeStats[row, "interaction_type"], "\t", mergeStats[row, "values"]), "\n")
}
EOF


get_sample_directory() {
local input_sample_name="$1"
local base_sample_name
local sample_dir
if [[ "$input_sample_name" =~ _version_[0-9]+$ ]]; then
    base_sample_name=$(echo "$input_sample_name" | grep -oP '.*(?=_version_[0-9]+$)' )
    sample_dir="$data_dir/$G/$base_sample_name/$input_sample_name"
else
    base_sample_name="$input_sample_name"
     local sample_info=$(awk -F$'\t' -v sample="$base_sample_name" -v genome="$G" '$1 == sample && $6 == "1" && $2 == genome {print; exit}' "$sample_sheet")
    if [ ! -z "$sample_info" ]; then
        local version_number=$(echo "$sample_info" | cut -d$'\t' -f5)
        sample_dir="$data_dir/$G/$base_sample_name/${base_sample_name}_version_${version_number}"
    else
        return 1
    fi
fi
if [ -d "$sample_dir" ]; then
    echo "$sample_dir"
else
    return 1
fi
}



if [ "$two_sample" == "FALSE" ] ; then

    if [ "$H" == "blank" ]; then

        # Update sample A
        sample_dir_A=$(get_sample_directory "$A")
        if [ $? -eq 0 ]; then
            version_dir_A=$(basename "$sample_dir_A")
            base_dir_A=$(basename "$(dirname "$sample_dir_A")")
            A="${base_dir_A}/${version_dir_A}"
        else
            echo "Sample directory not found. Exiting."
            exit 1
        fi

        path_hic="$data_dir"/"$G"/"$A"/"${version_dir_A}".allValidPairs.hic
        path_mgs="$data_dir"/"$G"/"$A"/mergeStats.txt 

    else 

        path_hic=$H
        Rscript /tmp/calculate_stats.r "$path_hic" "sample" "$G" > /tmp/mergeStats.txt
        path_mgs=/tmp/mergeStats.txt

    fi

fi


if [ "$two_sample" == "TRUE" ] ; then
   # Handle first sample (A or H)
   if [ "$H" == "blank" ]; then
       # Using sample A
       sample_dir_A=$(get_sample_directory "$A")
       if [ $? -eq 0 ]; then
           version_dir_A=$(basename "$sample_dir_A")
           base_dir_A=$(basename "$(dirname "$sample_dir_A")")
           A="${base_dir_A}/${version_dir_A}"
           path_hic_A="$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic"
           path_mgs_A="$data_dir/$G/$A/mergeStats.txt"
       else
           echo "Sample A directory not found. Exiting."
           exit 1
       fi
   else
       # Using hic1
       path_hic_A=$H
       Rscript /tmp/calculate_stats.r "$path_hic_A" "sample_A" "$G" > /tmp/mergeStats_A.txt
       path_mgs_A=/tmp/mergeStats_A.txt
   fi

   # Handle second sample (B or I)
   if [ "$I" == "blank" ]; then
       # Using sample B
       sample_dir_B=$(get_sample_directory "$B")
       if [ $? -eq 0 ]; then
           version_dir_B=$(basename "$sample_dir_B")
           base_dir_B=$(basename "$(dirname "$sample_dir_B")")
           B="${base_dir_B}/${version_dir_B}"
           path_hic_B="$data_dir/$G/$B/${version_dir_B}.allValidPairs.hic"
           path_mgs_B="$data_dir/$G/$B/mergeStats.txt"
       else
           echo "Sample B directory not found. Exiting."
           exit 1
       fi
   else
       # Using hic2
       path_hic_B=$I
       Rscript /tmp/calculate_stats.r "$path_hic_B" "sample_B" "$G" > /tmp/mergeStats_B.txt
       path_mgs_B=/tmp/mergeStats_B.txt
   fi
fi



col1=$(head -n 1  $P | cut -f 1 )
col1="${col1:0:3}"

col4=$(head -n 1  $P | cut -f 4 )
col4="${col4:0:3}"

if [[ $col1 != $col4  ]]; then
  echo "Please provide a 6-col bedpe file without headers"
  exit
fi

#----------------------------------

num_loops=`cat "$P" | wc -l`


###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if [ "$two_sample" == "FALSE" ]
then
    Rscript \
    $aqua_dir/query_bedpe.r \
      $P \
      $R \
      $path_hic \
      $path_mgs \
      $Q \
      $G \
      $num_loops \
      $f $F $S $s $p $e $i $data_dir/$G/$A $m
fi

###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if [ "$two_sample" == "TRUE" ]
then
    Rscript \
    $aqua_dir/query_bedpe.r \
      $P \
      $R \
      $path_hic_A \
      $path_mgs_A \
      $path_hic_B \
      $path_mgs_B \
      $Q \
      $G \
      $num_loops \
      $f $F $S $s $p $e $data_dir/$G/$A $data_dir/$G/$B $m $i
fi


if [ -f "/tmp/calculate_stats.r" ]; then
    rm /tmp/calculate_stats.r
fi
if [ -f "/tmp/mergeStats.txt" ]; then
    rm /tmp/mergeStats.txt
fi
if [ -f "/tmp/mergeStats_A.txt" ]; then
    rm /tmp/mergeStats_A.txt
fi
if [ -f "/tmp/mergeStats_B.txt" ]; then
    rm /tmp/mergeStats_B.txt
fi
