#!/bin/bash

OPTIND=1

win_size=10

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

sample_sheet="$HOME/setup/sample_sheet.txt"


juicer_tools='java -jar /home/ubuntu/juicer_tools_1.19.02.jar'


function usage {
    echo -e "usage : plot_APA.sh -P PATH_TO_GENOMIC_PAIRS_FILE_YOU_WANT_TO_PLOT -A NAME_OF_FIRST_SAMPLE -G GENOME_BUILD -O OUTPUT_DIRECTORY [-B NAME_OF_SECOND_SAMPLE] [-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Create APA plot with CPM/AQuA normalized contact values"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -P|--bedpe                  : Path to the bedpe file, without headers"
    echo "    -A|--sample1                : Name of the sample as it appears on the Tinkerbox"
    echo "    -G|--genome                 : Genome build used for sample processing"
    echo " [     --hic1                ]  : If not using the tinkerbox, directly specify the full path to the .hic"
    echo " [  -O|--out-dir             ]  : Name of the directory to store the output APA plot in"
    echo " [  -B|--sample2             ]  : For two sample delta plots, name of the second sample"
    echo " [     --hic2                ]  : If not using the tinkerbox, provide the path to the second .hic"
    echo " [  -Q|--norm                ]  : Which normalization to use: none, cpm, or aqua in lower case"
    echo " [  -r|--resolution          ]  : Bin size you want to use for the APA plots. Default 5000"
    echo " [     --max_cap             ]  : Upper limit of the plot range. Defaults to max bin value"
    echo " [     --min_cap             ]  : Lower limit of the plot range. Default 0"
    echo " [     --max_cap_delta       ]  : Upper limit of delta plot range. Defaults to max bin value"
    echo " [     --loop_norm           ]  : If TRUE, normalizes APA values by loop count in bedpe. Default FALSE"
    echo " [  -h|--help                ]   Help message"
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
      "--bedpe")                set -- "$@" "-P" ;;
      "--sample1")              set -- "$@" "-A" ;;
      "--genome")               set -- "$@" "-G" ;;
      "--out-dir")              set -- "$@" "-O" ;;
      "--sample2")              set -- "$@" "-B" ;;
      "--norm")                 set -- "$@" "-Q" ;;
      "--resolution")           set -- "$@" "-r" ;;
      "--max_cap")              set -- "$@" "-c" ;;
      "--min_cap")              set -- "$@" "-e" ;;
      "--max_cap_delta")        set -- "$@" "-d" ;;
      "--loop_norm")            set -- "$@" "-l" ;;
      "--hic1")                 set -- "$@" "-H" ;;
      "--hic2")                 set -- "$@" "-I" ;;
      "--help")                 set -- "$@" "-h" ;;
       *)                       set -- "$@" "$arg"
  esac
done

r=5000
c="no_cap"
d="no_cap"
e="no_cap"
l=FALSE
H=blank  
I=blank

while getopts ":P:A:G:O:B:Q:r:c:e:d:f:l:H:I:h" OPT
do
    case $OPT in
  P) P=$OPTARG;;
  A) A=$OPTARG;;
  G) G=$OPTARG;;
  B) B=$OPTARG;;
  O) O=$OPTARG;;
  r) r=$OPTARG;;
  c) c=$OPTARG;;
  Q) Q=$OPTARG;;
  l) l=$OPTARG;;
  d) d=$OPTARG;;
  e) e=$OPTARG;;
  H) H=$OPTARG;;
  I) I=$OPTARG;;
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

#----------------------------------

if [[ -z $P ]];
then
    usage
    exit
fi

if [[ ! -f "$P" ]]; then
    echo "\n'$P' is not a valid file\n"
    exit
fi

#----------------------------------

if [[ -z $G ]];
then
  usage
  exit
fi

#----------------------------------

# Validate sample1/hic1 inputs
if [[ -z "$A" && "$H" == "blank" ]]; then
    echo -e "\nNo sample1 or hic1 provided. Exiting.\n"
    usage
    exit 1
fi

if [[ ! -z "$A" && "$H" != "blank" ]]; then
    echo -e "\nNeed either sample1 name or hic1 path, not both. Exiting.\n"
    usage
    exit 1
fi

if [[ ! -z "$A" && "$A" == *.hic ]]; then
    echo -e "\n-A is for sample folders on Tinker. Did you mean to use the --hic1 option instead?\n"
    exit 1
fi

if [[ "$H" != "blank" && "$H" != *.hic ]]; then
    echo -e "\nFile '$H' is not a .hic file type. Please provide a valid .hic file.\n"
    exit 1
fi

# Initialize as single sample analysis
two_sample=FALSE

# Check sample2/hic2 inputs
if [[ ! -z "$B" || "$I" != "blank" ]]; then
    if [[ ! -z "$B" && "$I" != "blank" ]]; then
        echo "Need either sample2 name or hic2 path, not both. Exiting."
        usage
        exit 1
    fi
    two_sample=TRUE
fi

#----------------------------------

col1=$(head -n 1  "$P" | cut -f 1 )
col1="${col1:0:3}"

col4=$(head -n 1  "$P" | cut -f 4 )
col4="${col4:0:3}"


if [[ $col1 != "$col4"  ]]; then
echo "Please provide a 6-col bedpe file without headers"
exit
fi

first_line=$(head -n 1 "$P")
chr1=$(echo "$first_line" | awk '{print $1}')
chr2=$(echo "$first_line" | awk '{print $4}')

if [ "$chr1" != "$chr2" ]; then
  echo "APA statistics are not designed for interchromosomal regions, for more information see https://github.com/aidenlab/juicer/wiki/APA"
  exit 1
fi

#----------------------------------



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

#----------------------------------

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
        hic_basename=$(basename "$path_hic")

    else 

        path_hic=$H
        hic_basename=$(basename "$H")
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
           hic_A_basename=$(basename "$path_hic_A")
       else
           echo "Sample A directory not found. Exiting."
           exit 1
       fi
   else
       # Using hic1
       path_hic_A=$H
       hic_A_basename=$(basename "$H")
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
           hic_B_basename=$(basename "$path_hic_B")
       else
           echo "Sample B directory not found. Exiting."
           exit 1
       fi
   else
       # Using hic2
       path_hic_B=$I
       hic_B_basename=$(basename "$I")
       Rscript /tmp/calculate_stats.r "$path_hic_B" "sample_B" "$G" > /tmp/mergeStats_B.txt
       path_mgs_B=/tmp/mergeStats_B.txt
   fi
fi

#----------------------------------

# Determine the generated output name based on whether -O is provided
if [[ -z $O ]]; then
    random_num=$RANDOM
    
    if [ "$two_sample" == "TRUE" ] ; then
        # For a two-sample test
        O="${hic_A_basename}_v_${hic_B_basename}_APA_${random_num}"
    else
        # For a one-sample test
        O="${hic_basename}_APA_${random_num}"
    fi
    echo -e "\nNo output directory name provided. Generating output directory...\n"
fi

#----------------------------------

if [ -n "$Q" ]; then
    if [ "$Q" != "aqua" ] && [ "$Q" != "cpm" ] && [ "$Q" != "none" ]; then
        echo -e "\nInvalid value for -Q|--norm. Please use 'aqua', 'cpm', or 'none' in lower case."
        exit 1
    fi
else
    Q="blank"
fi

#---------------Binify-------------------

# Create temporary files in /tmp
tmp_trimmed=$(mktemp /tmp/trimmed_bedpe.XXXXXX)
tmp_binified=$(mktemp /tmp/binified_bedpe.XXXXXX)

# Create 6 col bedpe
cut -f1-6 "$P" > "$tmp_trimmed"

awk -v r="$r" 'BEGIN { OFS="\t" }
{
  if ((($3 - $2) < r) || (($6 - $5) < r)) {
    # If either interval is too short to bin, output the original coordinates
    print $1, $2, $3, $4, $5, $6;
  } else {
    end1_adj = $3 - r;
    end2_adj = $6 - r;
    for (i = $2; i <= end1_adj; i += r) {
      for (j = $5; j <= end2_adj; j += r) {
        print $1, i, i + r, $4, j, j + r;
      }
    }
  }
}' "$tmp_trimmed" > "$tmp_binified"


# Convert tmp_binified to absolute path
tmp_binified=$(readlink -f "$tmp_binified")


###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if [ "$two_sample" == "FALSE" ] ; then

#----------------------------------

    if [[ $d != "no_cap" ]]; then
    echo "Cannot have a delta cap in single-sample APA analysis"
    exit
    fi

#----------------------------------

    num_loops=`cat "$P" | wc -l`

    out_dir="$O"
    mkdir -p "$out_dir"

#----------------------------------

    pair1=$path_hic
    if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
  
    stat1=$path_mgs
    if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi


    # Call Juicer APA https://github.com/aidenlab/juicer/wiki/APA
    echo -e "\nBeginning APA...\n"
    $juicer_tools apa --threads 1 -k NONE -n 0 -r "$r" -w $win_size "$pair1" "$P" "$out_dir" &> /dev/null

    out_mat=$out_dir/$r/gw/APA.txt
    
    # Check expected file exists before attempting to move
    if [ -f "$out_mat" ]; then
        mv "$out_mat" "$out_dir/APA_raw.txt"
    else
        echo -e "\nAPA failed (no output file produced). Please check input files and parameters.\n"
        exit 1
    fi


    awk -F'[][,[]' '{ for(i=2; i<NF-1; i++) printf "%s\t", $i; printf "%s\n", $(NF-1) }' "$out_dir"/APA_raw.txt > "$out_dir"/APA_formatted.txt

    matrix_sum=$(awk '{ for(i=1;i<=NF;i++) sum+=$i } END { print sum }' "$out_dir/APA_formatted.txt")
    if [ "$matrix_sum" -eq 0 ]; then
        echo "APA matrix is all zeros - no valid contacts were detected."
        exit 1
    fi
    
    Rscript $aqua_dir/plot_APA.r "$out_dir/APA_formatted.txt" "${hic_basename}" "${c}" "${d}" "$out_dir/${hic_basename}_APA.pdf" ${win_size} "${r}" "${tmp_binified}" "${Q}" "$stat1" "${l}" "$num_loops" "$out_dir" "${e}"

    # Delete the specified directories and files after R script is called
    rm -r "$out_dir/$r"
    rm "$out_dir"/APA_raw.txt
    rm "$out_dir"/APA_formatted.txt
    rm "$tmp_trimmed"
    rm "$tmp_binified"

fi




###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if [ "$two_sample" == "TRUE" ] ; then

    num_loops=`cat "$P" | wc -l`

    out_dir=$O
    mkdir -p "$out_dir"

    #----------------------------------

    pair1=$path_hic_A
    pair2=$path_hic_B
    if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
    if [[ ! -f $pair2 ]]; then echo "cannot find $pair2"; exit 1; fi

    stat1=$path_mgs_A
    stat2=$path_mgs_B
    if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi
    if [[ ! -f $stat2 ]]; then echo "cannot find $stat2"; exit 1; fi

    #----------------------------------

    # https://github.com/aidenlab/juicer/wiki/APA

    echo -e "\nBeginning APA for sample A...\n"

    $juicer_tools apa --threads 1 -k NONE -n 0 -r "$r" -w $win_size "$pair1" "$tmp_binified" "$out_dir" &> /dev/null
    out_mat_A="$out_dir/$r/gw/APA.txt"
    if [ -f "$out_mat_A" ]; then
        mv "$out_mat_A" "$out_dir/APA_raw_${hic_A_basename}.txt"
    else
        echo -e "\nAPA failed for sample A (no output file produced). Please check input files and parameters.\n"
        exit 1
    fi

    echo -e "\nBeginning APA for sample B...\n"

    $juicer_tools apa --threads 1 -k NONE -n 0 -r "$r" -w $win_size "$pair2" "$tmp_binified" "$out_dir" &> /dev/null
    out_mat_B="$out_dir/$r/gw/APA.txt"
    if [ -f "$out_mat_B" ]; then
        mv "$out_mat_B" "$out_dir/APA_raw_${hic_B_basename}.txt"
    else
        echo -e "\nAPA failed for sample B (no output file produced). Please check input files and parameters\n"
        exit 1
    fi

    #----------------------------------


    awk -F'[][,[]' '{ for(i=2; i<NF-1; i++) printf "%s\t", $i; printf "%s\n", $(NF-1) }' "$out_dir"/APA_raw_$hic_A_basename.txt > "$out_dir"/APA_formatted_$hic_A_basename.txt

    awk -F'[][,[]' '{ for(i=2; i<NF-1; i++) printf "%s\t", $i; printf "%s\n", $(NF-1) }' "$out_dir"/APA_raw_$hic_B_basename.txt > "$out_dir"/APA_formatted_$hic_B_basename.txt

    # Sum values in each formatted matrix
    matrix_sum_A=$(awk '{ for(i=1;i<=NF;i++) sum+=$i } END { print sum }' "$out_dir/APA_formatted_${hic_A_basename}.txt")
    matrix_sum_B=$(awk '{ for(i=1;i<=NF;i++) sum+=$i } END { print sum }' "$out_dir/APA_formatted_${hic_B_basename}.txt")

    # Check if both matrices are all zeros
    if [ "$matrix_sum_A" -eq 0 ] && [ "$matrix_sum_B" -eq 0 ]; then
        echo "Both APA matrices are all zeros - no valid contacts were detected."
        exit 1
    fi

    Rscript $aqua_dir/plot_APA.r "$out_dir/APA_formatted_$hic_A_basename.txt" "$out_dir/APA_formatted_$hic_B_basename.txt" "${hic_A_basename}" "${hic_B_basename}" "${c}" "${d}" $"$out_dir/${hic_A_basename}_v_${hic_B_basename}_APA.pdf" ${win_size} "${r}" "${tmp_binified}" "${Q}" "$stat1" "$stat2" "${l}" "$num_loops" "$out_dir" "${e}"

    # Delete the specified directories and files after R script is called
    rm -r "$out_dir/$r"
    rm "$out_dir"/APA_raw_$hic_A_basename.txt
    rm "$out_dir"/APA_raw_$hic_B_basename.txt
    rm "$out_dir"/APA_formatted_$hic_A_basename.txt
    rm "$out_dir"/APA_formatted_$hic_B_basename.txt
    rm "$tmp_binified"
    rm "$tmp_trimmed"

fi
