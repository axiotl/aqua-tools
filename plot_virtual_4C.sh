#!/bin/bash


OPTIND=1
norm_method=NONE
unit=BP

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools


sample_sheet="/home/ubuntu/setup/sample_sheet.txt"


juicer_tools='java -jar ~/juicer_tools_1.19.02.jar'


function usage {
    echo -e "usage : plot_virtual_4C.sh -A NAME_OF_FIRST_SAMPLE -G GENOME_BUILD -R RANGE -V VIEWPOINT [-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Plot Virtual 4C"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -A|--sample1            : Name of the sample as it appears on the Tinkerbox"
    echo "    -G|--genome             : Genome build used for sample processing"
    echo "    -R|--range              : The genomic range to plot in chr:start:end format; -R chr1:40280000:40530000"
    echo "    -V|--viewpoint          : The viewpoint to consider in chr:start format; -V chr1:40400000"
    echo "  [    --hic1           ]   : If not using the tinkerbox, directly specify the full path to the .hic"
    echo "  [ -B|--sample2        ]   : For two sample delta plots, name of the second sample"
    echo "  [    --hic2           ]   : If not using the tinkerbox, directly specify the second .hic"
    echo "  [ -Q|--norm           ]   : Which normalization to use: none, cpm, or aqua in lower case"
    echo "  [ -r|--resolution     ]   : Resolution of sample in basepairs. Default 5000"
    echo "  [ -O|--output_name    ]   : Optional name for the plot"
    echo "  [    --quant_cut      ]   : Cap matrix values at a given percentile (0.00-1.00). Default 1.00"
    echo "  [    --max_cap        ]   : Cap matrix values by adjusting the maximum value, modifying color intensity"
    echo "  [    --width          ]   : Number of bins up/downstream of viewpoint included in profile. Default 0"
    echo "  [    --height         ]   : Numeric factor to control the height of the profile. Default 0.5"
    echo "  [ -i|--inherent       ]   : If TRUE, normalize using inherent normalization. Default FALSE"
    echo "  [    --inh_col_floor  ]   : Contact color for inherent values < 0, in RGB hexadecimal. Default FFFFFF"
    echo "  [    --inh_col_off    ]   : Contact color for inherent values ~ 0, in RGB hexadecimal. Default D4E4FB"
    echo "  [    --inh_col_on     ]   : Contact color for inherent values ~ 1, in RGB hexadecimal. Default FF0000"
    echo "  [    --inh_col_ceil   ]   : Contact color for inherent values > 1, in RGB hexadecimal. Default FF8D4A"
    echo "  [ -h|--help           ]     Help message"
    echo
    echo

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
      "--sample1")      set -- "$@" "-A" ;;
      "--genome")       set -- "$@" "-G" ;;
      "--range")        set -- "$@" "-R" ;;
      "--viewpoint")    set -- "$@" "-V" ;;
      "--hic1")         set -- "$@" "-H" ;;
      "--hic2")         set -- "$@" "-I" ;;
      "--sample2")      set -- "$@" "-B" ;;
      "--norm")         set -- "$@" "-Q" ;;
      "--resolution")   set -- "$@" "-r" ;;
      "--output_name")  set -- "$@" "-O" ;;
      "--quant_cut")    set -- "$@" "-q" ;;
      "--max_cap")      set -- "$@" "-m" ;;
      "--width")        set -- "$@" "-w" ;;
      "--height")       set -- "$@" "-t" ;;
      "--inherent")     set -- "$@" "-i" ;;
      "--inh_col_floor") set -- "$@" "-K" ;;
      "--inh_col_off")   set -- "$@" "-L" ;;
      "--inh_col_on")    set -- "$@" "-M" ;;
      "--inh_col_ceil")  set -- "$@" "-N" ;;
      "--help")          set -- "$@" "-h" ;;
       *)                set -- "$@" "$arg"
  esac
done


Q="blank"
r=5000
q=1
m="none"
w=0
t=0.5
i="FALSE"
K=ffffff
L=d4e4fb
M=ff0000
N=ff8d4a
H=blank  
I=blank  


while getopts ":A:G:R:V:H:I:B:Q:r:O:q:m:w:t:i:K:L:M:N:h" OPT
do
  case $OPT in
  A) A=$OPTARG;;
  G) G=$OPTARG;;
  R) R=$OPTARG;;
  V) V=$OPTARG;;
  H) H=$OPTARG;;
  I) I=$OPTARG;;
  B) B=$OPTARG;;
  Q) Q=$OPTARG;;
  r) r=$OPTARG;;
  O) O=$OPTARG;;
  q) q=$OPTARG;;
  m) m=$OPTARG;;
  w) w=$OPTARG;;
  t) t=$OPTARG;;
  i) i=$OPTARG;;
  K) K=$OPTARG;;
  L) L=$OPTARG;;
  M) M=$OPTARG;;
  N) N=$OPTARG;;
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

if [[ -z $R ]];
then
  usage
  exit
fi

#-----------------------------------

if [[ -z $G ]];
then
  usage
  exit
fi

#----------------------------------

if [[ -z $V ]];
then
  usage
  exit
fi

#----------------------------------



# Validate sample1/hic1 inputs
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

# Check sample2/hic2 inputs
if [[ ! -z "$B" || "$I" != "blank" ]]; then
    if [[ ! -z "$B" && "$I" != "blank" ]]; then
        echo "Need either sample2 name or hic2 path, not both. Exiting."
        usage
        exit 1
    fi
    two_sample=TRUE
fi
#-----------------------------------

if [[ "$Q" != "blank" && "$i" != "FALSE" ]]; then
  echo -e "\nInherent normalization is not compatible with other --norm methods. Continuing with --inherent TRUE ..."
  Q="blank"
  i="TRUE"
fi


# Generate output name if not provided
if [[ -z $O ]]; then
  IFS=":" read -ra RANGE_ARRAY <<< "$R"
  interval_chr="${RANGE_ARRAY[0]}"
  interval_start="${RANGE_ARRAY[1]}"
  interval_end="${RANGE_ARRAY[2]}"
  O="${interval_chr}_${interval_start}_${interval_end}_4C.pdf"
else
  # Add .pdf if it's missing
  if [[ $O != *.pdf ]]; then
    O="${O}.pdf"
  fi
fi



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



if [ "$two_sample" == "FALSE" ] ; then

    if [ "$H" == "blank" ]; then

        # Update sample A
        sample_dir_A=$(get_sample_directory "$A")
        if [ $? -eq 0 ]; then
            version_dir_A=$(basename "$sample_dir_A")
            base_dir_A=$(basename "$(dirname "$sample_dir_A")")
            A="${base_dir_A}/${version_dir_A}"
            sample_dir="$data_dir/$G/$A"
        else
            echo "Sample directory not found. Exiting."
            exit 1
        fi

        path_hic="$data_dir"/"$G"/"$A"/"${version_dir_A}".allValidPairs.hic
        path_mgs="$data_dir"/"$G"/"$A"/mergeStats.txt 

    else 

        path_hic=$H
        sample_dir=$(basename "$H")
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
           sample_dir="$data_dir/$G/$A"
       else
           echo "Sample A directory not found. Exiting."
           exit 1
       fi
   else
       # Using hic1
       path_hic_A=$H
       sample_dir_A=$(basename "$H")
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
           sample_dir_A="$data_dir/$G/$A"
           sample_dir_B="$data_dir/$G/$B"
       else
           echo "Sample B directory not found. Exiting."
           exit 1
       fi
   else
       # Using hic2
       path_hic_B=$I
       sample_dir_B=$(basename "$I")
       Rscript /tmp/calculate_stats.r "$path_hic_B" "sample_B" "$G" > /tmp/mergeStats_B.txt
       path_mgs_B=/tmp/mergeStats_B.txt
   fi
fi




range_text=`echo "$R" | sed 's/:/_/g'`

chr=`echo "$R" | cut -d":" -f1`
start=`echo "$R" | cut -d":" -f2`
end=`echo "$R" | cut -d":" -f3`

range_length=$((end - start))

if [[ $range_length -gt 3000000 && $r -le 5000 ]]; then
  echo -e "\nThe range length ($range_length bp) is too large for the specified resolution ($r). Try reducing the range or increasing the resolution."
  exit
fi


range=`printf "%s:%s:%s" $chr $start $end`

v_chr=`  echo "$V" | cut -d":" -f1`
v_start=`echo "$V" | cut -d":" -f2`




###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if [ "$two_sample" == "FALSE" ]; then

  analysis_type="single_sample"
  echo -e "\nCommencing $analysis_type analysis"

  pair1=$path_hic
  if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
  
  stat1=$path_mgs
  if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi

  
  hg_total1=`head -3 $stat1 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
  mm_total1=`head -3 $stat1 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`

  echo -e "\ntotal sample reads for $(basename "$path_hic"): $hg_total1"
  echo -e "total spike-in reads for $(basename "$path_hic"): $mm_total1\n"

  has_aqua=true
  total1=$hg_total1
  if [[ ! -z $mm_total1  ]];then
    echo "We have spike-In"
    total1=`echo "$hg_total1+$mm_total1" | bc`
    aqua_factor1=`echo "scale=4;$hg_total1/$mm_total1" | bc`
  fi
  norm_factor1=`echo "scale=4;1000000/$total1" | bc`
  
  echo -e "switching to R...\n"

  Rscript $aqua_dir/plot_virtual_4C.r \
    $analysis_type \
    $pair1 \
    $stat1 \
    $norm_method \
    $unit \
    $r \
    $O \
    $G \
    $Q \
    $chr $start $end \
    $q \
    $m \
    $w \
    $t \
    $v_chr \
    $v_start \
    $sample_dir \
    $i \
    $K \
    $L \
    $M \
    $N

  exit_status=$?

  # Report successful plot generation location
  if [ $exit_status -eq 0 ] && [ -f "$O" ]; then
      echo -e "\nDone! Created $O\n"
  fi


fi


###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if [ "$two_sample" == "TRUE" ]; then

  analysis_type="two_sample"
  echo -e "\nCommencing $analysis_type analysis\n"

  if [[ $i == "TRUE" ]];
  then
    echo -e "--inherent only applies to one-sample plots for now, exiting...\n"
    exit
  fi

  
  pair1=$path_hic_A
  pair2=$path_hic_B
  if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
  if [[ ! -f $pair2 ]]; then echo "cannot find $pair2"; exit 1; fi

  stat1=$path_mgs_A
  stat2=$path_mgs_B
  if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi
  if [[ ! -f $stat2 ]]; then echo "cannot find $stat2"; exit 1; fi
  
  
  
  hg_total1=`head -3 $stat1 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
  hg_total2=`head -3 $stat2 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
  mm_total1=`head -3 $stat1 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`
  mm_total2=`head -3 $stat2 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`
  
  echo -e "\ntotal sample reads for $(basename "$path_hic_A"): $hg_total1"
  echo -e "total spike-in reads for $(basename "$path_hic_A"): $mm_total1"
  echo -e "total sample reads for $(basename "$path_hic_B"): $hg_total2"
  echo -e "total spike-in reads for $(basename "$path_hic_B"): $mm_total2\n"

  has_aqua=true
  total1=$hg_total1
  total2=$hg_total2
  if [[ ! -z $mm_total1 && ! -z $mm_total2 ]];then
    echo "We have spike-In"
    total1=`echo "$hg_total1+$mm_total1" | bc`
    total2=`echo "$hg_total2+$mm_total2" | bc`
    aqua_factor1=`echo "scale=4;$hg_total1/$mm_total1" | bc`
    aqua_factor2=`echo "scale=4;$hg_total2/$mm_total2" | bc`

  fi
  norm_factor1=`echo "scale=4;1000000/$total1" | bc`
  norm_factor2=`echo "scale=4;1000000/$total2" | bc`

  echo -e "\nswitching to R...\n"

  Rscript $aqua_dir/plot_virtual_4C.r \
    $analysis_type \
    $pair1 \
    $stat1 \
    $pair2 \
    $stat2 \
    $norm_method \
    $unit \
    $r \
    $O \
    $G \
    $Q \
    $chr $start $end \
    $tag \
    $q \
    $m \
    $w \
    $t \
    $v_chr \
    $v_start \
    $v_tag \
    $sample_dir_A $sample_dir_B \
    $i \
    $K \
    $L \
    $M \
    $N
  

  exit_status=$?

  # Report successful plot generation location
  if [ $exit_status -eq 0 ] && [ -f "$O" ]; then
      echo -e "\nDone! Created $O\n"
  fi
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
