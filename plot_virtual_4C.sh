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
    echo "  [ -B|--sample2        ]   : For two sample analyses, name of the second sample"
    echo "  [    --sample1_dir    ]   : If not using the tinkerbox, specify the full path to the directory containing sample data"
    echo "  [    --sample2_dir    ]   : If not using the tinkerbox, full path to the second sample directory"
    echo "  [ -Q|--norm           ]   : Which normalization to use: none, cpm, or aqua in lower case"
    echo "  [ -r|--resolution     ]   : Resolution of sample in basepairs. Default 5000"
    echo "  [ -O|--output_name    ]   : Optional name for the plot"
    echo "  [    --quant_cut      ]   : Cap matrix values at a given percentile (0.00-1.00). Default 1.00"
    echo "  [    --max_cap        ]   : Cap matrix values by adjusting the maximum value, modifying color intensity"
    echo "  [    --width          ]   : Number of bins up/downstream of viewpoint included in profile. Default 0"
    echo "  [    --height         ]   : Value to control the height of the profile. Default calculated automatically"
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
      "--sample1_dir")  set -- "$@" "-H" ;;
      "--sample2_dir")  set -- "$@" "-I" ;;
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
t="blank"
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

if [[ -n "$A" && -z "$G" ]]; then
  usage
  exit 1
fi

#----------------------------------

if [[ -z $V ]];
then
  usage
  exit
fi

#----------------------------------


# Validate sample1 inputs
if [[ -z "$A" && "$H" == "blank" ]]; then
    echo "No sample provided. Exiting."
    usage
    exit 1
fi

if [[ ! -z "$A" && "$H" != "blank" ]]; then
    echo "Need either sample1 name or sample1 directory path, not both. Exiting."
    usage
    exit 1
fi

# Check if a directory was given to A
if [[ "$A" == *"/"* ]]; then
  echo "The -A parameter expects a sample name, not a path. Did you mean to use --sample1_dir?"
  exit 1
fi

# If A is provided, then G must be one of the allowed values: hg38, hg19, or mm10
if [[ -n "$A" ]]; then
  if [[ "$G" != "hg38" && "$G" != "hg19" && "$G" != "mm10" ]]; then
    echo "-G must be one of: hg38, hg19, or mm10."
    exit 1
  fi
fi

# Check for inherent normalization attempt with local samples
if [[ -n "$H" && "$i" == "TRUE" ]]; then
  echo -e "\n--inherent normalization is not compatible with local installation. Try tinker-public/Docker for access to inherent normalized public samples, or sign up for tinker-private to use inherent normalization with any sample. Continuing without --inherent TRUE..."
  i="FALSE"
fi


# Initialize as single sample analysis
two_sample=FALSE

# Check sample2/hic2 inputs
if [[ ! -z "$B" || "$I" != "blank" ]]; then
    if [[ ! -z "$B" && "$I" != "blank" ]]; then
        echo "Need either sample2 name or sample1 directory path, not both. Exiting."
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


validate_samples_and_stats() {
  local sample_dir_input="$1"  # either -H or -I value
  local sample_label="$2"      # sample1 or sample2 label
  local hic_path mgs_path sample_name
  
  # Validate that the directory exists
  if [[ ! -d "$sample_dir_input" ]]; then
      echo "Directory '$sample_dir_input' does not exist. Exiting."
      exit 1
  fi

  # Extract sample name from the directory name
  sample_name=$(basename "$sample_dir_input")
  
  # Define expected file paths
  hic_path="${sample_dir_input}/${sample_name}.hic"
  mgs_path="${sample_dir_input}/${sample_name}.mergeStats.txt"
  
  # Check if the .hic file exists
  if [[ ! -f "$hic_path" ]]; then
      echo "Missing .hic file: $hic_path"
      exit 1
  fi
  
  # Check if the mergeStats.txt file exists
  if [[ ! -f "$mgs_path" ]]; then
      echo "Missing mergeStats.txt file: $mgs_path"
      exit 1
  fi

  # Check for the required row and count data columns
  local valid_row num_fields data_columns
  valid_row=$(grep "^valid_interaction_rmdup" "$mgs_path")
  if [ -z "$valid_row" ]; then
      echo "The mergeStats file ($mgs_path) does not contain the required row 'valid_interaction_rmdup'. Please refer to the GitHub repository for formatting instructions: https://github.com/axiotl/aqua-tools"
      exit 1
  fi
  num_fields=$(echo "$valid_row" | awk -F'\t' '{print NF}')
  data_columns=$(( num_fields - 1 ))
  
  # Check for at least one data column (human reads)
  if [[ "$data_columns" -lt 1 ]]; then
      echo "$mgs_path is not properly formatted. Please refer to the GitHub repository for formatting instructions: https://github.com/axiotl/aqua-tools"
      exit 1
  fi

  # If aqua normalization is requested, check for a second data column (mouse reads)
  if [[ "$Q" == "aqua" ]]; then
      if [[ "$data_columns" -lt 2 ]]; then
          echo "Aqua normalization requested, but the mergeStats file ($mgs_path) does not contain a column for mouse reads. Please refer to the GitHub repository for formatting instructions: https://github.com/axiotl/aqua-tools"
          exit 1
      fi
  fi

  # Return the values as _A or _B
  if [ "$sample_label" == "sample1" ]; then
      path_hic_A="$hic_path"
      path_mgs_A="$mgs_path"
      sample_dir_A="$sample_dir_input"
  elif [ "$sample_label" == "sample2" ]; then
      path_hic_B="$hic_path"
      path_mgs_B="$mgs_path"
      sample_dir_B="$sample_dir_input"
  fi
}


if [ "$two_sample" == "FALSE" ]; then
    if [ "$H" == "blank" ]; then
        # Using sample A via sample name
        sample_dir_A=$(get_sample_directory "$A")
        if [ $? -eq 0 ]; then
            version_dir_A=$(basename "$sample_dir_A")
            base_dir_A=$(basename "$(dirname "$sample_dir_A")")
            A="${base_dir_A}/${version_dir_A}"
            sample_dir="$data_dir/$G/$A"
            path_hic_A="$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic"
            path_mgs_A="$data_dir/$G/$A/mergeStats.txt"
        else
            echo "Sample directory not found. Exiting."
            exit 1
        fi
    else
        # Use G if provided or set to placeholder value
        G=${G:-"G"}
        # Process -H as sample1 directory
        validate_samples_and_stats "$H" "sample1"
    fi
fi

if [ "$two_sample" == "TRUE" ]; then
    # Handle first sample (A or H)
    if [ "$H" == "blank" ]; then
        sample_dir_A=$(get_sample_directory "$A")
        if [ $? -eq 0 ]; then
            version_dir_A=$(basename "$sample_dir_A")
            base_dir_A=$(basename "$(dirname "$sample_dir_A")")
            A="${base_dir_A}/${version_dir_A}"
            path_hic_A="$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic"
            path_mgs_A="$data_dir/$G/$A/mergeStats.txt"
            sample_dir_A="$data_dir/$G/$A"
        else
            echo "Sample A directory not found. Exiting."
            exit 1
        fi
    else
        G=${G:-"G"}
        validate_samples_and_stats "$H" "sample1"
    fi

    # Handle second sample (B or I)
    if [ "$I" == "blank" ]; then
        sample_dir_B=$(get_sample_directory "$B")
        if [ $? -eq 0 ]; then
            version_dir_B=$(basename "$sample_dir_B")
            base_dir_B=$(basename "$(dirname "$sample_dir_B")")
            B="${base_dir_B}/${version_dir_B}"
            path_hic_B="$data_dir/$G/$B/${version_dir_B}.allValidPairs.hic"
            path_mgs_B="$data_dir/$G/$B/mergeStats.txt"
            sample_dir_B="$data_dir/$G/$B"
        else
            echo "Sample B directory not found. Exiting."
            exit 1
        fi
    else
        validate_samples_and_stats "$I" "sample2"
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

  pair1=$path_hic_A
  if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
  
  stat1=$path_mgs_A
  if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi

  
  hg_total1=`head -3 $stat1 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
  mm_total1=`head -3 $stat1 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`

  echo -e "\ntotal sample reads for $(basename "$path_hic"): $hg_total1"
  echo -e "total spike-in reads for $(basename "$path_hic"): $mm_total1\n"

  has_aqua=true
  total1=$hg_total1
  if [[ ! -z $mm_total1  ]];then
    echo "spike-in values detected"
    total1=`echo "$hg_total1+$mm_total1" | bc`
    aqua_factor1=`echo "scale=4;$hg_total1/$mm_total1" | bc`
  fi
  norm_factor1=`echo "scale=4;1000000/$total1" | bc`
  
  echo -e "\nswitching to R...\n"

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
    $sample_dir_A \
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
    echo "spike-in values detected"
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
