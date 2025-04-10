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
    echo "  [ -B|--sample2       ] : For two sample delta contacts, name of the second sample"
    echo "  [    --sample1_dir   ] : If not using the tinkerbox, specify the full path to the directory containing sample data"
    echo "  [    --sample2_dir   ] : If not using the tinkerbox, full path to the second sample directory"
    echo "  [ -Q|--norm          ] : Which normalization to use: none, cpm, aqua, or abc in lower case"
    echo "  [ -R|--resolution    ] : Resolution in bp for contact value calculation. Default 5000"
    echo "  [ -f|--formula       ] : Arithmetic for contact values: center, max, average, or sum. Default center"
    echo "  [ -F|--fix           ] : If FALSE, reports new coordinates using center or max arithmetic. Default TRUE"
    echo "  [    --expand        ] : Expands bedpe feet in both directions by supplied value in base pairs. Default 0"
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
      "--sample1_dir")   set -- "$@" "-H" ;;
      "--sample2_dir")   set -- "$@" "-I" ;;
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

if [[ -n "$A" && -z "$G" ]]; then
  usage
  exit 1
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

# Check if a path was given to A
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
# Checks for $P

if [[ -z $P ]];
then
    usage
    exit
fi

if [[ ! -f "$P" ]]; then
    echo "\n'$P' is not a valid file\n"
    exit
fi

first_line=$(head -n 1 "$P")

col1=$(echo "$first_line" | cut -f1)
col4=$(echo "$first_line" | cut -f4)
if [[ "${col1:0:3}" != "${col4:0:3}" ]]; then
    echo "Please provide a 6-col BEDPE file without headers."
    exit 1
fi

awk_output=$(awk '
        BEGIN {FS="\t"; errors=0; reversed=0}
        NF < 6 {
            print "BEDPE error: Line " NR " has fewer than 6 columns"
            errors++
        }
        $3 < $2 || $6 < $5 {
            print "BEDPE error: Reversed feet in line " NR
            reversed++
        }
        END {
            if (errors > 0) {
                print "BEDPE error: file should contain at least 6 columns in all rows."
                exit 1
            }
            if (reversed > 0) {
                print "Please ensure that column 3 > column 2 and column 6 > column 5 for all entries."
                exit 1
            }
            exit 0
        }' "$P")


if [[ $? -ne 0 ]]; then
    echo "$awk_output"
    exit 1
fi

#----------------------------------

if [[ "$i" == "TRUE" && ( "$Q" != "blank" && "$Q" != "none" ) ]];
then
    echo -e "\n\n--inherent is not compatible with aqua, cpm, or abc normalization. Choose either --norm $Q or --inherent TRUE\n"
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
            sample_dir_A="$data_dir/$G/$A"
            path_hic_A="$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic"
            path_mgs_A="$data_dir/$G/$A/mergeStats.txt"
        else
            echo "Sample directory not found. Exiting."
            exit 1
        fi
    else
        if [ "$i" == "TRUE" ]; then
           echo -e "\nInherent normalization is not compatible with this sample."
           exit 1
        fi 
        # Use G if provided or set to placeholder value
        G=${G:-"G"}
        A="blank"
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
        if [ "$i" == "TRUE" ]; then
            echo -e "\nInherent normalization is not compatible with this sample."
            exit 1
        fi
        G=${G:-"G"}
        A="blank"
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
        B="blank"
    fi
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
      $path_hic_A \
      $path_mgs_A \
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
