#!/bin/bash

OPTIND=1

win_size=10

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

sample_sheet="$HOME/setup/sample_sheet.txt"

juicer_tools="java -jar $HOME/juicer_tools_1.19.02.jar"


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
    echo " [  -O|--out-dir             ]  : Name of the directory to store the output APA plot in"
    echo " [  -B|--sample2             ]  : For two sample delta plots, name of the second sample"
    echo " [     --sample1_dir         ]  : If not using the tinkerbox, specify the full path to the directory containing sample data"
    echo " [     --sample2_dir         ]  : If not using the tinkerbox, full path to the second sample directory"
    echo " [  -Q|--norm                ]  : Which normalization to use: none, cpm, or aqua in lower case"
    echo " [  -r|--resolution          ]  : Bin size you want to use for the APA plots. Default 5000"
    echo " [     --max_cap             ]  : Upper limit of the plot range. Defaults to max bin value"
    echo " [     --min_cap             ]  : Lower limit of the plot range. Default 0"
    echo " [     --max_cap_delta       ]  : Upper limit of delta plot range. Defaults to max bin value"
    echo " [     --loop_norm           ]  : If TRUE, normalizes APA values by loop count in bedpe. Default FALSE"
    echo " [  -s|--scores              ]  : If TRUE, plots peak bin bin value and ratios of peak:corner means. Default TRUE"
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
      "--sample1_dir")          set -- "$@" "-H" ;;
      "--sample2_dir")          set -- "$@" "-I" ;;
      "--scores")               set -- "$@" "-s" ;;
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
s=TRUE

while getopts ":P:A:G:O:B:Q:r:c:e:d:f:l:H:I:s:h" OPT
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
  s) s=$OPTARG;;
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
    echo -e "\n'$P' is not a valid file\n"
    exit
fi

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

first_line=$(head -n 1 "$P")

# Check for tab-separated bedpe with at least 6 columns
tab_check=$(echo "$first_line" | awk -F'\t' '{print NF}')
if [[ "$tab_check" -lt 6 ]]; then
    echo ".bedpe file must have at least 6 TAB-delimited columns"
    exit 1
fi

# Check for trailing whitespaces
if [[ "$first_line" =~ [[:space:]]+$ ]]; then
  sed -i '1s/[[:space:]]\+$//' "$P"
fi

# Check for a header with column 2 as numeric
col2=$(echo "$first_line" | cut -f2)
if ! [[ "$col2" =~ ^[0-9]+$ ]]; then
    echo "Please provide a 6-col BEDPE file without headers."
    exit 1
fi

# Check for intrachromosomal interactions.
chr1=$(echo "$first_line" | awk '{print $1}')
chr2=$(echo "$first_line" | awk '{print $4}')
if [[ "$chr1" != "$chr2" ]]; then
    echo "APA statistics are not designed for interchromosomal regions. See https://github.com/aidenlab/juicer/wiki/APA"
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
            basename=$version_dir_A
        else
            echo "Sample directory not found. Exiting."
            exit 1
        fi
    else
        # Use G if provided or set to placeholder value
        G=${G:-"G"}
        # Process -H as sample1 directory
        validate_samples_and_stats "$H" "sample1"
        basename=$(basename "$H")
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
            basename_A=$version_dir_A
        else
            echo "Sample A directory not found. Exiting."
            exit 1
        fi
    else
        G=${G:-"G"}
        validate_samples_and_stats "$H" "sample1"
        basename_A=$(basename "$H")
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
            basename_B=$version_dir_B
        else
            echo "Sample B directory not found. Exiting."
            exit 1
        fi
    else
        validate_samples_and_stats "$I" "sample2"
        basename_B=$(basename "$I")
    fi
fi

#----------------------------------

# Determine the generated output name based on whether -O is provided
if [[ -z $O ]]; then
    random_num=$RANDOM
    
    if [ "$two_sample" == "TRUE" ] ; then
        # For a two-sample test
        O="${basename_A}_v_${basename_B}_APA_${random_num}"
    else
        # For a one-sample test
        O="${basename}_APA_${random_num}"
    fi
    echo -e "\nNo output directory name provided. Generating output directory..."
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

#----------------------------------
case "$r" in
    1000|5000|10000|25000|50000|100000|250000|500000|1000000|2500000)
    # Value is accepted
    ;;
    *)
    echo -e "\n--resolution '$r' is not accepted.\nAccepted resolutions are: 1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000\n"
    exit 1
    ;;
esac

#----------------------------------

if [[ "$s" != "TRUE" && "$s" != "FALSE" ]]; then
  echo "Variable \$s must be either TRUE or FALSE. Got: $s"
  exit 1
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

    pair1=$path_hic_A
    if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
  
    stat1=$path_mgs_A
    if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi

    if [ "$r" -le 5000 ]; then
        q=3
    else
        q=6
    fi

    # Call Juicer APA https://github.com/aidenlab/juicer/wiki/APA
    echo -e "\nBeginning APA...\n"

    $juicer_tools apa -k NONE -n 0 -r "$r" -q "$q" -w $win_size "$pair1" "$tmp_binified" "$out_dir" &> /dev/null

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
    
    Rscript $aqua_dir/plot_APA.r \
      "$out_dir/APA_formatted.txt" \
      "${basename}" \
      "${c}" \
      "${d}" \
      "$out_dir/${basename}_APA.pdf" \
      ${win_size} \
      "${r}" \
      "${tmp_binified}" \
      "${Q}" \
      "$stat1" \
      "${l}" \
      "$num_loops" \
      "$out_dir" \
      "${e}" \
      "${s}"

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

    out_dir_A="$out_dir/$basename_A"
    mkdir -p "$out_dir_A"

    out_dir_B="$out_dir/$basename_B"
    mkdir -p "$out_dir_B"

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

    if [ "$r" -le 5000 ]; then
        q=3
    else
        q=6
    fi

    # https://github.com/aidenlab/juicer/wiki/APA

    echo -e "\nBeginning APA...\n"

    $juicer_tools apa -k NONE -n 0 -r "$r" -q "$q" -w $win_size "$pair1" "$tmp_binified" "$out_dir_A" &> /dev/null
    out_mat_A="$out_dir_A/$r/gw/APA.txt"
    if [ -f "$out_mat_A" ]; then
        mv "$out_mat_A" "$out_dir/APA_raw_${basename_A}.txt"
    else
        echo -e "\nAPA failed for sample A (no output file produced). Please check input files and parameters.\n"
        exit 1
    fi

    $juicer_tools apa -k NONE -n 0 -r "$r" -q "$q" -w $win_size "$pair2" "$tmp_binified" "$out_dir_B" &> /dev/null
    out_mat_B="$out_dir_B/$r/gw/APA.txt"
    if [ -f "$out_mat_B" ]; then
        mv "$out_mat_B" "$out_dir/APA_raw_${basename_B}.txt"
    else
        echo -e "\nAPA failed for sample B (no output file produced). Please check input files and parameters\n"
        exit 1
    fi

    #----------------------------------

    awk -F'[][,[]' '{ for(i=2; i<NF-1; i++) printf "%s\t", $i; printf "%s\n", $(NF-1) }' "$out_dir"/APA_raw_$basename_A.txt > "$out_dir"/APA_formatted_$basename_A.txt

    awk -F'[][,[]' '{ for(i=2; i<NF-1; i++) printf "%s\t", $i; printf "%s\n", $(NF-1) }' "$out_dir"/APA_raw_$basename_B.txt > "$out_dir"/APA_formatted_$basename_B.txt

    # Sum values in each formatted matrix
    matrix_sum_A=$(awk '{ for(i=1;i<=NF;i++) sum+=$i } END { print sum }' "$out_dir/APA_formatted_${basename_A}.txt")
    matrix_sum_B=$(awk '{ for(i=1;i<=NF;i++) sum+=$i } END { print sum }' "$out_dir/APA_formatted_${basename_B}.txt")

    # Check if both matrices are all zeros
    if [ "$matrix_sum_A" -eq 0 ] && [ "$matrix_sum_B" -eq 0 ]; then
        echo "Both APA matrices are all zeros - no valid contacts were detected."
        exit 1
    fi

    Rscript $aqua_dir/plot_APA.r \
      "$out_dir/APA_formatted_$basename_A.txt" \
      "$out_dir/APA_formatted_$basename_B.txt" \
      "${basename_A}" "${basename_B}" \
      "${c}" \
      "${d}" \
      $"$out_dir/${basename_A}_v_${basename_B}_APA.pdf" \
      ${win_size} \
      "${r}" \
      "${tmp_binified}" \
      "${Q}" \
      "$stat1" \
      "$stat2" \
      "${l}" \
      "$num_loops" \
      "$out_dir" \
      "${e}" \
      "${s}"

    # Delete the specified directories and files after R script is called
    rm -r "$out_dir_A"
    rm -r "$out_dir_B"
    rm "$out_dir"/APA_raw_$basename_A.txt
    rm "$out_dir"/APA_raw_$basename_B.txt
    rm "$out_dir"/APA_formatted_$basename_A.txt
    rm "$out_dir"/APA_formatted_$basename_B.txt
    rm "$tmp_binified"
    rm "$tmp_trimmed"

fi
