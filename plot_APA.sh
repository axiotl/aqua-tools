#!/bin/bash

OPTIND=1
win_size=10
data_dir="${LAB_DATA_DIR:-$HOME/lab-data}"
aqua_dir="${AQUA_TOOLS_DIR:-$HOME/aqua_tools}"

sample_sheet="$HOME/setup/sample_sheet.txt"


function usage {
    echo -e "usage : plot_APA.sh -P PATH_TO_GENOMIC_PAIRS_FILE_YOU_WANT_TO_PLOT -A NAME_OF_FIRST_SAMPLE -G GENOME_BUILD -O OUTPUT_DIRECTORY [-B NAME_OF_SECOND_SAMPLE] [-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Create APA plot with CPM/AQuA/inherent normalized contact values"
    echo
    echo "Multiple samples can be passed to -A or -B to average as a cohort, e.g. --sample1 s1 s2 s3"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -P|--bedpe                  : Path to the bedpe file, without headers"
    echo "    -A|--sample1                : Name of the sample as it appears on the Tinkerbox."
    echo "    -G|--genome                 : Genome build used for sample processing"
    echo " [  -O|--out-dir             ]  : Name of the directory to store the output APA plot in"
    echo " [  -B|--sample2             ]  : For two sample delta plots, name of the second sample."
    echo " [  -Q|--norm                ]  : Which normalization to use: none, cpm, or aqua in lower case"
    echo " [  -r|--resolution          ]  : Bin size you want to use for the APA plots. Default 5000"
    echo " [     --max_cap             ]  : Upper limit of the plot range. Defaults to max bin value"
    echo " [     --min_cap             ]  : Lower limit of the plot range. Default 0"
    echo " [     --max_cap_delta       ]  : Upper limit of delta plot range. Defaults to max bin value"
    echo " [     --loop_norm           ]  : If TRUE, normalizes APA values by loop count in bedpe. Default FALSE if --inherent FALSE, else TRUE"
    echo " [  -s|--scores              ]  : If TRUE, plots peak value and ratios of peak:corner means. Default TRUE"
    echo " [  -i|--inherent            ]  : If TRUE, normalize contacts using inherent normalization. Default FALSE"
    echo " [     --sample1_dir         ]  : If not using the tinkerbox, specify the full path to the directory containing sample data"
    echo " [     --sample2_dir         ]  : If not using the tinkerbox, full path to the second sample directory"
    echo " [  -h|--help                ]   Help message"
    exit;
}

if [ $# -lt 1 ]
    then
    usage
    exit
fi

COHORT_SAMPLES_1=()
COHORT_SAMPLES_2=()

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
      "--scores")               set -- "$@" "-s" ;;
      "--inherent")             set -- "$@" "-i" ;;
      "--sample1_dir")          set -- "$@" "-H" ;;
      "--sample2_dir")          set -- "$@" "-I" ;;
      "--help")                 set -- "$@" "-h" ;;
       *)                       set -- "$@" "$arg"
  esac
done

r=5000
c="no_cap"
d="no_cap"
e="no_cap"
l="blank"
s=TRUE
i=FALSE
H=blank
I=blank

while getopts ":P:A:G:O:B:Q:r:c:e:d:l:s:i:H:I:h" OPT
do
    case $OPT in
  P) P=$OPTARG;;
  A) 
    # Accept multiple samples
    COHORT_SAMPLES_1+=("$OPTARG")
    while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
        COHORT_SAMPLES_1+=("${!OPTIND}")
        ((OPTIND++))
    done
    ;;
  G) G=$OPTARG;;
  B) 
    # Accept multiple samples
    COHORT_SAMPLES_2+=("$OPTARG")
    while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
        COHORT_SAMPLES_2+=("${!OPTIND}")
        ((OPTIND++))
    done
    ;;
  O) O=$OPTARG;;
  r) r=$OPTARG;;
  c) c=$OPTARG;;
  Q) Q=$OPTARG;;
  l) l=$OPTARG;;
  d) d=$OPTARG;;
  e) e=$OPTARG;;
  s) s=$OPTARG;;
  i) i=$OPTARG;;
  H) H=$OPTARG;;
  I) I=$OPTARG;;
  h) help ;;
  \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
  :)
      case "$OPTARG" in
        l) echo -e "\n--loop_norm requires a value of TRUE or FALSE (--loop_norm TRUE)\n" ;;
        s) echo -e "\n--scores requires a value of TRUE or FALSE (--scores FALSE)\n" ;;
        i) echo -e "\n--inherent requires a value of TRUE or FALSE (--inherent TRUE)\n" ;;
        *) echo -e "\nOption -$OPTARG requires an argument.\n" ;;
      esac
      usage
      exit 1
      ;;
  esac
done

# check for missing true/false flags
validate_bool_param() {
    local name="$1"
    local value="$2"
    if [[ -z "$value" || "$value" == -* ]]; then
        echo -e "\n--${name} requires a value of TRUE or FALSE (--${name} TRUE)\n"
        exit 1
    fi
    if [[ "$value" != "TRUE" && "$value" != "FALSE" ]]; then
        echo -e "\n--${name} must be TRUE or FALSE. Got: '$value'\n"
        exit 1
    fi
}

# Resolve loop_norm default based on inherent flag
if [[ "$l" == "blank" ]]; then
    if [[ "$i" == "TRUE" ]]; then
        l="TRUE"
    else
        l="FALSE"
    fi
fi

validate_bool_param "loop_norm" "$l"
validate_bool_param "scores" "$s"
validate_bool_param "inherent" "$i"

#----------------------------------

if [[ -z $G ]];
then
    usage
    exit
fi

#----------------------------------

if [ ${#COHORT_SAMPLES_1[@]} -eq 0 ] && [ "$H" == "blank" ];
then
    usage
    exit
fi

if [ ${#COHORT_SAMPLES_1[@]} -gt 0 ] && [ "$H" != "blank" ]; then
    echo "Need either sample1 name or sample1 directory path, not both. Exiting."
    usage
    exit 1
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
  local sample_dir_input="$1"
  local sample_label="$2"
  local hic_path mgs_path sample_name

  if [[ ! -d "$sample_dir_input" ]]; then
      echo "Directory '$sample_dir_input' does not exist. Exiting."
      exit 1
  fi

  sample_name=$(basename "$sample_dir_input")

  hic_path="${sample_dir_input}/${sample_name}.hic"
  mgs_path="${sample_dir_input}/${sample_name}.mergeStats.txt"

  if [[ ! -f "$hic_path" ]]; then
      echo "Missing .hic file: $hic_path"
      exit 1
  fi

  if [[ ! -f "$mgs_path" ]]; then
      echo "Missing mergeStats.txt file: $mgs_path"
      exit 1
  fi
}

# Resolve all cohort A samples
hic_paths_A=()
stat_paths_A=()
version_dirs_A=()

if [ "$H" != "blank" ]; then
    # Using --sample1_dir
    validate_samples_and_stats "$H" "sample1"
    local_name=$(basename "$H")
    version_dirs_A+=("$local_name")
    hic_paths_A+=("$H/${local_name}.hic")
    stat_paths_A+=("$H/${local_name}.mergeStats.txt")
else
    for sample_name in "${COHORT_SAMPLES_1[@]}"; do
        sample_dir=$(get_sample_directory "$sample_name")
        if [ $? -ne 0 ]; then
            echo "Sample directory not found for '$sample_name'. Exiting."
            exit 1
        fi
        vdir=$(basename "$sample_dir")
        version_dirs_A+=("$vdir")
        hic_paths_A+=("$sample_dir/${vdir}.allValidPairs.hic")
        stat_paths_A+=("$sample_dir/mergeStats.txt")
    done
fi

# Same for cohort B if present
if [ ${#COHORT_SAMPLES_2[@]} -gt 0 ] || [ "$I" != "blank" ]; then
    hic_paths_B=()
    stat_paths_B=()
    version_dirs_B=()
    if [ "$I" != "blank" ]; then
        # Using --sample2_dir
        validate_samples_and_stats "$I" "sample2"
        local_name=$(basename "$I")
        version_dirs_B+=("$local_name")
        hic_paths_B+=("$I/${local_name}.hic")
        stat_paths_B+=("$I/${local_name}.mergeStats.txt")
    else
        for sample_name in "${COHORT_SAMPLES_2[@]}"; do
            sample_dir=$(get_sample_directory "$sample_name")
            if [ $? -ne 0 ]; then
                echo "Sample directory not found for '$sample_name'. Exiting."
                exit 1
            fi
            vdir=$(basename "$sample_dir")
            version_dirs_B+=("$vdir")
            hic_paths_B+=("$sample_dir/${vdir}.allValidPairs.hic")
            stat_paths_B+=("$sample_dir/mergeStats.txt")
        done
    fi
fi

#----------------------------------

# Checks for $P
if [[ -z $P ]];
then
    usage
    exit
fi

if [[ ! -f "$P" ]]; then
    echo -e "\n'$P' is not a valid file\n"
    exit
fi

first_line=$(head -n 1 "$P")

# Check for tab-separated bedpe with at least 6 columns
tab_check=$(echo "$first_line" | awk -F'\t' '{print NF}')
if [[ "$tab_check" -lt 6 ]]; then
    echo ".bedpe file must have at least 6 TAB-delimited columns"
    exit 1
fi

# Check for a header with column 2 as numeric
col2=$(echo "$first_line" | cut -f2)
if ! [[ "$col2" =~ ^[0-9]+$ ]]; then
    echo "Please provide a 6-col BEDPE file without headers"
    exit 1
fi

awk_output=$(awk '
        BEGIN {FS="\t"; errors=0; reversed=0; inter=0; max_print=5}
        NF < 6 {
            if (errors < max_print) print "BEDPE error: Line " NR " has fewer than 6 columns"
            errors++
        }
        $3 < $2 || $6 < $5 {
            if (reversed < max_print) print "BEDPE error: Reversed feet in line " NR
            reversed++
        }
        END {
            if (errors > 0) {
                printf "BEDPE error: %d line(s) have fewer than 6 columns.\n", errors
                exit 1
            }
            if (reversed > 0) {
                printf "BEDPE error: %d line(s) have reversed feet (require col3 > col2 and col6 > col5).\n", reversed
                exit 1
            }
            exit 0
        }' "$P")

if [[ $? -ne 0 ]]; then
    echo "$awk_output"
    exit 1
fi

#------ Check for mixed cis/trans bedpe ------
 
bedpe_type=$(awk '
    BEGIN {FS="\t"; cis=0; trans=0}
    {
        if ($1 == $4) cis++
        else trans++
    }
    END {
        if (cis > 0 && trans > 0) {
            printf "BEDPE error: File contains both cis (%d) and trans (%d) pairs.\nPlease provide a bedpe with only cis or only trans pairs.\n", cis, trans
            exit 1
        }
        if (trans > 0) print "trans"
        else print "cis"
    }' "$P")
 
if [[ $? -ne 0 ]]; then
    echo "$bedpe_type"
    exit 1
fi

if [[ "$bedpe_type" == "trans" && "$i" == "TRUE" ]]; then
    echo -e "\nInherent normalization is not compatible with interchromosomal bedpe pairs."
    exit 1
fi

#----------------------------------

if [ ${#version_dirs_A[@]} -eq 1 ]; then
    label_A="${version_dirs_A[0]}"
else
    label_A="${version_dirs_A[0]}_cohort"
fi

if [ ${#COHORT_SAMPLES_2[@]} -gt 0 ] || [ "$I" != "blank" ]; then
    if [ ${#version_dirs_B[@]} -eq 1 ]; then
        label_B="${version_dirs_B[0]}"
    else
        label_B="${version_dirs_B[0]}_cohort"
    fi
fi

if [[ -z $O ]]; then
    random_num=$RANDOM
    if [ ${#COHORT_SAMPLES_2[@]} -gt 0 ] || [ "$I" != "blank" ]; then
        O="${label_A}_v_${label_B}_APA_${random_num}"
    else
        O="${label_A}_APA_${random_num}"
    fi
    echo -e "\nNo output directory name provided. Generating output directory..."
fi

#----------------------------------
# Validate normalization

if [ -n "$Q" ]; then
    if [ "$Q" != "aqua" ] && [ "$Q" != "cpm" ] && [ "$Q" != "none" ]; then
        echo -e "\nInvalid value for -Q|--norm. Please use 'aqua', 'cpm', or 'none' in lower case."
        exit 1
    fi
else
    Q="blank"
fi

# Check for aqua/cpm/inh normalization clash

if [[ "$Q" != "blank" && "$i" != "FALSE" ]]; then
  echo -e "\nInherent normalization is not compatible with other --norm methods. Continuing with --inherent TRUE ..."
  Q="blank"
  i="TRUE"
fi

#----------------------------------

# Inherent requires loop_norm TRUE
if [[ "$i" == "TRUE" && "$l" == "FALSE" ]]; then
    echo -e "\n--inherent TRUE requires --loop_norm TRUE. Setting --loop_norm TRUE ..."
    l="TRUE"
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

#---------------Enforce cis anchor ordering-------------------

# Enforce anchor1 start <= anchor2 start for cis pairs to
# match straw's convention for intrachromosomal data
# trans anchor ordering happens in R

tmp_ordered=$(mktemp /tmp/ordered_bedpe.XXXXXX)

awk '
BEGIN {
    OFS="\t"
}
{
    # Remove carriage return from all fields if present
    for (k=1; k<=NF; k++) gsub(/\r$/, "", $k);

    # For cis pairs only: swap so anchor1 start <= anchor2 start
    if ($1 == $4 && ($2 > $5 || ($2 == $5 && $3 > $6))) {
        temp1 = $1; temp2 = $2; temp3 = $3;
        $1 = $4; $2 = $5; $3 = $6;
        $4 = temp1; $5 = temp2; $6 = temp3;
    }
    print
}' "$P" > "$tmp_ordered"

awk_status=$?
if [ $awk_status -ne 0 ]; then
    echo "Error during bedpe anchor ordering. Exiting."
    rm -f "$tmp_ordered"
    exit $awk_status
fi

#---------------Binify-------------------

# Create temporary files in /tmp
tmp_trimmed=$(mktemp /tmp/trimmed_bedpe.XXXXXX)
tmp_binified=$(mktemp /tmp/binified_bedpe.XXXXXX)

# Create 6 col bedpe from the ordered file
cut -f1-6 "$tmp_ordered" > "$tmp_trimmed"

awk -v r="$r" 'BEGIN { OFS="\t" }
{
  len1 = $3 - $2
  len2 = $6 - $5

  if (len1 <= r && len2 <= r) {
    # too small to binify, print as is
    print $1, $2, $3, $4, $5, $6

  } else if (len1 > r && len2 > r) {
    # binify both anchors
    for (i = $2; i <= ($3 - r); i += r) {
      for (j = $5; j <= ($6 - r); j += r) {
        print $1, i, i + r, $4, j, j + r
      }
    }

  } else if (len1 > r) {
    # only binify anchor1
    for (i = $2; i <= ($3 - r); i += r) {
      print $1, i, i + r, $4, $5, $6
    }

  } else if (len2 > r) {
    # only binify anchor2
    for (j = $5; j <= ($6 - r); j += r) {
      print $1, $2, $3, $4, j, j + r
    }
  }
}' "$tmp_trimmed" > "$tmp_binified"

# Convert tmp_binified to absolute path
tmp_binified=$(readlink -f "$tmp_binified")

# Clean up the ordered temp file
rm -f "$tmp_ordered"

join_array() { local IFS=","; echo "$*"; }

###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if [ ${#COHORT_SAMPLES_2[@]} -eq 0 ] && [ "$I" == "blank" ]
then 

    #----------------------------------

    if [[ $d != "no_cap" ]]; then
        echo "Cannot have a delta cap in single-sample APA analysis"
        exit
    fi

    #----------------------------------

    num_loops=$(wc -l < "$P")

    out_dir="$(readlink -f "$O")"
    mkdir -p "$out_dir"

    #----------------------------------

    for f in "${hic_paths_A[@]}"; do
        if [[ ! -f "$f" ]]; then echo "cannot find $f"; exit 1; fi
    done

    for f in "${stat_paths_A[@]}"; do
        if [[ ! -f "$f" ]]; then echo "cannot find $f"; exit 1; fi
    done

    #----------------------------------

    hic_str_A=$(join_array "${hic_paths_A[@]}")
    stat_str_A=$(join_array "${stat_paths_A[@]}")
    label_str_A=$(join_array "${version_dirs_A[@]}")

    Rscript $aqua_dir/plot_APA.r \
      "${hic_str_A}" \
      "${label_str_A}" \
      "${c}" \
      "${d}" \
      "$out_dir/${label_A}_APA.pdf" \
      "${win_size}" \
      "${r}" \
      "${tmp_binified}" \
      "${Q}" \
      "${stat_str_A}" \
      "${l}" \
      "$num_loops" \
      "$out_dir" \
      "${e}" \
      "${s}" \
      "${i}" \
      "${tmp_trimmed}"

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


if [ ${#COHORT_SAMPLES_2[@]} -gt 0 ] || [ "$I" != "blank" ]
then

    num_loops=$(wc -l < "$P")

    out_dir="$(readlink -f "$O")"
    mkdir -p "$out_dir"

    for f in "${hic_paths_A[@]}"; do
        if [[ ! -f "$f" ]]; then echo "cannot find $f"; exit 1; fi
    done

    for f in "${stat_paths_A[@]}"; do
        if [[ ! -f "$f" ]]; then echo "cannot find $f"; exit 1; fi
    done

    for f in "${hic_paths_B[@]}"; do
        if [[ ! -f "$f" ]]; then echo "cannot find $f"; exit 1; fi
    done

    for f in "${stat_paths_B[@]}"; do
        if [[ ! -f "$f" ]]; then echo "cannot find $f"; exit 1; fi
    done

    #----------------------------------

    hic_str_A=$(join_array "${hic_paths_A[@]}")
    stat_str_A=$(join_array "${stat_paths_A[@]}")
    label_str_A=$(join_array "${version_dirs_A[@]}")

    hic_str_B=$(join_array "${hic_paths_B[@]}")
    stat_str_B=$(join_array "${stat_paths_B[@]}")
    label_str_B=$(join_array "${version_dirs_B[@]}")

    #----------------------------------

    Rscript $aqua_dir/plot_APA.r \
      "${hic_str_A}" \
      "${hic_str_B}" \
      "${label_str_A}" "${label_str_B}" \
      "${c}" \
      "${d}" \
      "$out_dir/${label_A}_v_${label_B}_APA.pdf" \
      "${win_size}" \
      "${r}" \
      "${tmp_binified}" \
      "${Q}" \
      "${stat_str_A}" \
      "${stat_str_B}" \
      "${l}" \
      "$num_loops" \
      "$out_dir" \
      "${e}" \
      "${s}" \
      "${i}" \
      "${tmp_trimmed}"

    # Clean up temp files
    rm "$tmp_binified"
    rm "$tmp_trimmed"

fi
