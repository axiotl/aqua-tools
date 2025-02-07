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
    echo " [  -O|--out-dir             ]  : Name of the directory to store the output APA plot in"
    echo " [  -B|--sample2             ]  : For two sample delta plots, name of the second sample"
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
      "--help")                 set -- "$@" "-h" ;;
       *)                       set -- "$@" "$arg"
  esac
done

r=5000
c="no_cap"
d="no_cap"
e="no_cap"
l=FALSE

while getopts ":P:A:G:O:B:Q:r:c:e:d:f:l:h" OPT
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



# If no sample B, the variable is empty:
if [ -z "$B" ]; then
  case $B in
  B) B="";;
  esac
fi

#----------------------------------

if [[ -z $P ]];
then
    usage
    exit
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

#----------------------------------

if [[ -z $A ]];
then
    usage
    exit
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

#----------------------------------

if [[ -z $G ]];
then
    usage
    exit
fi

#----------------------------------

# Check if there is a parameter provided for -B
if [[ -n $B ]]; then
    # Get sample directory B
    sample_dir_B=$(get_sample_directory "$B")
    if [ $? -eq 0 ]; then
        version_dir_B=$(basename "$sample_dir_B")
        base_dir_B=$(basename "$(dirname "$sample_dir_B")")
        B="${base_dir_B}/${version_dir_B}"
    else
        echo "Sample directory not found. Exiting."
        exit 1
    fi
fi

# Determine the generated output name based on whether -O is provided
if [[ -z $O ]]; then
    random_num=$RANDOM
    
    if [[ -n $B ]]; then
        # For a two-sample test
        O="${version_dir_A}_v_${version_dir_B}_APA_${random_num}"
    else
        # For a one-sample test
        O="${version_dir_A}_APA_${random_num}"
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



###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if [ -z "$B" ]
then 

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

pair1=$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic
if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi

#----------------------------------

stat1=$data_dir/$G/$A/mergeStats.txt
if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi

# Call Juicer APA https://github.com/aidenlab/juicer/wiki/APA
echo -e "\nBeginning APA...\n"
$juicer_tools apa --threads 1 -k NONE -n 0 -r "$r" -w $win_size "$pair1" "$P" "$out_dir" &> /dev/null

# Copying the generated APA.txt to APA_raw.txt
out_mat=$out_dir/$r/gw/APA.txt
mv "$out_mat" "$out_dir"/APA_raw.txt


awk -F'[][,[]' '{ for(i=2; i<NF-1; i++) printf "%s\t", $i; printf "%s\n", $(NF-1) }' "$out_dir"/APA_raw.txt > "$out_dir"/APA_formatted.txt

Rscript $aqua_dir/plot_APA.r "$out_dir/APA_formatted.txt" "${version_dir_A}" "${c}" "${d}" "$out_dir/${version_dir_A}_APA.pdf" ${win_size} "${r}" "${P}" "${Q}" "$stat1" "${l}" "$num_loops" "$out_dir" "${e}"

# Delete the specified directories and files after R script is called
rm -r "$out_dir/$r"
rm "$out_dir"/APA_raw.txt
rm "$out_dir"/APA_formatted.txt

fi




###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if [ -n "$B" ]
then 

    # Do all necessary parameter checks

    num_loops=`cat "$P" | wc -l`

    out_dir=$O
    mkdir -p "$out_dir"

    #----------------------------------

    pair1=$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic
    pair2=$data_dir/$G/$B/${version_dir_B}.allValidPairs.hic
    if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
    if [[ ! -f $pair2 ]]; then echo "cannot find $pair2"; exit 1; fi

    #----------------------------------

    stat1=$data_dir/$G/$A/mergeStats.txt
    stat2=$data_dir/$G/$B/mergeStats.txt
    if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi
    if [[ ! -f $stat2 ]]; then echo "cannot find $stat2"; exit 1; fi

    #----------------------------------

    # https://github.com/aidenlab/juicer/wiki/APA

    echo -e "\nBeginning APA...\n"

    $juicer_tools apa --threads 1 -k NONE -n 0 -r "$r" -w $win_size "$pair1" "$P" "$out_dir" &> /dev/null
    out_mat_A=$out_dir/$r/gw/APA.txt
    mv "$out_mat_A" "$out_dir"/APA_raw_${version_dir_A}.txt

    $juicer_tools apa --threads 1 -k NONE -n 0 -r "$r" -w $win_size "$pair2" "$P" "$out_dir" &> /dev/null
    out_mat_B=$out_dir/$r/gw/APA.txt
    mv "$out_mat_B" "$out_dir"/APA_raw_${version_dir_B}.txt

    #----------------------------------


    awk -F'[][,[]' '{ for(i=2; i<NF-1; i++) printf "%s\t", $i; printf "%s\n", $(NF-1) }' "$out_dir"/APA_raw_$version_dir_A.txt > "$out_dir"/APA_formatted_$version_dir_A.txt

    awk -F'[][,[]' '{ for(i=2; i<NF-1; i++) printf "%s\t", $i; printf "%s\n", $(NF-1) }' "$out_dir"/APA_raw_$version_dir_B.txt > "$out_dir"/APA_formatted_$version_dir_B.txt

    Rscript $aqua_dir/plot_APA.r "$out_dir/APA_formatted_$version_dir_A.txt" "$out_dir/APA_formatted_$version_dir_B.txt" "${version_dir_A}" "${version_dir_B}" "${c}" "${d}" $"$out_dir/${version_dir_A}_v_${version_dir_B}_APA.pdf" ${win_size} "${r}" "${P}" "${Q}" "$stat1" "$stat2" "${l}" "$num_loops" "$out_dir" "${e}"

    # Delete the specified directories and files after R script is called
    rm -r "$out_dir/$r"
    rm "$out_dir"/APA_raw_$version_dir_A.txt
    rm "$out_dir"/APA_raw_$version_dir_B.txt
    rm "$out_dir"/APA_formatted_$version_dir_A.txt
    rm "$out_dir"/APA_formatted_$version_dir_B.txt

fi
