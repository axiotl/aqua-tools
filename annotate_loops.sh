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
    echo "    -P|--bedpe           : Full path to the bedpe file you want to annotate, without headers! "
    echo "    -A|--sample1         : Name of the sample"
    echo "    -G|--genome          : The genome build the sample(s) has been processed using. Strictly hg19, hg38, or mm10"
    echo "  [ -Q|--norm          ] : Which normalization to use. Strictly 'none', 'cpm', 'aqua', or 'abc' in lower case. Non-spike-in samples default to cpm. Spike-in samples default to aqua. "
    echo "  [ -B|--sample2       ] : For two sample delta contacts, name of the second sample"
    echo "  [ -R|--resolution    ] : Resolution of sample in basepairs, using which the contact values should be calculated. Default 5000. Accepted resolutions- 1000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000"
    echo "  [ -f|--formula       ] : Arithmetic to use to report contact calues. Options- center, max, average, sum. Default = center"
    echo "  [ -F|--fix           ] : If FALSE, reports new coordinates based on arithmetic center or max. Default = TRUE"
    echo "  [    --expand        ] : Expands 1D bedpe feet in both directions based on supplied value (in bin units). Default = 0"
    echo "  [ -i|--inherent      ] : If TRUE, hic values transformed to inherent units. Default = FALSE"
    echo "  [ -m|--preserve_meta ] : If TRUE, bedpe metadata columns will be preserved. Default = TRUE"
    echo "  [ -h|--help          ] : Help message. Primer can be found at https://rb.gy/zyfjxc"
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
      "--genome")        set -- "$@" "-G" ;;
      "--norm")          set -- "$@" "-Q" ;;
      "--sample2")       set -- "$@" "-B" ;;
      "--resolution")    set -- "$@" "-R" ;;
      "--formula")       set -- "$@" "-f" ;;
      "--fix")           set -- "$@" "-F" ;;
      "--split")         set -- "$@" "-S" ;; 
      "--shrink_wrap")   set -- "$@" "-s" ;;
      "--padding")       set -- "$@" "-p" ;;
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

while getopts ":P:A:G:Q:B:R:f:F:S:s:p:e:I:m:h" OPT
do
    case $OPT in
  P) P=$OPTARG;;
  A) A=$OPTARG;;
  G) G=$OPTARG;;
  Q) Q=$OPTARG;;
  B) B=$OPTARG;;
  R) R=$OPTARG;;
  f) f=$OPTARG;;
  F) F=$OPTARG;;
  S) S=$OPTARG;;
  s) s=$OPTARG;;
  p) p=$OPTARG;;
  e) e=$OPTARG;;
  I) I=$OPTARG;;
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


if [ -z "$A" ]; then
    echo "No sample name provided. Exiting."
    usage
    exit 1
fi

if [ -z "$G" ]; then
    echo "No genome build provided. Exiting."
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


# If no sample B, the variable is empty:
if [ -z "$B" ]; then
  case $B in
  B) B="";;
  esac
fi


# Do all necessary parameter checks
#----------------------------------

if [[ -z $P || ! -f $P ]]; then
    echo "Invalid file specified in '$P'"
    usage
    exit 1
fi



if [[ "$i" == "TRUE" && ( "$Q" != "blank" && "$Q" != "none" ) ]];
then
    echo -e "\n\n--inherent is not compatible with aqua, cpm, or abc normalization. Choose either --norm $Q or --inherent TRUE\n"
    exit 1
fi


#----------------------------------


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

if [ -z "$B" ]
then
        Rscript \
        $aqua_dir/annotate_loops.r \
          $P \
          $R \
          $data_dir/$G/$A/${version_dir_A}.allValidPairs.hic \
          $data_dir/$G/$A/mergeStats.txt \
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

if [ -n "$B" ]
then

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

 
        Rscript \
        $aqua_dir/annotate_loops.r \
          $P \
          $R \
          $data_dir/$G/$A/${version_dir_A}.allValidPairs.hic \
          $data_dir/$G/$A/mergeStats.txt \
          $data_dir/$G/$B/${version_dir_B}.allValidPairs.hic \
          $data_dir/$G/$B/mergeStats.txt \
          $Q \
          $G \
          $num_loops \
          $f $F $S $s $p $e $data_dir/$G/$A $data_dir/$G/$B $m $i
    
fi
