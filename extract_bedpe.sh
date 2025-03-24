#!/bin/bash

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

sample_sheet="$HOME/setup/sample_sheet.txt"


function usage {
    echo "usage : extract_bedpe.sh \\"
    echo "  -A SAMPLE_NAME \\"
    echo "  -G GENOME \\"
    echo "  -R RANGE \\"
    echo " [-h]"
    echo "Use option -h|--help for more information"
}


function help {
    echo
    echo "Obtain clusters of interacting loops based on score threshold"
    echo "Scores are calculated using inherent normalization (hello@axiotl.com)"
    echo "Prints in bedpe format to standard out."
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -A|--sample1        : Name of the sample as it appears on the Tinkerbox"
    echo "   -G|--genome         : Genome build used for sample processing"
    echo "   -R|--range          : Range to obtain clusters from, in chr:start:end format"
    echo "  [-r|--resolution   ] : Resolution in base pairs. Only 5000 and 1000 supported. Default 5000"
    echo "  [-T|--TAD          ] : Full path to TAD file, the boundaries of which will be used to obtain clusters"
    echo "  [-S|--score        ] : Inherent score to seed cluster formation. Default 0.7"
    echo "  [-m|--mode         ] : Shape of bedpe to be called: loop, flare, minimal or glob. Default glob"
    echo "  [   --radius       ] : Bin distance units to search for neighbours. Default 1"
    echo "  [   --min_dist     ] : Distance in basepairs to filter out extracted elements. Default 0"
    echo "  [-h|--help         ]   Help message"
    exit;
}

# no arguments
if [ $# -lt 1 ]
    then
    usage
    exit
fi

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--sample1")    set -- "$@" "-A" ;;
      "--genome")      set -- "$@" "-G" ;;
      "--range")       set -- "$@" "-R" ;;
      "--resolution")  set -- "$@" "-r" ;;
      "--TAD")         set -- "$@" "-T" ;;
      "--score")       set -- "$@" "-S" ;;
      "--radius")      set -- "$@" "-d" ;;
      "--min_dist")    set -- "$@" "-x" ;;
      "--mode")        set -- "$@" "-m" ;;
      "--help")        set -- "$@" "-h" ;;
       *)              set -- "$@" "$arg"
  esac
done


r=5000
R="range"
T="NULL"
S=0.7
d=1
m="glob"
x=0
while getopts ":A:G:R:r:T:S:d:x:m:h" OPT
do
    case $OPT in
  A) A=$OPTARG;;
  G) G=$OPTARG;;
  R) R=$OPTARG;;
  r) r=$OPTARG;;
  T) T=$OPTARG;;
  S) S=$OPTARG;;
  d) d=$OPTARG;;
  x) x=$OPTARG;;
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


if [[ -z $G ]];
then
    usage
    exit
fi



if [ "$R" == "range"    ] && [ "$T" == "NULL"   ] ; then
    echo "Please supply at least a range or TAD file"
    exit
fi


if [ "$m" != "loop" ] && [ "$m" != "flare" ] && [ "$m" != "glob" ] && [ "$m" != "minimal" ]; then
    echo "mode stritcly either loop, flare, minimal or glob"
    exit 1
fi


###########################################################################
###########################################################################
###                                                                     ###
###                                RANGE                                ###
###                                                                     ###
###########################################################################
###########################################################################


if [ "$R" != "range"    ] && [ "$T" == "NULL"   ] ; then

    analysis_type="range"

    range_text=`echo "$R" | sed 's/:/_/g'`
    chr=`echo "$R" | cut -d":" -f1`
    start=`echo "$R" | cut -d":" -f2`
    end=`echo "$R" | cut -d":" -f3`


    path_hic=$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic
    if [[ ! -f $path_hic ]]; then echo "Cannot find .hic file"; exit 1; fi

    sample_dir=$data_dir/$G/$A

    Rscript $aqua_dir/extract_bedpe.r $analysis_type $path_hic $chr $start $end $r $S $sample_dir $d $m $x
fi



###########################################################################
###########################################################################
###                                                                     ###
###                                 TAD                                 ###
###                                                                     ###
###########################################################################
###########################################################################


if [ "$R" == "range"    ] && [ "$T" != "NULL"   ] ; then

    analysis_type="TAD"

    path_hic=$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic
    if [[ ! -f $path_hic ]]; then echo "Cannot find .hic file"; exit 1; fi

    sample_dir=$data_dir/$G/$A

    Rscript $aqua_dir/extract_bedpe.r $analysis_type $path_hic $T $r $S $sample_dir $d $m $x
fi
