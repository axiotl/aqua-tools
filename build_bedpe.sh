#!/bin/bash

data_dir=~/lab-data
aqua_dir=~/aqua_tools

function usage {
    echo "usage : build_bedpe.sh \\"
    echo "  -A PATH_TO_FIRST_BED_FILE \\"
    echo "  -B PATH_TO_SECOND_BED_FILE \\"
    echo " [-T PATH_TO_TAD_FILE] \\"
    echo " [-d MIN_DROP_DISTANCE] \\"
    echo " [-D MAX_DROP_DISTANCE] \\"
    echo " [-h]"
    echo "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Builds pairs between elements in two bed files."
    echo "Pairs can be constrained by a third bed file (usually TADs)"
    echo "or by a minimum distance between them."
    echo "Prints pairs in bedpe format to standard out."
    echo
    echo "Pairs can be use to query .hic files with annotate_loops.sh"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -A|--bed_A      PATH_TO_FIRST_BED_FILE    : Path of the first bed file you want to build the bedpe of"
    echo "   -B|--bed_B      PATH_TO_SECOND_BED_FILE   : Path of the second bed file you want to build the bedpe of"
    echo "  [-T|--TAD      ] PATH_TO_TAD_FILE          : Path of the TAD file that checks if genomic regions fall within them"
    echo "  [-d|--min_dist ] MIN_DROP_DISTANCE         : minimum distance between pairs used to drop results. Default 0 bp"
    echo "  [-D|--max_dist ] MAX_DROP_DISTANCE         : maximum distance between pairs used to drop results. Default 5 Mb"
    echo "  [-h|--help     ] Help message"
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
      "--bed_A")     set -- "$@" "-A" ;;
      "--bed_B")     set -- "$@" "-B" ;;
      "--TAD")       set -- "$@" "-T" ;;
      "--min_dist")  set -- "$@" "-d" ;;
      "--max_dist")  set -- "$@" "-D" ;;
      "--help")      set -- "$@" "-h" ;;
       *)            set -- "$@" "$arg"
  esac
done



while getopts ":A:B:T:d:D:h" OPT
do
    case $OPT in
  A) A=$OPTARG;;
  B) B=$OPTARG;;
  T) T=$OPTARG;;
  d) d=$OPTARG;;
  D) D=$OPTARG;;
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

if [ -z "$T" ]; then T="NULL";  fi
if [ -z "$d" ]; then d=0;       fi
if [ -z "$D" ]; then D=5000000; fi



Rscript $aqua_dir/build_bedpe.r $A $B $T $d $D
