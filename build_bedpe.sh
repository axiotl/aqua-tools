#!/bin/bash

data_dir=~/lab-data
aqua_dir=~/aqua_tools

function usage {
    echo "usage : build_bedpe.sh \\"
    echo "  -A PATH_TO_FIRST_BED_FILE \\"
    echo "  -B PATH_TO_SECOND_BED_FILE \\"
    echo " [-T PATH_TO_TAD_FILE] \\"
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
    echo "  [-D|--drop_dist] DROP_DISTANCE             : minimum distance between centers of pairs used to drop results"
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
      "--drop_dist") set -- "$@" "-D" ;;
      "--help")      set -- "$@" "-h" ;;
       *)            set -- "$@" "$arg"
  esac
done



while getopts ":A:B:T:D:h" OPT
do
    case $OPT in
  A) A=$OPTARG;;
  B) B=$OPTARG;;
  T) T=$OPTARG;;
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

if [ -z "$T" ]; then T="NULL"; fi
if [ -z "$D" ]; then D=0;      fi

# printf "A: %s\nB: %s\nT: %s\nD: %s\n" $A $B $T $D

Rscript $aqua_dir/build_bedpe.r $A $B $T $D
