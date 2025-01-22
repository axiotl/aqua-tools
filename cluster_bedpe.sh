#!/bin/bash

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

temp_dir=$(mktemp -d cluster_bedpe-XXXXXXXXXX)


ctrlc_count=0

function no_sigint {

    let ctrlc_count++
    if [[ $ctrlc_count == 1 ]]; then
        sudo rm -rf $temp_dir
    else
        :
    fi

}
trap no_sigint EXIT



function usage {
    echo -e "usage : cluster_bedpe.sh -P PATH_TO_BEDPE_FILE [-F FLANK_IN_BASE_PAIR] [-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Given a bedpe file, clusters bedpe rows based on basepair intersections"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -P|--bedpe        : Path of the bedpe file you want to use"
    echo "  [-F|--flank     ]  : Genome distance in bp to increase the size of bedpe rows. Default is 0"
    echo "  [-h|--help      ]    Help message"
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
      "--bedpe")     set -- "$@" "-P" ;;
      "--flank")     set -- "$@" "-F" ;;
      "--help")      set -- "$@" "-h" ;;
       *)            set -- "$@" "$arg"
  esac
done

F=0
while getopts ":P:F:h" OPT
do
    case $OPT in
  P) P=$OPTARG;;
  F) F=$OPTARG;;
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




Rscript $aqua_dir/cluster_bedpe.r $P $F $temp_dir




trap "sudo rm -rf $temp_dir" EXIT
