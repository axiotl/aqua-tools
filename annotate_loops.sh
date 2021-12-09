#!/bin/bash

data_dir=~/lab-data
aqua_dir=~/aqua_tools

function usage {
    echo -e "usage: "
    echo -e "  annotate_pairs.sh \\"
    echo -e "    -G GENOME_BUILD \\"
    echo -e "    -P PATH_TO_BEDPE_YOU_WANT_TO_ANNOTATE \\"
    echo -e "    -A NAME_OF_FIRST_SAMPLE \\"
    echo -e "   [-Q USE_AQUA_FACTORS] \\"
    echo -e "   [-B NAME_OF_SECOND_SAMPLE] \\"
    echo -e "   [-R RESOLUTION] \\"
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
    echo "   -P|--pair           PATH_TO_GENOMIC_PAIRS_FILE  : Full path to the bedpe (pairs) file you want to annotate. If bedpe contains headers, see manual on how to remove them!"
    echo "   -A|--sample1        NAME_OF_FIRST_SAMPLE        : Name of the sample you want to use the AQuA normalized values to annotate it with, as it appears on the Tinkerbox"
    echo "   -G|--genome         GENOME_BUILD                : The genome build the sample(s) has been processed using. Strictly hg19 or hg38"
    echo "  [-Q|--AQuA       ]   USE_AQUA_FACTORS            : TRUE or FALSE"
    echo "  [-B|--sample2    ]   NAME_OF_SECOND_SAMPLE       : The name of the second sample. If triggered, calculates the delta AQuA normalized values from both samples for that pair. Useful in case vs control."
    echo "  [-R|--resolution ]   RESOLUTION                  : Resolution of sample in basepairs, using which the contact values should be calculated. Default 5000. Accepted resolutions- 1000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000"
    echo "  [-h|--help       ]   Help message"
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
      "--pair")         set -- "$@" "-P" ;;
      "--sample1")      set -- "$@" "-A" ;;
      "--genome")       set -- "$@" "-G" ;;
      "--AQuA")         set -- "$@" "-Q" ;;
      "--sample2")      set -- "$@" "-B" ;;
      "--resolution")   set -- "$@" "-R" ;;  
      "--help")         set -- "$@" "-h" ;;
       *)               set -- "$@" "$arg"
  esac
done

R=5000

while getopts ":P:A:G:Q:B:R:h" OPT
do
    case $OPT in
  P) P=$OPTARG;;
  A) A=$OPTARG;;
  G) G=$OPTARG;;
  Q) Q=$OPTARG;;
  B) B=$OPTARG;;
  R) R=$OPTARG;;
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


if [ -z "$B" ]; then
  case $B in
  B) B="";;
  esac
fi


if [[ -z $P ]];
then
    usage
    exit
fi

if [[ -z $Q ]]; then Q="TRUE"; fi



col1=$(head -n 1  $P | cut -f 1 )
col1="${col1:0:3}"

col4=$(head -n 1  $P | cut -f 4 )
col4="${col4:0:3}"

if [[ $col1 != $col4  ]]; then
  echo "Please provide a 6-col bedpe file without headers"
  exit
fi

if [[ -z $A ]];
then
    usage
    exit
fi

if [[ -z $G ]];
then
    usage
    exit
fi


if [ -z "$B" ]
then

#echo "Commencing one sample analysis..."
Rscript \
$aqua_dir/annotate_loops.r \
  $P \
  $R \
  $data_dir/$G/$A/$A.allValidPairs.hic \
  $data_dir/$G/$A/mergeStats.txt \
  $Q
fi

if [ -n "$B" ]
then 

#echo "Commencing two sample analysis..."  

Rscript \
$aqua_dir/annotate_loops.r \
  $P \
  $R \
  $data_dir/$G/$A/$A.allValidPairs.hic \
  $data_dir/$G/$A/mergeStats.txt \
  $data_dir/$G/$B/$B.allValidPairs.hic \
  $data_dir/$G/$B/mergeStats.txt \
  $Q
fi


