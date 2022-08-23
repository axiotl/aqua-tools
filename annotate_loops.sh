#!/bin/bash

data_dir=~/lab-data
aqua_dir=~/aqua_tools
ram_disk=/media/ram_disk


ctrlc_count=0

function no_sigint {

    let ctrlc_count++
    if [[ $ctrlc_count == 1 ]]; then
        echo "Aborting process"
        if [ -d "$ram_disk" ]; then
        sudo umount "$ram_disk"
        cd /media
        sudo rm -r ram_disk
        else
        :
        fi
    fi

}

trap no_sigint SIGINT

function usage {
    echo -e "usage: "
    echo -e "  annotate_loops.sh \\"
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
    echo "   -P|--pair            PATH_TO_GENOMIC_PAIRS_FILE  : Full path to the bedpe (pairs) file you want to annotate, without headers! "
    echo "   -A|--sample1         NAME_OF_FIRST_SAMPLE        : Name of the sample you want to use as it appears on the Tinkerbox"
    echo "   -G|--genome          GENOME_BUILD                : The genome build the sample(s) has been processed using. Strictly hg19 or hg38"
    echo "  [-Q|--norm        ]   NORMALIZATION               : Which normalization to use. Strictly 'none', 'cpm' or 'aqua' in lower case"
    echo "  [-B|--sample2     ]   NAME_OF_SECOND_SAMPLE       : The name of the second sample. If triggered, calculates the delta contact values for that pair. Useful in case vs control"
    echo "  [-R|--resolution  ]   RESOLUTION                  : Resolution of sample in basepairs, using which the contact values should be calculated. Default 5000. Accepted resolutions- 1000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000"
    echo "  [-f|--formula     ]   FORMULA                     : Arithmetic to use to report contact calues. Options- center, max, average, sum. Default = center"
    echo "  [-F|--fix         ]   FIXED_COORDINATES           : If FALSE, reports new coordinates based on arithmetic center or max. Default = TRUE"
    echo "  [   --shrink_wrap ]   SHRINK_WRAP                 : Squeezes a 2D bedpe interval until supplied value(in raw read count units) is reached. Default = FALSE"
    echo "  [   --split       ]   SPLIT                       : Splits a 2D bedpe interval into multiple sub-intervals greater than supplied value(in raw read count units). Default = FALSE"
    echo "  [   --padding     ]   PADDING                     : Joins sub-intervals in 2D space reported by --split, based on supplied value(in bin units). Default = 2"
    echo "  [   --expand      ]   EXPAND                      : Expands 1D bedpe feet in both directions based on supplied value(in bin units). Default = 0"
    echo "  [-h|--help        ]   Help message"
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
      "--norm")         set -- "$@" "-Q" ;;
      "--sample2")      set -- "$@" "-B" ;;
      "--resolution")   set -- "$@" "-R" ;;
      "--formula")      set -- "$@" "-f" ;;
      "--fix")          set -- "$@" "-F" ;;
      "--split")        set -- "$@" "-S" ;; 
      "--shrink_wrap")  set -- "$@" "-s" ;;
      "--padding")      set -- "$@" "-p" ;;
      "--expand")       set -- "$@" "-e" ;;  
      "--help")         set -- "$@" "-h" ;;
       *)               set -- "$@" "$arg"
  esac
done


Q=none
R=5000
f=center
F=TRUE
S=FALSE
s=FALSE
p=2
e=0

while getopts ":P:A:G:Q:B:R:f:F:S:s:p:e:h" OPT
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


# Do all necessary parameter checks
#----------------------------------

if [[ -z $P ]];
then
    usage
    exit
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

if [[ -z $A ]]
then
    usage
    exit
fi

#----------------------------------

if [[ -z $G ]]
then
    usage
    exit
fi


num_loops=`cat "$P" | wc -l`



# RAM-disk update
if [ ! -d "$ram_disk" ]; then
 ## Proceeding with RAM_disk annotate_loops
 :
else
  echo "RAM disk currently in use on this machine. Please wait for the other process to stop."
  exit
fi

sudo mkdir -p "$ram_disk"

buffer=100000000

###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if [ -z "$B" ]
then
    
    size_file=`du -b $data_dir/$G/$A/$A.allValidPairs.hic | cut -f 1`
    size_disk=`echo "$size_file+$buffer" | bc`

    sudo mount -t tmpfs -o size="$size_disk" tmpfs "$ram_disk"

    ## echo "RAM disk successfully created and mounted"

    sudo cp $data_dir/$G/$A/$A.allValidPairs.hic /media/ram_disk

    ## echo ".hic copied to RAM disk"

    Rscript \
    $aqua_dir/annotate_loops.r \
      $P \
      $R \
      $ram_disk/$A.allValidPairs.hic \
      $data_dir/$G/$A/mergeStats.txt \
      $Q \
      $G \
      $num_loops \
      $f $F $S $s $p $e
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

    size_file1=`du -b $data_dir/$G/$A/$A.allValidPairs.hic | cut -f 1`
    size_file2=`du -b $data_dir/$G/$B/$B.allValidPairs.hic | cut -f 1`
    size_disk=`echo "$size_file1+$size_file2+$buffer" | bc`

    sudo mount -t tmpfs -o size="$size_disk" tmpfs "$ram_disk"

    ## echo "RAM disk successfully created and mounted"

    sudo cp $data_dir/$G/$A/$A.allValidPairs.hic /media/ram_disk
    sudo cp $data_dir/$G/$B/$B.allValidPairs.hic /media/ram_disk

    ## echo ".hics copied to RAM disk"

    Rscript \
    $aqua_dir/annotate_loops.r \
      $P \
      $R \
      $ram_disk/$A.allValidPairs.hic \
      $data_dir/$G/$A/mergeStats.txt \
      $ram_disk/$B.allValidPairs.hic \
      $data_dir/$G/$B/mergeStats.txt \
      $Q \
      $G \
      $num_loops \
      $f $F $S $s $p $e
fi


if [[ $ctrlc_count == 0 ]]; then
    sudo umount "$ram_disk"
    cd /media
    sudo rm -r ram_disk
fi

## echo "RAM disk unmounted and removed"
