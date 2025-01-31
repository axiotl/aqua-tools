#!/bin/bash

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools



temp_dir=$(mktemp -d /tmp/call_clusters-XXXXXXXXXX)



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
    echo "usage : call_clusters.sh \\"
    echo "  -A SAMPLE_NAME \\"
    echo "  -G GENOME \\"
    echo " [-h]"
    echo "Use option -h|--help for more information"
}


function help {
    echo
    echo "Call large-scale clusters of interacting regions from HiChIP using a sliding window"
    echo
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -A|--sample1        : Name of the sample as it appears on the Tinkerbox"
    echo "   -G|--genome         : Genome build used for sample processing"
    echo "  [-O|--output_dir    ]: Optional: provide a path for the final clusters file"
    echo "  [-S|--score         ]: Inherent score to seed cluster formation. Default 0.7"
    echo "  [-r|--resolution    ]: Resolution in base pairs. Only 5000 and 1000 supported. Default 5000"
    echo "  [   --radius        ]: Bin distance units to search for neighbours. Default 1"
    echo "  [   --min_dist      ]: Distance in basepairs to filter out extracted elements. Default 0"
    echo "  [   --window_size   ]: Size of sliding window in base pairs. Default 3000000"
    echo "  [   --overlap_size  ]: Size of overlap between two consecutive windows in base pairs. Default 1200000"
    echo "  [-h|--help          ]  Help message"
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
      "--sample1")         set -- "$@" "-A" ;;
      "--genome")          set -- "$@" "-G" ;;
      "--output_dir")      set -- "$@" "-O" ;;
      "--score")           set -- "$@" "-S" ;;
      "--resolution")      set -- "$@" "-r" ;;
      "--radius")          set -- "$@" "-d" ;;
      "--min_dist")        set -- "$@" "-x" ;;
      "--window_size")     set -- "$@" "-w" ;;
      "--overlap_size")    set -- "$@" "-o" ;;
      "--help")            set -- "$@" "-h" ;;
       *)                  set -- "$@" "$arg"
  esac
done


S=0.7
r=5000
d=1
x=0
O="blank"
w=3000000
o=1200000
while getopts ":A:G:O:S:r:d:x:w:o:h" OPT
do
    case $OPT in
  A) A=$OPTARG;;
  G) G=$OPTARG;;
  O) O=$OPTARG;;
  S) S=$OPTARG;;
  r) r=$OPTARG;;
  d) d=$OPTARG;;
  x) x=$OPTARG;;
  w) w=$OPTARG;;
  o) o=$OPTARG;;
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


if [[ -z $G ]];
then
    usage
    exit
fi


if [ $O == "blank" ];
then
    O=$(pwd)
else
    if [ ! -e $O ]; then
      mkdir -p $O
    fi
fi





# 1. Create sliding windows

echo
echo "Step1: Creating sliding windows...."


ref_dir="$data_dir"/"$G"/reference

chrom_sizes_file="$ref_dir"/"$G"_filtered.sizes.txt
centromeres_file="$ref_dir"/centromeres_"$G".bed


awk -v win_size="$w" -v stride="$o" '{
    chrom=$1;
    chrom_size=$2;
    start=0;
    end=start+win_size;
    while (end <= chrom_size) {
        print chrom"\t"start"\t"end;
        start=start+stride;
        end=start+win_size;
    }
}' $chrom_sizes_file > $temp_dir/sliding_windows.bed

if [ "$G" == "hg38" ] || [ "$G" == "hg19" ]; then

    # 2. Remove centromeres

    echo
    echo "Step2: Removing centromeric regions...."

    bedtools intersect \
     -a $temp_dir/sliding_windows.bed \
     -b $centromeres_file -v > $temp_dir/sliding_windows_without-centromeres.bed


    echo
    echo "Step3: Extracting bedpe from $A using sliding windows.... this takes a while"

    extract_bedpe.sh \
     --sample1    $A \
     --genome     $G \
     --TAD        $temp_dir/sliding_windows_without-centromeres.bed \
     --score      $S \
     --radius     $d \
     --min_dist   $x \
     --resolution $r > $temp_dir/"$A"_extract_tmp.bedpe


    if [ -s $temp_dir/"$A"_extract_tmp.bedpe ]
    then

        cut -f 1-6 $temp_dir/"$A"_extract_tmp.bedpe | sort  |  uniq > $temp_dir/"$A"_extract_whole-genome.bedpe

        echo
        echo "Step4: Clustering extracted bedpe..."

        cluster_bedpe.sh \
            -P $temp_dir/"$A"_extract_whole-genome.bedpe > $O/"$A"_"clusters_score-"$S"_radius-"$d".bedpe"
    else
        echo "No extracted bedpe found, exiting.."
        exit 1
    fi

    echo
    echo "Done! Clusters are in $O/"$A"_clusters_score-"$S"_radius-"$d".bedpe"

elif [ "$G" == "mm10" ]; then
    
    # 2. Remove centromeres

    echo
    echo "Step2: Extracting bedpe from $A using sliding windows.... this takes a while"

    extract_bedpe.sh \
     --sample1    $A \
     --genome     $G \
     --TAD        $temp_dir/sliding_windows.bed \
     --score      $S \
     --radius     $d \
     --min_dist   $x \
     --resolution $r > $temp_dir/"$A"_extract_tmp.bedpe


    if [ -s $temp_dir/"$A"_extract_tmp.bedpe ]
    then

        cut -f 1-6 $temp_dir/"$A"_extract_tmp.bedpe | sort  |  uniq > $temp_dir/"$A"_extract_whole-genome.bedpe

        echo
        echo "Step3: Clustering extracted bedpe..."

        cluster_bedpe.sh \
            -P $temp_dir/"$A"_extract_whole-genome.bedpe > $O/"$A"_"clusters_score-"$S"_radius-"$d".bedpe"
    else
        echo "No extracted bedpe found, exiting.."
        exit 1
    fi


    echo
    echo "Done! Clusters are in $O/"$A"_clusters_score-"$S"_radius-"$d".bedpe"

fi



trap "sudo rm -rf $temp_dir" EXIT
