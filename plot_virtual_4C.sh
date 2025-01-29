#!/bin/bash


OPTIND=1
norm_method=NONE
unit=BP

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools


sample_sheet="~/setup/sample_sheet.txt"


juicer_tools='java -jar ~/juicer_tools_1.19.02.jar'


function usage {
    echo -e "usage : plot_virtual_4C.sh -A NAME_OF_FIRST_SAMPLE -G GENOME_BUILD -R RANGE -V VIEWPOINT [-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Plot Virtual 4C"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -A|--sample1            : Name of the first sample you want to use as it appears on the Tinkerbox"
    echo "    -G|--genome             : The genome build the sample(s) has been processed using. Strictly hg19, hg38, or mm10"
    echo "    -R|--range              : The genomic range that is to be plotted, in chr:start:end format. For example: -R chr1:40280000:40530000"
    echo "    -V|--viewpoint          : The viewpoint that is to be considered, in chr:start format. For example: -V chr1:40400000"
    echo "  [ -B|--sample2        ]   : For two sample delta plots, name of the second sample"
    echo "  [ -Q|--norm           ]   : Which normalization to use. Strictly 'none', 'cpm' or 'aqua' in lower case. Non-spike-in samples default to cpm. Spike-in samples default to aqua."
    echo "  [ -r|--resolution     ]   : Resolution of sample in basepairs. Default 5000. Accepted resolutions- 1000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000"
    echo "  [ -O|--output_name    ]   : Optional: provide a name for the plot"
    echo "  [    --quant_cut      ]   : Cap matrix values at a given percentile (0.00-1.00). Default = 1.00"
    echo "  [    --max_cap        ]   : Cap matrix values by adjusting the maximum value, modifying color intensity"
    echo "  [    --width          ]   : Number of bins up and downstream of viewpoint locus to be considered for drawing profiles. Default 0"
    echo "  [    --height         ]   : Numeric factor to control the height of the Virtual 4C profile in the plot. Default 0.5"
    echo "  [ -i|--inherent       ]   : TRUE or FALSE. If TRUE, normalize the contact plot using inherent normalization"
    echo "  [    --inh_col_floor  ]   : Contact color for inherent values < 0, in RGB hexadecimal. Default = FFFFFF"
    echo "  [    --inh_col_off    ]   : Contact color for inherent values ~ 0, in RGB hexadecimal. Default = D4E4FB"
    echo "  [    --inh_col_on     ]   : Contact color for inherent values ~ 1, in RGB hexadecimal. Default = FF0000"
    echo "  [    --inh_col_ceil   ]   : Contact color for inherent values > 1, in RGB hexadecimal. Default = FF8D4A"
    echo "  [ -h|--help           ]     Help message"
    echo
    echo

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
      "--sample1")      set -- "$@" "-A" ;;
      "--genome")       set -- "$@" "-G" ;;
      "--range")        set -- "$@" "-R" ;;
      "--viewpoint")    set -- "$@" "-V" ;;
      "--sample2")      set -- "$@" "-B" ;;
      "--norm")         set -- "$@" "-Q" ;;
      "--resolution")   set -- "$@" "-r" ;;
      "--output_name")  set -- "$@" "-O" ;;
      "--quant_cut")    set -- "$@" "-q" ;;
      "--max_cap")      set -- "$@" "-m" ;;
      "--width")        set -- "$@" "-w" ;;
      "--height")       set -- "$@" "-H" ;;
      "--inherent")     set -- "$@" "-i" ;;
      "--inh_col_floor") set -- "$@" "-K" ;;
      "--inh_col_off")   set -- "$@" "-L" ;;
      "--inh_col_on")    set -- "$@" "-M" ;;
      "--inh_col_ceil")  set -- "$@" "-N" ;;
      "--help")          set -- "$@" "-h" ;;
       *)                set -- "$@" "$arg"
  esac
done


Q="blank"
r=5000
q=1
m="none"
w=0
H=0.5
i="FALSE"
K=ffffff
L=d4e4fb
M=ff0000
N=ff8d4a


while getopts ":A:G:R:V:B:Q:r:O:q:m:w:H:i:K:L:M:N:h" OPT
do
  case $OPT in
  A) A=$OPTARG;;
  G) G=$OPTARG;;
  R) R=$OPTARG;;
  V) V=$OPTARG;;
  B) B=$OPTARG;;
  Q) Q=$OPTARG;;
  r) r=$OPTARG;;
  O) O=$OPTARG;;
  q) q=$OPTARG;;
  m) m=$OPTARG;;
  w) w=$OPTARG;;
  H) H=$OPTARG;;
  i) i=$OPTARG;;
  K) K=$OPTARG;;
  L) L=$OPTARG;;
  M) M=$OPTARG;;
  N) N=$OPTARG;;
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

if [[ -z $R ]];
then
  usage
  exit
fi

#----------------------------------

if [[ -z $A ]];
then
  usage
  exit
fi

#-----------------------------------


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

if [[ -z $V ]];
then
  usage
  exit
fi

#----------------------------------

if [[ "$Q" != "blank" && "$i" != "FALSE" ]]; then
  echo -e "\nInherent normalization is not compatible with other --norm methods. Continuing with --inherent TRUE ..."
  Q="blank"
  i="TRUE"
fi


# Generate output name if not provided
if [[ -z $O ]]; then
  IFS=":" read -ra RANGE_ARRAY <<< "$R"
  interval_chr="${RANGE_ARRAY[0]}"
  interval_start="${RANGE_ARRAY[1]}"
  interval_end="${RANGE_ARRAY[2]}"
  O="${interval_chr}_${interval_start}_${interval_end}_4C.pdf"
else
  # Add .pdf if it's missing
  if [[ $O != *.pdf ]]; then
    O="${O}.pdf"
  fi
fi




range_text=`echo "$R" | sed 's/:/_/g'`

chr=`echo "$R" | cut -d":" -f1`
start=`echo "$R" | cut -d":" -f2`
end=`echo "$R" | cut -d":" -f3`

range_length=$((end - start))

if [[ $range_length -gt 3000000 && $r -le 5000 ]]; then
  echo -e "\nThe range length ($range_length bp) is too large for the specified resolution ($r). Try reducing the range or increasing the resolution."
  exit
fi


range=`printf "%s:%s:%s" $chr $start $end`

v_chr=`  echo "$V" | cut -d":" -f1`
v_start=`echo "$V" | cut -d":" -f2`




###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################

if [ -z "$B" ]
then

  analysis_type="single_sample"
  echo -e "\nCommencing $analysis_type analysis"

  
  pair1=$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic
  if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
  
  stat1=$data_dir/$G/$A/mergeStats.txt
  if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi

  
  
  hg_total1=`head -3 $stat1 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
  mm_total1=`head -3 $stat1 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`
  
  
  echo -e "\ntotal sample reads for $version_dir_A: $hg_total1"
  echo -e "total spike-in reads for $version_dir_A: $mm_total1\n"
  
  has_aqua=true
  total1=$hg_total1
  if [[ ! -z $mm_total1  ]];then
    echo "We have spike-In"
    total1=`echo "$hg_total1+$mm_total1" | bc`
    aqua_factor1=`echo "scale=4;$hg_total1/$mm_total1" | bc`
  fi
  norm_factor1=`echo "scale=4;1000000/$total1" | bc`

  sample_dir=$data_dir/$G/$A
  
  echo -e "switching to R...\n"

  Rscript $aqua_dir/plot_virtual_4C.r \
    $analysis_type \
    $pair1 \
    $stat1 \
    $norm_method \
    $unit \
    $r \
    $O \
    $G \
    $Q \
    $chr $start $end \
    $q \
    $m \
    $w \
    $H \
    $v_chr \
    $v_start \
    $sample_dir \
    $i \
    $K \
    $L \
    $M \
    $N

  exit_status=$?

  # Report successful plot generation location
  if [ $exit_status -eq 0 ] && [ -f "$O" ]; then
      echo -e "\nDone! Created $O\n"
  fi


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

  analysis_type="two_sample"
  echo -e "\nCommencing $analysis_type analysis\n"

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

  if [[ $i == "TRUE" ]];
  then
    echo -e "--inherent only applies to one-sample plots for now, exiting...\n"
    exit
  fi

  
  pair1=$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic
  pair2=$data_dir/$G/$B/${version_dir_B}.allValidPairs.hic
  if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
  if [[ ! -f $pair2 ]]; then echo "cannot find $pair2"; exit 1; fi

  stat1=$data_dir/$G/$A/mergeStats.txt
  stat2=$data_dir/$G/$B/mergeStats.txt
  if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi
  if [[ ! -f $stat2 ]]; then echo "cannot find $stat2"; exit 1; fi
  
  
  
  hg_total1=`head -3 $stat1 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
  hg_total2=`head -3 $stat2 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
  mm_total1=`head -3 $stat1 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`
  mm_total2=`head -3 $stat2 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`
  
  
  echo "total sample reads for $version_dir_A: $hg_total1"
  echo "total sample reads for $version_dir_B: $hg_total2"
  echo "total spike-in reads for $version_dir_A: $mm_total1"
  echo "total spike-in reads for $version_dir_B: $mm_total2"
  
  has_aqua=true
  total1=$hg_total1
  total2=$hg_total2
  if [[ ! -z $mm_total1 && ! -z $mm_total2 ]];then
    echo "We have spike-In"
    total1=`echo "$hg_total1+$mm_total1" | bc`
    total2=`echo "$hg_total2+$mm_total2" | bc`
    aqua_factor1=`echo "scale=4;$hg_total1/$mm_total1" | bc`
    aqua_factor2=`echo "scale=4;$hg_total2/$mm_total2" | bc`

  fi
  norm_factor1=`echo "scale=4;1000000/$total1" | bc`
  norm_factor2=`echo "scale=4;1000000/$total2" | bc`
  
  sample_dirA=$data_dir/$G/$A
  sample_dirB=$data_dir/$G/$B

  echo -e "\nswitching to R...\n"

  Rscript $aqua_dir/plot_virtual_4C.r \
    $analysis_type \
    $pair1 \
    $stat1 \
    $pair2 \
    $stat2 \
    $norm_method \
    $unit \
    $r \
    $O \
    $G \
    $Q \
    $chr $start $end \
    $tag \
    $q \
    $m \
    $w \
    $H \
    $v_chr \
    $v_start \
    $v_tag \
    $sample_dir_A $sample_dir_B \
    $i \
    $K \
    $L \
    $M \
    $N
  

  exit_status=$?

  # Report successful plot generation location
  if [ $exit_status -eq 0 ] && [ -f "$O" ]; then
      echo -e "\nDone! Created $O\n"
  fi
fi


