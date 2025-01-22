#!/bin/bash

OPTIND=1
norm_method=NONE
unit=BP
aqua_genome=mm10

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

if [[ -z "${SAMPLE_SHEET}" ]]; then
    sample_sheet="$HOME/setup/sample_sheet.txt"    
else
    sample_sheet="${SAMPLE_SHEET}"
fi


juicer_tools='java -jar $HOME/juicer_tools_1.19.02.jar'

# Check if the current directory is the shared-files folder
current_directory=$(pwd)
shared_files_path="/home/ubuntu/shared-files"
if [ "$current_directory" == "$shared_files_path" ]; then
  echo -e "\nError: You cannot generate plots from the ~/shared-files folder. Please navigate to a different directory and generate the plot there.\n"
  exit 1
fi

function usage {
    echo -e "usage : plot_contacts.sh -A NAME_OF_FIRST_SAMPLE -R RANGE -G GENOME_BUILD [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    echo 
    echo "Create contact plots with CPM/AQuA normalized contact values"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -A|--sample1                         : Name of sample you want to use to create the contact plot, name it as it appears on the Tinkerbox"
    echo "    -R|--range                           : The genomic range that is to be plotted, in chr:start:end format. Example: -R chr1:40280000:40530000. For interchromosomal pairs use chr1:start1:end1-chr2:start2:end2 format. Example: chr2:219000000:225000000-chr13:37500000:43500000"    
    echo "    -G|--genome                          : The genome build the sample has been processed using. Strictly hg19 or hg38"
    echo " [  -O|--output_name                 ]   : Optional: provide a name for the plot"
    echo " [  -B|--sample2                     ]   : For two sample delta plots, name of the second sample"
    echo " [  -Q|--norm                        ]   : Which normalization to use. Strictly 'none', 'cpm' or 'aqua' in lower case. Non-spike-in samples default to cpm. Spike-in samples default to aqua."
    echo " [  -r|--resolution                  ]   : Resolution of sample in basepairs. Default 5000. Accepted resolutions- 1000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000"
    echo " [  -p|--profiles                    ]   : If contact profiles should be drawn along the diagonal, x axis and y axis. Default = FALSE"
    echo " [  -o|--color_one_sample            ]   : Color for contacts for single sample plots in RGB hexadecimal, ex: red = FF0000 (RRGGBB). Default = FF0000"
    echo " [  -t|--color_two_sample            ]   : Color for contacts for two sample plots (delta) in RGB hexadecimal separated by '-', ex: 1E90FF-C71585"
    echo " [     --annotations_default         ]   : Draw bed annotations; TSSs, ENCODE 3 enhancers, CpG islands. Default = TRUE"
    echo " [     --annotations_custom          ]   : Path to bed file to draw custom annotations. For two bed files provide space separated paths. ex: --annotations_custom bed_A.bed bed_B.bed"
    echo " [     --annotations_custom_color    ]   : Color for supplied bed in RGB hexadecimal, ex: C71585 for one bed file or C71585 1E90FF for two bed files "
    echo " [     --quant_cut                   ]   : Between 0.00-1.00. Rather than using the max value of the matrix as the highest color, cap the values at a given percentile. Default 1.00"
    echo " [     --max_cap                     ]   : Adjusts the cap value; either increases color intensity by reducing higher values to cap or decreases color intensity if set above maximum"
    echo " [     --get_matrix                  ]   : TRUE or FALSE. Obtain raw contact matrices instead of contact plot. Default FALSE"
    echo " [     --bedpe                       ]   : Supply path to a bedpe file to highlight tiles of interacting bedpe feet"
    echo " [     --bedpe_color                 ]   : Color for supplied bedpe in RGB hexadecimal. ex: C71585"
    echo " [  -i|--inherent                    ]   : TRUE or FALSE. If TRUE, normalize the contact plot using inherent normalization"
    echo " [     --inh_col_floor               ]   : Contact color for inherent values < 0, in RGB hexadecimal. Default = FFFFFF"
    echo " [     --inh_col_off                 ]   : Contact color for inherent values ~ 0, in RGB hexadecimal. Default = D4E4FB"
    echo " [     --inh_col_on                  ]   : Contact color for inherent values ~ 1, in RGB hexadecimal. Default = FF0000"
    echo " [     --inh_col_ceil                ]   : Contact color for inherent values > 1, in RGB hexadecimal. Default = FF8D4A"
    echo " [  -w|--width                       ]   : Manually set width of printed bin between 0-1. Default width calculated automatically."
    echo " [  -g|--gene                        ]   : Provide a gene name instead of an interval range (-R). Example: --gene sox8. For interchromosomal pairs use gene1,gene2. Example: pax3,foxo1"
    echo " [  -f|--flank                       ]   : Change interval range by flank value in bp. Example: --flank 5000"
    echo " [  -h|--help                        ]   : Help message. Primer can be found at https://rb.gy/fjkwkr"
    exit;
}


if [ $# -lt 1 ]
    then
    usage
    exit
fi


# Initialize arrays for annotations
ANNOTATION_FILES=()
ANNOTATION_COLORS=()

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--sample1")                      set -- "$@" "-A" ;;
      "--range")                        set -- "$@" "-R" ;;
      "--genome")                       set -- "$@" "-G" ;;
      "--output_name")                  set -- "$@" "-O" ;;
      "--sample2")                      set -- "$@" "-B" ;;
      "--norm")                         set -- "$@" "-Q" ;;
      "--resolution")                   set -- "$@" "-r" ;;
      "--profiles")                     set -- "$@" "-p" ;;
      "--color_one_sample")             set -- "$@" "-o" ;;
      "--color_two_sample")             set -- "$@" "-t" ;;
      "--annotations_default")          set -- "$@" "-d" ;;
      "--annotations_custom")           set -- "$@" "-x" ;;
      "--annotations_custom_color")     set -- "$@" "-c" ;;
      "--quant_cut")                    set -- "$@" "-q" ;;
      "--max_cap")                      set -- "$@" "-m" ;;
      "--get_matrix")                   set -- "$@" "-D" ;;
      "--bedpe")                        set -- "$@" "-b" ;;
      "--bedpe_color")                  set -- "$@" "-y" ;;
      "--inherent")                     set -- "$@" "-i" ;;
      "--inh_col_floor")                set -- "$@" "-K" ;;
      "--inh_col_off")                  set -- "$@" "-L" ;;
      "--inh_col_on")                   set -- "$@" "-M" ;;
      "--inh_col_ceil")                 set -- "$@" "-N" ;;
      "--width")                        set -- "$@" "-w" ;;
      "--gene")                         set -- "$@" "-g" ;;
      "--flank")                        set -- "$@" "-f" ;;
      "--help")                         set -- "$@" "-h" ;;
      *)                                set -- "$@" "$arg"
  esac
done


Q="blank"
r=5000
p=FALSE
d=TRUE
x=NONE
c=NONE
q=1
m=none
D=FALSE
b=FALSE
y=000000
i="FALSE"
w="blank"
K=ffffff
L=d4e4fb
M=ff0000
N=ff8d4a
T=FALSE
f=0 

prefix="_"

# Process all arguments
while getopts ":A:R:G:O:B:Q:r:p:o:t:d:x:c:q:m:D:b:y:i:K:L:M:N:w:g:T:f:h" OPT
do
  case $OPT in
    A) A=$OPTARG;;
    R) R=$OPTARG;;
    G) G=$OPTARG;;
    O) O=$OPTARG;;
    B) B=$OPTARG;;
    Q) Q=$OPTARG;;
    r) r=$OPTARG;;
    p) p=$OPTARG;;
    o) o=$OPTARG;;
    t) t=$OPTARG;;
    d) d=$OPTARG;;
    x) 
        # Accept two bed annotation files
        ANNOTATION_FILES+=("$OPTARG")
        while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
            ANNOTATION_FILES+=("${!OPTIND}")
            ((OPTIND++))
        done
        ;;
    c)
        # Accept two bed annotation colors
        ANNOTATION_COLORS+=("$OPTARG")
        while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
            ANNOTATION_COLORS+=("${!OPTIND}")
            ((OPTIND++))
        done
        ;;
    q) q=$OPTARG;;
    m) m=$OPTARG;;
    D) D=$OPTARG;;
    b) b=$OPTARG;;
    y) y=$OPTARG;;
    i) i=$OPTARG;;
    K) K=$OPTARG;;
    L) L=$OPTARG;;
    M) M=$OPTARG;;
    N) N=$OPTARG;;
    w) w=$OPTARG;;
    g) g=$OPTARG;;
    T) T=$OPTARG;;
    f) f=$OPTARG;;
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

# Set custom annotation files and colors or return to default
if [ ${#ANNOTATION_FILES[@]} -eq 0 ]; then
    x="NONE"
else
    x="${ANNOTATION_FILES[*]}"
fi

if [ ${#ANNOTATION_COLORS[@]} -eq 0 ]; then
    c="NONE"
else
    c="${ANNOTATION_COLORS[*]}"
fi

#----------------------------------

# If no sample B, the variable is empty:
if [ -z "$B" ]; then
  case $B in
  B) B="";;
  esac
fi

#----------------------------------

if [[ -z $A ]];
then
  usage
  exit
fi

#---------------------------------

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

if [[ -z $G ]]; then
  usage
  exit
elif [[ $G != "hg38" && $G != "hg19" && $G != "mm10" ]]; then
  echo -e "\nInvalid value for G. It must be 'hg38', 'hg19', or 'mm10'.\n"
  usage
  exit 1
fi

#----------------------------------

# Check if both -R and -g are provided
if [[ -n $R && -n $g ]]; then
    echo -e "\nOptions -R and --gene cannot be used together. Please provide only one of them.\n"
    usage  
    exit 1  
fi

if [[ -z $R && -z $g ]]; then
    usage  
    exit 1
fi


# Inform -R with --gene
if [[ ! -z $g ]]; then

    # Set files based on the value of $G
    bed_tss=""
    if [[ $G == "hg19" ]]; then
        bed_tss="$data_dir/$G/reference/TSS.3000.bed"
    elif [[ $G == "hg38" ]]; then
        bed_tss="$data_dir/$G/reference/GENCODE_TSSs_hg38.bed"
    fi

    # Check if two genes are provided for interchromosomal range
    if [[ $g == *","* ]]; then
        # Logic for two genes (interchromosomal)
        IFS=',' read -ra GENES <<< "$g"
        chr_array=()
        coordinates_array=()

        for gene in "${GENES[@]}"; do
            gene_uppercase=$(echo "$gene" | tr '[:lower:]' '[:upper:]' | xargs)
            gene_coordinates=$(grep -P "\t$gene_uppercase\t" "$bed_tss")

            if [[ -z $gene_coordinates ]]; then
                echo "Gene '$gene' not found in '$bed_tss'." >&2
                exit 1
            fi

            chr=$(echo "$gene_coordinates" | awk '{print $1}')
            expanded_coordinates=$(echo "$gene_coordinates" | awk -v OFS="\t" '{print $1, $2-500000, $3+500000}' | tr '[:blank:]' ':')
            
            chr_array+=("$chr")
            coordinates_array+=("$expanded_coordinates")
        done

        # Check if genes are on different chromosomes
        if [[ "${chr_array[0]}" == "${chr_array[1]}" ]]; then
            echo "Error: Both genes are on the same chromosome. Please provide genes on different chromosomes for interchromosomal analysis." >&2
            exit 1
        fi

        R="${coordinates_array[0]}-${coordinates_array[1]}"
        
    else
        # Logic for one gene (intrachromosomal)
        gene_uppercase=$(echo "$g" | tr '[:lower:]' '[:upper:]')
        gene_coordinates=$(grep -P "\t$gene_uppercase\t" "$bed_tss")
        if [[ -z $gene_coordinates ]]; then
            echo "Gene '$g' not found in '$bed_tss'." >&2
            exit 1
        fi
        expanded_coordinates=$(echo "$gene_coordinates" | awk -v OFS="\t" '{print $1, $2-500000, $3+500000}' | tr '[:blank:]' ':')
        R="$expanded_coordinates"
        
    fi
fi



# Function to extract chromosome number for ordering later 
extract_chr_number() {
    echo "$1" | sed -e 's/chr\([0-9XY]*\).*/\1/'
}

# The correct_format regex matches strings in the format: 
# chr<number>:<number>:<number>-chr<number>:<number>:<number>
# likely_attempt is a looser pattern to catch common formatting errors in the input.

# Regular expressions for correctly formatted intra- and interchromosomal input
intra_format="^chr[0-9XY]+:[0-9]+:[0-9]+(:[a-zA-Z0-9_]+)?$"
inter_format="^chr[0-9XY]+:[0-9]+:[0-9]+-chr[0-9XY]+:[0-9]+:[0-9]+$"

# Enhanced regular expression for a likely attempt at interchromosomal input
likely_attempt="chr[0-9XY]*[: -]*[0-9]*[: -]*[0-9]*[: -]*chr[0-9XY]*"

# Check if input is intrachromosomal
if [[ $R =~ $intra_format ]]; then
  
    # Check and notify about the tag, if present
    if [[ $R == *:*:*:* ]]; then
      IFS=":" read -ra RANGE_ARRAY <<< "$R"
      tag="${RANGE_ARRAY[3]}"
      echo -e "\nNotice: The tag ($tag) is no longer required for -R ($R)\n"
      R="${RANGE_ARRAY[0]}:${RANGE_ARRAY[1]}:${RANGE_ARRAY[2]}"
fi

# Check if input is interchromosomal
elif [[ $R =~ $inter_format ]]; then
    T=TRUE
    echo -e "\nInterchromosomal input detected, continuing...\n"

    # Ensure that inherent normalization is FALSE for interchromosomal ranges
    if [[ $i == TRUE ]]; then
        echo -e "\nInherent normalization is not compatible with interchromosomal ranges.\n"
        exit 1
    fi
    
    # Split the input into two parts for interchromosomal range
    IFS='-' read -ra ADDR <<< "$R"
    chr1=${ADDR[0]}
    chr2=${ADDR[1]}

    # Extract chromosome numbers
    num1=$(extract_chr_number "$chr1")
    num2=$(extract_chr_number "$chr2")

    # Compare and rearrange if necessary (for strawr later)
    if [[ $num1 -gt $num2 ]]; then
        R="$chr2-$chr1"
        echo "Range rearranged to: $R"
    fi

# Check if input is a likely attempt at interchromosomal format
elif [[ $R =~ $likely_attempt ]]; then
    echo -e "\nThe range format for interchromosomal input seems incorrect. Please ensure it follows the chr1:start1:end1-chr2:start2:end2 format with no extra spaces or misplaced characters.\n"
    exit 1

# Catch-all for invalid formats
else
    echo -e "\nInvalid range. Exiting...\n"
    exit 1
fi

bin_size=$r


adjust_coordinates() {
    local coordinate=$1
    local bin_size=$2
    echo $(( coordinate / bin_size * bin_size ))
}

flank=$f

# Check if flank is a multiple of resolution
if (( flank % bin_size != 0 )); then
    echo -e "\nFlank value must be a multiple of the resolution, continuing with nearest multiple of resolution...\n"
    # Adjust flank to nearest multiple of resolution
    flank=$(( (flank / bin_size + 1) * bin_size ))
fi


# Adjust range for bin size and add flank
IFS='-' read -ra ADDR <<< "$R"
for addr in "${ADDR[@]}"; do
    IFS=':' read -ra COORDS <<< "$addr"
    start=${COORDS[1]}
    end=${COORDS[2]}

    # Adjust start and end coordinates to align with bin size
    adjusted_start=$(adjust_coordinates $start $bin_size)
    adjusted_end=$(adjust_coordinates $end $bin_size)

    # Add flank to the adjusted range
    # Ensure flank does not make the start negative
    adjusted_start_with_flank=$(( adjusted_start - flank ))
    if [ $adjusted_start_with_flank -lt 0 ]; then
        adjusted_start_with_flank=0
    fi

    # Add flank to the end
    adjusted_end_with_flank=$(( adjusted_end + flank ))
    adjusted_ranges+=("${COORDS[0]}:$adjusted_start_with_flank:$adjusted_end_with_flank")
done
R=$(IFS='-'; echo "${adjusted_ranges[*]}")



# Generate output name if not provided
if [[ -z $O ]]; then
  IFS=":" read -ra RANGE_ARRAY <<< "$R"
  interval_chr="${RANGE_ARRAY[0]}"
  interval_start="${RANGE_ARRAY[1]}"
  interval_end="${RANGE_ARRAY[2]}"
  O="${interval_chr}_${interval_start}_${interval_end}.pdf"
else
  # Add .pdf if it's missing
  if [[ $O != *.pdf ]]; then
    O="${O}.pdf"
  fi
fi



# Check bedpe
if [ "$b" = "FALSE" ]; then
    echo -e "\nNo BEDPE file supplied to highlight contact plot\n"
else
    # Check if file exists
    if [ ! -f "$b" ]; then
        echo "BEDPE file $b does not exist"
        echo "Continuing without --bedpe ..."
        echo
        b="FALSE"
    else
    # Check number of columns and coordinate ordering
        awk_output=$(awk '
        BEGIN {FS="\t"; errors=0; reversed=0}
        NF < 6 {
            print "BEDPE error: Line " NR " has fewer than 6 columns"
            errors++
        }
        $3 < $2 || $6 < $5 {
            print "BEDPE error: Reversed feet in line " NR
            reversed++
        }
        END {
            if (errors > 0) {
                print "BEDPE error: file should contain at least 6 columns in all rows."
                exit 1
            }
            if (reversed > 0) {
                print "Please ensure that column 3 > column 2 and column 6 > column 5 for all entries."
                exit 1
            }
            exit 0
        }' "$b")
        
        # Check the exit status of awk
        if [ $? -ne 0 ]; then
            echo "$awk_output"
            echo "Continuing without --bedpe $b ..."
            echo
            b="FALSE"
        fi
    fi
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

  analysis_type="single_sample"

  echo "Commencing $analysis_type analysis"

  #----------------------------------

  if [[ -z $o ]]; then
    # Set default color based on the value of 'inter'
    if [ "$T" = "TRUE" ]; then
        o="004554"
    else
        o="FF0000"
    fi
fi

  if [[ -n $t ]]; 
  then 
    echo "Oops! Looks like you provided the wrong parameter for contact color!"
    usage
    exit
  fi

  #----------------------------------

  echo "resolution: $r"

  if [[ $T == "FALSE" ]]; then
     range_text=`echo "$R" | sed 's/:/_/g'`

     chr=`echo "$R" | cut -d":" -f1`
     start=`echo "$R" | cut -d":" -f2`
     end=`echo "$R" | cut -d":" -f3`

  fi

  if [[ $T == "TRUE" ]]; then
    range_text=`echo "$R" | sed 's/[:-]/_/g'`

    # interchromosomal ranges calculated in R
    chr=NA
    start=NA
    end=NA
    tag=NA
fi

  pair1=$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic

  if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi

  stat1=$data_dir/$G/$A/mergeStats.txt
  if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi

  
  sample_dir=$data_dir/$G/$A
  out_dir=`pwd`

  echo "switch to R"
  Rscript $aqua_dir/plot_contacts.r \
    $analysis_type \
    $T \
    $pair1 \
    $stat1 \
    $norm_method \
    $unit \
    $bin_size \
    $prefix \
    "${range_text}_${bin_size}" \
    $out_dir \
    $chr $start $end \
    $G \
    $p \
    $o \
    $d \
    "$x" \
    "$c" \
    $m \
    $Q \
    $q \
    $O \
    $b \
    $y \
    $i \
    $sample_dir \
    $w \
    $K \
    $L \
    $M \
    $N \
    $D

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

  if [[ $i == "TRUE" ]];
  then
    echo -e "\n--inherent only applies to one-sample plots for now, exiting...\n"
    exit
  fi

  analysis_type="two_sample"

  echo "Commencing $analysis_type analysis"

  #----------------------------------
  

  if [[ -z $t ]]; then
    if [ "$T" = "TRUE" ]; then
        t="837900-a63200"
    else
        t="1E90FF-C71585"
    fi
  fi


  if [[ -n $o ]]; 
  then 
    echo "\nOops! Looks like you provided the wrong parameter for contact color\n!"
    usage
    exit
  fi

  #----------------------------------

  echo "resolution: $r"

  if [[ $T == "FALSE" ]]; then
     range_text=`echo "$R" | sed 's/:/_/g'`

     chr=`echo "$R" | cut -d":" -f1`
     start=`echo "$R" | cut -d":" -f2`
     end=`echo "$R" | cut -d":" -f3`

  fi

  if [[ $T == "TRUE" ]]; then
    range_text=`echo "$R" | sed 's/[:-]/_/g'`

    chr=NA
    start=NA
    end=NA
    tag=NA
  fi

  pair1=$data_dir/$G/$A/${version_dir_A}.allValidPairs.hic
  pair2=$data_dir/$G/$B/${version_dir_B}.allValidPairs.hic
  if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
  if [[ ! -f $pair2 ]]; then echo "cannot find $pair2"; exit 1; fi

  stat1=$data_dir/$G/$A/mergeStats.txt
  stat2=$data_dir/$G/$B/mergeStats.txt
  if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi
  if [[ ! -f $stat2 ]]; then echo "cannot find $stat2"; exit 1; fi


  out_dir=`pwd`
  
  sample_dirA=$data_dir/$G/$A
  sample_dirB=$data_dir/$G/$B

  echo "switch to R"
  Rscript $aqua_dir/plot_contacts.r \
    $analysis_type \
    $T \
    $pair1 \
    $pair2 \
    $stat1 \
    $stat2 \
    $norm_method \
    $unit \
    $bin_size \
    $prefix \
    "${range_text}_${bin_size}" \
    $out_dir \
    $chr $start $end \
    $G \
    $p \
    $t \
    $d \
    "$x" \
    "$c" \
    $m \
    $Q \
    $q \
    $O \
    $b \
    $y \
    $w \
    $sample_dirA $sample_dirB \
    $D \
    $K \
    $L \
    $M \
    $N \
    $i

fi

exit_status=$?

