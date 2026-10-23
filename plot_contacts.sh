#!/bin/bash

OPTIND=1
norm_method=NONE
unit=BP
aqua_genome=mm10

data_dir="${LAB_DATA_DIR:-$HOME/lab-data}"
aqua_dir="${AQUA_TOOLS_DIR:-$HOME/aqua_tools}"

if [[ -z "${SAMPLE_SHEET}" ]]; then
    sample_sheet="$HOME/setup/sample_sheet.txt"
else
    sample_sheet="${SAMPLE_SHEET}"
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
    echo "    -A|--sample1                         : Name of the sample as it appears on the Tinker box"
    echo "    -R|--range                           : Genomic range to plot in chr:start:end or interchromosomal chr1:start1:end1-chr2:start2:end2 format"
    echo "    -G|--genome                          : The genome build the sample(s) has been processed using: hg19, hg38, or mm10"
    echo " [     --sample1_dir                 ]   : If not using the Tinker box, full path to the directory containing sample data (<dir>/<dir>.hic and <dir>/<dir>.mergeStats.txt)"
    echo " [     --sample2_dir                 ]   : If not using the Tinker box, full path to the second sample directory"
    echo " [  -O|--output_name                 ]   : Name for the plot or matrix output"
    echo " [  -B|--sample2                     ]   : For two sample delta plots, name of the second sample"
    echo " [  -Q|--norm                        ]   : Which normalization to use: none, cpm, or aqua in lower case"
    echo " [  -r|--resolution                  ]   : Resolution of sample in basepairs. Default 5000"
    echo " [  -o|--color_one_sample            ]   : Color for single sample plots in RGB hexadecimal. Default FF0000"
    echo " [  -t|--color_two_sample            ]   : Color for delta plots in RGB hexadecimal separated by '-'; 1E90FF-C71585"
    echo " [     --annotations_default         ]   : Annotation style options: standard, axiotl, or none. Default standard"
    echo " [     --annotations_custom          ]   : Path to bed file(s) for custom annotations; --annotations_custom 1.bed 2.bed"
    echo " [     --annotations_custom_color    ]   : Color for supplied bed in RGB hexadecimal; C71585 1E90FF for two bed files"
    echo " [     --annotations_custom_label    ]   : Draw names/strands from supplied BED6-format bed(s). TRUE or FALSE per file. Default FALSE"
    echo " [     --quant_cut                   ]   : Cap matrix values at a given percentile (0.00-1.00). Default 1.00"
    echo " [     --max_cap                     ]   : Cap matrix values by adjusting the maximum value, modifying color intensity"
    echo " [     --get_matrix                  ]   : Obtain raw contact matrices instead of contact plot. Default FALSE"
    echo " [  -P|--bedpe                       ]   : Path to bedpe file(s) to highlight tiles of interacting bedpe feet; --bedpe 1.bedpe 2.bedpe"
    echo " [     --bedpe_color                 ]   : Color for supplied bedpe(s) in RGB hexadecimal; --bedpe_color C71585 1E90FF"
    echo " [     --bedpe_lwd                   ]   : Line width for supplied bedpe(s); --bedpe_lwd 0.5 1.0. Default 0.8"
    echo " [  -i|--inherent                    ]   : If TRUE, normalize contacts using inherent normalization. Default FALSE"
    echo " [     --inh_col_floor               ]   : Contact color for inherent values < 0, in RGB hexadecimal. Default FFFFFF"
    echo " [     --inh_col_off                 ]   : Contact color for inherent values ~ 0, in RGB hexadecimal. Default D4E4FB"
    echo " [     --inh_col_on                  ]   : Contact color for inherent values ~ 1, in RGB hexadecimal. Default FF0000"
    echo " [     --inh_col_ceil                ]   : Contact color for inherent values > 1, in RGB hexadecimal. Default FF8D4A"
    echo " [  -w|--width                       ]   : Manually set width of printed bin between 0-1. Default calculated automatically"
    echo " [  -g|--gene                        ]   : Provide a gene name instead of an interval range; --gene sox8 or pax3,foxo1"
    echo " [  -f|--flank                       ]   : Change interval range by flank value in bp; --flank 5000"
    echo " [  -v|--svg                         ]   : If TRUE, write an SVG instead of a PDF. Default FALSE"
    echo " [  -s|--expand_viewport             ]   : Expand the plot view to include rectangular viewport. Default TRUE"
    echo " [  -h|--help                        ]   : Help message"
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
ANNOTATION_LABELS=()

# Initialize arrays for bedpe
BEDPE_FILES=()
BEDPE_COLORS=()
BEDPE_LWDS=()

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--sample1")                      set -- "$@" "-A" ;;
      "--sample1_dir")                  set -- "$@" "-H" ;;
      "--sample2_dir")                  set -- "$@" "-I" ;;
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
      "--annotations_custom_label")     set -- "$@" "-X" ;;
      "--quant_cut")                    set -- "$@" "-q" ;;
      "--max_cap")                      set -- "$@" "-m" ;;
      "--get_matrix")                   set -- "$@" "-D" ;;
      "--bedpe")                        set -- "$@" "-P" ;;
      "--bedpe_color")                  set -- "$@" "-y" ;;
      "--bedpe_lwd")                    set -- "$@" "-l" ;;
      "--inherent")                     set -- "$@" "-i" ;;
      "--inh_col_floor")                set -- "$@" "-K" ;;
      "--inh_col_off")                  set -- "$@" "-L" ;;
      "--inh_col_on")                   set -- "$@" "-M" ;;
      "--inh_col_ceil")                 set -- "$@" "-N" ;;
      "--width")                        set -- "$@" "-w" ;;
      "--gene")                         set -- "$@" "-g" ;;
      "--flank")                        set -- "$@" "-f" ;;
      "--svg")                          set -- "$@" "-v" ;;
      "--expand_viewport")              set -- "$@" "-s" ;;
      "--help")                         set -- "$@" "-h" ;;
      *)                                set -- "$@" "$arg"
  esac
done

Q="blank"
r=5000
p=FALSE
d="standard"
x=NONE
c=NONE
n=NONE
q=1
m=none
D=FALSE
P=FALSE
y=000000
l=0.8
i="FALSE"
w="blank"
K=ffffff
L=d4e4fb
M=ff0000
N=ff8d4a
T=FALSE
f=0
v=FALSE
s=TRUE
H=blank
I=blank

prefix="_"

# Process all arguments
while getopts ":A:H:I:R:G:O:B:Q:r:p:o:t:d:x:c:X:q:m:D:P:y:l:i:K:L:M:N:w:g:T:f:v:s:h" OPT
do
  case $OPT in
    A) A=$OPTARG;;
    H) H=$OPTARG;;
    I) I=$OPTARG;;
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
        # Accept multiple bed annotation files
        ANNOTATION_FILES+=("$OPTARG")
        while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
            ANNOTATION_FILES+=("${!OPTIND}")
            ((OPTIND++))
        done
        ;;
    c)
        # Accept multiple bed annotation colors
        ANNOTATION_COLORS+=("$OPTARG")
        while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
            ANNOTATION_COLORS+=("${!OPTIND}")
            ((OPTIND++))
        done
        ;;
    X)
        # Accept multiple bed annotation label flags
        ANNOTATION_LABELS+=("$OPTARG")
        while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
            ANNOTATION_LABELS+=("${!OPTIND}")
            ((OPTIND++))
        done
        ;;
    q) q=$OPTARG;;
    m) m=$OPTARG;;
    D) D=$OPTARG;;
    P)
        # Accept multiple bedpe files
        BEDPE_FILES+=("$OPTARG")
        while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
            BEDPE_FILES+=("${!OPTIND}")
            ((OPTIND++))
        done
        ;;
    y)
        # Accept multiple bedpe colors
        BEDPE_COLORS+=("$OPTARG")
        while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
            BEDPE_COLORS+=("${!OPTIND}")
            ((OPTIND++))
        done
        ;;
    l)
        # Accept multiple bedpe lwds
        BEDPE_LWDS+=("$OPTARG")
        while [[ $OPTIND -le $# && ! "${!OPTIND}" =~ ^- ]]; do
            BEDPE_LWDS+=("${!OPTIND}")
            ((OPTIND++))
        done
        ;;
    i) i=$OPTARG;;
    K) K=$OPTARG;;
    L) L=$OPTARG;;
    M) M=$OPTARG;;
    N) N=$OPTARG;;
    w) w=$OPTARG;;
    g) g=$OPTARG;;
    T) T=$OPTARG;;
    f) f=$OPTARG;;
    s) s=$OPTARG;;
    v) v=$OPTARG;;
    h) help ;;
  \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
  :)
    case "$OPTARG" in
      i) echo -e "\n--inherent requires a value of TRUE or FALSE (--inherent TRUE)\n" ;;
      D) echo -e "\n--get_matrix requires a value of TRUE or FALSE (--get_matrix TRUE)\n" ;;
      s) echo -e "\n--expand_viewport requires a value of TRUE or FALSE (--expand_viewport TRUE)\n" ;;
      v) echo -e "\n--svg requires a value of TRUE or FALSE (--svg TRUE)\n" ;;
      X) echo -e "\n--annotations_custom_label requires TRUE or FALSE per bed file (--annotations_custom_label TRUE FALSE)\n" ;;
      *) echo -e "\nOption -$OPTARG requires an argument.\n" ;;
    esac
    usage
    exit 1
    ;;
  esac
done

echo

#----------------------------------
# Sample inputs: a sample name (-A/-B) or a local directory (--sample1_dir/--sample2_dir)

if [[ -z $A && "$H" == "blank" ]];
then
  usage
  exit
fi

if [[ -n $A && "$H" != "blank" ]]; then
  echo -e "\nNeed either sample1 name or sample1 directory path, not both.\n"
  usage
  exit 1
fi

if [[ "$A" == *"/"* ]]; then
  echo -e "\nThe -A parameter expects a sample name, not a path. Did you mean to use --sample1_dir?\n"
  exit 1
fi

if [[ -n $B && "$I" != "blank" ]]; then
  echo -e "\nNeed either sample2 name or sample2 directory path, not both.\n"
  usage
  exit 1
fi

if [[ "$B" == *"/"* ]]; then
  echo -e "\nThe -B parameter expects a sample name, not a path. Did you mean to use --sample2_dir?\n"
  exit 1
fi

two_sample=FALSE
if [[ -n $B || "$I" != "blank" ]]; then
  two_sample=TRUE
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

# Resolve a sample name or a local sample directory into
# hic path, mergeStats path, sample dir and a display name
# Sets: res_hic res_mgs res_dir res_name
resolve_sample() {
    local name="$1"
    local local_dir="$2"
    local label="$3"

    if [[ "$local_dir" == "blank" ]]; then
        local sdir
        sdir=$(get_sample_directory "$name")
        if [ $? -ne 0 ]; then
            echo "Sample '$name' not found for genome $G. Exiting." >&2
            exit 1
        fi
        local vdir
        vdir=$(basename "$sdir")
        res_hic="$sdir/${vdir}.allValidPairs.hic"
        res_mgs="$sdir/mergeStats.txt"
        res_dir="$sdir"
        res_name=$(basename "$(dirname "$sdir")")
    else
        local_dir="${local_dir%/}"
        if [[ ! -d "$local_dir" ]]; then
            echo "Directory '$local_dir' does not exist. Exiting." >&2
            exit 1
        fi
        local sname
        sname=$(basename "$local_dir")
        res_hic="$local_dir/${sname}.hic"
        res_mgs="$local_dir/${sname}.mergeStats.txt"
        res_dir="$local_dir"
        res_name="$sname"

        if ! grep -q "^valid_interaction_rmdup" "$res_mgs" 2>/dev/null; then
            if [[ -f "$res_mgs" ]]; then
                echo "The mergeStats file ($res_mgs) does not contain the required row 'valid_interaction_rmdup'. Please refer to the GitHub repository for formatting instructions: https://github.com/axiotl/aqua-tools" >&2
                exit 1
            fi
        fi
    fi

    if [[ ! -f "$res_hic" ]]; then echo "cannot find $res_hic"; exit 1; fi
    if [[ ! -f "$res_mgs" ]]; then echo "cannot find $res_mgs"; exit 1; fi
}

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

# check for missing true/false flags
validate_bool_param() {
    local name="$1"
    local value="$2"
    if [[ -z "$value" || "$value" == -* ]]; then
        echo -e "\n--${name} requires a value of TRUE or FALSE (--${name} TRUE)\n"
        exit 1
    fi
    if [[ "$value" != "TRUE" && "$value" != "FALSE" ]]; then
        echo -e "\n--${name} must be TRUE or FALSE. Got: '$value'\n"
        exit 1
    fi
}

validate_bool_param "inherent" "$i"
validate_bool_param "get_matrix" "$D"
validate_bool_param "expand_viewport" "$s"
validate_bool_param "svg" "$v"

#----------------------------------

if [[ "$i" == "TRUE" && ( "$H" != "blank" || "$I" != "blank" ) ]]; then
  echo -e "\nInherent normalization is not compatible with this sample.\n"
  exit 1
fi

# Resolve sample A
resolve_sample "$A" "$H" "sample1"
pair1="$res_hic"
stat1="$res_mgs"
sample_dir_A="$res_dir"
base_dir_A="$res_name"

#----------------------------------

case "$r" in
  1000|5000|10000|25000|50000|100000|250000|500000|1000000|2500000)
    # Value is accepted
    ;;
  *)
    echo -e "\n--resolution '$r' is not accepted.\nAccepted resolutions are: 1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000\n"
    exit 1
    ;;
esac

#----------------------------------
# check for multiple normalization methods

if [[ "$Q" != "blank" && "$i" != "FALSE" ]]; then
  echo -e "\nInherent normalization is not compatible with other --norm methods. Continuing with --inherent TRUE ..."
  Q="blank"
  i="TRUE"
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
        bed_tss="$data_dir/$G/reference/GENCODE_TSSs_hg19.bed"
    elif [[ $G == "hg38" ]]; then
        bed_tss="$data_dir/$G/reference/GENCODE_TSSs_hg38.bed"
    elif [[ $G == "mm10" ]]; then
        bed_tss="$data_dir/$G/reference/GENCODE_TSSs_200bp_mm10.bed"
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
        if [[ $r == 1000 ]]; then
            expanded_coordinates=$(echo "$gene_coordinates" | awk -v OFS="\t" '{print $1, $2-250000, $3+250000}' | tr '[:blank:]' ':')
        else
            expanded_coordinates=$(echo "$gene_coordinates" | awk -v OFS="\t" '{print $1, $2-500000, $3+500000}' | tr '[:blank:]' ':')
        fi
        R="$expanded_coordinates"

    fi
fi

#----------------------------------

# Extract chromosome number for ordering later
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

    intra_start=$(echo "$R" | cut -d':' -f2)
    intra_end=$(echo "$R" | cut -d':' -f3)
    if [[ $intra_end -le $intra_start ]]; then
        echo -e "\nEnd coordinate must be greater than start in '$R'\n"
        exit 1
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
    num_chr1=$(extract_chr_number "$chr1")
    num_chr2=$(extract_chr_number "$chr2")

    # Capture the original orientation
    orig_range1="$(echo "$chr1" | cut -d':' -f1-2)"
    orig_range2="$(echo "$chr2" | cut -d':' -f1-2)"

    # Pull out starts and ends for the ordering swap and overlap guard below
    orig_start1="$(echo "$chr1" | cut -d':' -f2)"
    orig_start2="$(echo "$chr2" | cut -d':' -f2)"
    orig_end1="$(echo "$chr1" | cut -d':' -f3)"
    orig_end2="$(echo "$chr2" | cut -d':' -f3)"

    # Each window must have end > start
    if [[ $orig_end1 -le $orig_start1 ]]; then
        echo -e "\nEnd coordinate must be greater than start in '$chr1'.\n"
        exit 1
    fi
    if [[ $orig_end2 -le $orig_start2 ]]; then
        echo -e "\nEnd coordinate must be greater than start in '$chr2'.\n"
        exit 1
    fi

    # Same-chromosome windows must not overlap
    if [[ "$num_chr1" == "$num_chr2" ]]; then
        if [[ $orig_start1 -lt $orig_end2 && $orig_start2 -lt $orig_end1 ]]; then
            echo -e "Same-chromosome interchromosomal plots require non-overlapping windows."
            exit 1
        fi
    fi

    # Canonicalize to lower-first order for straw (upper triangle)
    # Different chr: by chr number; Same chr: by start coordinate
    swapped=FALSE
    if [[ "$num_chr1" != "$num_chr2" ]]; then
        lower=$(printf '%s\n%s\n' "$num_chr1" "$num_chr2" | sort -V | head -n1)
        if [[ "$lower" != "$num_chr1" ]]; then
            R="$chr2-$chr1"; swapped=TRUE
        fi
    elif [[ $orig_start1 -gt $orig_start2 ]]; then
        R="$chr2-$chr1"; swapped=TRUE
    fi

    if [[ "$v" == "TRUE" ]]; then
        echo -e "\nSVG output is not yet supported for interchromosomal plots. Using PDF instead...\n"
        v=FALSE
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


#----------------------------------

bin_size=$r

# interval_to_bins() in aqua_functions.r does the binning in r
# Binning here only builds the default output name

flank=$f

# Check if flank is a multiple of resolution
if (( flank % bin_size != 0 )); then
    echo -e "\nFlank value must be a multiple of the resolution, continuing with next multiple of resolution...\n"
    # Adjust flank to nearest multiple of resolution
    flank=$(( (flank / bin_size + 1) * bin_size ))
fi


# Adjust range for bin size and add flank
IFS='-' read -ra ADDR <<< "$R"
for addr in "${ADDR[@]}"; do
    IFS=':' read -ra COORDS <<< "$addr"
    start=${COORDS[1]}
    end=${COORDS[2]}

    # Apply flank to the raw coordinates
    start_with_flank=$(( start - flank ))
    if [ $start_with_flank -lt 0 ]; then
        start_with_flank=0
    fi
    end_with_flank=$(( end + flank ))

    # Range that becomes the contact matrix in r
    adjusted_ranges+=("${COORDS[0]}:$start_with_flank:$end_with_flank")

    # Label: the half-open span the plot will cover, [floor(start), ceil(end))
    # Reads back as what the user typed whenever the range is bin-aligned
    name_start=$(( start_with_flank / bin_size * bin_size ))
    name_end=$(( (end_with_flank + bin_size - 1) / bin_size * bin_size ))
    name_ranges+=("${COORDS[0]}:$name_start:$name_end")
done
R=$(IFS='-'; echo "${adjusted_ranges[*]}")
R_name=$(IFS='-'; echo "${name_ranges[*]}")

if [[ $T == TRUE && "$num_chr1" == "$num_chr2" ]]; then
    IFS=':' read -r _ s1 e1 <<< "${adjusted_ranges[0]}"
    IFS=':' read -r _ s2 e2 <<< "${adjusted_ranges[1]}"
    last_bin1=$(( (e1 - 1) / bin_size ))
    first_bin2=$(( s2 / bin_size ))
    if (( last_bin1 > first_bin2 )); then
        echo -e "\nSame-chromosome windows overlap after flanking (flank $flank). Reduce the flank or move the windows apart.\n"
        exit 1
    fi
fi


#----------------------------------

# Output name / extension checks
if [ "$D" == TRUE ]; then
  # Matrix mode: write a .txt file, ignore --svg
  ext="txt"
else
  case "$v" in
    TRUE)  ext="svg" ;;
    FALSE) ext="pdf" ;;
    *)     echo "--svg must be TRUE or FALSE (got '$v')." >&2; exit 1 ;;
  esac
fi


#----------------------------------

# Validate per-file label flags before pairing them to files
for val in "${ANNOTATION_LABELS[@]}"; do
  if [[ "$val" != "TRUE" && "$val" != "FALSE" ]]; then
    echo -e "\n--annotations_custom_label takes one TRUE or FALSE per bed file, in the same order as --annotations_custom. Got: '$val'\n  e.g. --annotations_custom a.bed b.bed --annotations_custom_label TRUE FALSE\n" >&2
    exit 1
  fi
done

# Filter annotation files: remove files that don't exist or don't have at least 3 tab-separated columns
VALID_ANNOTATION_FILES=()
VALID_ANNOTATION_COLORS=()
VALID_ANNOTATION_LABELS=()
for idx in "${!ANNOTATION_FILES[@]}"; do
  file="${ANNOTATION_FILES[$idx]}"

  if [ ! -f "$file" ]; then
    echo -e "Annotation file '$file' does not exist. It will be skipped.\n" >&2
    continue
  fi

  if ! head -n 5 "$file" | awk -F'[ \t]+' '{ if (NF < 3) { exit 1 } }'; then
    echo -e "Annotation file '$file' does not have at least 3 tab or space separated columns. It will be skipped.\n" >&2
    continue
  fi

  VALID_ANNOTATION_FILES+=("$file")
  # keep this file's color, if the user supplied one for this position
  if [ "$idx" -lt "${#ANNOTATION_COLORS[@]}" ]; then
    VALID_ANNOTATION_COLORS+=("${ANNOTATION_COLORS[$idx]}")
  fi
  # keep this file's label flag; unspecified positions default to FALSE
  if [ "$idx" -lt "${#ANNOTATION_LABELS[@]}" ]; then
    VALID_ANNOTATION_LABELS+=("${ANNOTATION_LABELS[$idx]}")
  else
    VALID_ANNOTATION_LABELS+=("FALSE")
  fi
done

ANNOTATION_FILES=("${VALID_ANNOTATION_FILES[@]}")
ANNOTATION_COLORS=("${VALID_ANNOTATION_COLORS[@]}")
ANNOTATION_LABELS=("${VALID_ANNOTATION_LABELS[@]}")

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

if [ ${#ANNOTATION_LABELS[@]} -eq 0 ]; then
    n="NONE"
else
    n="${ANNOTATION_LABELS[*]}"
fi

#----------------------------------
# Bedpe file checks

VALIDATED_BEDPE_FILES=()
VALIDATED_BEDPE_COLORS=()
VALIDATED_BEDPE_LWDS=()
for idx in "${!BEDPE_FILES[@]}"; do
  file="${BEDPE_FILES[$idx]}"

  if [[ ! -f "$file" ]]; then
    echo -e "Bedpe file '$file' does not exist. It will be skipped.\n" >&2
    continue
  fi

  awk_output=$(
    awk -F'[ \t]+' '
      BEGIN { code=0; violations=0 }
      { sub(/\r$/, "", $0) }                   # strip CR (Windows endings)
      /^[[:space:]]*$/ { next }                # ignore blank lines

      NF < 6 {
        printf "BEDPE error: line %d has fewer than 6 columns\n", NR
        violations=1; next
      }

      # Header: if column 2 is not an integer, treat as header (fatal)
      $2 !~ /^[0-9]+$/ {
        printf "BEDPE error: headers are not allowed (line %d)\n", NR
        code=1; next
      }

      # Other coordinate columns must be integers
      ($3 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ || $6 !~ /^[0-9]+$/) {
        printf "BEDPE error: non-integer coordinates at line %d\n", NR
        violations=1; next
      }

      # Feet must be ordered
      ($3 <= $2 || $6 <= $5) {
        printf "BEDPE error: reversed/zero-length foot at line %d (require col3>col2 and col6>col5)\n", NR
        violations=1; next
      }

      END {
        if (code)       exit code   # 1 = header found (exit)
        if (violations) exit 2      # 2 = validation issues (skip)
        exit 0                      # 0 = OK
      }
    ' "$file"
  )
  status=$?

  if (( status == 1 )); then
    echo "$awk_output" >&2
    echo "Fix '$file' and retry." >&2
    exit 1
  elif (( status == 2 )); then
    echo "$awk_output" >&2
    echo "Continuing without --bedpe $file ..." >&2
    echo
    continue
  fi

  VALIDATED_BEDPE_FILES+=("$file")
  if [ "$idx" -lt "${#BEDPE_COLORS[@]}" ]; then
    VALIDATED_BEDPE_COLORS+=("${BEDPE_COLORS[$idx]}")
  fi
  if [ "$idx" -lt "${#BEDPE_LWDS[@]}" ]; then
    VALIDATED_BEDPE_LWDS+=("${BEDPE_LWDS[$idx]}")
  fi
done

# Keep only fully validated files
BEDPE_FILES=("${VALIDATED_BEDPE_FILES[@]}")
BEDPE_COLORS=("${VALIDATED_BEDPE_COLORS[@]}")
BEDPE_LWDS=("${VALIDATED_BEDPE_LWDS[@]}")

if [ ${#BEDPE_FILES[@]} -eq 0 ]; then
  P=FALSE
else
  P="${BEDPE_FILES[*]}"
fi

if [ ${#BEDPE_COLORS[@]} -eq 0 ]; then
  y=000000
else
  y="${BEDPE_COLORS[*]}"
fi

for v_col in "${BEDPE_COLORS[@]}"; do
  if ! [[ "$v_col" =~ ^[0-9A-Fa-f]{6}$ ]]; then
    echo -e "\n--bedpe_color must be 6 hex digits (e.g. C71585). Got: '$v_col'\n"
    exit 1
  fi
done

if [ ${#BEDPE_LWDS[@]} -eq 0 ]; then
  l=0.8
else
  l="${BEDPE_LWDS[*]}"
fi

for v_lwd in "${BEDPE_LWDS[@]}"; do
  if ! [[ "$v_lwd" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo -e "\n--bedpe_lwd must be a positive number (e.g. 0.8 or 1.0). Got: '$v_lwd'\n"
    exit 1
  fi
done

#----------------------------------
# Prepare output directory

prepare_output_target() {
  local requested="$1"

  # Expand leading ~
  if [[ "$requested" == "~/"* ]]; then
    requested="$HOME/${requested:2}"
  elif [[ "$requested" == "~" ]]; then
    requested="$HOME"
  fi

  local dirpart namepart
  dirpart="$(dirname "$requested")"
  namepart="$(basename "$requested")"

  # If user didn't include a path, write to current working directory
  if [[ "$requested" != */* || "$dirpart" == "." ]]; then
    out_dir="$PWD"
    O="$namepart"
  else
    # If path is relative, anchor it to current working directory
    if [[ "$dirpart" == /* ]]; then
      out_dir="$dirpart"
    else
      out_dir="$PWD/$dirpart"
    fi
    O="$namepart"
  fi

  # Create output directory if needed
  if [[ ! -d "$out_dir" ]]; then
    mkdir -p "$out_dir" || {
      echo -e "\nError: Failed to create output directory: $out_dir\n" >&2
      exit 1
    }
  fi
}

#----------------------------------


if [[ "$p" == "TRUE" ]]; then
  echo -e "--profiles has been deprecated. Contact Axiotl if needed. Continuing...\n"
  p=FALSE
fi

if [[ "$m" != "none" && "$q" != "1" ]]; then
  echo -e "\nPlease provide either max_cap or quant_cut, not both\n"
  exit 1
fi

if [[ "$Q" != "blank" && "$Q" != "none" && "$Q" != "cpm" && "$Q" != "aqua" ]]; then
  echo -e "\n--norm should strictly be none, cpm, or aqua in lower case\n"
  exit 1
fi

if [[ "$d" != "standard" && "$d" != "axiotl" && "$d" != "none" ]]; then
  echo -e "\n--annotations_default should strictly be standard, axiotl, or none. Got $d\n"
  exit 1
fi

###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if [ "$two_sample" == "FALSE" ]
then

  analysis_type="single_sample"

  # Validate contact color
  if [[ -z $o ]]; then
    # Set default color based on the value of 'inter'
    if [ "$T" = "TRUE" ]; then
        o="004554"
    else
        o="FF0000"
    fi
  fi

  if ! [[ "$o" =~ ^[0-9A-Fa-f]{6}$ ]]; then
    echo -e "\n--color_one_sample must be exactly 6 hex digits (e.g. FF0000). Got: $o\n"
    usage
    exit 1
  fi

  if [[ -n $t ]];
  then
    echo -e "\nOops! Looks like you provided the wrong parameter for contact color!\n"
    echo -e "Use --color_one_sample for single sample analyses\n"
    usage
    exit 1
  fi

  #----------------------------------
  echo "Commencing $analysis_type analysis"

  echo "resolution: $r"

  if [[ $T == "FALSE" ]]; then
     range_text=`echo "$R" | sed 's/:/_/g'`
     name_range=`echo "$R_name" | sed 's/:/_/g'` # for default output name

     chr=`echo "$R" | cut -d":" -f1`
     start=`echo "$R" | cut -d":" -f2`
     end=`echo "$R" | cut -d":" -f3`

  fi

  if [[ $T == "TRUE" ]]; then
    range_text=`echo "$R" | sed 's/[:-]/_/g'`
    if [[ "$swapped" == "TRUE" ]]; then
        IFS='-' read -ra _BR <<< "$R_name"
        name_range=$(echo "${_BR[1]}-${_BR[0]}" | sed 's/[:-]/_/g')
    else
        name_range=$(echo "$R_name" | sed 's/[:-]/_/g')
    fi

    # interchromosomal ranges calculated in R
    chr=NA
    start=NA
    end=NA
    tag=NA
  fi

  sample_dir="$sample_dir_A"
  out_dir=`pwd`

  # Build output name
  if [[ -z "${O:-}" ]]; then
    O="${base_dir_A}_${name_range}.${ext}"
  else
    case "$O" in
      *.pdf|*.PDF|*.svg|*.SVG|*.txt|*.TXT)
        O="${O%.*}.${ext}"
        ;;
      *)
        O="${O}.${ext}"
        ;;
    esac
  fi

  prepare_output_target "$O"

  # Build argument list
  args=(
  "$analysis_type"
  "$T"
  "$pair1" "$stat1"
  "$norm_method" "$unit" "$bin_size" "$prefix"
  "${range_text}_${bin_size}"
  "$out_dir"
  "$chr" "$start" "$end"
  "$G" "$p" "$o" "$d" "$x" "$c"
  "$m" "$Q" "$q" "$O" "$P" "$y" "$l" "$i"
  "$sample_dir"
  "$w" "$K" "$L" "$M" "$N" "$D" "$v" "$s"
  "$n"
  )

# Keep user supplied inter chr order
  if [[ $T == "TRUE" ]]; then
    args+=("$orig_range1" "$orig_range2")
  fi

  Rscript "$aqua_dir/plot_contacts.r" "${args[@]}" 2>&1 | grep -v "Fontconfig warning"
  exit_status=${PIPESTATUS[0]}

fi


###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if [ "$two_sample" == "TRUE" ]
then

  # Resolve sample B
  resolve_sample "$B" "$I" "sample2"
  pair2="$res_hic"
  stat2="$res_mgs"
  sample_dir_B="$res_dir"
  base_dir_B="$res_name"

  if [[ $i == "TRUE" ]];
  then
    echo -e "\n--inherent only applies to one-sample plots for now, exiting...\n"
    exit
  fi

  analysis_type="two_sample"

  #----------------------------------


  if [[ -z $t ]]; then
    if [ "$T" = "TRUE" ]; then
        t="837900-a63200"
    else
        t="1E90FF-C71585"
    fi
  fi

  # Validate two-sample color format: RRGGBB-RRGGBB
  if ! [[ $t =~ ^[0-9A-Fa-f]{6}-[0-9A-Fa-f]{6}$ ]]; then
    echo -e "\n--color_two_sample (-t) must be two 6-digit hex colors separated by '-'."
    echo -e "For example: 1E90FF-C71585. Got: $t\n"
    usage
    exit 1
  fi

  if [[ -n $o ]];
  then
    echo -e "\nOops! Looks like you provided the wrong parameter for contact color!\n"
    echo -e "Use --color_two_sample for two sample analyses\n"
    usage
    exit
  fi

  #----------------------------------
  echo "Commencing $analysis_type analysis"

  echo "resolution: $r"

  if [[ $T == "FALSE" ]]; then
     range_text=`echo "$R" | sed 's/:/_/g'`
     name_range=`echo "$R_name" | sed 's/:/_/g'`

     chr=`echo "$R" | cut -d":" -f1`
     start=`echo "$R" | cut -d":" -f2`
     end=`echo "$R" | cut -d":" -f3`

  fi

  if [[ $T == "TRUE" ]]; then
    range_text=`echo "$R" | sed 's/[:-]/_/g'`
    if [[ "$swapped" == "TRUE" ]]; then
        IFS='-' read -ra _BR <<< "$R_name"
        name_range=$(echo "${_BR[1]}-${_BR[0]}" | sed 's/[:-]/_/g')
    else
        name_range=$(echo "$R_name" | sed 's/[:-]/_/g')
    fi

    # interchromosomal ranges calculated in R
    chr=NA
    start=NA
    end=NA
    tag=NA
  fi

  out_dir=`pwd`

  sample_dirA="$sample_dir_A"
  sample_dirB="$sample_dir_B"

  # Build output name
  if [[ -z "${O:-}" ]]; then
    O="${base_dir_A}_${base_dir_B}_${name_range}.${ext}"
  else
    case "$O" in
      *.pdf|*.PDF|*.svg|*.SVG|*.txt|*.TXT)
        O="${O%.*}.${ext}"
        ;;
      *)
        O="${O}.${ext}"
        ;;
    esac
  fi

  prepare_output_target "$O"

  # Build argument list
  args=(
    "$analysis_type"
    "$T"
    "$pair1" "$pair2" "$stat1" "$stat2"
    "$norm_method" "$unit" "$bin_size" "$prefix"
    "${range_text}_${bin_size}"
    "$out_dir"
    "$chr" "$start" "$end"
    "$G" "$p" "$t" "$d" "$x" "$c"
    "$m" "$Q" "$q" "$O" "$P" "$y" "$l" "$w"
    "$sample_dirA" "$sample_dirB"
    "$D" "$K" "$L" "$M" "$N" "$i" "$v" "$s"
    "$n"
    )

  # Keep user supplied inter chr order
  if [[ $T == "TRUE" ]]; then
    args+=("$orig_range1" "$orig_range2")
  fi

  Rscript "$aqua_dir/plot_contacts.r" "${args[@]}" 2>&1 | grep -v "Fontconfig warning"
  exit_status=${PIPESTATUS[0]}

fi

if [ "$exit_status" -eq 0 ]; then
  out_to_copy="${out_dir}/${O}"

  if [ -f "$out_to_copy" ]; then
    echo -e "\nDone! Created $out_to_copy\n"
  fi
fi
