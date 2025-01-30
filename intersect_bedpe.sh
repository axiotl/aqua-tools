#!/bin/bash

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

temp_dir=$(mktemp -d /tmp/intersect_bedpe-XXXXXXXXXX)

function no_sigint {

    let ctrlc_count++
    if [[ $ctrlc_count == 1 ]]; then
        rm -rf $temp_dir
    else
        :
    fi

}
trap no_sigint EXIT

function usage {
    echo -e "usage : intersect_bedpe.sh -A PATH_TO_FIRST_BED_FILE -P PATH_TO_BEDPE_FILE [-F FLANK_THRESHOLD_IN_BASE_PAIR] [-B PATH_TO_SECOND_BED_FILE] [-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Given a bedpe file, prints those rows of the bedpe in standard out that intersect with rows of given bed file(s) on either foot of the pair"
    echo
    echo "Useful in extracting biological subsets from the bedpe"
    echo
    echo "----------------------------------"
    echo "OPTIONS"
    echo
    echo "   -A|--bed_A            : Path of the first bed file"
    echo "   -P|--bedpe            : Path of the bedpe file"
    echo "  [-B|--bed_B         ]  : Path of the second bed file"
    echo "  [-f|--flank         ]  : Genome distance (bp) the BED should be near either foot. Default is 0"
    echo "  [-v|--absence       ]  : If TRUE, reports bedpe rows that do not intersect bed rows. Default FALSE"
    echo "  [   --print_bed     ]  : If TRUE, reports intersecting rows of bed instead of bedpe. Default FALSE"
    echo "  [   --print_bool    ]  : If TRUE, keeps all bedpe rows and adds columns marking intersections and absences. Default FALSE"
    echo "  [-h|--help          ]     Help message"
    echo
    echo
    echo "----------------------------------"
    echo "EXPECTED OUTPUT"
    echo
    echo "<chr_up>    <start_up>    <end_up>    <chr_down>    <start_down>    <end_down>  |  <intersection_up>    <intersection_down>"
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
      "--bed_A")     set -- "$@" "-A" ;;
      "--bedpe")     set -- "$@" "-P" ;;
      "--flank")     set -- "$@" "-f" ;;
      "--bed_B")     set -- "$@" "-B" ;;
      "--absence")   set -- "$@" "-v" ;;
      "--print_bed") set -- "$@" "-b" ;;
      "--print_bool") set -- "$@" "-l" ;;
      "--display_mode") set -- "$@" "-m" ;;
      "--help")      set -- "$@" "-h" ;;
       *)            set -- "$@" "$arg"
  esac
done

f=0
b=FALSE
l=FALSE
v=FALSE
m="blank"

while getopts ":A:B:P:f:b:v:l:m:h" OPT
do
    case $OPT in
  A) A=$OPTARG;;
  P) P=$OPTARG;;
  f) f=$OPTARG;;
  B) B=$OPTARG;;
  v) v=$OPTARG;;
  b) b=$OPTARG;;
  l) l=$OPTARG;;
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

flank=$f

# Check if both $P and $A are provided
if [ -z "$P" ] || [ -z "$A" ]; then
    echo "Both --bedpe (-P) and --bed_A (-A) are required."
    usage
    exit 1
fi

# Check absence is true with display mode options
if [ "$v" = "TRUE" ]; then
    if [ "$m" != "blank" ]; then
        echo "--absence cannot be combined with --display_mode options."
        usage
        exit 1
    fi
fi

# Check print_bool is true and absence is also true
if [ "$v" = "TRUE" ]; then
    if [ "$l" = "TRUE" ]; then
        echo "Use --print_bool or --absence, not both."
        usage
        exit 1
    fi
fi

# Check print_bool is true and print_bed is also true
if [ "$b" != "FALSE" ] && [ "$l" != "FALSE" ]; then
    echo -e "\nChoose either --print_bool or --print_bed.\n"
    usage
    exit 1
fi

# Check if mode is specified with print_bed or print_bool
if [ "$m" != "blank" ]; then
    if [ "$b" = "TRUE" ] || [ "$l" = "TRUE" ]; then
    echo -e "\n--display_mode cannot be combined with print_bed or print_bool. \n"
    usage
    exit 1
    fi
fi


awk '
BEGIN {
    OFS="\t"
}
{
    # Check for tab separator
    if ($0 !~ /\t/) {
        print "The BEDPE file does not contain tab-separated values on line " NR | "cat 1>&2";
        exit 1;
    }

    # Remove carriage return (\r) from all fields if present
    for (i=1; i<=NF; i++) gsub(/\r$/, "", $i);

    # Validate the format of the BEDPE file
    if (NR == 1) {
        if (!($1 ~ /^chr/ && $4 ~ /^chr/ && $2+0 == $2 && $3+0 == $3 && $5+0 == $5 && $6+0 == $6)) {
            print "The BEDPE file does not have the expected format in line " NR ". Expected chr1 start1 end1 chr2 start2 end2" | "cat 1>&2";
            exit 1;
        }
    }

    # Check if the start is less than the end for both intervals
    if ($2 > $3 || $5 > $6) {
        print "Start coordinate is greater than end coordinate in the BEDPE file on line " NR | "cat 1>&2";
        exit 1;
    }

    # Compare and swap columns if necessary
    if ($2 > $5 || ($2 == $5 && $3 > $6)) {
        # Swap columns 1-3 with 4-6
        temp1 = $1; temp2 = $2; temp3 = $3;
        $1 = $4; $2 = $5; $3 = $6;
        $4 = temp1; $5 = temp2; $6 = temp3;
    }

    # Print the modified line
    print $0 > "'"$temp_dir/temp.bedpe"'"

    # Print the first 6 columns with the line number appended
    print $1, $2, $3, $4, $5, $6, NR > "'"$temp_dir/annotated.bedpe"'"
}' "$P"


# Capture the exit status of awk
awk_status=$?

# Check if awk exited with an error
if [ $awk_status -ne 0 ]; then
    echo "Exiting."
    exit $awk_status
fi


awk '{print $1"\t"$2"\t"$3"\t"$7}' "$temp_dir/annotated.bedpe" > "$temp_dir/foot1.bed"
awk '{print $4"\t"$5"\t"$6"\t"$7}' "$temp_dir/annotated.bedpe" > "$temp_dir/foot2.bed"

# Adjust the start and end positions by the flank distance
if [ "$flank" -ne 0 ]; then
    # Adjust foot1.bed
    awk -v flank="$flank" 'BEGIN {OFS="\t"} {print $1, ($2 - flank < 0 ? 0 : $2 - flank), $3 + flank, $4}' "$temp_dir/foot1.bed" > "$temp_dir/temp_foot1.bed"
    mv "$temp_dir/temp_foot1.bed" "$temp_dir/foot1.bed"

    # Adjust foot2.bed
    awk -v flank="$flank" 'BEGIN {OFS="\t"} {print $1, ($2 - flank < 0 ? 0 : $2 - flank), $3 + flank, $4}' "$temp_dir/foot2.bed" > "$temp_dir/temp_foot2.bed"
    mv "$temp_dir/temp_foot2.bed" "$temp_dir/foot2.bed"
fi

if [ -z "$B" ]; then
    two_bed="FALSE"
  if [ "$v" = "TRUE" ]; then
    if [ "$b" = "FALSE" ]; then
  
        # -v option for non-intersecting rows
        bedtools intersect -a "$temp_dir/foot1.bed" -b $A -v > $temp_dir/non_intersect_foot1.bed
        bedtools intersect -a "$temp_dir/foot2.bed" -b $A -v > $temp_dir/non_intersect_foot2.bed
        Rscript $aqua_dir/intersect_bedpe.r "$temp_dir/non_intersect_foot1.bed" "$temp_dir/non_intersect_foot2.bed" "$temp_dir/temp.bedpe" $b $l $v $m $two_bed

     elif [ "$b" = "TRUE" ]; then
        # switch -a and -b for non intersecting bed regions
        bedtools intersect -a $A -b "$temp_dir/foot1.bed" -v > "$temp_dir/non_intersect_A_with_foot1.bed"
        bedtools intersect -a $A -b "$temp_dir/foot2.bed" -v > "$temp_dir/non_intersect_A_with_foot2.bed"
        Rscript $aqua_dir/intersect_bedpe.r "$temp_dir/non_intersect_A_with_foot1.bed" "$temp_dir/non_intersect_A_with_foot2.bed" "$temp_dir/temp.bedpe" $b $l $v $m $two_bed
    fi
  else
    # get intersections
    bedtools intersect -a "$temp_dir/foot1.bed" -b $A -wa -wb > $temp_dir/intersect_foot1.bed
    bedtools intersect -a "$temp_dir/foot2.bed" -b $A -wa -wb > $temp_dir/intersect_foot2.bed
    Rscript $aqua_dir/intersect_bedpe.r $temp_dir/intersect_foot1.bed $temp_dir/intersect_foot2.bed "$temp_dir/temp.bedpe" $b $l $v $m $two_bed
  fi
fi

if [ -n "$B" ]; then
    two_bed="TRUE"
    if [ "$v" = "TRUE" ]; then
        if [ "$b" = "FALSE" ]; then
            # -v option for non-intersecting rows
            bedtools intersect -a "$temp_dir/foot1.bed" -b $A $B -v > $temp_dir/non_intersect_foot1.bed
            bedtools intersect -a "$temp_dir/foot2.bed" -b $A $B -v > $temp_dir/non_intersect_foot2.bed
            Rscript $aqua_dir/intersect_bedpe.r $temp_dir/non_intersect_foot1.bed $temp_dir/non_intersect_foot2.bed "$temp_dir/temp.bedpe" $b $l $v $m $two_bed 
        elif [ "$b" = "TRUE" ]; then
            # switch -a and -b for non intersecting bed regions
            bedtools intersect -a $A -b "$temp_dir/foot1.bed" -v > "$temp_dir/non_intersect_A_with_foot1.bed"
            bedtools intersect -a $A -b "$temp_dir/foot2.bed" -v > "$temp_dir/non_intersect_A_with_foot2.bed"
            bedtools intersect -a $B -b "$temp_dir/foot1.bed" -v > "$temp_dir/non_intersect_B_with_foot1.bed"
            bedtools intersect -a $B -b "$temp_dir/foot2.bed" -v > "$temp_dir/non_intersect_B_with_foot2.bed"

            # Concatenate non-intersecting regions 
            cat "$temp_dir/non_intersect_A_with_foot1.bed" "$temp_dir/non_intersect_B_with_foot1.bed" > "$temp_dir/non_intersect_combined_with_foot1.bed"
            cat "$temp_dir/non_intersect_A_with_foot2.bed" "$temp_dir/non_intersect_B_with_foot2.bed" > "$temp_dir/non_intersect_combined_with_foot2.bed"

            Rscript $aqua_dir/intersect_bedpe.r "$temp_dir/non_intersect_combined_with_foot1.bed" "$temp_dir/non_intersect_combined_with_foot2.bed" "$temp_dir/temp.bedpe" $b $l $v $m $two_bed 
        fi
    else
        # get intersecting rows
        bedtools intersect -a "$temp_dir/foot1.bed" -b $A $B -wa -wb -names bed_A bed_B > $temp_dir/intersect_foot1.bed
        bedtools intersect -a "$temp_dir/foot2.bed" -b $A $B -wa -wb -names bed_A bed_B > $temp_dir/intersect_foot2.bed
        Rscript $aqua_dir/intersect_bedpe.r $temp_dir/intersect_foot1.bed $temp_dir/intersect_foot2.bed "$temp_dir/temp.bedpe" $b $l $v $m $two_bed
    fi
fi

rm -r $temp_dir
