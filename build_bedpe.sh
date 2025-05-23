#!/bin/bash

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

temp_dir=$(mktemp -d build_bedpe-XXXXXXXXXX)

ctrlc_count=0

function no_sigint {

    let ctrlc_count++
    if [[ $ctrlc_count == 1 ]]; then
        rm -rf $temp_dir
    else
        :
    fi

}
trap no_sigint EXIT

filter_duplicate_feet() {
    cat $1 | perl -ne 'chomp;
        my @d = split /\t/, $_;
        
        my $original = join("\t", @d[0..5]);
        my $swapped = join("\t", @d[3..5], @d[0..2]);
        
        unless ($seen{$original} || $seen{$swapped}) {
            print "$_\n";
            $seen{$original} = 1;
            $seen{$swapped} = 1;
        }
    '
}

order_feet() {
    awk '
    BEGIN { OFS = "\t" }
    {
        if ($1 == $4 && $2 <= $5) {
            print
        } else {
            n = NF
            meta = ""
            if (n > 6) {
                for (i = 7; i <= n; i++) meta = meta OFS $i
            }
            print $4, $5, $6, $1, $2, $3 meta
        }
    }'
}

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
    echo "or by user-defined minimum and maximum distances."
    echo "Prints pairs in bedpe format to standard out."
    echo
    echo "Pairs can be use to query .hic files with query_bedpe"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -A|--bed_A          : Path to the first bed file"
    echo "   -B|--bed_B          : Path to the second bed file"
    echo "  [-T|--TAD          ] : Path to the TAD file (to restrict pairings outside of the TAD)"
    echo "  [-d|--min_dist     ] : Minimum distance between pairs used to drop results. Default 0"
    echo "  [-D|--max_dist     ] : Maximum distance between pairs used to drop results. Default 5000000"
    echo "  [-t|--get_trans    ] : If TRUE, prints trans pairs only. Default FALSE"
    echo "  [-h|--help         ]   Help message"
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
      "--bed_A")           set -- "$@" "-A" ;;
      "--bed_B")           set -- "$@" "-B" ;;
      "--TAD")             set -- "$@" "-T" ;;
      "--min_dist")        set -- "$@" "-d" ;;
      "--max_dist")        set -- "$@" "-D" ;;
      "--get_trans")       set -- "$@" "-t" ;;
      "--help")            set -- "$@" "-h" ;;
       *)                  set -- "$@" "$arg"
  esac
done



T="NULL"
d=0
D=5000000
t="FALSE"

while getopts ":A:B:T:d:D:t:h" OPT
do
    case $OPT in
  A) A=$OPTARG;;
  B) B=$OPTARG;;
  T) T=$OPTARG;;
  d) d=$OPTARG;;
  D) D=$OPTARG;;
  t) t=$OPTARG;;
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

# -----------------------------------------
# Parameter checks

if [ -z "$A" ]; then
  echo -e "\nMissing required argument -A (path to first BED file)." >&2
  exit 1
elif [ ! -f "$A" ]; then
  echo -e "\nFile '$A' does not exist." >&2
  exit 1
fi

if [ -z "$B" ]; then
  echo -e "\nMissing required argument -B (path to second BED file)." >&2
  exit 1
elif [ ! -f "$B" ]; then
  echo -e "\nFile '$B' does not exist." >&2
  exit 1
fi

if [[ "$T" != "NULL" ]] && [ ! -f "$T" ]; then
  echo -e "\nTAD file '$T' does not exist." >&2
  exit 1
fi

if ! [[ "$d" =~ ^[0-9]+$ ]]; then
  echo -e "\nMinimum distance (-d) must be numeric only (e.g., 0, 100, 100000)" >&2
  exit 1
fi

if ! [[ "$D" =~ ^[0-9]+$ ]]; then
  echo -e "\nMaximum distance (-D) must be numeric only (e.g., 0, 100, 100000)" >&2
  exit 1
fi

if [[ "$d" -ge "$D" ]]; then
  echo -e "\nMinimum distance (d=$d) must be less than maximum distance (D=$D)" >&2
  exit 1
fi

if [[ "$d" -lt 0 ]]; then 
  echo "min dist cannot be less than 0"
  exit 1
fi

if [[ "$t" != "TRUE" && "$t" != "FALSE" ]]; then 
  echo "--get_trans can either be TRUE or FALSE"
  exit 1
fi

if [[ "$T" != "NULL" ]] && [[ "$t" == "TRUE" ]]; then
  echo -e "\nCannot use a TAD (-T) file for trans pairs (-t TRUE)" >&2
  exit 1
fi

# -----------------------------------------

if [[ $t == "FALSE" ]]; then
    if [[ $T == "NULL" ]]; then
        awk -v D="$D" -v d="$d" '
        BEGIN { OFS=FS="\t" }
        {
            original_2 = $2;
            original_3 = $3;
            $2 = $2 - D; if ($2 < 1) $2 = 1;
            $3 = original_2 - d; if ($3 < 1) $3 = 1;
            print $0, $1, original_2, original_3;
        }' "$B" > "$temp_dir/file_b1.bed"

        awk -v D="$D" -v d="$d" '
        BEGIN { OFS=FS="\t" }
        {
            original_2 = $2;
            original_3 = $3;
            $2 = $2 + d; if ($2 < 1) $2 = 1;
            $3 = original_2 + D; if ($3 < 1) $3 = 1;
            print $0, $1, original_2, original_3;
        }' "$B" > "$temp_dir/file_b2.bed"

        bedtools intersect -a "$A" -b "$temp_dir/file_b1.bed" -wa -wb > "$temp_dir/intersect_left.bedpe"
        bedtools intersect -a "$A" -b "$temp_dir/file_b2.bed" -wa -wb > "$temp_dir/intersect_right.bedpe"

        num_columns_A=$(awk -F'\t' 'NR==1 {print NF; exit}' "$A")
        num_columns_B=$(awk -F'\t' 'NR==1 {print NF; exit}' "$B")
        ml=$((num_columns_A - 3))
        mr=$((num_columns_B - 3))

        cat "$temp_dir/intersect_left.bedpe" "$temp_dir/intersect_right.bedpe" > "$temp_dir/intersected_file.bedpe"

        cat "$temp_dir/intersected_file.bedpe" | 
        
        awk -v ml="$ml" -v mr="$mr" 'BEGIN {OFS="\t"} {
            # Print A columns (1–3) without tab before 1st value, rest with leading tab
            printf $1
            for (i = 2; i <= 3; i++) {
                printf OFS $i
            }

            # Print B columns (7+ml+mr – 9+ml+mr)
            for (i = 7+ml+mr; i <= 9+ml+mr; i++) {
                printf OFS $i
            }

            # A metadata (4–3+ml)
            if (ml > 0) {
                for (i = 4; i <= 3+ml; i++) {
                    printf OFS $i
                }
            }

            # B metadata (7+ml – 6+ml+mr)
            if (mr > 0) {
                for (i = 7+ml; i <= 6+ml+mr; i++) {
                    printf OFS $i
                }
            }

            printf "\n"
        }' "$temp_dir/intersected_file.bedpe" > "$temp_dir/desired_file.bedpe"


        # Process the desired_file.bedpe to enforce the distance constraints
        awk -v d="$d" -v D="$D" 'BEGIN { OFS="\t" } {
            diff = ($5 > $2 ? $5 - $2 : $2 - $5);
            if (diff >= d && diff <= D)
                print $0;
        }' "$temp_dir/desired_file.bedpe" | uniq | filter_duplicate_feet | order_feet

    elif [[ $T != "NULL" ]]; then

        num_columns_A=$(awk -F'\t' 'NR==1 {print NF; exit}' "$A")
        num_columns_B=$(awk -F'\t' 'NR==1 {print NF; exit}' "$B")
        num_columns_T=$(awk -F'\t' 'NR==1 {print NF; exit}' "$T")

        ma=$((num_columns_A - 3))
        mb=$((num_columns_B - 3))
        mt=$((num_columns_T - 3))

        bedtools intersect -a "$T" -b "$A" -wa -wb -F 1 | \
        awk -v n="$num_columns_T" 'BEGIN {OFS="\t"} {printf $1"_"$2"_"$3 OFS; for (i=n+1; i<=NF; i++) printf $i OFS; print ""}' | \
        sort -k1,1 > "$temp_dir/temp_bed_A.bed"

        bedtools intersect -a "$T" -b "$B" -wa -wb -F 1 | \
        awk -v n="$num_columns_T" 'BEGIN {OFS="\t"} {printf $1"_"$2"_"$3 OFS; for (i=n+1; i<=NF; i++) printf $i OFS; print ""}' | \
        sort -k1,1 > "$temp_dir/temp_bed_B.bed"

        join -t $'\t' -j 1 "$temp_dir/temp_bed_A.bed" "$temp_dir/temp_bed_B.bed" | \
        awk 'BEGIN {OFS="\t"} {$1=$1; print}' | \
        cut -f 2- | \
        awk -v ml="$ma" -v mr="$mb" 'BEGIN {OFS="\t"} {
            # Print A columns (1–3)
            printf $1
            for (i=2; i<=3; i++) {
                printf OFS $i
            }

            # Print B columns (4+ml to 6+ml)
            for (i=4+ml; i<=6+ml; i++) {
                printf OFS $i
            }

            # Print A-meta columns (4 to 3+ml)
            if (ml > 0) {
                for (i=4; i<=3+ml; i++) {
                    printf OFS $i
                }
            }

            # Print B-meta columns (7+ml to 6+ml+mr)
            if (mr > 0) {
                for (i=7+ml; i<=6+ml+mr; i++) {
                    printf OFS $i
                }
            }

            printf "\n"
        }' | uniq | filter_duplicate_feet | order_feet > "$temp_dir/tad_output.bedpe"


        # If non-default distance constraints are specified, apply additional filtering
        if [ "$d" -ne 0 ] || [ "$D" -ne 5000000 ]; then
            awk -v d="$d" -v D="$D" 'BEGIN { OFS="\t" } {
                diff = ($5 > $2 ? $5 - $2 : $2 - $5);
                if (diff >= d && diff <= D)
                    print $0;
            }' "$temp_dir/tad_output.bedpe"
        else
            cat "$temp_dir/tad_output.bedpe"
        fi

    fi

elif [[ $t == "TRUE" ]]; then

    num_columns_A=$(awk -F'\t' 'NR==1 {print NF; exit}' "$A")
    num_columns_B=$(awk -F'\t' 'NR==1 {print NF; exit}' "$B")
    ma=$((num_columns_A - 3))
    mb=$((num_columns_B - 3))
    
    chroms1=$(awk '{print $1}' "$A" | sort -u)
    chroms2=$(awk '{print $1}' "$B" | sort -u)

    for chr1 in $chroms1; do
        for chr2 in $chroms2; do
            if [ "$chr1" != "$chr2" ]; then
                awk -v chr="$chr2" 'BEGIN {OFS="\t"} $1 == chr {print}' "$B" | \
                while read -r lineB; do
                    awk -v lineB="$lineB" -v chr="$chr1" -v ma="$ma" -v mb="$mb" 'BEGIN {OFS="\t"} $1 == chr {
                        split(lineB, arrB, OFS);

                        # Print As 3 columns
                        printf $1
                        printf OFS $2
                        printf OFS $3

                        # Print Bs 3 columns
                        printf OFS arrB[1]
                        printf OFS arrB[2]
                        printf OFS arrB[3]

                        # Print A meta
                        for (i=4; i<=3+ma; i++) {
                            printf OFS $i
                        }

                        # Print B meta
                        for (i=4; i<=3+mb; i++) {
                            printf OFS arrB[i]
                        }

                        printf "\n"
                    }' "$A"
                done
            fi
        done
    done
fi


trap "rm -rf $temp_dir" EXIT
