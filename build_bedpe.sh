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
    echo "  [-d|--min_dist     ] : Minimum distance between pairs used to drop results. Default 0 bp"
    echo "  [-D|--max_dist     ] : Maximum distance between pairs used to drop results. Default 5 Mb"
    echo "  [-t|--get_trans    ] : If TRUE, prints trans pairs only. Default = FALSE"
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


if [[ "$d" -lt 0 ]]; then 
    echo "min dist cannot be less than 0"
    exit 1
fi

if [[ "$t" != "TRUE" && "$t" != "FALSE" ]]; then 
    echo "--get_trans can either be TRUE or FALSE"
    exit 1
fi





if [[ $t == "FALSE" ]]; then
   if [[ $T == "NULL" ]]; then

       awk -v D="$D" -v d="$d" '
       BEGIN {
          OFS=FS="\t";
       }
       {
           original_2 = $2;
           original_3 = $3;

           $2 = $2 - D;
           if ($2 < 1) $2 = 1;


           $3 = original_2 - d;
           if ($3 < 1) $3 = 1;

           print $0, $1, original_2, original_3;
       }' $B > $temp_dir/file_b1.bed


       awk -v D="$D" -v d="$d" '
       BEGIN {
           OFS=FS="\t";
       }
       {
           original_2 = $2;
           original_3 = $3;

           $2 = $2 + d;
           if ($2 < 1) $2 = 1;

           $3 = original_2 + D;
           if ($3 < 1) $3 = 1;

           print $0, $1, original_2, original_3;
       }' $B > $temp_dir/file_b2.bed


       bedtools intersect -a $A -b $temp_dir/file_b1.bed -wa -wb > $temp_dir/intersect_left.bedpe
       bedtools intersect -a $A -b $temp_dir/file_b2.bed -wa -wb > $temp_dir/intersect_right.bedpe

       num_columns_A=$(awk -F'\t' 'NR==1 {print NF; exit}' "$A")
       num_columns_B=$(awk -F'\t' 'NR==1 {print NF; exit}' "$B")

       ml=$((num_columns_A - 3))
       mr=$((num_columns_B - 3))


       cat $temp_dir/intersect_left.bedpe $temp_dir/intersect_right.bedpe > $temp_dir/intersected_file.bedpe

       cat $temp_dir/intersected_file.bedpe | 
       awk -v ml=$ml -v mr=$mr 'BEGIN {OFS="\t"} {
           # Print A columns
           for (i=1; i<=3; i++) {
               printf $i OFS
           }
        
           # Print C columns
           for (i=7+ml+mr; i<=9+ml+mr; i++) {
               printf $i OFS
           }
        
           # Print A-meta columns
           if (ml > 0) {
               for (i=4; i<=3+ml; i++) {
                   printf $i OFS
               }
           }
        
           # Print B-meta columns
           if (mr > 0) {
               for (i=7+ml; i<=6+ml+mr; i++) {
                   printf $i OFS
               }
           }
        
           printf "\n"
       }' $temp_dir/intersected_file.bedpe > $temp_dir/desired_file.bedpe
       
       # Process the desired_file.bedpe to enforce the distance constraints
       awk -v d=$d -v D=$D 'BEGIN { OFS="\t" } {
          # Calculate the absolute difference between the start positions in columns 2 and 5
          diff = ($5 > $2 ? $5 - $2 : $2 - $5);

          # Print the line only if the difference is within the specified constraints
          if (diff >= d && diff <= D) {
          print $0;
          }
    }' $temp_dir/desired_file.bedpe | \
    uniq | filter_duplicate_feet

fi


   if [[ $T != "NULL" ]]; then
    
    
       num_columns_A=$(awk -F'\t' 'NR==1 {print NF; exit}' "$A")
       num_columns_B=$(awk -F'\t' 'NR==1 {print NF; exit}' "$B")
       num_columns_T=$(awk -F'\t' 'NR==1 {print NF; exit}' "$T")

       ma=$((num_columns_A - 3))
       mb=$((num_columns_B - 3))
       mt=$((num_columns_T - 3))
    

       bedtools intersect -a $T -b $A -wa -wb -F 1 | \
       awk -v n=$num_columns_T 'BEGIN {OFS="\t"} {printf $1"_"$2"_"$3 OFS; for (i=n+1; i<=NF; i++) printf $i OFS; print ""}' | \
       sort -k1,1 > $temp_dir/temp_bed_A.bed

       bedtools intersect -a $T -b $B -wa -wb -F 1 | \
       awk -v n=$num_columns_T 'BEGIN {OFS="\t"} {printf $1"_"$2"_"$3 OFS; for (i=n+1; i<=NF; i++) printf $i OFS; print ""}' | \
       sort -k1,1 > $temp_dir/temp_bed_B.bed

       join -t $'\t' -j 1 $temp_dir/temp_bed_A.bed $temp_dir/temp_bed_B.bed |

       awk 'BEGIN {OFS="\t"} {$1=$1; print}' | \
       cut -f 2- | \
       awk -v ml=$ma -v mr=$mb 'BEGIN {OFS="\t"} {
           # Print A columns
           for (i=1; i<=3; i++) {
               printf $i OFS
           }
        
           # Print C columns
           for (i=4+ml; i<=6+ml; i++) {
               printf $i OFS
           }
        
           # Print A-meta columns
           if (ml > 0) {
               for (i=4; i<=3+ml; i++) {
                   printf $i OFS
               }
           }
        
           # Print B-meta columns
           if (mr > 0) {
               for (i=7+ml; i<=6+ml+mr; i++) {
                   printf $i OFS
               }
           }
        
           printf "\n"
       }' | \
       uniq | filter_duplicate_feet
   fi

fi




if [[ $t == "TRUE" ]]; then
    
    
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
                        printf $1 OFS $2 OFS $3 OFS;
                        printf arrB[1] OFS arrB[2] OFS arrB[3] OFS;
                        for (i=4; i<4+ma; i++) {
                            printf $i OFS;
                        }
                        for (i=4; i<4+mb; i++) {
                            printf arrB[i] OFS;
                        }
                        printf "\n";
                    }' "$A"
                done
            fi
        done
    done


fi

trap "rm -rf $temp_dir" EXIT
