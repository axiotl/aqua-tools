#!/bin/bash

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

temp_dir=$(mktemp -d /tmp/union_bedpe-XXXXXXXXXX)


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



function usage {
    echo -e "usage : union_bedpe.sh -P PATHS_TO_BEDPE_FILES[-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Given n bedpe files, print union of intervals to standard out"
    echo
    echo "Multiple bedpe files can be supplied as -P <bedpe1> -P <bedpe2> -P <bedpe3> etc"
    echo 
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -P|--bedpe     : Paths of the bedpe files you want to use"
    echo "  [ -h|--help   ]    Help message"
    exit;
}


if [ $# -lt 1 ]
    then
    usage
    exit
fi



for arg in "$@"; do
  shift
  case "$arg" in
      "--bedpe")     set -- "$@" "-P" ;;
       *)            set -- "$@" "$arg"
  esac
done



BEDPE_FILES=()
while getopts ":P:h" OPT; do
    case $OPT in
        P) BEDPE_FILES+=("$OPTARG");;  # Add each -P argument to the array
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


if [ ${#BEDPE_FILES[@]} -eq 0 ]; then
    echo "Error: No .bedpe files provided"
    usage
    exit 1
fi


for bedpe in "${BEDPE_FILES[@]}"; do
    if [ ! -f "$bedpe" ]; then
        echo "Error: File '$bedpe' not found"
        exit 1
    fi

    awk -v fname="$bedpe" '
    BEGIN {
        OFS = "\t"
    }
    {
        # Check for tab separator
        if ($0 !~ /\t/) {
            print "BEDPE error in file " fname ": does not contain tab-separated values on line " NR | "cat 1>&2";
            exit 1;
        }

        # Remove carriage return (\r) from all fields if present
        for (i = 1; i <= NF; i++) gsub(/\r$/, "", $i);

        # Validate the format of the BEDPE file
        if (NR == 1) {
            if (!($1 ~ /^chr/ && $4 ~ /^chr/ && $2+0 == $2 && $3+0 == $3 && $5+0 == $5 && $6+0 == $6)) {
                print "BEDPE error in file " fname ": expected format not found on line " NR ". Expected chr1 start1 end1 chr2 start2 end2" | "cat 1>&2";
                exit 1;
            }
        }

        # Check if the start is less than the end for both intervals
        if ($2 > $3 || $5 > $6) {
            print "BEDPE error in file " fname ": start coordinate is greater than end coordinate on line " NR | "cat 1>&2";
            exit 1;
        }

        # Canonical ordering: swap if needed
        if ($2 > $5 || ($2 == $5 && $3 > $6)) {
            temp1 = $1; temp2 = $2; temp3 = $3;
            $1 = $4; $2 = $5; $3 = $6;
            $4 = temp1; $5 = temp2; $6 = temp3;
        }

        print $1, $2, $3, $4, $5, $6
    }' "$bedpe" > "$temp_dir/$(basename "$bedpe")"

    if [ $? -ne 0 ]; then
        echo "Validation failed for '$bedpe'"
        exit 1
    fi
done



temp_bedpe_files=($temp_dir/*.bedpe)

if [ ${#temp_bedpe_files[@]} -lt 2 ]; then
    echo "Error: Need at least two .bedpe files to perform union"
    exit 1
fi



output_file="$temp_dir/union_output.bedpe"
Rscript "$aqua_dir/union_bedpe.r" "${temp_bedpe_files[0]}" "${temp_bedpe_files[1]}" > "$output_file"


for (( i=2; i<${#temp_bedpe_files[@]}; i++ )); do
    input_file="$output_file"
    output_file="$temp_dir/union_output_$i.bedpe"
    Rscript "$aqua_dir/union_bedpe.r" "$input_file" "${temp_bedpe_files[$i]}" >  "$output_file"
done



final_output="${output_file}"
cat "$final_output" | uniq


trap "rm -rf $temp_dir" EXIT
