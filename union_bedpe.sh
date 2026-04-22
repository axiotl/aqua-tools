#!/bin/bash

data_dir="${LAB_DATA_DIR:-$HOME/lab-data}"
aqua_dir="${AQUA_TOOLS_DIR:-$HOME/aqua_tools}"

temp_dir=$(mktemp -d /tmp/union_bedpe-XXXXXXXXXX)

ctrlc_count=0
function no_sigint {
    let ctrlc_count++
    if [[ $ctrlc_count == 1 ]]; then
        rm -rf "$temp_dir"
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
    echo "A single bedpe file can be supplied to merge overlapping pairs within that file"
    echo "Multiple bedpe files can be supplied as -P <bedpe1> -P <bedpe2> -P <bedpe3> etc"
    echo 
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -P|--bedpe     : Paths of the bedpe files you want to use"
    echo "  [ -h|--help   ]    Help message"
    exit;
}

if [ $# -lt 1 ]; then
    usage
    exit
fi

for arg in "$@"; do
  shift
  case "$arg" in
      "--bedpe")     set -- "$@" "-P" ;;
       *)            set -- "$@" "$arg" ;;
  esac
done

BEDPE_FILES=()
while getopts ":P:h" OPT; do
    case $OPT in
        P) BEDPE_FILES+=("$OPTARG");;
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
    echo "Error: No input files provided"
    usage
    exit 1
fi

# Record the cleaned temp outputs here (forced .bedpe extension)
temp_bedpe_files=()
idx=0

for bedpe in "${BEDPE_FILES[@]}"; do
    if [ ! -f "$bedpe" ]; then
        echo "Error: File '$bedpe' not found"
        exit 1
    fi

    base="$(basename "$bedpe")"
    # strip any extension so input can be .bedpe/.bed/.txt/etc.
    base="${base%.*}"

    out="$temp_dir/${base}_${idx}.bedpe"
    idx=$((idx+1))

    awk -v fname="$bedpe" '
    BEGIN { OFS = "\t" }
    {
        # Check for tab separator
        if ($0 !~ /\t/) {
            print "BEDPE error in file " fname ": does not contain tab-separated values on line " NR | "cat 1>&2";
            exit 1;
        }

        # Remove carriage return (\r) from all fields if present
        for (i = 1; i <= NF; i++) gsub(/\r$/, "", $i);

        # Validate chromosome and numeric fields on every line
        if (!($1 ~ /^chr/ && $4 ~ /^chr/)) {
            print "BEDPE error in file " fname ": non-chromosome value in column 1 or 4 on line " NR | "cat 1>&2";
            exit 1;
        }
        if ($2+0 != $2 || $3+0 != $3 || $5+0 != $5 || $6+0 != $6) {
            print "BEDPE error in file " fname ": non-numeric coordinate on line " NR | "cat 1>&2";
            exit 1;
        }

        # Check if the start is less than the end for both intervals
        if ($2 > $3 || $5 > $6) {
            print "BEDPE error in file " fname ": start coordinate is greater than end coordinate on line " NR | "cat 1>&2";
            exit 1;
        }
        # canonicalize cis and trans
        if ($1 > $4 || ($1 == $4 && ($2 > $5 || ($2 == $5 && $3 > $6)))) {
            temp1 = $1; temp2 = $2; temp3 = $3;
            $1 = $4; $2 = $5; $3 = $6;
            $4 = temp1; $5 = temp2; $6 = temp3;
        }

        print $1, $2, $3, $4, $5, $6
    }' "$bedpe" > "$out" || exit 1

    temp_bedpe_files+=("$out")
done


output_file="$temp_dir/union_output.bedpe"

if [ ${#temp_bedpe_files[@]} -eq 1 ]; then
    # Self-merge: compare the file against itself to merge overlapping pairs
    Rscript "$aqua_dir/union_bedpe.r" "${temp_bedpe_files[0]}" "${temp_bedpe_files[0]}" > "$output_file"
else
    # Multiple files: merge first two, then iteratively add the rest
    Rscript "$aqua_dir/union_bedpe.r" "${temp_bedpe_files[0]}" "${temp_bedpe_files[1]}" > "$output_file"

    for (( i=2; i<${#temp_bedpe_files[@]}; i++ )); do
        input_file="$output_file"
        output_file="$temp_dir/union_output_$i.bedpe"
        Rscript "$aqua_dir/union_bedpe.r" "$input_file" "${temp_bedpe_files[$i]}" >  "$output_file"
    done
fi

final_output="${output_file}"

sort -u "$final_output"

