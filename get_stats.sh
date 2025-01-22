#!/bin/bash

data_dir=$HOME/lab-data
sample_sheet="$HOME/setup/sample_sheet.txt"

# Ensure the sample sheet exists
#if [ ! -f "$sample_sheet" ]; then
#    python3 /home/ubuntu/aqua_tools/restore_sample_sheet.py > /dev/null 2>&1
#    if [ ! -f "$sample_sheet" ]; then
#        echo "Failed to restore sample sheet. Exiting."
#        exit 1
#    fi
#fi

function usage {
    echo "usage : get_stats.sh --sample SAMPLE_NAME --genome GENOME"
    echo -e "Use option -h|--help for more information"
}

function help {
    echo 
    echo "Retrieve statistics for the specified sample and genome build"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    --sample     : Name of the sample (or comma-separated list of sample names such as sample1,sample2)"
    echo "    --genome     : Genome build - strictly hg19, hg38, or mm10"
    echo "  -h|--help      : Help message"
    exit;
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --sample)
            SAMPLE_NAMES="$2"
            shift
            shift
            ;;
        --genome)
            G="$2"
            shift
            shift
            ;;
        -h|--help)
            help
            exit
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done

# Validate arguments
if [ -z "$SAMPLE_NAMES" ] || [ -z "$G" ]; then
    usage
    exit 1
fi

# Check genome parameter
if [[ "$G" != "hg19" && "$G" != "hg38" && "$G" != "mm10" ]]; then
    echo "Error: Invalid genome parameter. Allowed values are hg19, hg38, or mm10."
    exit 1
fi

get_sample_directory() {

    local input_sample_name="$1"
    local base_sample_name
    local sample_dir

    if [[ "$input_sample_name" =~ _version_[0-9]+$ ]]; then
        base_sample_name=$(echo "$input_sample_name" | grep -oP '.*(?=_version_[0-9]+$)' )
        sample_dir="$data_dir/$G/$base_sample_name/$input_sample_name"
    else
        base_sample_name="$input_sample_name"
        local sample_info=$(awk -F$'\t' -v sample="$base_sample_name" '$1 == sample && $6 == "1" {print; exit}' "$sample_sheet")
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

# Prepare header for output with tab-separated columns
header="Sample Name\tTotal_pairs\tvalid_int\tvalid_int_rmdup\ttrans_int\tcis_int\tcis_short\tcis_long"
echo -e "$header"

# Process each sample name
IFS=',' read -ra SAMPLES <<< "$SAMPLE_NAMES"
for A in "${SAMPLES[@]}"; do
    sample_dir_A=$(get_sample_directory "$A")
    if [ $? -eq 0 ]; then
        version_dir_A=$(basename "$sample_dir_A")
        base_dir_A=$(basename "$(dirname "$sample_dir_A")")
        full_sample_name="${base_dir_A}/${version_dir_A}"

        # Get Total_pairs_processed value
        pairstats_file="$sample_dir_A/${version_dir_A}.pairStats.txt"
        if [ -f "$pairstats_file" ]; then
            total_pairs_processed=$(awk '$1 == "Total_pairs_processed" {print $2}' "$pairstats_file")
            total_pairs_processed="${total_pairs_processed:-.}"
        else
            total_pairs_processed="."
        fi

        # Initialize array for mergestats values
        mergestats_values=()

        # Collect mergestats values
        mergestats_file="$sample_dir_A/mergeStats.txt"
        if [ -f "$mergestats_file" ]; then
            for header in "valid_interaction" "valid_interaction_rmdup" "trans_interaction" "cis_interaction" "cis_shortRange" "cis_longRange"; do
                value=$(awk -v key="$header" '$1 == key {print $2}' "$mergestats_file")
                mergestats_values+=("${value:-.}")
            done
        else
            for _ in {1..6}; do
                mergestats_values+=(".")
            done
        fi

        # Replace any dash values with dots
        for i in "${!mergestats_values[@]}"; do
            if [[ "${mergestats_values[$i]}" == "-" ]]; then
                mergestats_values[$i]="."
            fi
        done

        # Combine values with sample name as the first value
        combined_values=("$A" "$total_pairs_processed" "${mergestats_values[@]}")
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${combined_values[@]}"
    else
        echo "Sample directory not found for $A"
    fi
done
