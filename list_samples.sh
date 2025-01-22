#!/bin/bash

sample_sheet="$HOME/setup/sample_sheet.txt"


print_samples_for_genome() {
    local G="$1"  # Genome build
    local has_samples=0  # Flag to check if samples exist
    local printed_header=0  # Flag to check if header is printed

    # Process the sample sheet
    awk -v genome="$G" -F '\t' '
     function add_commas(num) {
        if (num == "." || num == "" || num == "0") return "."
        result = ""
        while (length(num) > 3) {
            result = "," substr(num, length(num)-2) result
            num = substr(num, 1, length(num)-3)
        }
        return num result
    }

    # Process each line of the sample sheet
    NR>1 && $2 == genome {
        default_val = ($6 == "1" ? "Y" : "N");  # Determine default value
        sample_dir="/home/ubuntu/lab-data/" genome "/" $1  # Sample directory path
        if (system("[ -d \"" sample_dir "\" ]") == 0) {  # Check if sample directory exists
            has_samples = 1  # Set flag if samples exist
            versioned_sample = $1 "_version_" $5  # Construct versioned sample name
            cmd = "bash get_stats.sh --genome " genome " --sample " versioned_sample # Call get_stats.sh
            while ((cmd | getline stats_output) > 0) {
                if (stats_output !~ /^Sample Name/) {  # Skip header line
                    split(stats_output, stats_array, "\t")  # Split stats output into array
                    valid_int_rmdup = add_commas(stats_array[4])  # Extract and format valid_int_rmdup
                    cis_short = add_commas(stats_array[7])  # Extract and format cis_short
                    cis_long = add_commas(stats_array[8])  # Extract and format cis_long
                    if (printed_header == 0) {  # Print header if not already printed
                        printf "\n"
                        printf "--------------------------------\n"
                        printf "|                              |\n"
                        printf "|             %s             |\n", genome
                        printf "|                              |\n"
                        printf "--------------------------------\n"
                        printf "\n"
                        printf "%-30s|%-6s|%-9s|%-7s|%-7s|%-15s|%-15s|%-15s\n", "Sample Name", "AQuA", "Inherent", "Version", "Default", "valid_int_rmdup", "cis_short", "cis_long"
                        printf "------------------------------+------+---------+-------+-------+---------------+---------------+---------------\n"
                        printed_header = 1
                    }
                    # Print sample details
                    printf "%-30s|%-6s|%-9s|%-7s|%-7s|%-15s|%-15s|%-15s\n", $1, $3, $4, $5, default_val, valid_int_rmdup, cis_short, cis_long
                }
            }
            close(cmd)
        }
    }
    ' "$sample_sheet"
    
    # Print a newline if samples were found
    if [[ $has_samples -eq 1 ]]; then
        echo
    fi
}

# Call the function for each genome build
print_samples_for_genome "hg19"
print_samples_for_genome "hg38"
print_samples_for_genome "mm10"
