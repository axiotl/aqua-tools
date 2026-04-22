#!/bin/bash

data_dir="${LAB_DATA_DIR:-$HOME/lab-data}"
aqua_dir="${AQUA_TOOLS_DIR:-$HOME/aqua_tools}"

temp_dir=$(mktemp -d cluster_bedpe-XXXXXXXXXX)

cleanup() {
    rm -rf "$temp_dir"
}
trap cleanup EXIT

usage() {
    echo -e "usage : cluster_bedpe.sh -P PATH_TO_BEDPE_FILE [-f FLANK_IN_BASE_PAIR] [-h]"
    echo -e "Use option -h|--help for more information"
}

help() {
    echo
    echo "Given a bedpe file, clusters bedpe rows based on basepair intersections"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -P|--bedpe       : Bedpe file path"
    echo "  [-f|--flank    ]  : Genome distance in bp to increase bedpe pairs. Default 0"
    echo "  [-h|--help     ]  : Help message"
    exit
}

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

# Transform long options to short ones
for arg in "$@"; do
    shift
    case "$arg" in
        --bedpe) set -- "$@" "-P" ;;
        --flank) set -- "$@" "-f" ;;
        --help)  set -- "$@" "-h" ;;
        *)       set -- "$@" "$arg" ;;
    esac
done

f=0
while getopts ":P:f:h" OPT; do
    case $OPT in
        P) P=$OPTARG ;;
        f) f=$OPTARG ;;
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

if [ -z "$P" ]; then
    echo "Error: -P|--bedpe is required." >&2
    usage
    exit 1
fi

Rscript "$aqua_dir/cluster_bedpe.r" "$P" "$f" "$temp_dir"
