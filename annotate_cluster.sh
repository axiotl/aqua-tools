#!/bin/bash 

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools
ref_dir=$data_dir/hg38/reference

temp_dir=$(mktemp -d /tmp/annotate_cluster-XXXXXXXXXX)

ctrlc_count=0

function no_sigint {

    let ctrlc_count++
    if [[ $ctrlc_count == 1 ]]; then
        sudo rm -rf $temp_dir
    else
        :
    fi

}
trap no_sigint EXIT



function usage {
    echo "usage : annotate_cluster.sh \\"
    echo "  -P PATH_TO_BEDPE_FILE \\"
    echo "  -A SAMPLE_NAME \\"
    echo "  -G GENOME \\"
    echo "  -B PATH_TO_PEAKS_BED_FILE \\"
    echo " [-h]"
    echo "Use option -h|--help for more information"
}


function help {
    echo
    echo "Annotate a cluster_bedpe object with biological signals"
    echo
    echo
    echo "Input bedpe should be 7-columns only, with column 7 indicating cluster membership"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -P|--bedpe          : Path to 7-col .bedpe"
    echo "   -A|--sample1        : Name of sample"
    echo "   -G|--genome         : Genome build used for sample processing"
    echo "   -B|--bed            : Path to the sample's H3K27ac/ATAC/DNAse peaks .bed file"
    echo "  [-U|--user_bed      ]: Paths to any other .bed files"
    echo "  [-h|--help          ]  Help message"
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
      "--bedpe")           set -- "$@" "-P" ;;
      "--sample1")         set -- "$@" "-A" ;;
      "--genome")          set -- "$@" "-G" ;;
      "--bed")             set -- "$@" "-B" ;;
      "--user_bed")        set -- "$@" "-U" ;;
      "--help")            set -- "$@" "-h" ;;
       *)                  set -- "$@" "$arg"
  esac
done


USER_BED_FILES=()
while getopts ":P:A:G:B:U:h" OPT; do
    case $OPT in
        P) P=$OPTARG;;
        A) A=$OPTARG;;
        G) G=$OPTARG;;
        B) B=$OPTARG;;
        U) USER_BED_FILES+=("$OPTARG");;  # Add each -U argument to the array
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





if [[ -z $P ]];
then
    usage
    exit
fi

if [[ -z $A ]];
then
    usage
    exit
fi


if [[ -z $G ]];
then
    usage
    exit
fi

if [ $G != "hg38" ];
then
    echo "Currently only hg38 build is supported"
    exit
fi

if [[ -z $B ]];
then
    usage
    exit
fi



column_count=$(head -n 1 $P | awk '{print NF}')
if [ "$column_count" -ne 7 ]; then
  echo "Error: bedpe does not have exactly 7 columns."
  exit 1
fi



TSS=$ref_dir/FANTOM_all-alternate-TSSs_hg38.bed # this is made using FANTOM, GENCODE and ENCODE-3
bedtools intersect \
 -a $TSS \
 -b $B > $temp_dir/sample_tss.bed


housekeeping=$ref_dir/FANTOM_housekeeping-alternate-TSSs_human-proteome-atlas_hg38.bed
bedtools intersect \
 -a $housekeeping \
 -b $B > $temp_dir/sample_housekeeping-tss.bed

ENH=$ref_dir/ENCODE3_cCRE-enhancers_hg38.bed
bedtools intersect \
 -a $ENH \
 -b $B > $temp_dir/sample_enh.bed




clusters=$(cut -f 7 $P | sort -V | uniq)


header="Cluster\tchr\tstart\tend\tDegree\tTotal_sum\tWedge_short\tWedge_long\tRange_span\tBin_span\tPeak_span\t\
Num_-B_peaks\tNum_Alternate_TSSs\tNum_lncRNA\tNum_Housekeeping_Genes\tNum_All_Genes(protein_coding)\tNum_ENCODE-3_Enh\tNum_ENCODE-3_CTCF\tNum_UCSC_CpG\t\
LncRNAs\tHousekeeping_Genes\tAll_Genes(protein_coding)\tNum_Enh-Enh_loops\tNum_Enh-Pro_loops\tNum_Pro-Pro_loops"


for user_bed in "${USER_BED_FILES[@]}"; do
    header="${header}\t$(basename $user_bed)"
done

echo -e "$header"


for cluster in $clusters; do

    # subset cluster
    awk -v clust="$cluster" '$7 == clust' $P > $temp_dir/"${cluster}_subset.bedpe"


    # 1. coordinate range
    chr=$(awk '{print $1; print $4}' $temp_dir/"${cluster}_subset.bedpe" | sort | uniq )
    start=$(awk '{print $2; print $3; print $5; print $6}' $temp_dir/"${cluster}_subset.bedpe" | sort -k1,1n | uniq | head -1 )
    end=$(awk '{print $2; print $3; print $5; print $6}' $temp_dir/"${cluster}_subset.bedpe"   | sort -k1,1n | uniq | tail -1 )
    range=$(echo "$chr:$start:$end")


    # 2. degree
    degree=$(wc -l < $temp_dir/"${cluster}_subset.bedpe")



    # 3. query_bedpe
    $aqua_dir/annotate_loops.sh \
     --pair $temp_dir/"${cluster}_subset.bedpe" \
     --sample1 $A \
     --genome $G \
     --resolution 5000 > $temp_dir/"${cluster}_subset_C-annotated.bedpe"

    total_count=$(awk '{sum += $7} END {print sum}' $temp_dir/"${cluster}_subset_C-annotated.bedpe")


    # 4. summarise_interval
    echo -e "$chr\t$start\t$end" > $temp_dir/"${cluster}_subset_range.bed"

    $aqua_dir/summarize_interval.sh \
     --input $temp_dir/"${cluster}_subset_range.bed" \
     --sample1 $A \
     --genome $G \
     --resolution 5000 > $temp_dir/"${cluster}_subset_summary.bed"

    short_count=$(tail -1 $temp_dir/"${cluster}_subset_summary.bed" | cut -f 4)
    long_count=$(tail -1 $temp_dir/"${cluster}_subset_summary.bed"  | cut -f 5)


    # 5. range_span
    range_span=$((end - start))


    # 6. bin_span
    ## cluster bed file
    awk '{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' $temp_dir/"${cluster}_subset.bedpe" | \
    sort -k1,1 -k2,2n | \
    uniq | \
    bedtools merge -i - > $temp_dir/"${cluster}_subset_merged.bed"

    bin_span=$(awk '{print $0"\t"$3-$2}' $temp_dir/"${cluster}_subset_merged.bed" | awk '{sum += $4} END {print sum}' )


    # 7. peak_span
    peak_span=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $B | \
        awk '{print $0"\t"$3-$2}' | \
        awk '{sum += $4} END {if (NR == 0) print 0; else print sum}')



    # 8. Fuctional intersections


    ## 8a. TSSs, lncRNAs & genes
    num_TSS=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $temp_dir/sample_tss.bed \
         -wa -wb | \
        awk '$8=="protein_coding"' | \
        cut -f 4-9 | \
        sort | \
        uniq | \
        wc -l )

    num_genes=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $temp_dir/sample_tss.bed \
         -wa -wb | \
        awk '$8=="protein_coding"' | \
        cut -f 7 | \
        sort | \
        uniq | \
        wc -l )

    num_lncRNA=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $temp_dir/sample_tss.bed \
         -wa -wb | \
        awk '$8=="lncRNA"' | \
        cut -f 7 | \
        sort | \
        uniq | \
        wc -l )

    genes=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $temp_dir/sample_tss.bed \
         -wa -wb | \
        awk '$8=="protein_coding"' | \
        cut -f 7 | \
        sort | \
        uniq | \
        paste -sd, - )

    lncRNAs=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $temp_dir/sample_tss.bed \
         -wa -wb | \
        awk '$8=="lncRNA"' | \
        cut -f 7 | \
        sort | \
        uniq | \
        paste -sd, - )


    ## 8b. CpG islands
    CpG=$ref_dir/UCSC_CpG-islands_hg38.bed

    num_CpGs=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $CpG \
         -wa -wb | \
        cut -f 4-9 | \
        sort | \
        uniq | wc -l )


    ## 8c. Enhancers and peaks
    num_enhs=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $temp_dir/sample_enh.bed -wa -wb | \
        cut -f 4-9 | \
        sort | \
        uniq | wc -l )

    num_peaks=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $B \
         -wa -wb | \
        cut -f 4-6 | \
        sort | \
        uniq | wc -l )


    ## 8d. CTCF sites
    CTCF=$ref_dir/ENCODE3_cCRE-CTCF_hg38.bed

    num_CTCF=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $CTCF -wa -wb | \
        cut -f 4-9 | \
        sort | \
        uniq | wc -l )

    ##8e. Housekeeping genes
    num_housekeeping_genes=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $temp_dir/sample_housekeeping-tss.bed \
         -wa -wb | \
        cut -f 7 | \
        sort | \
        uniq | \
        wc -l )

    housekeeping_genes=$(\
        bedtools intersect \
         -a $temp_dir/"${cluster}_subset_merged.bed" \
         -b $temp_dir/sample_housekeeping-tss.bed \
         -wa -wb | \
        cut -f 7 | \
        sort | \
        uniq | \
        paste -sd, - )


    #9. Loop breakdown
    E_E=0
    E_P=0
    P_P=0

    $aqua_dir/intersect_bedpe.sh \
     --bedpe $temp_dir/"${cluster}_subset.bedpe" \
     --bed_A $temp_dir/sample_tss.bed \
     --bed_B $temp_dir/sample_enh.bed > $temp_dir/"${cluster}_subset-intersect.bedpe"

    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $NF}' $temp_dir/"${cluster}_subset-intersect.bedpe" | \
    sort | uniq | \
    cut -f 7 | \
    sort | uniq -c > $temp_dir/"${cluster}_loop-breakdown.txt"

    while read -r count pair; do
        case $pair in
            "B-B")
                E_E=$((E_E + count))
                ;;
            "A-B"|"B-A")
                E_P=$((E_P + count))
                ;;
            "A-A")
                P_P=$((P_P + count))
                ;;
        esac
    done < <(awk '{print $1, $2}' $temp_dir/"${cluster}_loop-breakdown.txt")



    output="${cluster}\t${chr}\t${start}\t${end}\t${degree}\t${total_count}\t${short_count}\t${long_count}\t${range_span}\t${bin_span}\t${peak_span}\t"\
"${num_peaks}\t${num_TSS}\t${num_lncRNA}\t${num_housekeeping_genes}\t${num_genes}\t${num_enhs}\t${num_CTCF}\t${num_CpGs}\t"\
"${lncRNAs}\t${housekeeping_genes}\t${genes}\t${E_E}\t${E_P}\t${P_P}"


    for user_bed in "${USER_BED_FILES[@]}"; do
            num_user=$(\
                bedtools intersect \
                 -a $temp_dir/"${cluster}_subset_merged.bed" \
                 -b $user_bed \
                 -wa -wb | \
                cut -f 4-6 | \
                sort | \
                uniq | wc -l )
            output="${output}\t${num_user}"
    done


    echo -e "${output}"

done




trap "sudo rm -rf $temp_dir" EXIT
