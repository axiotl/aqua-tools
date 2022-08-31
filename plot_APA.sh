#!/bin/bash

OPTIND=1

win_size=10


aqua_dir=~/aqua_tools

juicer_tools='java -jar /home/ubuntu/juicer_tools_1.19.02.jar'


function usage {
    echo -e "usage : plot_APA.sh -P PATH_TO_GENOMIC_PAIRS_FILE_YOU_WANT_TO_PLOT -A NAME_OF_FIRST_SAMPLE -G GENOME_BUILD -O OUTPUT_DIRECTORY [-B NAME_OF_SECOND_SAMPLE] [-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo 
    echo "Create APA plot with AQuA normalized contact values"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -P|--pair                    : Path to the bedpe (pairs) file you want to use."
    echo "    -A|--sample1                 : Full path to folder that houses the .hic and mergestats files"
    echo "    -O|--out-dir                 : Full path of the directory you want to store the output plots in"
    echo " [  -B|--sample2             ]   : Full path to folder that houses the .hic and mergestats files. If triggered, plots the delta AQuA normalized values from both samples for that pair. Useful in case vs control"
    echo " [     --cpml                ]   : No input required. If --cpml is specified, CPM and AQuA APA values get normalised by the number of loops in the bedpe"
    echo " [     --bin_size            ]   : Bin size you want to use for the APA plots. If not specified, default 5000 will be used"
    echo " [     --hard_cap_cpm        ]   : Upper limit of the CPM plot range. If not specified, upper limit will be calcualted using max bin value"
    echo " [     --hard_cap_cpm_delta  ]   : Upper limit of the CPM delta plot range. Only for two sample analysis. If not specified, upper limit will be calcualted using max delta value"
    echo " [     --hard_cap_aqua       ]   : Upper limit of the AQuA plot range. If not specified, upper limit will be calcualted using max bin value"
    echo " [     --hard_cap_aqua_delta ]   : Upper limit of the AQuA delta plot range. Only for two sample analysis. If not specified, upper limit will be calcualted using max delta value"
    echo " [  -h|--help                ]   : Help message"
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
      "--pair")                 set -- "$@" "-P" ;;
      "--sample1")              set -- "$@" "-A" ;;
      "--out-dir")              set -- "$@" "-O" ;;
      "--help")                 set -- "$@" "-h" ;;
      "--sample2")              set -- "$@" "-B" ;;
      "--cpml")                 set -- "$@" "-N" ;;
      "--bin_size")             set -- "$@" "-b" ;;
      "--hard_cap_cpm")         set -- "$@" "-c" ;;
      "--hard_cap_cpm_delta")   set -- "$@" "-e" ;;
      "--hard_cap_aqua")        set -- "$@" "-a" ;;
      "--hard_cap_aqua_delta")  set -- "$@" "-d" ;;
       *)                       set -- "$@" "$arg"
  esac
done

b=5000
c="no_cap"
a="no_cap"
e="no_cap"
d="no_cap"
N=0

while getopts ":P:A:B:O:b:c:a:e:d:Nh" OPT
do
    case $OPT in
  P) P=$OPTARG;;
  A) A=$OPTARG;;
  B) B=$OPTARG;;
  O) O=$OPTARG;;
  b) b=$OPTARG;;
  c) c=$OPTARG;;
  a) a=$OPTARG;;
  e) e=$OPTARG;;
  d) d=$OPTARG;;
  N) N=1;;
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



# If no sample B, the variable is empty:
if [ -z "$B" ]; then
  case $B in
  B) B="";;
  esac
fi


# Do all necessary paramter checks
#----------------------------------

if [[ $b == 5000 ]]; then
 echo "Using default bin size 5000"
else
 echo "Using bin size: ${b}" 
fi

#----------------------------------

if [[ -z $P ]];
then
    usage
    exit
fi

#----------------------------------

col1=$(head -n 1  $P | cut -f 1 )
col1="${col1:0:3}"

col4=$(head -n 1  $P | cut -f 4 )
col4="${col4:0:3}"


if [[ $col1 != $col4  ]]; then
echo "Please provide a 6-col bedpe file without headers"
exit
fi

#----------------------------------

if [[ -z $A ]];
then
    usage
    exit
fi

#----------------------------------

if [[ -z $O ]];
then
    usage
    exit
fi

#----------------------------------





###########################################################################
###########################################################################
###                                                                     ###
###                         ONE-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if [ -z "$B" ]
then 



# Do all necessary paramter checks
#----------------------------------

if [[ $c == "no_cap" ]]; then
 echo "No CPM hard cap for plotting specified"
else
 echo "CPM hard cap: ${c}" 
fi

#----------------------------------

if [[ $a == "no_cap" ]]; then
 echo "No AQuA hard cap for plotting specified"
else
 echo "AQuA hard cap: ${a}" 
fi

#----------------------------------

if [[ $e != "no_cap" ]]; then
 echo "Cannot have a delta cap in single-sample APA analysis"
 exit
fi

#----------------------------------

if [[ $d != "no_cap" ]]; then
 echo "Cannot have a delta cap in single-sample APA analysis"
 exit
fi

#----------------------------------

num_loops=`cat "$P" | wc -l`


id=`basename $P`
id_create=$(echo "$id" | cut -f 1 -d '.')

if [[ $id_create == $id ]]; then
  echo "Please add an extension to the input pairs file (a simple .txt would suffice), or consider changing output directory"
  exit
fi


sample1=`basename $A`


out_dir=$O/$id_create
mkdir -p $out_dir
out_dir=$out_dir/$sample1
mkdir -p $out_dir
mkdir -p $out_dir/$sample1



#----------------------------------

pair1=$A/$sample1.hic
if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi

#----------------------------------

stat1=$A/$sample1.mergeStats.txt
if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi

#----------------------------------



hg_total1=`head -3 $stat1 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
mm_total1=`head -3 $stat1 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`


echo "total human reads for $sample1: $hg_total1"


has_aqua=true
total1=$hg_total1


# AQuA Factor
if [[ ! -z $mm_total1 ]];then
  echo "We have spike-In"
  total1=`echo "$hg_total1+$mm_total1" | bc`
  aqua_factor1=`echo "scale=4;$hg_total1/$mm_total1" | bc`
  echo "aqua_factor: $aqua_factor1"
fi


# Normalization factor
norm_factor1=`echo "scale=10;1000000/$total1" | bc`
echo "norm_factor: $norm_factor1"


# Call Juicer APA https://github.com/aidenlab/juicer/wiki/APA
echo "juicer_tools"
$juicer_tools apa --threads 1 -k NONE -n 0 -r $b -w $win_size $pair1 $P $out_dir/$sample1 #&> /dev/null


out_cpm_matrix1=$out_dir/$sample1/$b/gw/APA.cpm.txt
out_aqua_matrix1=$out_dir/$sample1/$b/gw/APA.aqua.txt

# If no --cpml
if [[ $N == 0 ]]; then

echo "norm_factor"
cat $out_dir/${sample1}/${b}/gw/APA.txt | tr -d '[]' | awk -v factor="${norm_factor1}" -F',' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_cpm_matrix1
  
echo "aqua_factor"
cat $out_cpm_matrix1 | awk -v factor="${aqua_factor1}" -F$'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_aqua_matrix1


Rscript $aqua_dir/plot_APA.r $out_cpm_matrix1  ${sample1} ${c} "HiChIP" $out_dir/APA_cpm.pdf   ${win_size} ${b} ${P}
Rscript $aqua_dir/plot_APA.r $out_aqua_matrix1 ${sample1} ${a} "AQuA"   $out_dir/APA_aqua.pdf  ${win_size} ${b} ${P}

fi

# If --cpml
if [[ $N == 1 ]]; then



out_aqua_cpml_matrix1=$out_dir/$A/$b/gw/APA.aqua.cpml.txt

echo "norm_factor"
cat $out_dir/${sample1}/${b}/gw/APA.txt | tr -d '[]' | awk -v factor="${norm_factor1}" -F',' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_cpm_matrix1
  
echo "aqua_factor"
cat $out_cpm_matrix1 | awk -v factor="${aqua_factor1}" -F$'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_aqua_matrix1

echo "cpml"
cat $out_aqua_matrix1 | awk -v factor="${num_loops}" -F$'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)/=factor; print}' > $out_aqua_cpml_matrix1


Rscript $aqua_dir/plot_APA.r $out_cpm_matrix1  ${sample1} ${c} "HiChIP"  $out_dir/APA_cpm.pdf  ${win_size} ${b} ${P}
Rscript $aqua_dir/plot_APA.r $out_aqua_cpml_matrix1 ${sample1} ${c} "AQuA"  $out_dir/APA_aqua.pdf  ${win_size} ${b} ${P} 

fi



# end one-sample loop
fi





###########################################################################
###########################################################################
###                                                                     ###
###                         TWO-SAMPLE ANALYSIS                         ###
###                                                                     ###
###########################################################################
###########################################################################


if [ -n "$B" ]
then 



# Do all necessary paramter checks
#----------------------------------

if [[ $c == "no_cap" ]]; then
 echo "No CPM hard cap for plotting specified"
else
 echo "CPM hard cap: ${c}" 
fi

#----------------------------------

if [[ $a == "no_cap" ]]; then
 echo "No AQuA hard cap for plotting specified"
else
 echo "AQuA hard cap: ${a}" 
fi

#----------------------------------

if [[ $e == "no_cap" ]]; then
 echo "No CPM delta hard cap for plotting specified"
else
 echo "CPM delta hard cap: ${e}" 
fi

#----------------------------------

if [[ $d == "no_cap" ]]; then
 echo "No AQuA delta hard cap for plotting specified"
else
 echo "AQuA delta hard cap: ${d}" 
fi

#----------------------------------

num_loops=`cat "$P" | wc -l`


id=`basename $P`
id_create=$(echo "$id" | cut -f 1 -d '.')

if [[ $id_create == $id ]]; then
  echo "Please add an extension to the input pairs file, or consider changing output directory"
  exit
fi


sample1=`basename $A`
sample2=`basename $B`

out_dir=$O/$id_create
mkdir -p $out_dir
out_dir=$out_dir/${sample2}__vs__${sample1}
mkdir -p $out_dir
mkdir -p $out_dir/${sample2}
mkdir -p $out_dir/${sample1}


#----------------------------------

pair1=$A/$sample1.hic
pair2=$B/$sample2.hic
if [[ ! -f $pair1 ]]; then echo "cannot find $pair1"; exit 1; fi
if [[ ! -f $pair2 ]]; then echo "cannot find $pair2"; exit 1; fi

#----------------------------------

stat1=$A/$sample1.mergeStats.txt
stat2=$B/$sample2.mergeStats.txt
if [[ ! -f $stat1 ]]; then echo "cannot find $stat1"; exit 1; fi
if [[ ! -f $stat2 ]]; then echo "cannot find $stat2"; exit 1; fi

#----------------------------------



hg_total1=`head -3 $stat1 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
hg_total2=`head -3 $stat2 | tail -1 | cut -f2 | perl -nle 's/\r//g; print;'`
mm_total1=`head -3 $stat1 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`
mm_total2=`head -3 $stat2 | tail -1 | cut -f3 | perl -nle 's/\r//g; print;'`

echo "total human reads for $A: $hg_total1"
echo "total human reads for $B: $hg_total2"
echo "total mouse reads for $A: $mm_total1"
echo "total mouse reads for $B: $mm_total2"

has_aqua=true
total1=$hg_total1
total2=$hg_total2

if [[ ! -z $mm_total1 && ! -z $mm_total2 ]];then
  echo "We have spike-In"
  total1=`echo "$hg_total1+$mm_total1" | bc`
  total2=`echo "$hg_total2+$mm_total2" | bc`
  aqua_factor1=`echo "scale=4;$hg_total1/$mm_total1" | bc`
  aqua_factor2=`echo "scale=4;$hg_total2/$mm_total2" | bc`
  echo "aqua_factor1: $aqua_factor1"
  echo "aqua_factor2: $aqua_factor2"
fi

norm_factor1=`echo "scale=10;1000000/$total1" | bc`
norm_factor2=`echo "scale=10;1000000/$total2" | bc`
echo "norm_factor1: $norm_factor1"
echo "norm_factor2: $norm_factor2"

# https://github.com/aidenlab/juicer/wiki/APA

echo "juicer_tools"
$juicer_tools apa --threads 1 -k NONE -n 0 -r $b -w $win_size $pair1 $P $out_dir/$sample1 #&> /dev/null
$juicer_tools apa --threads 1 -k NONE -n 0 -r $b -w $win_size $pair2 $P $out_dir/$sample2 #&> /dev/null


out_cpm_matrix1=$out_dir/$sample1/$b/gw/APA.cpm.txt
out_cpm_matrix2=$out_dir/$sample2/$b/gw/APA.cpm.txt
out_aqua_matrix1=$out_dir/$sample1/$b/gw/APA.aqua.txt
out_aqua_matrix2=$out_dir/$sample2/$b/gw/APA.aqua.txt

# If no --cpml
if [[ $N == 0 ]]; then

echo "norm_factor"
cat $out_dir/${sample1}/${b}/gw/APA.txt | tr -d '[]' | awk -v factor="${norm_factor1}" -F',' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_cpm_matrix1
cat $out_dir/${sample2}/${b}/gw/APA.txt | tr -d '[]' | awk -v factor="${norm_factor2}" -F',' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_cpm_matrix2
  
echo "aqua_factor"
cat $out_cpm_matrix1 | awk -v factor="${aqua_factor1}" -F$'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_aqua_matrix1
cat $out_cpm_matrix2 | awk -v factor="${aqua_factor2}" -F$'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_aqua_matrix2

Rscript $aqua_dir/plot_APA.r $out_cpm_matrix1  $out_cpm_matrix2  ${sample1} ${sample2} ${c} ${e} "HiChIP"  $out_dir/APA_cpm.pdf ${win_size} ${b} ${P}
Rscript $aqua_dir/plot_APA.r $out_aqua_matrix1 $out_aqua_matrix2 ${sample1} ${sample2} ${a} ${d} "AQuA"    $out_dir/APA_aqua.pdf ${win_size} ${b} ${P}

fi

# If --cpml
if [[ $N == 1 ]]; then

out_aqua_cpml_matrix1=$out_dir/$sample1/$b/gw/APA.aqua.cpml.txt
out_aqua_cpml_matrix2=$out_dir/$sample2/$b/gw/APA.aqua.cpml.txt

echo "norm_factor"
cat $out_dir/${sample1}/${b}/gw/APA.txt | tr -d '[]' | awk -v factor="${norm_factor1}" -F',' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_cpm_matrix1
cat $out_dir/${sample2}/${b}/gw/APA.txt | tr -d '[]' | awk -v factor="${norm_factor2}" -F',' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_cpm_matrix2
  
echo "aqua_factor"
cat $out_cpm_matrix1 | awk -v factor="${aqua_factor1}" -F$'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_aqua_matrix1
cat $out_cpm_matrix2 | awk -v factor="${aqua_factor2}" -F$'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)*=factor; print}' > $out_aqua_matrix2

echo "cpml"
cat $out_aqua_matrix1 | awk -v factor="${num_loops}" -F$'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)/=factor; print}' > $out_aqua_cpml_matrix1
cat $out_aqua_matrix2 | awk -v factor="${num_loops}" -F$'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $(i)/=factor; print}' > $out_aqua_cpml_matrix2


Rscript $aqua_dir/plot_APA.r $out_cpm_matrix1  $out_cpm_matrix2  ${sample1} ${sample2} ${c} ${e} "HiChIP"  $out_dir/APA_cpm.pdf ${win_size} ${b} ${P}
Rscript $aqua_dir/plot_APA.r $out_aqua_cpml_matrix1 $out_aqua_cpml_matrix2 ${sample1} ${sample2} ${a} ${d} "AQuA"  $out_dir/APA_aqua.pdf ${win_size} ${b} ${P}

fi




fi
