# AQuA-Tools

## About

3D genomics projects have rapidly expanded to include dozens of patients, treatment conditions, and multiple layers of information (H3K27ac, H3K4me3, PolII, CTCF, transcription factors…). [Axiotl Inc.](https://axiotl.com) and the [Gryder lab](https://gryderlab.com) at Case Western Reserve University are developing AQuA tools to facilitate the analysis and visualization of increasingly complex 3D data.

The name 'AQuA' derives from **A**bsolute **QU**antification of chromatin **A**rchitecture, an experimental design in which 3D contacts are counted from paired-end tags (PETs) from the human genome and are normalized to the total PETs from the mouse genome, enabling more quantitative insights into the topological determinants of biology in question. The publication can be found [here](https://www.nature.com/articles/s41596-019-0285-9).

>Gryder, B.E., Khan, J. & Stanton, B.Z. Measurement of differential chromatin interactions with absolute quantification of architecture (AQuA-HiChIP). Nat Protoc 15, 1209–1236 (2020). https://doi.org/10.1038/s41596-019-0285-9

For more information, please reach us at hello@axiotl.com.


## Installation and Dependencies


### Juicer: 

A working Java installation (version >= 1.8) on Windows, Linux, and Mac OSX is required. Quoting the Aiden lab-
>We recommend using the latest Java version available, but please do not use the Java Beta Version. Minimum system requirements for running Java can be found at https://java.com/en/download/help/sysreq.xml.
>To download and install the latest Java Runtime Environment (JRE), please go to https://www.java.com/download. 


1. Download Juicer Tools
```
cd
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.19.02.jar
```
2. Edit terminal start script (.bashrc/.zshrc) and add the following line
```
export juicer_tools='java -jar $HOME/juicer_tools_1.19.02.jar'
```
3. Restart terminal


### R: 

Install the following R packages- 

1. strawr
2. parallel
3. ggplot2
4. gridExtra
5. pheatmap


### AQuA Tools:

1. Download and unpack AQuA-Tools
```
cd
wget https://github.com/axiotl/aqua-tools/archive/refs/heads/main.zip
unzip main.zip
mv aqua-tools-main aqua_tools
cd aqua_tools
chmod +x *
```
2. Edit terminal start script (.bashrc/.zshrc) and add the following line
```
export PATH="$HOME/aqua_tools:$PATH"
```
3. Restart terminal


## Data and Organization

We use `.hic` files as the only format of data source. In addition, we need a `.mergeStats.txt` file containing valid interaction read counts for human and mouse genome alignments. This is a modified version of a typical `.mergestat` file (an output of the [HiC-Pro](https://github.com/nservant/HiC-Pro) pipeline) that contains QC information for human and mouse PETs in one file. A `.mergeStats.txt` file should look like the following-

|  | GM12878.hg38 | GM12878.mm10 |
| ---------- | -------- | ------ |
| valid_interaction       |  377103997    |  7165030  |
| valid_interaction_rmdup       |  274352694   | 5205894 |
| trans_interaction       |  80252048   |  708100   |
| cis_interaction       |   194100646   |  4497794  |
| cis-shortRange       |  170608942   |  3506516  |
| cis_longRange       |  23491704   |  991278  |

The two files above should be housed in a folder named after the sample, and should strictly be names as follows- 

```
   + ~/SAMPLE
       ++ SAMPLE.hic
       ++ SAMPLE.mergeStats.txt
```
If sample name is `GM12878`, the data structure should look like-
```
   + ~/GM12878
       ++ GM12878.hic
       ++ GM12878.mergeStats.txt
```

Please contact [Axiotl Inc.](https://axiotl.com) for any further questions.


## Usage

Once installed, all AQuA tools should be executable from anywhere and can be used by typing `-h` or `--help` followed by the name of the tool. For example-
```
build_bedpe.sh -h 


Builds pairs between elements in two bed files.
Pairs can be constrained by a third bed file (usually TADs)
or by a minimum distance between them.
Prints pairs in bedpe format to standard out.

Pairs can be use to query .hic files with annotate_loops.sh
---------------
OPTIONS

   -A|--bed_A      PATH_TO_FIRST_BED_FILE    : Path of the first bed file you want to build the bedpe of
   -B|--bed_B      PATH_TO_SECOND_BED_FILE   : Path of the second bed file you want to build the bedpe of
  [-T|--TAD      ] PATH_TO_TAD_FILE          : Path of the TAD file that checks if genomic regions fall within them
  [-d|--min_dist ] MIN_DROP_DISTANCE         : minimum distance between pairs used to drop results. Default 0 bp
  [-D|--max_dist ] MAX_DROP_DISTANCE         : maximum distance between pairs used to drop results. Default 5 Mb
  [-h|--help     ] Help message

```
