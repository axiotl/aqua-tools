# AQuA-Tools

## About

3D genomics projects have rapidly expanded to include dozens of patients, treatment conditions, and multiple layers of information (H3K27ac, H3K4me3, PolII, CTCF, transcription factors…). [Axiotl Inc.](https://axiotl.com) and the [Gryder lab](https://gryderlab.com) at Case Western Reserve University are developing AQuA tools to facilitate the analysis and visualization of increasingly complex 3D data.

The name 'AQuA' derives from **A**bsolute **QU**antification of chromatin **A**rchitecture, an experimental design in which 3D contacts are counted from paired-end tags (PETs) from the human genome and are normalized to the total PETs from the mouse genome, enabling more quantitative insights into the topological determinants of biology in question. The publication can be found [here](https://www.nature.com/articles/s41596-019-0285-9).

>Gryder, B.E., Khan, J. & Stanton, B.Z. Measurement of differential chromatin interactions with absolute quantification of architecture (AQuA-HiChIP). Nat Protoc 15, 1209–1236 (2020). https://doi.org/10.1038/s41596-019-0285-9


## Installation and Dependencies

### Docker Container: 
We highly recommend using Docker for an easy and cohesive experience to use AQuA tools.

Sample data and AQuA tools along with all their dependencies can be easily setup via Docker. Before getting started, ensure Docker Desktop is installed on your machine. 

Download the `Dockerfile` available in the Docker folder of this repo and build the image as follows- 


i. Set up working station
```
working_dir=$HOME/aqua_tools_container  # <-- change this path as per your convenience
mkdir -p $working_dir
cd $working_dir                         # <-- place the Dockerfile from GitHub in this directory
```

ii. Build the image (this can take upto 30 minutes)
```
sudo docker build -t aqua_tools .
```

iii. Run the container and get started!
(once inside the container, keep all outputs in `~/container_outputs` to access them after exiting)
```
sudo docker run -it -v $working_dir:/home/ubuntu/container_outputs aqua_tools
```


## Recipes

All AQuA tools are executable from anywhere inside the container. The container mimics [Tinker](https://tinker.axiotl.com/), a cloud platform built from ground up for 3D genomics analyses. 

Interested in getting your data to Tinker? Contact us at hello@axiotl.com

### Getting started:

The container comes preloaded with publically available H3K27ac HiChIP samples. 
This can be viewed using-
```
list_samples
```

All AQuA tools can be listed using-
```
list_tools
```

### I. Loop calling:

We use `extract_bedpe` to call loops using a sample of interest. 
The tool expects a genomic interval to call loops, or a TAD file to call loops genome wide.
All outputs are printed in standard out in `.bedpe` format

#### Loop calling for a range/interval:
```
sample=K562_H3K27ac
genome=hg38
range=chr8:127000000:130000000      # <-- MYC locus
output_dir=/home/ubuntu/container_outputs

extract_bedpe \
 --sample1 $sample \
 --genome $genome \
 --range $range > $output_dir/MYC.bedpe
```

#### Loop calling genome-wide:
```
sample=K562_H3K27ac
genome=hg38
output_dir=/home/ubuntu/container_outputs
TAD_file=/home/ubuntu/lab-data/hg38/reference/TAD_goldsorted_span_centromeres-removed_hg38.bed

extract_bedpe \
 --sample1 $sample \
 --genome $genome \
 --TAD $TAD_file > $output_dir/K562_genome-wide-loops.bedpe
```

### II. Visualization:

Seeing is believing!
We can easily visualize any locus of interest of a sample using `plot_contacts`, 
and overlay the called loops to assess our loop calls
```
sample=K562_H3K27ac
genome=hg38
range=chr8:127000000:130000000
output_dir=/home/ubuntu/container_outputs

plot_contacts \
 --sample1 $sample \
 --genome $genome \
 --range $range \
 --bedpe $output_dir/MYC.bedpe \
 --output_name $output_dir/MYC.pdf
```

The `.bedpe` and `.pdf` files will now be available in your local machine in `$working_dir`
