# AQuA-Tools

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

## Table of Contents
- [About](#about)
- [Use](#Use)
- [Installation](#installation)
- [Usage Recipes](#recipes)
  - [Getting Started](#getting-started)
  - [Loop Calling](#i-loop-calling)
  - [Visualization](#ii-visualization)
- [License](#license)


## About
3D genomics projects have rapidly expanded to include dozens of patients, treatment conditions, and multiple layers of information (H3K27ac, H3K4me3, PolII, CTCF, transcription factors…). [Axiotl Inc.](https://axiotl.com) and the [Gryder lab](https://gryderlab.com) at Case Western Reserve University are developing AQuA tools to facilitate the analysis and visualization of increasingly complex 3D data.

### Publication
>Gryder, B.E., Khan, J. & Stanton, B.Z. Measurement of differential chromatin interactions with absolute quantification of architecture (AQuA-HiChIP). Nat Protoc 15, 1209–1236 (2020). https://doi.org/10.1038/s41596-019-0285-9


## Use and Availability
AQuA-Tools can be used through any of the following entrypoints:
- [Tinker Cloud Platform](#tinker)
- [Docker Container](#docker)
- [Local Installation](#local)

AQuA tools `annotate_clusters`, `plot_contacts` and `extract_bedpe` will be incompatible in local installation, along with the use of the parameter `--inherent` for any tool.

<img src="https://github.com/user-attachments/assets/ba022609-b24b-4f9a-a121-3f99eef21bb3" width="750" height="550">

### Tinker
[Tinker](https://tinker.axiotl.com/public) is our cloud platform specifically designed for 3D genomics analyses. It offers:

**Key Features:**
- No installation required and immeditely accessible via browser
- Access to a collection of reference datasets and invariant genome annotations
- Access to all tools, including `annotate_clusters`, `extract_bedpe`, `plot_contacts` and parameter `--inherent` for `query_bedpe` and `plot_virtual_4C`
- Faster outputs for large `build_bedpe` calls

The free tier of Tinker is capped at 40 total hours of machine use time, following which users can request access to Tinker Pro which additionally comes with:
- **Collaborative Features**: Share files and analyses with team members
- **Scalable Computing**: Handle large datasets efficiently
- **Unlimited access**: No time cap for Tinker use time
- **Pre-built notebooks**: Ready to use notebooks shoehorned to 3D genomic research questions

To get started with Tinker, visit [tinker](https://tinker.axiotl.com/public)


### Docker
To use AQuA tools locally, we highly recommend using Docker for an easy and cohesive experience.

**Key Features:**
- Isolated environment with all dependencies
- Regular updates with addition of publically available HiChIPs

Similar to [Tinker](https://tinker.axiotl.com/public), the container built using the image will come with:
- Access to a collection of reference datasets and invariant genome annotations
- Access to all tools, including `annotate_clusters`, `extract_bedpe`, `plot_contacts` and parameter `--inherent` for `query_bedpe` and `plot_virtual_4C`

**Prerequisites:**
- **Docker Desktop**
- **System Requirements**:
  - CPU: 4+ cores recommended
  - RAM: 8GB or more recommended
  - Storage: 6GB free space for Docker image and data
- **Operating Systems**:
  - Linux (Ubuntu 18.04+, CentOS 7+)
  - macOS (10.15+)

#### Step 1: Set up working station
```bash
working_dir=$HOME/aqua_tools_container  # <-- change this path as per your convenience
mkdir -p $working_dir
cd $working_dir                         # <-- place the Dockerfile from GitHub in this directory
```

#### Step 2: Build the image
```bash
sudo docker build -t aqua_tools .       # This can take up to 30 minutes
```

#### Step 3: Run the container
```bash
sudo docker run -it -v $working_dir:/home/ubuntu/container_outputs aqua_tools
```

**Note**: Keep all outputs in `~/container_outputs` to access them after exiting the container.

### Local 
A subset of AQuA Tools can work directly on your local machine. 

AQuA tools require the following shell dependencies-
```bash
apt update && apt install -y \
    openjdk-8-jdk \
    openjdk-8-jre \
    software-properties-common \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    bedtools \
    bc \
    gnupg
```

AQuA tools require the following R dependencies-
```r
Rscript --slave --no-save --no-restore-history -e '
    install.packages(
        c(
            "data.table",
            "ggplot2",
            "dplyr",
            "tidyr",
            "remotes",
            "parallel",
            "pheatmap",
            "gridExtra",
            "RCurl",
            "dbscan",
            "igraph"
        ),
        repos = "http://www.freestatistics.org/cran/"
    )
'

Rscript --slave --no-save --no-restore-history -e '
    if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes", repos = "http://www.freestatistics.org/cran/")
    }
    remotes::install_github("aidenlab/straw/R")
'

Rscript --slave --no-save --no-restore-history -e '
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "http://www.freestatistics.org/cran/")
    }
    BiocManager::install(c(
        "GenomicRanges",
        "HiCcompare",
        "InteractionSet",
        "S4Vectors"
    ))
'
```

## Recipes
All AQuA tools are executable from anywhere inside the container. The container mimics [Tinker](https://tinker.axiotl.com/), a cloud platform built from ground up for 3D genomics analyses.
Detailed material on how to use the tools can be found at our [Docs](https://docs.axiotl.com/).

### Getting Started
View available samples:
```bash
list_samples
```

List available tools:
```bash
list_tools
```

### I. Loop Calling
We use `extract_bedpe` to call loops using a sample of interest. The tool supports both regional and genome-wide loop calling.

#### Loop calling for a range/interval:
```bash
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
```bash
sample=K562_H3K27ac
genome=hg38
output_dir=/home/ubuntu/container_outputs
TAD_file=/home/ubuntu/lab-data/hg38/reference/TAD_goldsorted_span_centromeres-removed_hg38.bed
extract_bedpe \
 --sample1 $sample \
 --genome $genome \
 --TAD $TAD_file > $output_dir/K562_genome-wide-loops.bedpe
```

### II. Visualization
Generate visual representations of genomic regions with loop overlay:

```bash
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

### III. Differential loops

#### Pre-loaded samples
```
# [type list_samples to view sample names]
sample1=RH4_DMSO_H3K27ac
sample2=RH4_HDACi_H3K27ac
genome_build=hg38
resolution=5000

# invariant reference TAD file 
TAD_file=~/lab-data/hg38/reference/TAD_goldsorted_span_centromeres-removed_hg38.bed
```
#### Call genome-wide loops for each sample
```
extract_bedpe \
 --sample1 $sample1 \
 --genome $genome_build \
 --TAD $TAD_file \
 --resolution $resolution > "$sample1"_gw-loops_"$genome_build".bedpe

extract_bedpe \
 --sample1 $sample2 \
 --genome $genome_build \
 --TAD $TAD_file \
 --resolution $resolution > "$sample2"_gw-loops_"$genome_build".bedpe
```
#### Merge the called genome-wide loops into one file
```
union_bedpe \
 --bedpe "$sample1"_gw-loops_"$genome_build".bedpe \
 --bedpe "$sample2"_gw-loops_"$genome_build".bedpe > union_scaffold.bedpe
```
#### Cluster the unioned scaffold to obtain bedpe membership
```
cluster_bedpe \
 --bedpe union_scaffold.bedpe > union_scaffold_clustered.bedpe
```
#### Annotate union scaffold bedpe with contact values
```
query_bedpe \
 --bedpe union_scaffold_clustered.bedpe \
 --sample1 $sample1 \
 --sample2 $sample2 \
 --genome $genome_build \
 --resolution $resolution \
 --formula max \
 --inherent TRUE > union_scaffold_inh-annotated.bedpe
```
#### Threshold with standard inherent score 1 
```
awk '$NF>= 1{print $0}' union_scaffold_inh-annotated.bedpe > loop-gains.bedpe
awk '$NF<=-1{print $0}' union_scaffold_inh-annotated.bedpe > loop-losses.bedpe
```
#### Visualise gains and losses using APA plotting
```
plot_APA \
 --bedpe loop-gains.bedpe \
 --sample1 $sample1 \
 --sample2 $sample2 \
 --out-dir $(pwd)

plot_APA \
 --bedpe loop-losses.bedpe \
 --sample1 $sample1 \
 --sample2 $sample2 \
 --out-dir $(pwd)
```

## License
AQuA tools is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License v3.0.
The full license text can be found in the [LICENSE.md](LICENSE.md) file.

### Note for Commercial Users
For commercial licensing options, please contact hello@axiotl.com.
