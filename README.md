# AQuA-Tools

[![License](https://img.shields.io/github/license/yourusername/aqua-tools)](LICENSE)

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

## Availability
tool chart and description

## Use
AQuA-Tools can be used through any of the following entrypoints:
- [Tinker Cloud Platform](#tinker)
- [Docker Container](#docker)
- [Local Installation](#local)

### Tinker
[Tinker](https://tinker.axiotl.com/public) is our cloud platform specifically designed for 3D genomics analyses. It offers:

- **Immediate Access**: No installation required and immeditely accessible via browser
- **Pre-loaded Data**: Access to a collection of reference datasets and invariant genome annotations
- **Access to all tools**: Access to all tools, including `extract_bedpe` and parameter `--inherent` for `plot_contacts`, `query_bedpe` and `plot_virtual_4C`
- **Optimized processing**: Faster outputs for large `build_bedpe` calls

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
- Easy updates and version management
- Portable analysis environment

Similar to [Tinker](https://tinker.axiotl.com/public), the container built using the image will come with:
- **Pre-loaded Data**: Access to common reference datasets and invariant genome annotations
- **Access to all tools**

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
For users who need maximum flexibility or have specific system requirements, AQuA-Tools can be installed directly on your local machine.

**Prerequisites:**
- Python 3.8 or higher
- R 4.0 or higher
- gcc/g++ compiler
- Development libraries for HDF5 and zlib

**Installation Steps:**

1. Install system dependencies:
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y \
    python3-dev \
    r-base \
    libhdf5-dev \
    zlib1g-dev \
    build-essential
```

2.

## Recipes
All AQuA tools are executable from anywhere inside the container. The container mimics [Tinker](https://tinker.axiotl.com/), a cloud platform built from ground up for 3D genomics analyses.

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
We use `extract_bedpe` to call loops using a sample of interest. The tool supports both targeted interval analysis and genome-wide loop calling.

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




## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

