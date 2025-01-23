# AQuA-Tools

## About

3D genomics projects have rapidly expanded to include dozens of patients, treatment conditions, and multiple layers of information (H3K27ac, H3K4me3, PolII, CTCF, transcription factors…). [Axiotl Inc.](https://axiotl.com) and the [Gryder lab](https://gryderlab.com) at Case Western Reserve University are developing AQuA tools to facilitate the analysis and visualization of increasingly complex 3D data.

The name 'AQuA' derives from **A**bsolute **QU**antification of chromatin **A**rchitecture, an experimental design in which 3D contacts are counted from paired-end tags (PETs) from the human genome and are normalized to the total PETs from the mouse genome, enabling more quantitative insights into the topological determinants of biology in question. The publication can be found [here](https://www.nature.com/articles/s41596-019-0285-9).

>Gryder, B.E., Khan, J. & Stanton, B.Z. Measurement of differential chromatin interactions with absolute quantification of architecture (AQuA-HiChIP). Nat Protoc 15, 1209–1236 (2020). https://doi.org/10.1038/s41596-019-0285-9

For more information, please reach us at hello@axiotl.com.


## Installation and Dependencies


### Docker: 

Sample data and AQuA tools along with all their dependencies can be easily setup via Docker. Before getting started, ensure Docker is intalled on your machine. 

Download the `Dockerfile` available in the Docker folder of this repo and build the image as follows- 

```
# set up working station

working_dir=$HOME/aqua_tools_container  # <-- change this path as per your convenience
mkdir -p $working_dir
cd $working_dir                         # <-- place the Dockerfile from GitHub in this directory


# build the image (this can take upto 30 minutes)
sudo docker build -t aqua_tools:1 .


# run the container and get started!
# once inside the container, please keep all outputs in 
# ~/container_outputs to access them after exiting
sudo docker run -it -v $working_dir:/home/ubuntu/container_outputs aqua_tools:1
```


## Recipes

All AQuA tools can now be used inside the container. Here are some tutorials-

### Getting started / Listing samples:
```

```

### Local loop calling:
```

```

### Genome-wide loop calling:
```

```