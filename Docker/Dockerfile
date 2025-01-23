FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

# Update Ubuntu packages
RUN apt update && apt install -y \
    software-properties-common \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    wget \
    curl \
    vim \
    nano \
    tar \
    bedtools \
    gnupg && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" > /etc/apt/sources.list.d/cran.list && \
    apt update && apt install -y r-base && \
    apt clean && rm -rf /var/lib/apt/lists/*


# Install R and CRAN packages
RUN Rscript --slave --no-save --no-restore-history -e 'install.packages(\
  c("data.table", "ggplot2", "dplyr", "tidyr", "remotes", "parallel", \
    "pheatmap", "gridExtra", "RCurl", "dbscan", "igraph"), \
  repos = "http://www.freestatistics.org/cran/"\
)'

# Install strawr 
RUN Rscript --slave --no-save --no-restore-history -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "http://www.freestatistics.org/cran/"); remotes::install_github("aidenlab/straw/R")'

# Install Bioconductor packages
RUN Rscript --slave --no-save --no-restore-history -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "http://www.freestatistics.org/cran/"); BiocManager::install(c("GenomicRanges", "HiCcompare", "InteractionSet", "S4Vectors"))'


RUN id -u ubuntu &>/dev/null || useradd -m -d /home/ubuntu -s /bin/bash ubuntu

WORKDIR /home/ubuntu

USER ubuntu

# Add environment variables to .bashrc
RUN echo "export juicer_tools='java -jar $HOME/juicer_tools_1.19.02.jar'" >> /home/ubuntu/.bashrc && \
    echo "export PATH=\"$HOME/aqua_tools:$PATH\"" >> /home/ubuntu/.bashrc && \
    chown ubuntu:ubuntu /home/ubuntu/.bashrc

# AQuA tools and sample data
RUN wget https://convergence-beta-public.s3.amazonaws.com/aqua_tools_publication.tar.gz && \
    tar -xzvf aqua_tools_publication.tar.gz && \
    mv aqua_tools_publication/* . && \
    rm -r aqua_tools_publication aqua_tools_publication.tar.gz

RUN echo "alias query_bedpe='/home/ubuntu/aqua_tools/annotate_loops.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias annotate_cluster='/home/ubuntu/aqua_tools/annotate_cluster.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias build_bedpe='/home/ubuntu/aqua_tools/build_bedpe.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias formula_bedpe='/home/ubuntu/aqua_tools/formula_bedpe.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias extract_bedpe='/home/ubuntu/aqua_tools/extract_bedpe.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias cluster_bedpe='/home/ubuntu/aqua_tools/cluster_bedpe.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias call_clusters='/home/ubuntu/aqua_tools/call_clusters.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias union_bedpe='/home/ubuntu/aqua_tools/union_bedpe.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias intersect_bedpe='/home/ubuntu/aqua_tools/intersect_bedpe.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias list_samples='/home/ubuntu/aqua_tools/list_samples.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias list_tools='/home/ubuntu/aqua_tools/list_tools.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias plot_APA='/home/ubuntu/aqua_tools/plot_APA.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias plot_contacts='/home/ubuntu/aqua_tools/plot_contacts.sh'" >> /home/ubuntu/.bashrc && \
    echo "alias plot_virtual_4C='/home/ubuntu/aqua_tools/plot_virtual_4C.sh'" >> /home/ubuntu/.bashrc


RUN mkdir -p container_outputs


# Set a default command to run an interactive shell
CMD ["/bin/bash"]
