# Use an official Ubuntu base image
FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && apt install -y \
    openjdk-8-jdk \
    openjdk-8-jre \
    software-properties-common \
    lsb-release \
    build-essential \
    pkg-config \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    wget \
    curl \
    vim \
    less \
    nano \
    tar \
    unzip \
    bedtools \
    bc \
    gnupg \
    git \
    ca-certificates && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" > /etc/apt/sources.list.d/cran.list && \
    apt update && apt install -y r-base r-base-dev && \
    apt clean && rm -rf /var/lib/apt/lists/*

# Install R and CRAN packages
RUN Rscript --slave --no-save --no-restore-history -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages( \
  c("data.table","ggplot2","dplyr","tidyr","remotes","parallel","pheatmap","gridExtra","RCurl","dbscan","igraph") \
)'

# Install strawr 
RUN Rscript --slave --no-save --no-restore-history -e 'options(repos=c(CRAN="https://cloud.r-project.org")); if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes"); remotes::install_github("aidenlab/straw/R")'

# Install Bioconductor packages
RUN Rscript --slave --no-save --no-restore-history -e 'options(repos=c(CRAN="https://cloud.r-project.org")); if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("GenomicRanges","HiCcompare","InteractionSet","S4Vectors"))'

# Create ubuntu user 
RUN id -u ubuntu &>/dev/null || useradd -m -d /home/ubuntu -s /bin/bash ubuntu
WORKDIR /home/ubuntu

# Download sample data (Always fetch the latest)
RUN wget -O latest_aqua_tools.tar.gz "https://convergence-beta-public.s3.amazonaws.com/aqua_tools_publication.tar.gz?cache-bust=$(date +%s)" && \
    tar -xzvf latest_aqua_tools.tar.gz && \
    mv aqua_tools_publication/* . && \
    rm -r aqua_tools_publication latest_aqua_tools.tar.gz

# Install AQuA tools from GitHub (Always fetch the latest)
RUN wget -O latest_aqua_tools.zip "https://github.com/axiotl/aqua-tools/archive/refs/heads/main.zip?cache-bust=$(date +%s)" && \
    unzip latest_aqua_tools.zip && \
    mv aqua-tools-main aqua_tools && \
    rm latest_aqua_tools.zip

# Make tool scripts executable and symlink them into PATH
WORKDIR /home/ubuntu/aqua_tools
RUN chmod +x *.sh && \
    for f in *.sh; do \
      # link without .sh (preferred command names)
      ln -sfn "/home/ubuntu/aqua_tools/$f" "/usr/local/bin/${f%.sh}"; \
      # link with .sh
      ln -sfn "/home/ubuntu/aqua_tools/$f" "/usr/local/bin/$f"; \
    done

# Switch to ubuntu user
USER ubuntu
WORKDIR /home/ubuntu

# Add environment variables to .bashrc
RUN echo 'export juicer_tools="java -jar $HOME/juicer_tools_1.19.02.jar"' >> /home/ubuntu/.bashrc

# Output directory
RUN mkdir -p container_outputs

# Default command
CMD ["/bin/bash"]
