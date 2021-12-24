FROM ubuntu:20.04
MAINTAINER Karthik G <gkarthik@scripps.edu>


RUN apt update
RUN apt install -y wget git parallel python3 python3-pip r-base

RUN DEBIAN_FRONTEND=noninteractive \
    apt install -y default-jre-headless &&\
    apt clean

RUN pip install --no-input numpy scipy matplotlib pandas biopython baltic tqdm unidecode

RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.490-1_amd64.deb && dpkg -i mafft_7.490-1_amd64.deb

RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.1.3/iqtree-2.1.3-Linux.tar.gz && \
    tar xf iqtree-2.1.3-Linux.tar.gz

ENV PATH="/iqtree-2.1.3-Linux/bin:${PATH}"

RUN mkdir nextflow &&\
    cd nextflow &&\
    wget -qO- https://get.nextflow.io | bash &&\
    chmod +x nextflow

ENV PATH="/nextflow:${PATH}"

RUN Rscript -e 'install.packages("ape", repos="https://cloud.r-project.org")'


# RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
#     /bin/bash Mambaforge-Linux-x86_64.sh -b &&\
#     eval "$(/root/mambaforge/bin/conda shell.bash hook)" &&\
#     conda init &&\
#     conda config --set auto_activate_base false &&\
#     mamba create --quiet -c conda-forge -c bioconda -n snakemake snakemake &&\
#     conda activate snakemake &&\
    