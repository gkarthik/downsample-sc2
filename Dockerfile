FROM ubuntu:20.04
MAINTAINER Karthik G <gkarthik@scripps.edu>


RUN apt update &&\
    DEBIAN_FRONTEND=noninteractive \
    apt install -y wget git parallel python3 python3-pip r-base default-jre-headless default-jdk ant git locales &&\
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

RUN wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.5pre_thorney_v0.1.1/BEASTGen_v0.3pre_thorney.tgz &&\
    tar xf BEASTGen_v0.3pre_thorney.tgz

ENV PATH="/BEASTGen_v0.3pre_thorney/bin:${PATH}"

RUN wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.5pre1/BEASTv1.10.5pre.tgz && \
tar xf BEASTv1.10.5pre.tgz && \
mv BEASTv1.10.5pre/bin/* /usr/local/bin && \
mv BEASTv1.10.5pre/lib/* /usr/local/lib && \
rm BEASTv1.10.5pre.tgz

# Set locale for ant
ENV LC_ALL "en_US.UTF-8"
ENV LC_CTYPE "en_US.UTF-8"

RUN locale-gen en_US.utf8

RUN git clone https://github.com/beast-dev/beast-mcmc.git &&\
    cd beast-mcmc &&\
    ant

# RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
#     /bin/bash Mambaforge-Linux-x86_64.sh -b &&\
#     eval "$(/root/mambaforge/bin/conda shell.bash hook)" &&\
#     conda init &&\
#     conda config --set auto_activate_base false &&\
#     mamba create --quiet -c conda-forge -c bioconda -n snakemake snakemake &&\
#     conda activate snakemake &&\
    