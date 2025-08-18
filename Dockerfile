FROM registry.servicemgr.gendow:5000/huql/base:v1
LABEL maintainer="huql"
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/usr/local/bin:/usr/bin:/usr/local/python3/bin:${PATH}"
ARG SETUPDIR=/tmp/toolbox/
RUN mkdir -p $SETUPDIR
WORKDIR $SETUPDIR

RUN apt-get update && apt-get -y install bedtools samtools \
    && mamba install -y snakemake racon minimap2 longshot flye filtlong trf \
    && mamba create -n nanoplot -c bioconda -c conda-forge nanoplot -y
    && pip install --upgrade pip && pip install --no-cache-dir -U quast numpy pandas argparse pysam biopython
RUN apt-get autoclean clean && rm -rf $SETUPDIR /var/lib/apt/lists/* /var/tmp/* /tmp/* /root/.cache/pip/*
WORKDIR /usr/bin
CMD ["/bin/bash"]