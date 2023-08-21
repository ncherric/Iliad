.. _projectinfo/container_images:

==================================================
Definition files for Singularity and Docker Images
==================================================




Singularity definition file
===========================

.. code-block:: none 

    Bootstrap: docker
    From: ubuntu:20.04

    %environment
        PATH=$PATH:/usr/bin:/usr/bin/bcftools-1.16/bin/:/usr/bin/samtools-1.16/bin/:/usr/bin/miniconda3/bin:$PATH

        export PATH
        export BCFTOOLS_PLUGINS=/usr/bin/bcftools-1.16/plugins/

    %post
        apt-get -y update
        apt -y update

        apt install software-properties-common -y

        # Install bwa
        apt-get install -y wget \
        bwa \
        curl \
        git \
        tar \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        unzip

        # build-essential for bcftools1.16
        # zlib1g-dev for bcftools1.16
        # libbz2-dev for bcftools1.16
        # liblzma-dev for bcftools1.16
        # libncurses5-dev for samtools1.16
        # unzip for illumina product files
        # https://www.devmanuals.net/install/ubuntu/ubuntu-20-04-focal-fossa/installing-bwa-on-ubuntu20-04.html

        mkdir -p /usr/share/man/man1/
        apt-get install -y openjdk-8-jdk
        apt-get install -y openjdk-8-jre
        update-alternatives --config java
        update-alternatives --config javac

        wget -P /usr/bin/ https://github.com/broadinstitute/picard/releases/download/2.27.2/picard.jar

        wget -P /usr/bin/ https://sourceforge.net/projects/samtools/files/samtools/1.16/bcftools-1.16.tar.bz2
        cd /usr/bin/
        tar xvfj /usr/bin/bcftools-1.16.tar.bz2
        cd /usr/bin/bcftools-1.16/
        ./configure --prefix=/usr/bin/bcftools-1.16
        make
        make install
        export PATH="/usr/bin/bcftools-1.16/bin:$PATH"
        rm /usr/bin/bcftools-1.16.tar.bz2

        # install samtools
        wget -P /usr/bin/ https://sourceforge.net/projects/samtools/files/samtools/1.16/samtools-1.16.tar.bz2
        cd /usr/bin/
        tar xvfj /usr/bin/samtools-1.16.tar.bz2
        cd /usr/bin/samtools-1.16/
        ./configure --prefix=/usr/bin/samtools-1.16
        make
        make install
        export PATH="/usr/bin/samtools-1.16/bin:$PATH"
        rm /usr/bin/samtools-1.16.tar.bz2

        wget -P /usr/bin/bcftools-1.16/plugins/ https://software.broadinstitute.org/software/gtc2vcf/gtc2vcf_1.16-20221221.zip
        cd /usr/bin/bcftools-1.16/plugins/
        unzip gtc2vcf_1.16-20221221.zip

        # Install miniconda 
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/bin/miniconda3/
        rm Miniconda3-latest-Linux-x86_64.sh

        . /usr/bin/miniconda3/etc/profile.d/conda.sh
        export PATH="/usr/bin/miniconda3/bin:$PATH"

        #Conda configuration of channels from .condarc file
        conda config --file /.condarc --add channels defaults
        conda config --file /.condarc --add channels conda-forge
        conda update conda
        #List installed environments
        conda list
        conda install mamba -n base -c conda-forge

    %labels
        Author noahherrick1@gmail.com
        Version 1.2.0
        Pipeline Iliad-Genomic-Data-Pipeline
        SoftwareToolsVersion 1.16


Docker 'definition' file
========================

.. code-block:: none 

    FROM ubuntu:20.04

    MAINTAINER	Noah Herrick	<noahherrick1@gmail.com>
    ARG BCFTOOLS_VERSION=1.16
    ARG SAMTOOLS_VERSION=1.16
    ARG GTC2VCF_VERSION=1.14
    ENV PATH /usr/bin:/usr/bin/bcftools-1.16/bin/:/usr/bin/samtools-1.16/bin/:/usr/bin/miniconda3/bin:${PATH}

    RUN apt-get -y update \
    && apt -y update \
    && apt install software-properties-common -y \
    && apt-get install -y \
    wget \
    bwa \
    curl \
    git \
    tar \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    unzip \
    screen \
    less \
    nano && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

    RUN export PATH && \
    export BCFTOOLS_PLUGINS=/usr/bin/bcftools-${BCFTOOLS_VERSION}/plugins/ 

    RUN wget -P /usr/bin/ https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd /usr/bin && \
    tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2
    WORKDIR /usr/bin/bcftools-${BCFTOOLS_VERSION}
    RUN ./configure --prefix=/usr/bin/bcftools-${BCFTOOLS_VERSION} && make && make install && \
    export PATH="/usr/bin/bcftools-${BCFTOOLS_VERSION}/bin:$PATH" && \
    rm /usr/bin/bcftools-${BCFTOOLS_VERSION}.tar.bz2

    RUN wget -P /usr/bin/ https://sourceforge.net/projects/samtools/files/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd /usr/bin/ && \
    tar xvfj /usr/bin/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd /usr/bin/samtools-${SAMTOOLS_VERSION}/ && \
    ./configure --prefix=/usr/bin/samtools-${SAMTOOLS_VERSION} && make && make install && \
    export PATH="/usr/bin/samtools-${SAMTOOLS_VERSION}/bin:$PATH" && \
    rm /usr/bin/samtools-${SAMTOOLS_VERSION}.tar.bz2

    RUN wget -P /usr/bin/bcftools-${BCFTOOLS_VERSION}/plugins/ https://software.broadinstitute.org/software/gtc2vcf/gtc2vcf_${GTC2VCF_VERSION}-20220112.zip && \
    cd /usr/bin/bcftools-${BCFTOOLS_VERSION}/plugins/ && \
    unzip gtc2vcf_${GTC2VCF_VERSION}-20220112.zip

    # Install miniconda
    ENV CONDA_DIR /opt/conda
    RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

    # Put conda in path so we can use conda activate
    ENV PATH=$CONDA_DIR/bin:$PATH

    RUN mkdir -p /usr/projects/ && cd /usr/projects/

    SHELL ["/bin/bash", "-c"]

    RUN conda init bash

    RUN source ~/.bashrc