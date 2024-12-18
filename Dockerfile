FROM golang:1.23.2-bookworm

ADD Dashboard ./Dashboard/

RUN apt-get update && apt-get install -y \
    build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools \
    libseccomp-dev wget pkg-config git glib-2.0 libfuse3-dev autoconf \
    libtool python3 python3-pip python3.11-venv unzip automake \
    default-jre bbmap fastqc trimmomatic r-base docker.io \
    && rm -rf /var/lib/apt/lists/*
    
RUN wget https://github.com/sylabs/singularity/releases/download/v4.2.1/singularity-ce-4.2.1.tar.gz
RUN tar -xzf singularity-ce-4.2.1.tar.gz
RUN cd singularity-ce-4.2.1 && ./mconfig && \
    cd builddir && make && \
    make install
RUN singularity pull --arch amd64 library://remiallio/default/mitofinder:v1.4.2

RUN python3 -m venv srodo \
    && . srodo/bin/activate \
    && pip install dash \
    dash-bootstrap-components \
    numpy \
    pandas \
    plotly \
    multiqc

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -u -p /opt/miniconda3 \
    && rm /tmp/miniconda.sh \
    && /opt/miniconda3/bin/conda init bash \
    && /opt/miniconda3/bin/conda init zsh

RUN wget -P programs https://raw.githubusercontent.com/chrishah/MITObim/master/misc_scripts/downsample.py \
    && wget -P programs/novo https://github.com/Edith1715/NOVOplasty/archive/refs/heads/master.zip \
    && unzip programs/novo/master.zip -d programs/ \
    && wget -P programs/mitfi https://github.com/RemiAllio/MitoFinder/archive/master.zip \
    && unzip programs/mitfi/master.zip -d programs/ \
    && mv programs/MitoFinder-master/install.sh programs/MitoFinder-master/install.sh.ok \
    && wget -P programs https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R

WORKDIR /Dashboard
EXPOSE 8050

CMD ["python3", "./main.py"]
