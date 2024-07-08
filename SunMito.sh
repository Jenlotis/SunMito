#!/bin/bash
sudo apt-get install --assume-yes python-pip python3 python3-pip
pip install dash
pip install numpy
pip install pandas
pip install plotly
pip install multiqc
pip install conda
sudo apt-get install automake autoconf
sudo apt install default-jre
wget -P programs https://raw.githubusercontent.com/chrishah/MITObim/master/misc_scripts/downsample.py
wget -P programs/novo https://github.com/Edith1715/NOVOplasty/archive/refs/heads/master.zip
unzip programs/novo/master.zip -d programs/ 
wget -P programs/mitfi https://github.com/RemiAllio/MitoFinder/archive/master.zip
unzip programs/master.zip -d programs/ 
conda create --name MitFi mitofinder -c bioconda
mv programs/MitoFinder-master/install.sh programs/MitoFinder-master/install.sh.ok
sudo apt-get install --assume-yes bbmap
sudo apt-get install --assume-yes fastqc
sudo apt-get install --assume-yes conda
conda install -c bioconda seqkit
conda install -c bioconda chopper
sudo apt-get install --assume-yes trimmomatic
sudo apt-get install --assume-yes r-base
wget -P programs https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R
sudo apt-get install --assume-yes docker.io
python3 main.py
