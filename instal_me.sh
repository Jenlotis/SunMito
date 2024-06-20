sudo apt-get install --assume-yes python-pip python3 python3-pip
#curl https://pyenv.run | bash
#echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
#echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
#echo 'eval "$(pyenv init -)"' >> ~/.bashrc
#exec "$SHELL"
#pyenv install 2.7.18
pip install dash
pip install numpy
pip install pandas
pip install plotly
pip install multiqc
sudo apt-get install automake autoconf
sudo apt install default-jre
wget -P programs https://raw.githubusercontent.com/chrishah/MITObim/master/misc_scripts/downsample.py
wget -P programs/novo https://github.com/Edith1715/NOVOplasty/archive/refs/heads/master.zip
unzip programs/novo/master.zip -d programs/ 
wget -P programs/mitfi https://github.com/RemiAllio/MitoFinder/archive/master.zip
unzip programs/master.zip -d programs/ 
sudo apt-get install --assume-yes bbmap
sudo apt-get install --assume-yes fastqc
sudo apt-get install --assume-yes conda
conda install -c bioconda seqkit
conda install -c bioconda chopper
sudo apt-get install --assume-yes trimmomatic
sudo apt-get install --assume-yes r-base
wget -P programs https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R
sudo apt-get install --assume-yes docker
