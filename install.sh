#!/bin/bash

#git clone https://github.com/adamtongji/rainbowRNA
# if error is not about conda, user can skip the step of conda create.
conda create -n rainbow_env --file trans_req.txt # -c R -c bioconda -y
# conda update --all -c bioconda -c R -c default
##
source activate rainbow_env
mkdir -p ~/bin
myPath=`pwd`
ln -s $myPath/bin/run.sh ~/bin/rainbowRNA
chmod 777 $myPath/bin/run.sh
chmod 777 ~/bin/rainbowRNA
echo $myPath > ~/bin/.rainbow_path_store
chmod 700 ~/bin/.rainbow_path_store
conda install java-jdk -c cyclus -y