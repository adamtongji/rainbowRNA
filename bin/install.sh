#!/bin/bash

#git clone https://github.com/adamtongji/rainbowRNA
conda install -c R -c bioconda -c default --file requirement.txt -n rainbow_env
source activate rainbow_env
mkdir -p ~/bin
ln -s rainbowRNA/__main__.py ~/bin/

