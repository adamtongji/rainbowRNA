#!/bin/bash
config=$1
source activate rainbow_env
softpath=`cat ~/bin/.rainbow_path_store`
python ${softpath}/__main__.py $config