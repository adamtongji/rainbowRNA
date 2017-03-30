#!/bin/bash
genes=$1
outfile=$2
softpath=$3
cd ${softpath}/lib/db/RBP_binding_sites
rbpfile=`ls`
for file in ${rbpfile}; do
    filepref=${file%.txt}
    overlap=$(intersectBed -b $genes -a $file -wa -s | cut -f 1-4 | uniq | wc -l)
    awk -v pref=$filepref -v num=$overlap 'BEGIN{printf "%s\t%s\n",pref,num}' >> $outfile

done