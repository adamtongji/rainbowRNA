#!/usr/bin/env python
#coding:utf-8
import os
import sys
import subprocess


def sh(args):
    return subprocess.call(args,shell=True)


def hisat2_single_run(hisat2_path, hisat2_index, outputdir, outprefix,files,
                       Annotationfile, feature_count_path):
    sh("{0} -p 8 -x {1} -1 {2} -t -p 8 -q -S {3}/mapping/{4}.sam > {3}/mapping/{4}_log.txt"\
       .format(hisat2_path, hisat2_index, files, outputdir,outprefix))
    sh("{0} -T 8 -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}.sam"\
       .format(feature_count_path, Annotationfile,outputdir,outprefix))
    awk_args = "{if (NR>2){print $0}}"
    sh("cut -f1,7 {0}/mapping/{1}_count.txt | awk '{2}' >{0}/expr/{1}.txt".format(outputdir,outprefix, awk_args))



def hisat2_pair_run(hisat2_path, hisat2_index, outputdir, outprefix,file1,file2,
                    Annotationfile, feature_count_path):
    sh("{} -p 8 -x {} -1 {} -2 {} -t -p 8 -q -S {}.sam"\
       .format(hisat2_path, hisat2_index, file1, file2, outprefix))
    sh("{0} -p -T 8 -t exon -g gene_id -a {1} -o {2}/expr/{3}_count.txt {2}/mapping/{3}.sam" \
       .format(feature_count_path, Annotationfile, outputdir, outprefix))
    sh("{0} -T 8 -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}.sam" \
       .format(feature_count_path, Annotationfile, outputdir, outprefix))
    awk_args = "{if (NR>2){print $0}}"
    sh("cut -f1,7 {0}/mapping/{1}_count.txt | awk '{2}' >{0}/expr/{1}.txt".format(outputdir, outprefix, awk_args))
