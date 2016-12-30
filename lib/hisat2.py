#!/usr/bin/env python
#coding:utf-8
import os
import sys
import subprocess


def sh(args):
    return subprocess.call(args,shell=True)


def hisat2_single_run(hisat2_path, hisat2_index, outputdir, outprefix,files,
                       Annotationfile, feature_count_path, libtype):

    sh("{0} -p 4 -x {1} -1 {2} -t -q -S {3}/mapping/{4}.sam 1>/dev/null 2>{3}/mapping/{4}_log.txt"\
           .format(hisat2_path, hisat2_index, files, outputdir,outprefix))
    if libtype == "strand-specific" or libtype == "strandspecific":
        sh("{0} -s 1 -T 4 -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}.sam" \
           .format(feature_count_path, Annotationfile, outputdir, outprefix))
    else:
        sh("{0} -T 4 -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}.sam"\
           .format(feature_count_path, Annotationfile,outputdir,outprefix))
    awk_args = "{if (NR>2){print $0}}"
    sh("cut -f1,7 {0}/mapping/{1}_count.txt | awk '{2}' >{0}/expr/{1}.txt".format(outputdir,outprefix, awk_args))



def hisat2_pair_run(hisat2_path, hisat2_index, outputdir, outprefix,file1,file2,
                    Annotationfile, feature_count_path, libtype):
    sh("{0} -p 4 -x {1} -1 {2} -2 {3} -t -q -S {4}/mapping/{5}.sam 1>/dev/null 2>{4}/mapping/{5}_log.txt"\
           .format(hisat2_path, hisat2_index, file1, file2,outputdir, outprefix))

    if libtype == "strand-specific" or libtype == "strandspecific":
        sh("{0} -T 4 -s 1  -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}.sam" \
           .format(feature_count_path, Annotationfile, outputdir, outprefix))
    else:
        sh("{0} -T 4 -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}.sam" \
           .format(feature_count_path, Annotationfile, outputdir, outprefix))
    awk_args = "{if (NR>2){print $0}}"
    sh("cut -f1,7 {0}/mapping/{1}_count.txt | awk '{2}' >{0}/expr/{1}.txt".format(outputdir, outprefix, awk_args))
