#!/usr/bin/env python
#coding:utf-8
import os
import sys
import subprocess


def sh(args):
    return subprocess.call(args,shell=True)


def starindex(starpath, star_index, genome_file, annotation_file, readlength):
    sh("{0} --runMode genomeGenerate --runThreadN 8 --genomeDir {1}\
     --genomeFastaFiles {2} --sjdbGTFfile {3} --sjdbOverhang {4}  &>>star.log"\
       .format(starpath, star_index,genome_file, annotation_file, str(int(readlength)-1)))
    print "STAR Index build successfully at {0}".format(star_index)


def star_pair_run(starpath, star_index,outputdir,outprefix, file1, file2,
    Annotationfile, feature_count_path, libtype):
    if file1.endswith(".gz"):
        if libtype=="strand-specific":
            sh("{0} --genomeDir {1} --runThreadN 8 --readFilesCommand zcat --outSAMattributes\
            Standard --outSAMtype BAM Unsorted  --readFilesIn {2} {4} --outFileNamePrefix {5}/mapping/{3} &>>star.log"\
           .format(starpath, star_index,file1, outprefix,file2,outputdir))
        else:
            sh("{0} --genomeDir {1} --runThreadN 8 --readFilesCommand zcat --outSAMattributes\
            Standard --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --readFilesIn\
             {2} {4} --outFileNamePrefix {5}/mapping/{3} &>>star.log".format(starpath, star_index,file1, outprefix,file2,outputdir))
    elif file1.endswith(".fastq") or file1.endswith(".fq"):
        if libtype=="strand-specific":
            sh("{0} --genomeDir {1} --runThreadN 8 --outSAMattributes\
            Standard --outSAMtype BAM Unsorted --readFilesIn {2} {4} --outFileNamePrefix {5}/mapping/{3} &>>star.log"\
           .format(starpath, star_index,file1, outprefix,file2,outputdir))
        else:
            sh("{0} --genomeDir {1} --runThreadN 8 --outSAMattributes\
            Standard --outSAMstrandField intronMotif --outSAMtype BAM Unsorted  --readFilesIn\
             {2} {4} --outFileNamePrefix {5}/mapping/{3} &>>star.log".format(starpath, star_index,file1, outprefix,file2,outputdir))
    else:
        print "please give fastq or gzipped fastq files"
        sys.exit(1)

    # sh("samtools sort -n {0}Aligned.sortedByCoord.out.bam {0}.sorted".format(outprefix))
    # sh("samtools index {0}.sorted.bam".format(outprefix))
    if libtype == "strand-specific" or libtype == "strandspecific":
        sh("{0} -p -s 1 -T 4 -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}Aligned.out.bam &>>count.log" \
           .format(feature_count_path, Annotationfile, outputdir, outprefix))
    else:
        sh("{0} -p -T 4 -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}Aligned.out.bam &>>count.log" \
           .format(feature_count_path, Annotationfile, outputdir, outprefix))
    awk_args = "{if (NR>2){print $0}}"
    sh("cut -f1,7 {0}/mapping/{1}_count.txt | awk '{2}' >{0}/expr/{1}.txt".format(outputdir, outprefix, awk_args))
    # sh("cut -f1,7 {0}/mapping/{1}_count.txt  >>{0}/expr/{1}.txt".format(outputdir, outprefix))


def star_single_run(starpath, star_index,outputdir,outprefix, file1,
    Annotationfile, feature_count_path, libtype):
    if file1.endswith(".gz"):
        if libtype=="strand-specific":
            sh("{0} --genomeDir {1} --runThreadN 8 --readFilesCommand zcat --outSAMattributes\
            Standard --outSAMtype BAM Unsorted  --readFilesIn {2} --outFileNamePrefix {4}/mapping/{3} &>>star.log"\
           .format(starpath, star_index,file1, outprefix,outputdir))
        else:
            sh("{0} --genomeDir {1} --runThreadN 8 --readFilesCommand zcat --outSAMattributes\
            Standard --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --readFilesIn\
             {2} --outFileNamePrefix {4}/mapping/{3} &>>star.log".format(starpath, star_index,file1, outprefix,outputdir))
    elif file1.endswith(".fastq") or file1.endswith(".fq"):
        if libtype=="strand-specific":
            sh("{0} --genomeDir {1} --runThreadN 8 --outSAMattributes\
            Standard --outSAMtype BAM Unsorted --readFilesIn {2} --outFileNamePrefix {4}/mapping/{3} &>>star.log"\
           .format(starpath, star_index,file1, outprefix,outputdir))
        else:
            sh("{0} --genomeDir {1} --runThreadN 8 --outSAMattributes\
            Standard --outSAMstrandField intronMotif --outSAMtype BAM Unsorted  --readFilesIn\
             {2} --outFileNamePrefix {4}/mapping/{3} &>>star.log".format(starpath, star_index,file1, outprefix,outputdir))
    else:
        print "please give fastq or gzipped fastq files"
        sys.exit(1)


    if libtype == "strand-specific" or libtype == "strandspecific":
        sh("{0} -s 1 -T 4 -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}Aligned.out.bam &>>count.log" \
           .format(feature_count_path, Annotationfile, outputdir, outprefix))
    else:
        sh("{0} -T 4 -t exon -g gene_id -a {1} -o {2}/mapping/{3}_count.txt {2}/mapping/{3}Aligned.out.bam &>>count.log" \
           .format(feature_count_path, Annotationfile, outputdir, outprefix))
    awk_args = "{if (NR>2){print $0}}"
    sh("cut -f1,7 {0}/mapping/{1}_count.txt | awk '{2}' >{0}/expr/{1}.txt".format(outputdir, outprefix, awk_args))
    # sh("cut -f1,7 {0}/mapping/{1}_count.txt  >>{0}/expr/{1}.txt".format(outputdir, outprefix))