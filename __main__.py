#!/usr/bin/env python
#coding:utf-8
"""
Licensed Materials - Property of Tongji University
(C) Copyright Tongji University LifeBio Dept. 2016, 2016 All Rights Reserved
--------------------------------------------------------------------------------------
File Name   :  __main__.py
Old File Name: rainbow_rnaseq.py
Description : RNAseq pipeline

Author: Shaliu Fu
Change activity:
v1.0.1
Support single end data and mouse data
Generate final files format for customers

v1.0.2
Add support to circRNA, need to add samtools path for different version of samtools

v1.0.3
Update with featurecounts, take exon regions for targetscan in circRNA mode
Version: v1.0.3
"""
import os,sys
from lib.decorator import time_dec,time_func
from lib.IO.reportIO import get_report, get_ciri_report
from lib.IO.inputIO import ConfigParser, input_check
from lib.hisat2 import hisat2_single_run, hisat2_pair_run
from lib.star import star_pair_run, star_single_run, starindex
from lib.ciri import ciri_run, ciri_process
from multiprocessing import Pool
from lib.generate_web import update_html
import subprocess


def sh(args):
    return subprocess.call(args,shell=True)


@time_func
def mapping_main(HISAT2_path, HISAT2_index, Treat, Control,Outputdir, Seqtype,
                 Inputcheck,Max_process, Annotationfile, feature_count_path):
    if Inputcheck:
        input_check(Treat,Control)
    if not os.path.exists(Outputdir):
        os.makedirs(Outputdir)
    if not os.path.exists('{0}/mapping/'.format(Outputdir)):
        os.makedirs('{0}/mapping/'.format(Outputdir))
    if not os.path.exists('{0}/expr/'.format(Outputdir)):
        os.makedirs('{0}/expr/'.format(Outputdir))

    data_type = Seqtype[0]
    pool = Pool(processes= int(Max_process))
    result = []
    if data_type.lower() == "pair_end" or data_type.lower() =="pairend":
        for index in range(0, len(Treat), 2):
            treatment = Treat[index:index + 2]
            outprefix = "treat{0}".format(str((index / 2) + 1))
            result.append(pool.apply_async(hisat2_pair_run,(HISAT2_path, HISAT2_index, Outputdir,
                                                            outprefix, treatment[0], treatment[1],
                                                            Annotationfile, feature_count_path, Seqtype[1],)))

        for index in range(0, len(Control), 2):
            controls = Control[index:index + 2]
            outprefix = "control{0}".format(str((index / 2) + 1))
            result.append(pool.apply_async(hisat2_pair_run,(HISAT2_path, HISAT2_index, Outputdir,
                                                            outprefix, controls[0], controls[1],
                                                            Annotationfile, feature_count_path, Seqtype[1],)))

    elif data_type.lower() == "single_end" or data_type.lower() =="singleend":
        for index in range(len(Treat)):
            treatment = Treat[index]
            outprefix = "treat{0}".format(str(index + 1))
            result.append(pool.apply_async(hisat2_single_run, (HISAT2_path, HISAT2_index,
                                                              Outputdir, outprefix, treatment,
                                                              Annotationfile, feature_count_path, Seqtype[1],)))
        for index in range(len(Control)):
            controls = Control[index]
            outprefix = "control{0}".format(str(index + 1))
            result.append(pool.apply_async(hisat2_single_run,(HISAT2_path, HISAT2_index,
                                                             Outputdir, outprefix, controls,
                                                             Annotationfile, feature_count_path, Seqtype[1],)))

    else:
        print "Unidentified library type!"
        sys.exit(1)
    pool.close()
    pool.join()


@time_func
def mapping_main_star(STARpath, STARindex, STARindexdir,Treat, Control,Outputdir, Seqtype,
                 Inputcheck,Max_process, Annotationfile, feature_count_path, Genomefile, ReadLength):
    if Inputcheck:
        input_check(Treat, Control)
    if not os.path.exists(Outputdir):
        os.makedirs(Outputdir)
    if not os.path.exists('{0}/mapping/'.format(Outputdir)):
        os.makedirs('{0}/mapping/'.format(Outputdir))
    if not os.path.exists('{0}/expr/'.format(Outputdir)):
        os.makedirs('{0}/expr/'.format(Outputdir))

    if STARindex.lower()=="false":
        starindex(STARpath, STARindexdir, Genomefile,Annotationfile,ReadLength)

    data_type = Seqtype[0]
    pool = Pool(processes=int(Max_process))
    result = []
    if data_type.lower() == "pair_end" or data_type.lower() == "pairend":
        for index in range(0, len(Treat), 2):

            treatment = Treat[index:index + 2]
            outprefix = "treat{0}".format(str((index / 2) + 1))
            result.append(pool.apply_async(star_pair_run, (STARpath, STARindexdir, Outputdir,
                                                             outprefix, treatment[0], treatment[1],
                                                             Annotationfile, feature_count_path, Seqtype[1],)))

        for index in range(0, len(Control), 2):
            controls = Control[index:index + 2]
            outprefix = "control{0}".format(str((index / 2) + 1))
            result.append(pool.apply_async(star_pair_run, (STARpath, STARindexdir, Outputdir,
                                                             outprefix, controls[0], controls[1],
                                                             Annotationfile, feature_count_path, Seqtype[1],)))

    elif data_type.lower() == "single_end" or data_type.lower() == "singleend":
        for index in range(len(Treat)):
            treatment = Treat[index]
            outprefix = "treat{0}".format(str(index + 1))
            result.append(pool.apply_async(star_single_run, (STARpath, STARindexdir,
                                                               Outputdir, outprefix, treatment,
                                                               Annotationfile, feature_count_path, Seqtype[1],)))
        for index in range(len(Control)):
            controls = Control[index]
            outprefix = "control{0}".format(str(index + 1))
            result.append(pool.apply_async(star_single_run, (STARpath, STARindexdir,
                                                               Outputdir, outprefix, controls,
                                                               Annotationfile, feature_count_path, Seqtype[1],)))

    else:
        print "Unidentified library type!"
        sys.exit(1)
    pool.close()
    pool.join()


@time_func
def deseq_main(Outputdir, Pair_rep, pvalue,Expr_dir=None):
    if Expr_dir is None:
        Expr_dir = "{}/expr".format(Outputdir)
    #if not os.path.exists('{0}/results'.format(Outputdir)):
    sh('mkdir -p {0}/results'.format(Outputdir))
    sh('mkdir -p {0}/results/up'.format(Outputdir))
    sh('mkdir -p {0}/results/down'.format(Outputdir))

    sh("Rscript {}/rscript/deseq2.r {} {} {} &>>deseq.log"\
       .format(SOFT_PATH, Expr_dir,  Pair_rep, pvalue))


@time_func
def downstream_main(Outputdir, Genome):
    if Genome.lower() == 'hg19' or Genome.lower() == 'hg38':
        sh("Rscript {1}/rscript/pathway.r {0}/results/ &>pathway.log".format(Outputdir,SOFT_PATH))
    if Genome.lower() == 'mm9' or Genome.lower() == 'mm10':
        sh("Rscript {1}/rscript/mouse_pathway.r {0}/results &>pathway.log".format(Outputdir,SOFT_PATH))
    if not os.path.exists('{0}/results/gsea/'.format(Outputdir)):
        sh('mkdir -p {0}/results/gsea/'.format(Outputdir))
    else:
        sh('rm -rf {0}/results/gsea/'.format(Outputdir))
        sh('mkdir -p {0}/results/gsea/'.format(Outputdir))
    gsea_file = [i.rstrip().split('\t') for i in open('{0}/results/Treat_vs_control_diff.txt'.format(Outputdir))]
    gsea_file = [i[:1] + i[4:5] for i in gsea_file]
    gsea_file = gsea_file[1:]
    gsea_file2 = []
    for _line in gsea_file:
        if _line[1] == 'NA':
            pass
        else:
            _line[0] = _line[0].upper()
            gsea_file2.append(_line)
    gsea_file = ['\t'.join(i) for i in gsea_file2]
    with open('{0}/results/gsea_input.txt'.format(Outputdir), 'w') as f:
        for line in gsea_file:
            print >> f, line

    sh("sort -k 2gr {0}/results/gsea_input.txt >{0}/results/gsea/gsea_input.rnk".format(Outputdir))
    sh("export LANG=en_US.UTF-8;java -cp {1}/bin/gsea2-2.2.4.jar -Xmx8g xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.bp.v5.1.symbols.gmt\
         -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk {0}/results/gsea/gsea_input.rnk -scoring_scheme weighted\
          -rpt_label bp -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500\
           -set_min 15 -zip_report false -out {0}/results/gsea/ -gui false &>gsea.log".format(Outputdir,SOFT_PATH))
    sh("export LANG=en_US.UTF-8; java -cp {1}/bin/gsea2-2.2.4.jar -Xmx8g xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v5.1.symbols.gmt\
         -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk {0}/results/gsea/gsea_input.rnk -scoring_scheme weighted\
          -rpt_label kegg -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500\
         -set_min 15 -zip_report false -out {0}/results/gsea/ -gui false &>>gsea.log".format(Outputdir,SOFT_PATH))
    if not os.path.exists('{0}/results/cytoscape/'.format(Outputdir)):
        sh('mkdir {0}/results/cytoscape'.format(Outputdir))
    if Genome.lower() == 'hg19' or Genome.lower() == 'hg38':
        sh('{1}/lib/DEG2network.py -p 0.05 -n 5 -k\
         {1}/lib/db/merged_KEGG.txt -i {0}/results/Treat_vs_control_diff.txt -d \
           {0}/results/cytoscape'.format(Outputdir,SOFT_PATH))
    elif Genome.lower() == 'mm9' or Genome.lower() == 'mm10':
        sh('{1}/lib/DEG2network_mouse.py -p 0.05 -n 5 -k\
         {1}/lib/db/mouse_merged_KEGG.txt -i {0}/results/Treat_vs_control_diff.txt -d \
           {0}/results/cytoscape'.format(Outputdir,SOFT_PATH))
    get_report(Outputdir)
    # sh("mkdir -p {}/html".format(Outputdir))
    if os.path.exists("{}/html".format(Outputdir)):
        sh("rm -rf {}/html".format(Outputdir))
    sh("cp -r {0}/lib/html_template {1}/html; mv {1}/html/index.html {1}/html/Web_Report.html"\
       .format(SOFT_PATH, Outputdir))
    update_html("{0}".format(Outputdir),"{0}/html/".format(Outputdir))


@time_func
def circ_mapping_main(Genomefile,Treat, Control,Outputdir, Seqtype,
                 Inputcheck,Max_process, Annotationfile):
    if Inputcheck:
        input_check(Treat, Control)
    if not os.path.exists(Outputdir):
        os.makedirs(Outputdir)
    if not os.path.exists('{0}/ciri_out/'.format(Outputdir)):
        os.makedirs('{0}/ciri_out/'.format(Outputdir))
    data_type = Seqtype[0]
    pool = Pool(processes=int(Max_process))
    result =[]
    if data_type.lower() =="pair_end" or data_type.lower() =="pairend":
        for index in range(0, len(Treat), 2):
            treatment = Treat[index:index + 2]
            outprefix = "{1}/ciri_out/treat{0}".format(str((index / 2) + 1),Outputdir)
            result.append(pool.apply_async(ciri_run, (outprefix, Genomefile, Annotationfile,
                                                      treatment[0],treatment[1],)))

        for index in range(0, len(Control), 2):
            controls = Control[index:index + 2]
            outprefix = "{1}/ciri_out/control{0}".format(str((index / 2) + 1),Outputdir)
            result.append(pool.apply_async(ciri_run, (outprefix, Genomefile, Annotationfile,
                                                      controls[0],controls[1],)))
    else:
        print "error in data type!"
        sys.exit(1)

    pool.close()
    pool.join()



@time_func
def circ_process_main(Outputdir, Pair_rep, pvalue,Annotationfile,Genome,Genomefile,Treat):
    Num_of_rep = int(len(Treat)/2)
    ciri_process(Outputdir,Num_of_rep, Pair_rep,pvalue,SOFT_PATH,Annotationfile,Genomefile,Genome)
    get_ciri_report(Outputdir)


@time_dec
def main():
    if len(sys.argv) == 2:
        config_file = sys.argv[1]
    else:
        config_file = False

    myParser = ConfigParser()
    if config_file:
        command, paramter = myParser.load_configs(config_file)
    else:
        try:
            f = open('config.txt')
        except IOError:
            print "Please give config files"
            sys.exit(1)
        command, paramter = myParser.load_configs()
    command_functions = {"mapping":mapping_main_star,
                         "deseq":deseq_main,
                         "downstream":downstream_main,
                         "circ_mapping":circ_mapping_main,
                         "circ_process":circ_process_main
                         }

    for subcom in command:
        command_functions[subcom](**paramter[subcom])


if __name__=="__main__":
    global SOFT_PATH
    usr_home=os.path.expanduser('~')
    SOFT_PATH=[i.rstrip() for i in open("{}/bin/.rainbow_path_store".format(usr_home))][0]
    main()


