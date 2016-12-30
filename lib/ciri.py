#!/usr/bin/env python
#coding:utf-8
import os
import sys
import subprocess

def sh(args):
    return subprocess.call(args,shell=True)


def ciri_run(outprefix,genomefile,annotationfile, *files):

    if len(files)!=2:
        print "error in ciri_run! exit"
        sys.exit(1)
    sh("bwa mem {0} {1} {2} -t 4 -T 19 -v 1 > {3}.sam".format(genomefile,files[0],files[1], outprefix))
    sh("CIRI2 -T 4 -I {0}.sam -O {0}.txt -F {1} -A {2}".format(outprefix, genomefile, annotationfile))


def ciri_process(outputdir,nrep, pairtype, pvalue,softpath,annotationfile,genomefile,genome="hg19"):
    if not os.path.isdir("{0}/results/".format(outputdir)):
        sh("mkdir {0}/results/".format(outputdir))
    if not os.path.isdir("{0}/expr/".format(outputdir)):
        sh("mkdir {0}/expr/".format(outputdir))
    # process ciri_expr files to merge_expr and merge_annotation files.
    sh('mkdir {0}/results/up;mkdir {0}/results/down {0}/results/targetscan'.format(outputdir))
    gene_list=[]
    for filetype in ["control","treat"]:
        for num in range(1,nrep+1):
            filename = "{0}/ciri_out/{1}{2}.txt".format(outputdir,filetype,str(num))
            temp = [i.rstrip().split("\t") for i in open(filename)]
            temp = temp[1:]
            for line in temp:
                genename = line[1]+'_'+line[2]+'_'+line[3]+'_'+line[10]
                if not genename in gene_list:
                    gene_list.append(genename)

    annot_table = {}
    for filetype in ["control","treat"]:
        for num in range(1,nrep+1):
            expr_table = [[i] for i in gene_list]
            filename = "{0}/ciri_out/{1}{2}.txt".format(outputdir,filetype,str(num))
            temp = [i.rstrip().split("\t") for i in open(filename)]
            temp = temp[1:]
            for line in temp:
                genename = line[1]+'_'+line[2]+'_'+line[3]+'_'+line[10]
                if not genename in gene_list:
                    print "loss keys!"
                if not annot_table.has_key(genename):
                    annot_table[genename]='\t'.join(line[1:4]+line[8:9]+[line[9].rstrip(",")]+line[10:11])

            genesummary=[[line[1]+'_'+line[2]+'_'+line[3]+'_'+line[10] for line in temp]]\
                        +[[line[4] for line in temp]]
            for gene_id in expr_table:
                if gene_id[0] in genesummary[0]:
                    for index, item in enumerate(genesummary[0]):
                        if gene_id[0]==item:
                            gene_id.append(genesummary[1][index])
                else:
                    gene_id.append("0")
            with open("{0}/expr/{1}{2}.txt".format(outputdir,filetype,str(num)),"w") as f:
                expr_table2 = ['\t'.join(i) for i in expr_table]
                for line in expr_table2:
                    print >>f, line

    with open("{0}/results/annotation_table.txt".format(outputdir),"w") as f:
        for mykey in annot_table.iterkeys():
            print >> f, "{1}\t{0}".format(mykey,annot_table[mykey])
    awk_arg='{printf "%s\\t%s\\t%s\\t%s\\n",$3,$2,$1,$4}'
    if genome.lower()=="mm9" or genome.lower()=="hg19":
        sh("intersectBed -a {0}/results/annotation_table.txt -b\
         {3}/lib/db/{2}_circRNA_database.txt -wa -wb -s -f 1.0 -r\
         | cut -f 4,5,7,11| awk '{1}' | uniq > {0}/results/annotation_table_circbase.txt"\
           .format(outputdir, awk_arg, genome, softpath))

    sh("Rscript {4}/rscript/circ_deseq.r {0}/expr/ {1} {2} {3}"\
       .format(outputdir, nrep, pairtype, pvalue, softpath))
    # make circ_down and up to bed format
    up_file=[i.rstrip() for i in open("{0}/results/circ_up.txt".format(outputdir))]
    with open("{0}/results/circ_up2.txt".format(outputdir),"w") as f:
        up_out=[i.split("_")[:3]+["."]+[i]+i.split("_")[3:] for i in up_file]
        up_out2 = ['\t'.join(i) for i in up_out]
        for line in up_out2:
            print >>f, line
    down_file=[i.rstrip() for i in open("{0}/results/circ_down.txt".format(outputdir))]
    with open("{0}/results/circ_down2.txt".format(outputdir),"w") as f:
        down_out=[i.split("_")[:3]+[i]+["."]+i.split("_")[3:] for i in down_file]
        down_out2 = ['\t'.join(i) for i in down_out]
        for line in down_out2:
            print >>f, line

    # hs is 9606
    if genome.lower() == 'hg19': #or genome.lower() == 'hg38':
        awk_arg = '{printf "%s\\t%s\\t%s\\n",$1,9606,$2}'
    elif genome.lower() == 'mm9':
        # need revise
        awk_arg = '{printf "%s\\t%s\\t%s\\n",$1,10090,$2}'
        # pass
    elif genome.lower()=="rn6":
        awk_arg = '{printf "%s\\t%s\\t%s\\n",$1,10116,$2}'
    else:
        print ""


    if awk_arg:
        sh("grep 'exon' {3} |intersectBed -a {1}/results/circ_up2.txt -b stdin\
        |bedtools getfasta -fi {0} -bed stdin -fo stdout -name -tab|\
          awk '{2}' > {1}/results/targetscan/circ_up_sequence.txt"\
           .format(genomefile,outputdir,awk_arg,annotationfile))
        sh("grep 'exon' {3} |intersectBed -a {1}/results/circ_down2.txt -b stdin\
            |bedtools getfasta -fi {0} -bed stdin -fo stdout -name -tab|\
              awk '{2}' > {1}/results/targetscan/circ_down_sequence.txt" \
           .format(genomefile, outputdir, awk_arg, annotationfile))

    sh("targetscan_60.pl /home/Public/software/targetscan/miR_Family_Info.txt {0}/results/targetscan/circ_up_sequence.txt\
       {0}/results/targetscan/circ_up_miRNA.txt ".format(outputdir))
    sh("targetscan_60.pl /home/Public/software/targetscan/miR_Family_Info.txt {0}/results/targetscan/circ_down_sequence.txt\
        {0}/results/targetscan/circ_down_miRNA.txt ".format(outputdir))
    sh("cp /home/Public/software/targetscan/TA_SPS_by_seed_region.txt ./")
    sh("targetscan_60_context_scores.pl /home/Public/software/targetscan/miR_for_context_scores.txt\
     {0}/results/targetscan/circ_down_miRNA.txt {0}/results/circ_down2.txt {0}/results/targetscan/circ_down_score.txt".format(outputdir))
    sh("targetscan_60_context_scores.pl /home/Public/software/targetscan/miR_for_context_scores.txt\
        {0}/results/targetscan/circ_up_miRNA.txt {0}/results/circ_up2.txt {0}/results/targetscan/circ_up_score.txt".format(outputdir))
    sh("rm ./TA_SPS_by_seed_region.txt")
    sh("cut -f1-13,17 {0}/results/targetscan/circ_up_score.txt >{0}/results/targetscan/circ_up_tmp.txt".format(outputdir))
    sh("cut -f1-13,17 {0}/results/targetscan/circ_down_score.txt >{0}/results/targetscan/circ_down_tmp.txt".format(outputdir))
    # filter high and medium targets
    #with open("{0}/results/targetscan/circ_down_full.txt",'w') as f:
    high_up = open("{0}/results/targetscan/circ_up_high.txt".format(outputdir),'w')
    pass_up = open("{0}/results/targetscan/circ_up_pass.txt".format(outputdir),'w')
    up_mi = [i.rstrip().split("\t") for i in open("{0}/results/targetscan/circ_up_tmp.txt".format(outputdir))]
    head = '\t'.join(up_mi[0])
    print >>high_up, head
    print >>pass_up, head
    for myline in range(1,len(up_mi),1):
        if int(up_mi[myline][12]) > 90:
            _mylines = "\t".join(up_mi[myline])
            print >>high_up, _mylines
        if int(up_mi[myline][12]) > 60:
            _mylines = "\t".join(up_mi[myline])
            print >> pass_up, _mylines

    high_down = open("{0}/results/targetscan/circ_down_high.txt".format(outputdir), 'w')
    pass_down = open("{0}/results/targetscan/circ_down_pass.txt".format(outputdir), 'w')
    down_mi = [i.rstrip().split("\t") for i in open("{0}/results/targetscan/circ_down_tmp.txt".format(outputdir))]
    head = '\t'.join(down_mi[0])
    print >> high_down, head
    print >> pass_down, head
    for myline in range(1, len(down_mi), 1):
        if int(down_mi[myline][12]) > 90:
            _mylines = "\t".join(down_mi[myline])
            print >> high_down, _mylines
        if int(down_mi[myline][12]) > 60:
            _mylines = "\t".join(down_mi[myline])
            print >> pass_down, _mylines
    high_up.close()
    high_down.close()
    pass_down.close()
    pass_up.close()

    # start with hg19 only first
    if genome.lower()=='hg19' or genome.lower()=='hg38':
        sh("Rscript {1}/rscript/pathway.r {0}/results/".format(outputdir,softpath))
    if genome.lower()=='mm9' or genome.lower()=='mm10':
        sh("Rscript {1}/rscript/mouse_pathway.r {0}/results".format(outputdir,softpath))

    ## RBP binding with circ and gene
    with open("{0}/results/circ_up.bed".format(outputdir),"w") as f:
        temp=[i.rstrip().split("_") for i in open("{0}/results/circ_up2.txt".format(outputdir))]
        temp2 = [i[:3]+["."]+["."]+[i[3]] for i in temp]
        temp3 = ['\t'.join(i) for i in temp2]
        for line in temp3:
            print >>f, line
    with open("{0}/results/circ_down.bed".format(outputdir), "w") as f:
        temp = [i.rstrip().split("_") for i in open("{0}/results/circ_down2.txt".format(outputdir))]
        temp2 = [i[:3] + ["."] + ["."] + [i[3]] for i in temp]
        temp3 = ['\t'.join(i) for i in temp2]
        for line in temp3:
            print >> f, line
    if genome.lower() == 'hg19' or genome.lower() == 'hg38':
        sh("bash {1}/lib/RBP_count.sh\
         {0}/results/circ_down.bed {0}/results/ciri_down_RBP.txt".format(outputdir,softpath))
        sh("bash {1}/lib/RBP_count.sh\
            {0}/results/circ_up.bed {0}/results/ciri_up_RBP.txt".format(outputdir,softpath))

    if not os.path.exists('{0}/results/gsea/'.format(outputdir)):
        sh('mkdir {0}/results/gsea/'.format(outputdir))
    else:
        sh('rm -rf {0}/results/gsea/'.format(outputdir))
        sh('mkdir {0}/results/gsea/'.format(outputdir))

    gsea_file = [i.rstrip().split('\t') for i in open('{0}/results/Treat_vs_control_diff.txt'.format(outputdir))]
    gsea_file = [i[:1]+i[4:5] for i in gsea_file]
    gsea_file = gsea_file[1:]
    gsea_file2 = []
    for _line in gsea_file:
        if _line[1]=='NA' or _line[1]=='n/a':
            pass
        else:
            _line[0]=_line[0].upper()
            gsea_file2.append(_line)
    gsea_file = ['\t'.join(i) for i in gsea_file2]
    with open('{0}/results/gsea_input.txt'.format(outputdir),'w') as f:
        for line in gsea_file:
            print >>f, line

    sh("sort -k 2gr {0}/results/gsea_input.txt >{0}/results/gsea/gsea_input.rnk".format(outputdir))
    sh("export LANG=en_US.UTF-8;java -cp /usr/local/src/gsea2-2.2.2.jar -Xmx2g xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.bp.v5.1.symbols.gmt\
     -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk {0}/results/gsea/gsea_input.rnk -scoring_scheme weighted\
      -rpt_label bp -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500\
       -set_min 15 -zip_report false -out {0}/results/gsea/ -gui false".format(outputdir))
    sh("export LANG=en_US.UTF-8; java -cp /usr/local/src/gsea2-2.2.2.jar -Xmx2g xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v5.1.symbols.gmt\
     -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk {0}/results/gsea/gsea_input.rnk -scoring_scheme weighted\
      -rpt_label kegg -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500\
     -set_min 15 -zip_report false -out {0}/results/gsea/ -gui false".format(outputdir))

    # step 4 cytoscape
    if not os.path.exists('{0}/results/cytoscape/'.format(outputdir)):
        sh('mkdir {0}/results/cytoscape'.format(outputdir))
    if genome.lower()=='hg19' or genome.lower()=='hg38':
        sh('{1}/lib/DEG2network.py -p 0.05 -n 5 -k\
         {1}/lib/db/merged_KEGG.txt -i {0}/results/Treat_vs_control_diff.txt -d \
           {0}/results/cytoscape'.format(outputdir,softpath))
    elif genome.lower()=='mm9' or genome.lower()=='mm10':
        sh('{1}/lib/DEG2network_mouse.py -p 0.05 -n 5 -k\
         {1}/lib/db/mouse_merged_KEGG.txt -i {0}/results/Treat_vs_control_diff.txt -d \
           {0}/results/cytoscape'.format(outputdir,softpath))

    ###! find is there other species


