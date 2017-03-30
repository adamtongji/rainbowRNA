#!/usr/bin/env python
#coding:utf-8
import os
import subprocess


def sh(args):
    return subprocess.call(args,shell=True)


def merge_mapping_rate(outputdir):
    rates=os.listdir("{0}/results/mapping_rate/".format(outputdir))
    temps = [i.rstrip().split('\t') for i in open("{0}/results/mapping_rate/{1}".format(outputdir, rates[0]))]
    prefix = rates[0].rstrip("Log.final.out")
    _total_mapping = [['',prefix]]+temps[5:]

    for _file in rates[1:]:
        temps = [i.rstrip().split('\t') for i in open("{0}/results/mapping_rate/{1}".format(outputdir, _file))]
        prefix = _file.rstrip("Log.final.out")
        for item in temps:
            if len(item)==1:
                item.append('')
        temp2 = [i[1] for i in temps]
        temp2 = temp2[5:]
        for i in range(len(_total_mapping)):
            if i == 0:
                _total_mapping[i].append(prefix)
            else:
                _total_mapping[i].append(temp2[i-1])
    _mapping_out = ['\t'.join(i) for i in _total_mapping]
    with open("{0}/results/mapping_rate/mapping_rate_summary.txt".format(outputdir),'w') as f:
        for line in _mapping_out:
            print >>f, line


def get_report(Outputdir):

    if not os.path.exists('{0}/final/'.format(Outputdir)):
        sh('mkdir {0}/final/ {0}/results/mapping_rate/'.format(Outputdir))
    else:
        sh('rm -rf {0}/final/ {0}/results/mapping_rate/'.format(Outputdir))
        sh('mkdir {0}/final/ {0}/results/mapping_rate/'.format(Outputdir))
    sh('mkdir {0}/final/phase1-AllExpGenes {0}/final/phase2-DiffExpGenes\
       {0}/final/phase3-GO_KEGG {0}/final/phase4-GSEA {0}/final/phase5-SignalNet'.format(Outputdir))

    # sh('cp {0}/mapping/*_log.txt {0}/results/mapping_rate/'.format(Outputdir))
    # STAR setting
    sh('cp {0}/mapping/*.final.out {0}/results/mapping_rate/'.format(Outputdir))
    merge_mapping_rate(Outputdir)
    ####

    sh('cp {0}/results/treat_and_control_expression.txt {0}/final/phase1-AllExpGenes/expression_level_all.xls'\
       .format(Outputdir))
    sh('cp {0}/results/genes_up.txt {0}/final/phase2-DiffExpGenes/genes_up.xls'\
       .format(Outputdir))
    sh('cp {0}/results/genes_down.txt {0}/final/phase2-DiffExpGenes/genes_down.xls'\
       .format(Outputdir))
    sh('cp {0}/results/Treat_vs_control_diff.txt {0}/final/phase2-DiffExpGenes/Treat_vs_control_diff.xls'\
       .format(Outputdir))
    sh('cp {0}/results/Treat_vs_control_diff.pdf {0}/final/phase2-DiffExpGenes/Treat_vs_control_diff.pdf'\
       .format(Outputdir))
    sh('cp -r {0}/results/down/ {0}/final/phase3-GO_KEGG/down/'.format(Outputdir)\
       .format(Outputdir))
    sh('cp -r {0}/results/up/ {0}/final/phase3-GO_KEGG/up/'.format(Outputdir)\
       .format(Outputdir))
    # sh('rm {0}/final/phase3-GO_KEGG/*/Rplots.pdf'.format(Outputdir))
    sh('cp -r {0}/results/gsea/bp*/ {0}/final/phase4-GSEA/BP'.format(Outputdir))
    sh('cp -r {0}/results/gsea/kegg*/ {0}/final/phase4-GSEA/KEGG'.format(Outputdir))
    sh('cp {0}/results/cytoscape/* {0}/final/phase5-SignalNet/'.format(Outputdir))


def get_ciri_report(Outputdir):
    if not os.path.exists('{0}/final/'.format(Outputdir)):
        sh('mkdir {0}/final/'.format(Outputdir))
    else:
        sh('rm -rf {0}/final/'.format(Outputdir))
        sh('mkdir {0}/final/'.format(Outputdir))
    sh('mkdir {0}/final/phase1-AllExpGenes {0}/final/phase2-DiffExpGenes\
       {0}/final/phase3-GO_KEGG {0}/final/phase4-GSEA {0}/final/phase5-SignalNet'.format(Outputdir))

    sh('cp {0}/results/treat_and_control_expression.txt {0}/final/phase1-AllExpGenes/expression_level_all.xls' \
       .format(Outputdir))
    sh('cp {0}/results/genes_up.txt {0}/final/phase2-DiffExpGenes/genes_up.xls' \
       .format(Outputdir))
    sh('cp {0}/results/genes_down.txt {0}/final/phase2-DiffExpGenes/genes_down.xls' \
       .format(Outputdir))
    sh('cp {0}/results/Treat_vs_control_diff.txt {0}/final/phase2-DiffExpGenes/Treat_vs_control_diff.xls' \
       .format(Outputdir))
    sh('cp {0}/results/Treat_vs_control_diff.pdf {0}/final/phase2-DiffExpGenes/Treat_vs_control_diff.pdf' \
       .format(Outputdir))
    sh('cp -r {0}/results/down/ {0}/final/phase3-GO_KEGG/down/'.format(Outputdir) \
       .format(Outputdir))
    sh('cp -r {0}/results/up/ {0}/final/phase3-GO_KEGG/up/'.format(Outputdir) \
       .format(Outputdir))
    sh('rm {0}/final/phase3-GO_KEGG/*/Rplots.pdf'.format(Outputdir))
    sh('cp -r {0}/results/gsea/bp*/ {0}/final/phase4-GSEA/BP'.format(Outputdir))
    sh('cp -r {0}/results/gsea/kegg*/ {0}/final/phase4-GSEA/KEGG'.format(Outputdir))
    sh('cp {0}/results/cytoscape/* {0}/final/phase5-SignalNet/'.format(Outputdir))
    sh("mkdir {0}/final/phase6-circRNA".format(Outputdir))
    sh("cp {0}/results/circ_up_RBP.txt {0}/final/phase6-circRNA/up_circRNA_with_RBP.txt".format(Outputdir))
    sh("cp {0}/results/circ_down_RBP.txt {0}/final/phase6-circRNA/down_circRNA_with_RBP.txt".format(Outputdir))
    sh("cp {0}/results/annotation_table.txt {0}/final/phase6-circRNA/circRNA_summary.xls".format(Outputdir))
    sh("cp {0}/results/annotation_table_circbase.txt {0}/final/phase6-circRNA/circRNA_circbase.xls".format(Outputdir))
    sh("cp {0}/results/targetscan/circ_up_high.txt {0}/final/phase6-circRNA/down_circRNA_miRNA_top.xls".format(Outputdir))
    sh("cp {0}/results/targetscan/circ_down_high.txt {0}/final/phase6-circRNA/up_circRNA_miRNA_top.xls".format(Outputdir))
    sh("cp {0}/results/targetscan/circ_up_pass.txt {0}/final/phase6-circRNA/down_circRNA_miRNA_okay.xls".format(
        Outputdir))
    sh("cp {0}/results/targetscan/circ_down_pass.txt {0}/final/phase6-circRNA/up_circRNA_miRNA_okay.xls".format(
        Outputdir))
    sh("cp {0}/results/targetscan/circ_up_miRNA.txt {0}/final/phase6-circRNA/down_circRNA_miRNA_all.xls".format(
        Outputdir))
    sh("cp {0}/results/targetscan/circ_down_miRNA.txt {0}/final/phase6-circRNA/up_circRNA_miRNA_all.xls".format(
        Outputdir))