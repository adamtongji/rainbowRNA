#!/usr/bin/env python
import os
import sys

resultdir = sys.argv[1]
webdir = sys.argv[2]
sh = os.system

def load_fig():
    sh("cp {}/results/pca.png {}/static/pic/images/figure3.5.5.png".format(resultdir, webdir))
    sh("cp {}/results/Treat_vs_control_diff.png {}/static/pic/images/figure3.5.7.png".format(resultdir, webdir))
    sh("cp {}/results/up/gene_up_go.png {}/static/pic/images/figure4.1.2.png".format(resultdir, webdir))
    sh("cp {}/results/down/gene_down_go.png {}/static/pic/images/figure4.1.1.png".format(resultdir, webdir))
    sh("cp {}/results/up/gene_up_kegg.png {}/static/pic/images/figure4.2.2.png".format(resultdir, webdir))
    sh("cp {}/results/down/gene_down_kegg.png {}/static/pic/images/figure4.2.1.png".format(resultdir, webdir))
    sh("cp {}/results/down/sample.png {}/static/pic/images/figure4.3.png".format(resultdir, webdir))
    sh("a=`ls {0}/final/phase4-GSEA/BP/enplot_*.png`;for i in $a; cp $i \
    {1}/static/pic/images/figure4.4.png;break;done ".format(resultdir, webdir))


def load_table():
    outf = [i.rstrip() for i in open("{}/templates/mapping_summary.html".format(webdir))]
    with open("{}/templates/mapping_summary.html".format(webdir), "w") as fo:
        table1 = [i.rstrip().split("\t") for i in open("{}/results/mapping_rate/mapping_rate_summary.txt" \
                                                       .format(resultdir))]
        outtext1 = process_head(table1[0]) + process_body(table1[1:])
        outf2 = outf[:30]+[outtext1]+outf[280:]
        for line in outf2:
            print >> fo, line

    outf=[i.rstrip() for i in open("{}/templates/expression.html".format(webdir))]
    with open("{}/templates/expression.html".format(webdir),"w") as fo:
        table1=[i.rstrip().split("\t") for i in open("{}/results/treat_and_control_expression.txt"\
                                                 .format(resultdir))]
        outtext1=process_head(table1[0])+process_body(table1[1:11])
        table2=[i.rstrip().split("\t")[:7] for i in open("{}/results/Treat_vs_control_diff.txt"\
                                                 .format(resultdir))]
        outtext2 = process_head(table2[0]) + process_body(table2[1:11])
        table3=[i.rstrip().split("\t")[:1]+i.rstrip().split("\t")[7:] for i in open("{}/results/Treat_vs_control_diff.txt"\
                                                 .format(resultdir))]
        outtext3 = process_head(table3[0]) + process_body(table3[1:11])
        outf2=outf[:26]+[outtext1]+outf[123:140]+[outtext2]+[outtext3]+outf[261:]
        for line in outf2:
            print >> fo, line

    outf = [i.rstrip() for i in open("{}/templates/function.html".format(webdir))]
    with open("{}/templates/expression.html".format(webdir), "w") as fo:
        table1 = [i.rstrip().split("\t")[:7] for i in open("{}/results/up/genes_up_go.xls" \
                                                       .format(resultdir))]
        outtext1 ="""<table class="table table-hover table-bordered">
          <thead style="background-color:#67a9cf; font-color:black;">
        <tr><th colspan="7" style="text-align:center">Up-regulated Genes Gene Ontology Analysis</th>
        </tr></thead>
        """+process_head2(table1[0])+process_body(table1[1:6])

        table2 = [i.rstrip().split("\t")[:7] for i in open("{}/results/down/genes_down_go.xls" \
                                                       .format(resultdir))]
        outtext2 = """<table class="table table-hover table-bordered">
            <thead style="background-color:#67a9cf; font-color:black;">
           <tr><th colspan="7" style="text-align:center">Down-regulated Genes Gene Ontology Analysis</th>
               </tr></thead>
               """ + process_head2(table2[0]) + process_body(table2[1:6])

        table3 = [i.rstrip().split("\t")[:7] for i in open("{}/results/up/genes_up_kegg.xls" \
                                                       .format(resultdir))]
        outtext3 = """<table class="table table-hover table-bordered">
         <thead style="background-color:#67a9cf; font-color:black;">
        <tr><th colspan="7" style="text-align:center">Up-regulated Genes KEGG Pathway Analysis</th>
            </tr></thead>
            """ + process_head2(table3[0]) + process_body(table3[1:6])

        table4 = [i.rstrip().split("\t")[:7] for i in open("{}/results/down/genes_down_kegg.xls" \
                                                           .format(resultdir))]
        outtext4 = """<table class="table table-hover table-bordered">
         <thead style="background-color:#67a9cf; font-color:black;">
        <tr><th colspan="7" style="text-align:center">Down-regulated Genes KEGG Pathway Analysis</th>
            </tr></thead>
            """ + process_head2(table4[0]) + process_body(table4[1:6])
        outf2 = outf[:30]+[outtext1]+[outtext2]+outf[141:166]+ [outtext3]+ [outtext4]+outf[273:]
        for line in outf2:
            print >> fo, line


def process_head(inheads):
    start="""
    <table class="table table-hover table-bordered">
        <thead style="background-color:#67a9cf; font-color:black;">
    <tr>
    """
    mid= "<th><b>" +"</b></th><th><b>".join(inheads)
    tail = "</b></th></tr></thead>"
    return start+mid+tail


def process_head2(inheads):
    start="""
        <thead style="background-color:#67a9cf; font-color:black;">
    <tr>
    """
    mid= "<th><b>" +"</b></th><th><b>".join(inheads)
    tail = "</b></th></tr></thead>"
    return start+mid+tail

def process_body(intab):
    start="<tbody>"
    tmp=["<tr><td>"+"</td><td>".join(i)+"</td></tr>" for i in intab]
    mid="".join(tmp)
    tail_tmp=["..." for i in range(len(intab[0]))]
    tail="<tr><td>"+"</td><td>".join(tail_tmp)+"</td></tr>"+"</tbody></table>"
    return start + mid + tail


def update_html():
    load_fig()
    load_table()


