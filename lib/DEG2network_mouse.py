#!/usr/bin/env python
#coding:utf-8
from optparse import OptionParser

parser = OptionParser(usage = """python /Users/wangcy/Desktop/signal_network/DEG2network.py -h \
python /Users/wangcy/Desktop/signal_network/DEG2network.py -i /Users/wangcy/Desktop/signal_network/CS_MF_L_diff.txt -p 0.05 -n 5 -k /Users/wangcy/Desktop/signal_network/merged_KEGG.txt -d /Users/wangcy/Desktop/signal_network""")
parser.add_option("-i","--input",dest="input_DEG",type="str",\
                  help="""-i directs to your differentially expressed file whose prototype is generated from DEseq2.\
                   The original differentally expressed genes file given by Deseq2 should be open with excel and then\
                    re-save in the format with delimiter "\\t". e.g. The input format should be as follows if you read\
                     it in python 'Gene\\tbaseMean\\tlog2FoldChange\\tlfcSE\\tstat\\tpvalue\\tpadj\\n"\
                     'SURF4\\t27756.50\\t1.004\\t0.182\\t5.514\\t3.49e-08\\t0.000108\\n""")#input_DEG

parser.add_option("-p","--p_value",dest="p_value",type="float",default=0.05,help="""-p refers to the thresfold of\
 p-value for differential expression. Note that it is for p-value rather than adj p-value. The default is 0.05.\
  You may need to adjust it when the output covers too many genes, e.g, the number of total lines of one of the\
   output files, network.txt, is not larger than 2000""")
parser.add_option("-n","--number",dest="number",type="int",default=5,help="""The maximum number of hub genes in network\
 file for cytoscape plotting""")
parser.add_option("-k","--kegg_file",dest="kegg_file",type="str",help="""-k refers to the merged network file from\
 KEGG, named by "merged_KEGG.txt". """)
parser.add_option("-d","--diretory",dest="output_directory",type="str",help="""the network.txt should be loaded into\
 cytoscape along with the -i file telling the colors of each gene node according to the extent of
##differental expression. Besides, The network.txt is complemented by network_degree.txt which counts the indegree\
 and outdegree for each gene node""")

(options,args)=parser.parse_args()

def DEG_to_network():
    DEG_file=open(options.input_DEG)
    DEG_file.readline()
    DEG_list=[]
    DEG_dic={}


    while True:
        line=DEG_file.readline()
        if not line:break
        line=line.strip().split()
        if line[5]=="NA":break
        if float(line[5])<options.p_value:  ## p-value less than 0.05
            # DEG_list.append(line[0].upper())
            # DEG_dic[line[0].upper()]=line[2]
            DEG_list.append(line[0])
            DEG_dic[line[0]]=line[2]
        else:
            DEG_dic[line[0]]=line[2]

    #print DEG_dic
    ###################################
    output_a=open(options.output_directory+'/'+'whole_network.txt','w')
    output_b=open(options.output_directory+'/'+'network_degree.txt','w')
    output_c=open(options.output_directory+'/'+'selected_network.txt','w')
    #output_c.writelines("Start	End	Relationship	KEGG_ID	Pathway\n")
    output_b.writelines("gene	logFC	outdegree	indegree	outdegree_genes	indegree_genes\n")

    my_kegg_file=open(options.kegg_file)  ##
    output_a.writelines(my_kegg_file.readline())
    degree_dic={} ##gene:[]
    while True:
        line=my_kegg_file.readline()
        if not line:break
        tt=line.strip().split()
        if tt[0] in DEG_list or tt[1] in DEG_list:
            if (not DEG_dic.has_key(tt[0])) or (not DEG_dic.has_key(tt[1])):

                continue
            output_a.writelines(line)
            try:
                degree_dic[tt[0]][1]+=1
                degree_dic[tt[0]][3].append(tt[1])
            except KeyError:
                if DEG_dic.has_key(tt[0]):
                    degree_dic[tt[0]]=[DEG_dic[tt[0]],1,0,[tt[1]],[]]
                else:
                    degree_dic[tt[0]]=['NA',1,0,[tt[1]],[]]

            try:
                degree_dic[tt[1]][2]+=1
                degree_dic[tt[1]][4].append(tt[0])
            except KeyError:
                if DEG_dic.has_key(tt[1]):
                    degree_dic[tt[1]]=[DEG_dic[tt[1]],0,1,[],[tt[0]]]
                else:
                    degree_dic[tt[1]]=['NA',0,1,[],[tt[0]]]

        else:
            pass

    degree_list=degree_dic.items()
    degree_list.sort(key=lambda x: int(x[1][1])+int(x[1][2]),reverse=True)
    top_list=[v[0] for v in degree_list[:options.number]]
    for value in degree_list:
        #print value
        gene=value[0]
        logFC=value[1][0]
        outdegree=str(value[1][1])
        indegree=str(value[1][2])
        outdegree_string=",".join(value[1][3])
        indegree_string=",".join(value[1][4])
        output_b.writelines("\t".join([gene,logFC,outdegree,indegree,outdegree_string,indegree_string])+'\n')


    output_a.close()
    output_b.close()
    my_kegg_file.close()

    file=open(options.output_directory+'/'+'whole_network.txt')
    output_c.writelines(file.readline())
    while True:
        line=file.readline()
        if not line:break
        line=line.strip().split()
        if line[0] in top_list or line[1] in top_list:
            output_c.writelines("\t".join(line)+'\n')
    output_c.close()
    file.close()

DEG_to_network()
