library(clusterProfiler)

args<-commandArgs(T)
outresultdir = args[1]
setwd(outresultdir)
filename="genes"
file_down=paste(filename,"down.txt",sep="_")
file_down_kegg=paste(filename,"down_kegg.xls",sep="_")
file_down_go=paste(filename,"down_go.xls",sep="_")

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
