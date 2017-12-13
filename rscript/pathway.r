library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(org.Hs.eg.db)
library(pathview)
rm(list = ls())
args<-commandArgs(T)
outresultdir = args[1]
setwd(outresultdir)
filedown<-read.table("genes_down.txt",sep='\t')
filedown<-as.matrix(filedown)

setwd("./down")
gene.down<-bitr(filedown, fromType='SYMBOL', toType=c('ENSEMBL',"ENTREZID"),OrgDb = org.Hs.eg.db)

down_go<- enrichGO(gene         = gene.down$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) #

write.table(down_go,"genes_down_go.xls",sep='\t',quote=F,col.names = T,row.names = F)
p<-dotplot(down_go)
ggsave(p,file="genes_down_go.pdf",width=8)
ggsave(p,file="genes_down_go.png",width=8)

keggs<-bitr_kegg(gene.down$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = "hsa")
down_kegg <- enrichKEGG(gene         = gene.down$ENTREZID,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
p<-dotplot(down_kegg)
ggsave(p,file="genes_down_kegg.pdf",width=8)
ggsave(p,file="genes_down_kegg.png",width=8)

write.table(down_kegg,"genes_down_kegg.xls",sep='\t',quote=F,col.names = T,row.names = F)


for (kegg_id in down_kegg$ID){
  pv <-pathview(gene.data = gene.down[,3], pathway.id = kegg_id,species = "hsa",kegg.native=T,
                low = list(gene = "green"),mid =list(gene = "gray"),  gene.idtype = "ENTREZID",
                high = list(gene = "green"),plot.col.key=FALSE)
}

setwd("../")
fileup<-read.table("genes_up.txt",sep='\t')
fileup<-as.matrix(fileup)
setwd("./up")
gene.up<-bitr(fileup, fromType='SYMBOL', toType=c('ENSEMBL',"ENTREZID"),OrgDb = org.Hs.eg.db)

up_go<- enrichGO(gene         = gene.up$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05) #  keytype       = 'SYMBOL',

write.table(up_go,"genes_up_go.xls",sep='\t',quote=F,col.names = T,row.names = F)
p<-dotplot(up_go)
ggsave(p,file="genes_up_go.pdf",width=8)
ggsave(p,file="genes_up_go.png",width=8)

keggs<-bitr_kegg(gene.up$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = "hsa")
up_kegg <- enrichKEGG(gene         = gene.up$ENTREZID,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)

write.table(up_kegg,"genes_up_kegg.xls",sep='\t',quote=F,col.names = T,row.names = F)
p<-dotplot(up_kegg)
ggsave(p,file="genes_up_kegg.pdf",width=8)
ggsave(p,file="genes_up_kegg.png",width=8)

fileup<-data.frame(fileup)
rownames(fileup)<-fileup[,1]
for (kegg_id in up_kegg$ID){
  pv <-pathview(gene.data = gene.up[,3], pathway.id = kegg_id,species = "hsa",kegg.native=T,
                low = list(gene = "red"),mid =list(gene = "gray"),  gene.idtype = "ENTREZID",
                high = list(gene = "red"),plot.col.key=FALSE)
}
args=paste("cp ",up_kegg$ID[1],".pathview.png sample.png",sep='')
system(args)


