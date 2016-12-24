rm(list=ls())
library(ggplot2)
library(reshape2)
library(pheatmap)
library(pathview)
library(Category)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GOstats)
library(KEGG.db)
library(GO.db)

args<-commandArgs(T)
outresultdir = args[1]
setwd(outresultdir)
filename="genes"
file_down=paste(filename,"down.txt",sep="_")
file_down_kegg=paste(filename,"down_kegg.xls",sep="_")
file_down_go=paste(filename,"down_go.xls",sep="_")


data=read.table(file_down,header=F)
setwd('./down')
a=t(data)
b=as.vector(a)
genes=b

#entrezIDs <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)


## gene ontology background

goAnn <- get("org.Mm.egGO")
universe <- Lkeys(goAnn)


params <- new("GOHyperGParams",
              geneIds=entrezIDs,
              universeGeneIds=universe,
              annotation="org.Mm.eg.db",
              ontology="BP",
              pvalueCutoff=0.05,
              testDirection="over")


over <- hyperGTest(params)
bp <- summary(over)

glist <- geneIdsByCategory(over)

glist <- sapply(glist, function(.ids) {
  .sym <- mget(.ids, envir=org.Mm.egSYMBOL, ifnotfound=NA)
  .sym[is.na(.sym)] <- .ids[is.na(.sym)]
  paste(.sym, collapse=";")
})

bp$Symbols <- glist[as.character(bp$GOBPID)]

write.table(bp,file_down_go,col.names=T,row.names=F,sep="\t")



## kegg enrichment background; do kegg,p<0.05
# keggAnn <- get("org.Mm.egPATH")
keggAnn <- get("org.Mm.egPATH")
universe <- Lkeys(keggAnn)

params <- new("KEGGHyperGParams",
              geneIds=entrezIDs,
              universeGeneIds=universe,
              annotation="org.Mm.eg.db",
              categoryName="KEGG",
              pvalueCutoff=0.05,
              testDirection="over")

over <- hyperGTest(params)
kegg <- summary(over)

glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {
  .sym <- mget(.ids, envir=org.Mm.egSYMBOL, ifnotfound=NA)
 	.sym[is.na(.sym)] <- .ids[is.na(.sym)]
 	paste(.sym, collapse=";")
})

kegg$Symbols <- glist[as.character(kegg$KEGGID)]

write.table(kegg,file_down_kegg,col.names=T,row.names=F,sep="\t")

gIds <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound=NA)
gEns <- unlist(gIds)
gene.data <- rep(1, length(gEns))
names(gene.data) <- gEns

for(i in 1:nrow(kegg)){
  pv.out <- pathview(gene.data, pathway.id=as.character(kegg$KEGGID)[i], species="mouse",
                     kegg.native=T,low = list(gene = "green"),
                     mid =list(gene = "gray"), high = list(gene = "green"),plot.col.key=FALSE)
}


library(ggplot2)
bp = read.table("genes_down_go.xls", header=T,row.names=1, sep="\t")


bp_sort=bp[order(bp$Pvalue),]
if (nrow(bp_sort)>=20){
    tmp=bp_sort[1:20,]
}else{
    tmp=bp_sort
}

tmp=tmp[order(tmp$Pvalue,decreasing = T),]
tmp$Term=factor(tmp$Term,levels = tmp$Term)
p<-ggplot(tmp)+geom_bar(aes(x=Term, y=-log(Pvalue,10)),width=0.6, stat="identity")+theme(axis.text.x=element_text(angle = 90,hjust=1))+coord_flip()+theme_bw()+theme(panel.grid.major=element_line(colour=NA))+xlab("GO Name")+ylab("-log10(p value)")

goplotname="GOplot_down.pdf"
ggsave(p,filename =goplotname,width=8)




kg = read.table("genes_down_kegg.xls", header=T,row.names=1, sep="\t")


kg_sort=kg[order(kg$Pvalue),]
if (nrow(kg_sort)>=20){
    tmp=kg_sort[1:20,]
}else{
    tmp=kg_sort
}
tmp=tmp[order(tmp$Pvalue,decreasing = T),]
tmp$Term=factor(tmp$Term,levels = tmp$Term)
p<-ggplot(tmp)+geom_bar(aes(x=Term, y=-log(Pvalue,10)),width=0.6, stat="identity")+theme(axis.text.x=element_text(angle = 90,hjust=1))+coord_flip()+theme_bw()+theme(panel.grid.major=element_line(colour=NA)) +xlab("KEGG Name")+ylab("-log10(p value)")

goplotname="KEGGplot_down.pdf"
ggsave(p,filename =goplotname,width=8)



setwd('../')


######## up genes #####
filename="genes"
file_up=paste(filename,"up.txt",sep="_")
file_up_kegg=paste(filename,"up_kegg.xls",sep="_")
file_up_go=paste(filename,"up_go.xls",sep="_")


data=read.table(file_up,header=T)
setwd('./up')
a=t(data)
b=as.vector(a)
genes=b

# entrezIDs <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)


## gene ontology background

# goAnn <- get("org.Mm.egGO")
goAnn <- get("org.Mm.egGO")
universe <- Lkeys(goAnn)

params <- new("GOHyperGParams",
              geneIds=entrezIDs,
              universeGeneIds=universe,
              annotation="org.Mm.eg.db",
              ontology="BP",
              pvalueCutoff=0.05,
              testDirection="over")

over <- hyperGTest(params)
bp <- summary(over)

glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {
  .sym <- mget(.ids, envir=org.Mm.egSYMBOL, ifnotfound=NA)
  .sym[is.na(.sym)] <- .ids[is.na(.sym)]
  paste(.sym, collapse=";")
})

bp$Symbols <- glist[as.character(bp$GOBPID)]

write.table(bp,file_up_go,col.names=T,row.names=F,sep="\t")



## kegg enrichment background; do kegg,p<0.05
keggAnn <- get("org.Mm.egPATH")
universe <- Lkeys(keggAnn)

params <- new("KEGGHyperGParams",
              geneIds=entrezIDs,
              universeGeneIds=universe,
              annotation="org.Mm.eg.db",
              categoryName="KEGG",
              pvalueCutoff=0.05,
              testDirection="over")

over <- hyperGTest(params)
kegg <- summary(over)

glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {
  .sym <- mget(.ids, envir=org.Mm.egSYMBOL, ifnotfound=NA)
 	.sym[is.na(.sym)] <- .ids[is.na(.sym)]
 	paste(.sym, collapse=";")
})

kegg$Symbols <- glist[as.character(kegg$KEGGID)]

write.table(kegg,file_up_kegg,col.names=T,row.names=F,sep="\t")

gIds <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound=NA)
gEns <- unlist(gIds)
gene.data <- rep(1, length(gEns))
names(gene.data) <- gEns

for(i in 1:nrow(kegg)){
  pv.out <- pathview(gene.data, pathway.id=as.character(kegg$KEGGID)[i], species="mouse",
                     kegg.native=T,low = list(gene = "red"),
                     mid =list(gene = "gray"), high = list(gene = "red"),plot.col.key=FALSE)
}


library(ggplot2)
bp = read.table("genes_up_go.xls", header=T,row.names=1, sep="\t")


bp_sort=bp[order(bp$Pvalue),]
if (nrow(bp_sort)>=20){
    tmp=bp_sort[1:20,]
}else{
    tmp=bp_sort
}
tmp=tmp[order(tmp$Pvalue,decreasing = T),]
tmp$Term=factor(tmp$Term,levels = tmp$Term)
p<-ggplot(tmp)+geom_bar(aes(x=Term, y=-log(Pvalue,10)),width=0.6, stat="identity")+theme(axis.text.x=element_text(angle = 90,hjust=1))+coord_flip()+theme_bw()+theme(panel.grid.major=element_line(colour=NA)) +xlab("GO Name")+ylab("-log10(p value)")

goplotname="GOplot_up.pdf"
ggsave(p,filename =goplotname,width=8)



kg = read.table("genes_up_kegg.xls", header=T,row.names=1, sep="\t")


kg_sort=kg[order(kg$Pvalue),]
if (nrow(kg_sort)>=20){
    tmp=kg_sort[1:20,]
}else{
    tmp=kg_sort
}
tmp=tmp[order(tmp$Pvalue,decreasing = T),]
tmp$Term=factor(tmp$Term,levels = tmp$Term)
p<-ggplot(tmp)+geom_bar(aes(x=Term, y=-log(Pvalue,10)),width=0.6, stat="identity")+theme(axis.text.x=element_text(angle = 90,hjust=1))+coord_flip()+theme_bw()+theme(panel.grid.major=element_line(colour=NA)) +xlab("KEGG Name")+ylab("-log10(p value)")

goplotname="KEGGplot_up.pdf"
ggsave(p,filename =goplotname,width=8)
