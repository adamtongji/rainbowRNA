rm(list=ls())
library(DESeq2)
library(reshape2)
library(pheatmap)

args<-commandArgs(T)
outexprdir = args[1]
pairtype = args[2]
p_value = args[3]
setwd(outexprdir)
pvalue<-as.numeric(p_value)
file_name = dir()
files = lapply(file_name, function(i) read.table(i,sep='\t', row.names=1))
control_len = as.integer(0.5*length(files))
control=files[[1]]
mycol=paste("control",1,sep='')
colnames(control) = mycol
if (control_len >1){
    for (i in c(2:control_len)){
        mycol=paste("control",i,sep='')
        temp = files[[i]]
        colnames(temp) = mycol
        control = cbind(control, temp)
    }
}
control = control[-c((nrow(control)-4):nrow(control)),]
treat = files[[control_len+1]]
treat_len = length(file_name)-control_len
mycol=paste("treat",1,sep='')
colnames(treat) = mycol
if (treat_len >1){
    for  (i in c((control_len+2):length(file_name))){
        mycol=paste("treat",i-control_len,sep='')
        temp = files[[i]]
        colnames(temp) = mycol
        treat = cbind(treat, temp)
    }
}
treat = treat[-c((nrow(treat)-4):nrow(treat)),]
countdata <-cbind(treat, control)
write.table(countdata, '../results/treat_and_control_expression.txt',row.names=T,col.names=T, quote=F, sep='\t')
condition = c(rep('treat',treat_len), rep('control', control_len))
#
#paired_end = pattern

if (pairtype == 'True'){
    pattern = factor(rep(c(1:control_len),2))
    print("pattern");
    print(pattern)
    print("condition")
    print(condition)
    print("countdata")
    print(head(countdata))
    coldata<- data.frame(row.names = colnames(countdata), condition, pattern)
    dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition+pattern)
}else{
    coldata<- data.frame(row.names = colnames(countdata), condition)
    dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
}


dds <- DESeq(dds)

pca_dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(pca_dds, normalized=TRUE)),colData=colData(pca_dds))
pdf('../results/pca.pdf')
plotPCA( DESeqTransform(se))
dev.off()
png('../results/pca.png')
plotPCA( DESeqTransform(se))
dev.off()
res <- results(dds,contrast = c("condition","treat", "control"))
res <- res[order(res$pvalue),]
head(res)

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

rownames(resdata)<-resdata[,1]
resdata<-resdata[,2:dim(resdata)[2]]
diff<-resdata[! is.na(resdata[,5]),]
diff<-diff[diff[,5]<p_value,]
up<-diff[diff[,4]>0,]
down<-diff[diff[,4]<0,]

write.table(rownames(down),file="../results/genes_down.txt",quote=F,row.names = F,col.names = F)
write.table(rownames(up),file="../results/genes_up.txt",quote=F,row.names = F,col.names = F)

newresdata=data.frame(gene_id=rownames(resdata),resdata)
write.table(newresdata, file="../results/Treat_vs_control_diff.txt",sep="\t",quote=F,row.names = FALSE)


pdata<-newresdata[newresdata$pvalue<p_value,]
pdata<-pdata[!is.na(pdata$pvalue),]
km <- kmeans(pdata[,5],2)
pdata$cl <- km$cluster
pdata <- pdata[order(pdata$cl),]
ann <- data.frame(cluster=pdata$cl)
rownames(ann) <- rownames(pdata)
ann_col=data.frame(sample=c(rep('treat',treat_len),rep('control',control_len)))
rownames(ann_col)=colnames(pdata)[8:(7+length(file_name))]
pdata<-t(scale(t(pdata[,8:(7+length(file_name))])))
pdf('../results/Treat_vs_control_diff.pdf')
pheatmap(pdata[,1:length(file_name)],col=colorRampPalette(c('blue','white','red'))(50),scale='none',cluster_rows=F,cluster_cols=F,show_rownames=F,cellwidth=30,annotation_row=ann,annotation_col=ann_col)
dev.off()
png('../results/Treat_vs_control_diff.png')
pheatmap(pdata[,1:length(file_name)],col=colorRampPalette(c('blue','white','red'))(50),scale='none',cluster_rows=F,cluster_cols=F,show_rownames=F,cellwidth=30,annotation_row=ann,annotation_col=ann_col)
dev.off()


####################
