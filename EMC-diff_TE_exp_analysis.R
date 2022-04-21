###     Author:     Xun Chen
###     Date:       2020/5/20
###     Contact:    xunchen85@gmail.com

##### load pacakges
library(fitdistrplus)
library(gplots)
library(RColorBrewer)
library(devtools)
library(MASS)
library(ggplot2)
library(grid)
library(ggrepel)
library(dplyr)
library(splitstackshape)
library(logspline)
library(DESeq2, quietly=T)

###########################
####### volc plot function
makeVolc<-function(inputData,outputFile){
  mutateddf <- mutate(inputData, sig=ifelse(inputData$padj<0.05, "padj<0.05", "Not Sig")) #Will have different colors depending on significance
  input <- cbind(gene=rownames(mutateddf), mutateddf ) #convert the rownames to a column
  input$gene<-rownames(inputData)
  p<-ggplot(input, aes(log2FoldChange, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
    geom_point(aes(col=sig),na.rm=TRUE) + #add points colored by significance
    scale_color_manual(values=c("black", "red")) + 
    ggtitle(outputFile) + #e.g. 'Volcanoplot DESeq2'
    geom_text_repel(data=head(input, 20), aes(label=head(input$gene,20))) + #adding text for the top 20 genes
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          
          axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
          axis.title=element_text(colour="black",size=rel(1.5)),
          #        legend.position=c(0.8,0.8),
          #        legend.position="bottom",
          legend.position="right",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1.2)))
  pdf(paste(outputFile,".pdf",sep=""),                      
      width = 6,
      height = 4.5,
      pointsize = 10)        # smaller font size
    grid.draw(p)
  dev.off()
}

###############################
####### Reading data
inputFile = "EMC_RNAseq_hg19.rename.CountM"

data<-read.delim2(paste("inputs/",inputFile,sep=""),header=T)
row.names(data)<-data[,1]
data<-data[,-1]
data2<-data[,c(colnames(data)[grepl("Flu",colnames(data))])]
data3<-data[,c(colnames(data)[!grepl("Flu",colnames(data))])]
TGroupLength<-ncol(data2)
CGroupLength<-ncol(data3)
data<-cbind(data2,data3)
groups <-factor(ifelse(grepl("_Flu",colnames(data)),"TGroup","CGroup"))  

####### DESEQ2 analyses, standadization
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
datamatrix <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
datamatrix$condition = relevel(datamatrix$groups,"CGroup")
datamatrix <- estimateSizeFactors(datamatrix)
datamatrix <- estimateDispersions(datamatrix)
datamatrix <- nbinomWaldTest(datamatrix)
save(datamatrix,file=paste(inputFile,"datamatrix_normalized",sep="_"))
tDataNorm <- as.data.frame(counts(datamatrix, normalized=TRUE))
tDataNorm2 <- data.frame(tDataNorm)

####### differential expression
res <- results(datamatrix,independentFiltering=F)
write.table(res, file=paste(inputFile,"pairedEnd_test_gene_TE_analysis.txt",sep="_"), sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) & (abs(res$log2FoldChange)> 0.000000)), ]
write.table(resSig, file=paste(inputFile,"pairedEnd_test_sigdiff_gene_TE.txt",sep="_"),sep="\t", quote=F)

