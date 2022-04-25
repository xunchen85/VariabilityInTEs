### Author:     Xun Chen
### Date:       2021/11/15
### Contact:    xunchen85@gmail.com

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
library(ggpmisc)
library(gridExtra)
library(NMF)
library(cowplot)
library(ggpubr)
library(DESeq2)
library(PCAtools)

################# step 1 load raw count
## 1. RNA raw counts
TE_exp_raw = read.delim("inputs/EMC_RNAseq_hg19.rename.CountM",header=T,sep="")

################# step 2 normalize
Conditions = c("RNA_gene","RNA_TE")
Condition = "RNA_gene"
Condition = "RNA_TE"

## 
for (Condition in Conditions){
  if (Condition == "RNA_gene"){
    TE_exp_raw_matrix = TE_exp_raw
    rownames(TE_exp_raw_matrix) = TE_exp_raw$gene.TE
    TE_exp_raw_matrix = TE_exp_raw_matrix[,-1]
    data = as.matrix(TE_exp_raw_matrix)
    groups <-factor(ifelse(grepl("_Flu",colnames(data)),"TGroup","CGroup"))
    min_read <- 1
    data <- data[apply(data,1,function(x){max(x)}) > min_read,]
    sampleInfo <- data.frame(groups,row.names=colnames(data))
    datamatrix <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
    datamatrix$condition = relevel(datamatrix$groups,"CGroup")
    datamatrix <- estimateSizeFactors(datamatrix)
    datamatrix <- estimateDispersions(datamatrix)
    datamatrix <- nbinomWaldTest(datamatrix)
    vst <- assay(vst(datamatrix))
    vst = vst[!grepl(":",rownames(vst)),]
    
  } else if (Condition == "RNA_TE"){
    TE_exp_raw_matrix = TE_exp_raw
    rownames(TE_exp_raw_matrix) = TE_exp_raw$gene.TE
    TE_exp_raw_matrix = TE_exp_raw_matrix[,-1]
    data = as.matrix(TE_exp_raw_matrix)
    groups <-factor(ifelse(grepl("_Flu",colnames(data)),"TGroup","CGroup"))
    min_read <- 1
    data <- data[apply(data,1,function(x){max(x)}) > min_read,]
    sampleInfo <- data.frame(groups,row.names=colnames(data))
    datamatrix <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
    datamatrix$condition = relevel(datamatrix$groups,"CGroup")
    datamatrix <- estimateSizeFactors(datamatrix)
    datamatrix <- estimateDispersions(datamatrix)
    datamatrix <- nbinomWaldTest(datamatrix)
    vst <- assay(vst(datamatrix))
    vst = vst[grepl(":",rownames(vst)),]
  }
  
  ## prepare info table
  sampleInfo = data.frame(colnames(data))
  colnames(sampleInfo)[1] = "sampleID"
  sampleInfo$sampleID2 = sampleInfo$sampleID
  sampleInfo = data.frame(cSplit(sampleInfo,"sampleID2",sep="_",type.convert = as.character))
  sampleInfo$Ethnicity = ifelse(grepl("AF",sampleInfo$sampleID),"AF","EU")
  ### correct sample ethnicity for AF22, AF36 and AF38
  sampleInfo$Ethnicity = ifelse(grepl("AF22",sampleInfo$sampleID) | grepl("AF38",sampleInfo$sampleID) | grepl("AF36",sampleInfo$sampleID),"EU",sampleInfo$Ethnicity)
  colnames(sampleInfo) = c("sampleID","individualID","Flu_status","Ethnicity")
  rownames(sampleInfo) = as.character(sampleInfo$sampleID)
  sampleInfo = sampleInfo[,-1]
  
  ## plot     
  p <- pca(vst,metadata = sampleInfo,removeVar = 0.1)
  pdf(paste("EMC-Figure1A-PCA_",Condition,".pdf",sep=""),                      
      width = 6,
      height = 5,
      pointsize = 10)        # smaller font size
    biplot(p,
         showLoadings = FALSE,
         colby = 'Flu_status', colkey = c('Flu'='#99000d', 'NI'='#4d4d4d'),
         shape = 'Ethnicity', shapekey = c('EU'=15, 'AF'=17),
         #hline = 0, 
         #vline = 0,
         lab = NULL,
         pointSize = 3,
         drawConnectors = FALSE,
         gridlines.major = FALSE, 
         gridlines.minor = FALSE,
         legendPosition = 'right')
  dev.off()
  
}

