###     Author:     Xun Chen
###     Date:       2022/3/27
###     Contact:    xunchen85@gmail.com

library(ggplot2)
library(grid)
library(dplyr)
library(ggpmisc)
library(splitstackshape)
library(logspline)
library(gridExtra)
library(ggrepel)
library("pROC")
library("MASS")
library("ggpubr")
library(FSA)
library(psych)
library(PerformanceAnalytics)

options(scipen = 999)
theme_1 = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(hjust = 0.5),
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

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

####################################
####### Set up enviroment;
colorset = c("Flu" = "#99000d","NI" = "#4d4d4d")
setwd(dir="/Users/xchen/OneDrive - Kyoto University/Projects_Bourque/EMC_project_2020_12_1/Final_Rscripts_2022_4_20/")

####### variables
my.formula = y~x

gene_foldChange_threadhold = 1
pValue = 0.05
################################################ section 1: 2021_10_3 start

## sample information
Viralload = read.csv("inputs/H1N1_Viral_Load_5_27.csv")
Viralload$individualID = gsub("_Flu|_NI","",Viralload$Sample_ID)

## TE list
TElist = read.delim("inputs/TElist.sep.27",header=T,sep="")

## log2FC
FoldChange = read.delim("inputs/EMC_RNAseq_hg19.rename.CountM_pairedEnd_test_gene_TE_analysis.txt",sep="",header=T)
FoldChange$is_sig = ifelse(FoldChange$log2FoldChange>=1 & FoldChange$padj<=pValue,"sig_up","not_sig")
FoldChange$is_sig = ifelse(FoldChange$log2FoldChange>=1 & FoldChange$padj<=pValue,"sig_down",FoldChange$is_sig)
Gene_foldChange = FoldChange[!grepl(":",rownames(FoldChange)),]
Gene_foldChange$geneID = rownames(Gene_foldChange)
TE_foldChange = FoldChange[grepl(":",rownames(FoldChange)),]
rm(FoldChange)

############################################### section 2: variables
################## TE expression (normalized counts, ratio)
TE_exp_normalized = read.delim("inputs/EMC_RNAseq_hg19.rename.CountM_Normalized_PvalueOutput.txt",header=T,sep="",row.names = NULL)
colnames(TE_exp_normalized)[1] = "gene.TE"
TE_exp_normalized = TE_exp_normalized[,1:79]
TE_exp = TE_exp_normalized
Sum1 = data.frame(apply(TE_exp[,2:ncol(TE_exp)],2,sum))
Sum1$sample=rownames(Sum1)
colnames(Sum1)[1] = "totalReads"
Sum2 = data.frame(apply(TE_exp[grepl(":",TE_exp$gene.TE),2:ncol(TE_exp)],2,sum))
Sum2$sample=rownames(Sum2)
colnames(Sum2)[1] = "TEReads"
SumTable = merge(Sum1,Sum2,by="sample")
SumTable$TE_transcript = SumTable$TEReads/SumTable$totalReads
SumTable$individualID = gsub("_Flu","",SumTable$sample)
SumTable$individualID = gsub("_NI","",SumTable$individualID)
SumTable = merge(SumTable,Viralload,by="individualID",all.x=T)
TE_transcript = SumTable[grepl("NI",SumTable$sample),]
TE_transcript$TE_transcript_perc = TE_transcript$TE_transcript*100
TE_transcript$ViralLoad_perc = TE_transcript$ViralLoad/100
rm(SumTable,Sum1,Sum2,TE_exp,TE_exp_normalized) 

################## 2. IFN genes
immuneGenes_IFN = read.csv("inputs/EMC-predictorModel_2021_5_10_Genes_IFN.csv")
immuneGenes_IFN = immuneGenes_IFN[,2:41]
colnames(immuneGenes_IFN)[2:40] = paste("IFN_",colnames(immuneGenes_IFN)[2:40],sep="")
immuneGenes_IFN$IFN_aggregated = apply(immuneGenes_IFN[,2:ncol(immuneGenes_IFN)],1,sum)
immuneGenes_IFN$IFN_median = apply(immuneGenes_IFN[,2:(ncol(immuneGenes_IFN)-1)],1,median)
immuneGenes_IFN$sampleID = paste(as.character(immuneGenes_IFN$individualID),"_NI",sep="")

################## All genes (TPM)
# Gene name
hg19_geneName = read.delim("inputs/hg19.ensGene.Genename",header=F,sep="")
colnames(hg19_geneName) = c("transcriptID","geneID","chr","direction","start","end","geneName")
hg19_geneName = hg19_geneName[!duplicated(hg19_geneName[,c("geneID","geneName")]),]

# individual level expression change (TPM)
Both_exp_tpm = read.csv("inputs/EMC_RNAseq_hg19.rename.CountM-TPM.csv")
rownames(Both_exp_tpm) = Both_exp_tpm$name
Both_exp_tpm = Both_exp_tpm[,c(-1,-2)]
Gene_exp_tpm_global = Both_exp_tpm[!grepl(":",rownames(Both_exp_tpm)),c(1:78)]
Gene_exp_tpm_global = Gene_exp_tpm_global[,colnames(Gene_exp_tpm_global) %in% c(paste(Viralload$individualID,"_NI",sep=""),paste(Viralload$individualID,"_Flu",sep=""))]
Gene_exp_tpm_global$meanFlu = apply(Gene_exp_tpm_global[grepl("_Flu",colnames(Gene_exp_tpm_global))],1,mean,na.rm=T)
Gene_exp_tpm_global$meanNI = apply(Gene_exp_tpm_global[grepl("_NI",colnames(Gene_exp_tpm_global))],1,mean,na.rm=T)
Gene_exp_tpm_global$geneID = rownames(Gene_exp_tpm_global)
Gene_exp_tpm_global = merge(Gene_exp_tpm_global,hg19_geneName,by.x="geneID",all.x=T)
Gene_exp_tpm_global$is_highExp = ifelse(Gene_exp_tpm_global$meanFlu>=1 | Gene_exp_tpm_global$meanNI >=1,"highExp","lowExp")
rm(Both_exp_tpm,hg19_geneName)

### TRIM28 and SETDB1
TRIM28_SETDB1 = Gene_exp_tpm_global[Gene_exp_tpm_global$geneName %in% c("TRIM28","SETDB1"),]

### All TFs
TFs = read.csv("inputs/gProfiler_inputList_2022_3_8.csv")
TFs$is_existed = ifelse(TFs$Motif %in% Gene_exp_tpm_global$geneName,"Matched","Missed")
colnames(TFs) = c("TE_cluster","geneName","is_existed")

############################################## Section 3: determine factors
### correlation of TFs and TRIM28 and SETDB1
Gene_correlation = Gene_exp_tpm_global[Gene_exp_tpm_global$geneName %in% c(TFs$geneName,"TRIM28","SETDB1"),]
Gene_correlation$r.squared_NI = NA
Gene_correlation$r.adj_NI = NA
Gene_correlation$pValue_NI = NA
Gene_correlation$direction_NI = NA
i = 1
for (i in 1:nrow(Gene_correlation)){
  Gene_correlation_sub = data.frame(t(Gene_correlation[i,grepl("_NI$",colnames(Gene_correlation))]))
  colnames(Gene_correlation_sub)[1] = "TPM"
  Gene_correlation_sub$individualID = gsub("_Flu","",rownames(Gene_correlation_sub))
  Gene_correlation_sub$individualID = gsub("_NI","",Gene_correlation_sub$individualID)
  
  Gene_correlation_sub = merge(Gene_correlation_sub,TE_transcript,by="individualID",all=T)
  Gene_correlation_sub=Gene_correlation_sub[!is.na(Gene_correlation_sub$ViralLoad_...),]
  Gene_correlation_sub$TPM = as.numeric(as.character(Gene_correlation_sub$TPM))
  Gene_correlation_sub$ViralLoad_... = as.numeric(as.character(Gene_correlation_sub$ViralLoad_...))
  
  if (!any(is.na(Gene_correlation_sub$TPM))){
    lm_tmp = summary(lm(Gene_correlation_sub$TPM ~Gene_correlation_sub$ViralLoad_...))
    Gene_correlation[i,]$r.squared_NI = lm_tmp$r.squared
    Gene_correlation[i,]$r.adj_NI = lm_tmp$adj.r.squared
    Gene_correlation[i,]$pValue_NI = lm_tmp$coefficients[,"Pr(>|t|)"][2]
    Gene_correlation[i,]$direction_NI = lm_tmp$coefficients[2,1]
  }
  rm(Gene_correlation_sub,lm_tmp)
}
Gene_correlation = merge(Gene_correlation,Gene_foldChange,by="geneID",all.x=T)

head(Gene_correlation)
Gene_correlation_correlated = Gene_correlation[Gene_correlation$is_highExp == "highExp" & Gene_correlation$r.squared_NI>=0.3,]
### barplots
GeneSets = c("highVar","lowVar","NI")
Gene_correlation_plot = Gene_correlation[Gene_correlation$geneName %in% TFs$geneName,]

### keep the high exp. genes only (TPM >=1)
Gene_correlation_plot = Gene_correlation_plot[order(-Gene_correlation_plot$r.squared_NI),]
Gene_correlation_plot = Gene_correlation_plot[Gene_correlation_plot$is_highExp == "highExp",]
Gene_correlation_plot$geneName = factor(Gene_correlation_plot$geneName,levels=unique(Gene_correlation_plot$geneName))

########################## add methylation
Methylation_all = read.delim("inputs/5mC_allhg19.bed_summary",header=F,sep="")
Methylation_lowVar = read.delim("inputs/Summary.5mC_lowVar.table.individual_family",header=T,sep="")
Methylation_lowVar = Methylation_lowVar[Methylation_lowVar$TE_group=="All",]
Methylation_lowVar$TE_group = "lowVar"

Methylation_highVar = read.delim("inputs/Summary.5mC_highVar.table.individual_family",header=T,sep="")
Methylation_highVar = Methylation_highVar[Methylation_highVar$TE_group=="All",]
Methylation_highVar$TE_group = "highVar"

Methylation_reduced = read.delim("inputs/Summary.5mC_reduced.table.individual_family",header=T,sep="")
Methylation_reduced = Methylation_reduced[Methylation_reduced$TE_group=="All",]
Methylation_reduced$TE_group = "reduced"

colnames(Methylation_all) = c("Name",colnames(Methylation_lowVar))
Methylation_all = Methylation_all[Methylation_all$TE_group=="All",]
Methylation_all[,2] = c("exon","gene","intron","gene_u1k","gene_u5k","TE_fourGroups","TE_all")
Methylation_all = Methylation_all[,-1]

Methylation_all = rbind(Methylation_all,Methylation_lowVar,Methylation_highVar,Methylation_reduced)
Methylation_all[,(ncol(Methylation_all)+1):(ncol(Methylation_all)+70)] = NA
colnames(Methylation_all)[(ncol(Methylation_all)-70+1):ncol(Methylation_all)] = gsub("_dup.Methyl|.Methyl","_methy",colnames(Methylation_all)[grepl("Methyl",colnames(Methylation_all))])
colnames(Methylation_all) = gsub("_dup.",".",colnames(Methylation_all))
individualIDs = gsub("_dup.Methyl|.Methyl","",colnames(Methylation_all)[grepl("Methyl",colnames(Methylation_all))])
individualID = "EU47_Flu"
for (individualID in individualIDs){
  Methylation_all[,which(colnames(Methylation_all)== paste(individualID,"_methy",sep=""))] = Methylation_all[,which(colnames(Methylation_all)== paste(individualID,".Methyl",sep=""))]/Methylation_all[,which(colnames(Methylation_all)== paste(individualID,".count",sep=""))]
}
rownames(Methylation_all) = Methylation_all$TE_group
Methylation_all = data.frame(t(Methylation_all[,grepl("methy",colnames(Methylation_all))]))
Methylation_all$sampleID = gsub("_methy","",rownames(Methylation_all))
Methylation_all$individualID = gsub("_Flu|_NI","",Methylation_all$sampleID)
Methylation_all_kept = data.frame(Methylation_all[grepl("_NI",Methylation_all$sampleID),])
colnames(Methylation_all_kept)[1:10] = paste(colnames(Methylation_all_kept)[1:10],"_methy",sep="")

################################# Integrate variables
# 2022_3_14
factorList = c("TRIM28","SETDB1","ZNF460","PLAGL1","ZKSCAN5","ELF3","ELF1","NFYC","MEF2D",
               "STAT2","IRF1","IRF9","STAT5A","IRF7","REL","RELA","STAT3","STAT5B")
factorList_immune = c("STAT2","IRF1","IRF9","STAT5A","IRF7","REL","RELA","STAT3","STAT5B")
factorList_TE = c("TRIM28","SETDB1","ZNF460","PLAGL1","TE_transcripts","TE_methyl")
factorList_IFN = c("TRIM28","SETDB1","IFN_median")
Variables = data.frame(t(Gene_correlation[Gene_correlation$geneName %in% factorList,grepl("_NI",colnames(Gene_correlation)) | grepl("geneName",colnames(Gene_correlation))]))
colnames(Variables) = as.character(c(Variables[rownames(Variables) == "geneName",]))
Variables$individualID = gsub("_NI","",rownames(Variables))
## TE transcripts, Age, ViralLoad
Variables = merge(Variables,TE_transcript[,c("Age","individualID","TE_transcript","ViralLoad_...")],by="individualID",all.y=T)
## Interferon
Variables = merge(Variables,immuneGenes_IFN[,c("individualID","IFN_median")],by="individualID",all.x=T)

## Methylation in TEs
## Interferon
Variables = merge(Variables,Methylation_all_kept[,c("individualID","TE_all_methy")],by="individualID",all.x=T)

################ modeling
## format
Variables_plot = Variables
colnames(Variables_plot)[which(colnames(Variables_plot)=="ViralLoad_...")] = "Viral_load"
colnames(Variables_plot)[which(colnames(Variables_plot)=="TE_transcript")] = "TE_transcripts"
colnames(Variables_plot)[which(colnames(Variables_plot)=="TE_all_methy")] = "TE_methyl"
rownames(Variables_plot) = Variables_plot$individualID
Variables_plot = Variables_plot[,-1]
for(i in 1:ncol(Variables_plot)) {
  Variables_plot[,i] = as.numeric(as.character(Variables_plot[,i]))
}

Variables_plot_35 = Variables_plot[!is.na(Variables_plot$TE_methyl),]
Variables_plot_35 = Variables_plot_35[,c("Viral_load","IFN_median","TE_transcripts","Age","TRIM28","SETDB1","TE_methyl","PLAGL1","ZNF460")]

############ plot 1: correlation plot
pdf("EMC-Figure_S7I_predictiveModel_2_2022_3_27.pdf",     
    width = 10,        # 5 x 300 pixels
    height = 10,
    pointsize = 10)
chart.Correlation(Variables_plot_35,
                  method="pearson",
                  histogram=TRUE,
                  pch=16)
dev.off()

############################################# Section 4: training 
Variables_plot_35 = Variables_plot[!is.na(Variables_plot$TE_methyl),]
######### final models
### IFN only
final.1_IFN1 = glm(Viral_load ~ IFN_median + Age, data = Variables_plot_35,family = gaussian())
final.1_IFN1.select = stepAIC(final.1_IFN1)
final.1_IFN1.select2 = lm(data = Variables_plot_35,Viral_load ~ IFN_median)
summary(final.1_IFN1.select2)
final.1_IFN2 = glm(Viral_load ~ IFN_median * Age, data = Variables_plot_35,family = gaussian())
final.1_IFN2.select = stepAIC(final.1_IFN2)
final.1_IFN2.select2 = lm(data = Variables_plot_35,Viral_load ~ IFN_median * Age)
summary(final.1_IFN2.select2)

### Immune genes only
final.1_immune1 = glm(Viral_load ~ STAT2 + IRF1 + IRF9 + STAT5A + IRF7 + REL + RELA + STAT3 + STAT5B + Age, data = Variables_plot_35,family = gaussian())
final.1_immune1.select = stepAIC(final.1_immune1)
final.1_immune1.select2 = lm(data = Variables_plot_35,Viral_load ~ IRF7 + IRF1)
summary(final.1_immune1.select2)
final.1_immune2 = glm(Viral_load ~ (STAT2 + IRF1 + IRF9 + STAT5A + IRF7 + REL + RELA + STAT3 + STAT5B)*Age, data = Variables_plot_35,family = gaussian())
final.1_immune2.select = stepAIC(final.1_immune2)
final.1_immune2.select2 = lm(data = Variables_plot_35,Viral_load ~ (STAT3 + REL + STAT5A + IRF7 + (RELA + IRF1 + IRF9 + STAT2)*Age))
summary(final.1_immune2.select2)

# ### testing
final.1_immune2.four = glm(Viral_load ~ (STAT2 + IRF1 + IRF7 + IRF9) * Age, data = Variables_plot_35,family = gaussian())
final.1_immune2.four.select = stepAIC(final.1_immune2.four)
final.1_immune2.four.select2 = lm(data = Variables_plot_35,Viral_load ~ IRF7 + (IRF1 + STAT2 + IRF9)*Age)
summary(final.1_immune2.four.select2)
final.1_immune1.five = glm(Viral_load ~ (STAT2 + IRF1 + IRF7 + IRF9 + STAT5A + Age), data = Variables_plot_35,family = gaussian())
final.1_immune1.five.select = stepAIC(final.1_immune1.five)
final.1_immune1.five.select2 = lm(data = Variables_plot_35,Viral_load ~ STAT2)
summary(final.1_immune1.five.select2)
final.1_immune2.five = glm(Viral_load ~ (STAT2 + IRF1 + IRF7 + IRF9 + STAT5A) * Age, data = Variables_plot_35,family = gaussian())
final.1_immune2.five.select = stepAIC(final.1_immune2.five)
final.1_immune2.five.select2 = lm(data = Variables_plot_35,Viral_load ~ IRF7 + (IRF1 + STAT2 + IRF9 + STAT5A)*Age)
summary(final.1_immune2.five.select2)
final.1_immune1.six = glm(Viral_load ~ (STAT2 + IRF1 + IRF7 + IRF9 + STAT5A + REL + Age), data = Variables_plot_35,family = gaussian())
final.1_immune1.six.select = stepAIC(final.1_immune1.six)
final.1_immune1.six.select2 = lm(data = Variables_plot_35,Viral_load ~ STAT2)
summary(final.1_immune1.six.select2)
final.1_immune2.six = glm(Viral_load ~ (STAT2 + IRF1 + IRF7 + IRF9 + STAT5A + REL) * Age, data = Variables_plot_35,family = gaussian())
final.1_immune2.six.select = stepAIC(final.1_immune2.six)
final.1_immune2.six.select2 = lm(data = Variables_plot_35,Viral_load ~ IRF7 + (IRF1 + STAT2 + IRF9 + STAT5A + REL)*Age)
summary(final.1_immune2.six.select2)

### TEs only
final.1_TE1 = glm(Viral_load ~ TRIM28 + SETDB1 + ZNF460 + PLAGL1 + TE_transcripts + TE_methyl + Age, data = Variables_plot_35,family = gaussian())
final.1_TE1.select = stepAIC(final.1_TE1)
final.1_TE1.select2 = lm(data = Variables_plot_35,Viral_load ~ SETDB1 + TRIM28 + TE_methyl + Age)
summary(final.1_TE1.select2)
final.1_TE2 = glm(Viral_load ~ (TRIM28 + SETDB1 + ZNF460 + PLAGL1 + TE_transcripts + TE_methyl) * Age, data = Variables_plot_35,family = gaussian())
final.1_TE2.select = stepAIC(final.1_TE2)
final.1_TE2.select2 = lm(data = Variables_plot_35,Viral_load ~ SETDB1 + TE_transcripts + TE_methyl + PLAGL1*Age)
summary(final.1_TE2.select2)

### testing 
final.1_TE2.five = glm(Viral_load ~ (SETDB1 + ZNF460 + PLAGL1 + TE_transcripts + TRIM28) * Age, data = Variables_plot_35,family = gaussian())
final.1_TE2.five.select = stepAIC(final.1_TE2.five)
final.1_TE2.five.select2 = lm(data = Variables_plot_35,Viral_load ~ TE_transcripts + PLAGL1*Age)
summary(final.1_TE2.four.select2)

### Combine IFN and TEs
#
final.1_IFNTE1 = glm(Viral_load ~ IFN_median + TRIM28 + SETDB1 + ZNF460 + PLAGL1 + TE_transcripts + TE_methyl + Age, data = Variables_plot_35,family = gaussian())
final.1_IFNTE1.select.IFN = stepAIC(final.1_IFNTE1)
final.1_IFNTE1.select2.IFN = lm(data = Variables_plot_35,Viral_load ~ TRIM28 + SETDB1 + TE_methyl + Age)
summary(final.1_IFNTE1.select2.IFN)
#
final.1_IFNTE2 = glm(Viral_load ~ (IFN_median + TRIM28 + SETDB1 + ZNF460 + PLAGL1 + TE_transcripts + TE_methyl) * Age, data = Variables_plot_35,family = gaussian())
final.1_IFNTE2.select.IFN = stepAIC(final.1_IFNTE2)
final.1_IFNTE2.select2.IFN = lm(data = Variables_plot_35,Viral_load ~ (IFN_median + TE_transcripts + TE_methyl + SETDB1)*Age)
summary(final.1_IFNTE2.select2.IFN)

# STAT2
final.1_IFNTE1 = glm(Viral_load ~ STAT2 + TRIM28 + SETDB1 + ZNF460 + PLAGL1 + TE_transcripts + TE_methyl + Age, data = Variables_plot_35,family = gaussian())
final.1_IFNTE1.select = stepAIC(final.1_IFNTE1)
final.1_IFNTE1.select2.STAT2 = lm(data = Variables_plot_35,Viral_load ~ (TE_methyl + STAT2))
summary(final.1_IFNTE1.select2.STAT2)
final.1_IFNTE2 = glm(Viral_load ~ (STAT2 + TRIM28 + SETDB1 + ZNF460 + PLAGL1 + TE_transcripts + TE_methyl) * Age, data = Variables_plot_35,family = gaussian())
final.1_IFNTE2.select = stepAIC(final.1_IFNTE2)
final.1_IFNTE2.select2.STAT2 = lm(data = Variables_plot_35,Viral_load ~ (ZNF460 + TE_transcripts + STAT2 + (TE_methyl + PLAGL1)*Age))
summary(final.1_IFNTE2.select2.STAT2)

final.1_IFNTE2.five = glm(Viral_load ~ (STAT2 + SETDB1 + ZNF460 + PLAGL1 + TE_transcripts + TRIM28) * Age, data = Variables_plot_35,family = gaussian())
final.1_IFNTE2.five.select = stepAIC(final.1_IFNTE2.five)
final.1_IFNTE2.five.select2 = lm(data = Variables_plot_35,Viral_load ~ (ZNF460 + TE_transcripts + (STAT2 + SETDB1 +PLAGL1) * Age))
summary(final.1_TE2.four.select2)

## final model
summary(final.1_IFN1.select2)
final.1_IFN2.select2
summary(final.1_immune1.select2)
summary(final.1_immune2.select2)
summary(final.1_immune2.four.select2)
summary(final.1_immune1.five.select2)
summary(final.1_immune2.five.select2)
summary(final.1_immune1.six.select2)
summary(final.1_immune2.six.select2)

Variables_plot_35$final.1_IFN1.select2 = predict(final.1_IFN1.select2, type = "response")
Variables_plot_35$final.1_IFN2.select2 = predict(final.1_IFN2.select2, type = "response")
Variables_plot_35$final.1_immune1.select2 = predict(final.1_immune1.select2, type = "response")
Variables_plot_35$final.1_immune2.select2 = predict(final.1_immune2.select2, type = "response")
Variables_plot_35$final.1_immune2.four.select2 = predict(final.1_immune2.four.select2, type = "response")
Variables_plot_35$final.1_immune1.five.select2 = predict(final.1_immune1.five.select2, type = "response")
Variables_plot_35$final.1_immune2.five.select2 = predict(final.1_immune2.five.select2, type = "response")
Variables_plot_35$final.1_immune1.six.select2 = predict(final.1_immune1.six.select2, type = "response")
Variables_plot_35$final.1_immune2.six.select2 = predict(final.1_immune2.six.select2, type = "response")
Variables_plot_35$final.1_TE1.select2 = predict(final.1_TE1.select2, type = "response")
Variables_plot_35$final.1_TE2.select2 = predict(final.1_TE2.select2, type = "response")
Variables_plot_35$final.1_TE2.four.select2 = predict(final.1_TE2.four.select2, type = "response")
Variables_plot_35$final.1_IFNTE1.select2.IFN = predict(final.1_IFNTE1.select2.IFN, type = "response")
Variables_plot_35$final.1_IFNTE2.select2.IFN = predict(final.1_IFNTE2.select2.IFN, type = "response")
Variables_plot_35$final.1_IFNTE1.select2.STAT2 = predict(final.1_IFNTE1.select2.STAT2, type = "response")
Variables_plot_35$final.1_IFNTE2.select2.STAT2 = predict(final.1_IFNTE2.select2.STAT2, type = "response")

########## plots 2022_4_2
glist = list()
Order = 1
colnames(Variables_plot_35)
colID = 32
for (colID in grep("final",colnames(Variables_plot_35))){
  print(colID)
  p1 = ggplot(Variables_plot_35,aes(x=Viral_load,y=Variables_plot_35[,colID])) + 
    geom_point(aes(),size=2)+
    xlab("Observed viral load (%)") + 
    ylab(colnames(Variables_plot_35)[colID]) +
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 p.digits = 3,
                 aes(label = paste(..adj.rr.label.., sep = "*plain(\",\")~")), 
                 parse = TRUE)+
    theme(
      #plot.title = element_text(hjust = 0.5, size = rel(1.5)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.text=element_text(colour="black",size=rel(1.1),angle = 0),
      axis.text.x=element_text(colour="black",size=rel(1.1)),
      axis.title=element_text(colour="black",size=rel(1.1)),
      legend.key = element_rect(colour = "transparent", fill = "white"),
      legend.position="right",
      legend.background = element_blank(),
      legend.text=element_text(size=rel(1.1)))
  glist[[Order]] <- ggplotGrob(p1)
  Order = Order + 1
}

pdf("EMC-Figure6HIJ_predictiveModel_3_2022_4_3.pdf",     
    width = 9,        # 5 x 300 pixels
    height = 18,
    pointsize = 10)
do.call("grid.arrange",c(glist,ncol=3))
dev.off()


################################################ Section 5: Each variable
#####
glist = list()
Order = 1
Variables_plot_tmp = Variables_plot     ### all samples
### IFN scores and viral loads
p1<-ggplot(Variables_plot_tmp, aes(Viral_load,IFN_median)) +
  geom_point(aes(color = Age),size=2) + 
  scale_color_gradient2(low="#f7fbff",high = "#08306b")+
  ylab("Basal IFN score")+
  xlab("Viral load post-infection (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

### TE transcripts and viral loads
p1<-ggplot(Variables_plot_tmp, aes(Viral_load,TE_transcripts)) +
  geom_point(aes(color = Age),size=2) + 
  scale_color_gradient2(low="#f7fbff",high = "#08306b")+
  ylab("Basal TE transcripts")+
  xlab("Viral load post-infection (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

### Viral loads and TRIM28
p1<-ggplot(Variables_plot_tmp, aes(Viral_load,TRIM28)) +
  geom_point(aes(color = Age),size=2) + 
  scale_color_gradient2(low="#f7fbff",high = "#08306b")+
  ylab("Basal TRIM28 exp. (TPM)")+
  xlab("Viral load post-infection (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

### Viral loads and SETDB1
p1<-ggplot(Variables_plot_tmp, aes(Viral_load,SETDB1)) +
  geom_point(aes(color = Age),size=2) + 
  scale_color_gradient2(low="#f7fbff",high = "#08306b")+
  ylab("Basal SETDB1 exp. (TPM)")+
  xlab("Viral load post-infection (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

### Viral loads and ZNF460
p1<-ggplot(Variables_plot_tmp, aes(Viral_load,ZNF460)) +
  geom_point(aes(color = Age),size=2) + 
  scale_color_gradient2(low="#f7fbff",high = "#08306b")+
  ylab("Basal ZNF460 exp. (TPM)")+
  xlab("Viral load post-infection (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

### Viral loads and PLAGL1
p1<-ggplot(Variables_plot_tmp, aes(Viral_load,PLAGL1)) +
  geom_point(aes(color = Age),size=2) + 
  scale_color_gradient2(low="#f7fbff",high = "#08306b")+
  ylab("Basal PLAGL1 exp. (TPM)")+
  xlab("Viral load post-infection (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))


glist[[Order]] = p1
Order = Order + 1

Gene_correlation_TFs_high[Gene_correlation_TFs_high$geneName == "PLAGL1",]

### Viral loads and Methylation
p1<-ggplot(Variables_plot_tmp, aes(Viral_load,TE_methyl)) +
  geom_point(aes(color = Age),size=2) + 
  scale_color_gradient2(low="#f7fbff",high = "#08306b")+
  ylab("Basal TE methyl. (%)")+
  xlab("Viral loads post-infection (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1


### Age and viral load
Variables_plot_tmp$TEexp = Variables_plot_tmp$TE_transcripts*100
p1<-ggplot(Variables_plot_tmp, aes(Age,Viral_load)) +
  geom_point(aes(color = TEexp),size=2) + 
  scale_color_gradient2(low="#f7fcfd",high = "#4d004b")+
  #guides(color=guide_legend(title="Basal TE transcripts"))+
  xlab("Age (yr)")+
  ylab("Viral load post-infection (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  #labs(title=type)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1


Variables_plot_tmp$VirL = Variables_plot_tmp$Viral_load
### Age and TE expression
p1<-ggplot(Variables_plot_tmp, aes(Age,TEexp)) +
  geom_point(aes(color = VirL),size=2) + 
  scale_color_gradient2(low="#f7fcfd",high = "#00441b")+
  #guides(color=guide_legend(title="Basal TE transcripts"))+
  xlab("Age (years old)")+
  ylab("Basal TE transcripts (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  #labs(title=type)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

### SETDB1 and TE expression
p1<-ggplot(Variables_plot_tmp, aes(SETDB1,TEexp)) +
  geom_point(aes(color = VirL),size=2) + 
  scale_color_gradient2(low="#f7fcfd",high = "#00441b")+
  #guides(color=guide_legend(title="Basal TE transcripts"))+
  xlab("Basal SETDB1 exp. (TPM)")+
  ylab("Basal TE transcripts (%)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  #labs(title=type)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

### SETDB1 and TE expression
p1<-ggplot(Variables_plot_tmp, aes(SETDB1,ZNF460)) +
  geom_point(aes(color = VirL),size=2) + 
  scale_color_gradient2(low="#f7fcfd",high = "#00441b")+
  #guides(color=guide_legend(title="Basal TE transcripts"))+
  xlab("Basal SETDB1 exp. (TPM)")+
  ylab("Basal ZNF460 exp. (TPM)")+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
  stat_poly_eq(formula = my.formula,
               eq.with.lhs = "italic(hat(y))~`=`~",
               p.digits = 10,
               aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
               parse = TRUE) +
  stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
  #labs(title=type)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

### plot
pdf("EMC-Figure_S7BCDEFG_predictiveModel4.pdf",     
    width = 12,        # 5 x 300 pixels
    height = 8,
    pointsize = 10)
do.call("grid.arrange",c(glist,ncol=4))
dev.off()

############################################## Section 6: Immune genes violin plot
####################### gene TPM global
## Immune gene list
ImmuneGeneList = read.delim("inputs/InnateDB_IRIS_ensemble_hg19.geneName_2021_10_12",header=F,sep="")
colnames(ImmuneGeneList) = "Immune_gene"

### TFs
TFs_all = read.delim("inputs/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.geneNames",header=F,sep="")
colnames(TFs_all) = "TF_name"
TFs_all$is_immuneRelated = ifelse(TFs_all$TF_name %in% ImmuneGeneList$Immune_gene,"Immune","Non-Immune")

### kept TF genes, achieve 685 TFs
# individual level expression change (TPM)

########################################### step3: determine TFs
### correlation of TFs and TRIM28 and SETDB1
Gene_correlation = Gene_exp_tpm_global[Gene_exp_tpm_global$geneName %in% TFs_all$TF_name,]
Gene_correlation$r.squared_NI = NA
Gene_correlation$r.adj_NI = NA
Gene_correlation$pValue_NI = NA
Gene_correlation$direction_NI = NA
i = 1
for (i in 1:nrow(Gene_correlation)){
  Gene_correlation_sub = data.frame(t(Gene_correlation[i,grepl("_NI$",colnames(Gene_correlation))]))
  colnames(Gene_correlation_sub)[1] = "TPM"
  Gene_correlation_sub$individualID = gsub("_Flu","",rownames(Gene_correlation_sub))
  Gene_correlation_sub$individualID = gsub("_NI","",Gene_correlation_sub$individualID)
  
  Gene_correlation_sub = merge(Gene_correlation_sub,TE_transcript,by="individualID",all=T)
  Gene_correlation_sub=Gene_correlation_sub[!is.na(Gene_correlation_sub$ViralLoad_...),]
  Gene_correlation_sub$TPM = as.numeric(as.character(Gene_correlation_sub$TPM))
  Gene_correlation_sub$ViralLoad_... = as.numeric(as.character(Gene_correlation_sub$ViralLoad_...))
  
  if (!any(is.na(Gene_correlation_sub$TPM))){
    lm_tmp = summary(lm(Gene_correlation_sub$TPM ~Gene_correlation_sub$ViralLoad_...))
    Gene_correlation[i,]$r.squared_NI = lm_tmp$r.squared
    Gene_correlation[i,]$r.adj_NI = lm_tmp$adj.r.squared
    Gene_correlation[i,]$pValue_NI = lm_tmp$coefficients[,"Pr(>|t|)"][2]
    Gene_correlation[i,]$direction_NI = lm_tmp$coefficients[2,1]
  }
  rm(Gene_correlation_sub,lm_tmp)
}
Gene_correlation = merge(Gene_correlation,Gene_foldChange,by="geneID",all.x=T)

### kept unique gene Name
Gene_correlation_TFs = merge(Gene_correlation,TFs_all,by.x="geneName",by.y="TF_name",all.x=T)
Gene_correlation_TFs = Gene_correlation_TFs[order(-Gene_correlation_TFs$baseMean),]
Gene_correlation_TFs = Gene_correlation_TFs[!duplicated(Gene_correlation_TFs$geneName),]

### TF groups, give high var. the highest priority
TFs_tmp = TFs[order(TFs$TE_cluster),]
TFs_tmp = TFs_tmp[!duplicated(TFs_tmp$geneName),]
Gene_correlation_TFs = merge(Gene_correlation_TFs,TFs_tmp,by.x="geneName",by.y="geneName",all.x=T)
Gene_correlation_TFs$TE_cluster = ifelse(is.na(Gene_correlation_TFs$TE_cluster),"Others",Gene_correlation_TFs$TE_cluster)
Gene_correlation_TFs$Label = ifelse(Gene_correlation_TFs$TE_cluster!="Others" & Gene_correlation_TFs$r.squared_NI>=0.3,Gene_correlation_TFs$geneName,NA)
Gene_correlation_TFs$TE_cluster_shape = factor(Gene_correlation_TFs$TE_cluster)
levels(Gene_correlation_TFs$TE_cluster_shape) = c(17,15,18,16)
Shapes = c("high_var"="17","low_var"="15","NI"="18","Others"="16")
colorSet = c("high_var" = "#8c510a","low_var" = "#4d4d4d","NI"="#542788","Others" = "#bdbdbd")

### plots
glist = list()
Order = 1

Gene_correlation_TFs_high = data.frame(Gene_correlation_TFs[Gene_correlation_TFs$is_highExp == "highExp",])
Gene_correlation_TFs_high_correlated = Gene_correlation_TFs_high[Gene_correlation_TFs_high$r.squared_NI>=0.3,]
Gene_correlation_TFs_high_correlated = Gene_correlation_TFs_high_correlated[Gene_correlation_TFs_high_correlated$is_immuneRelated == "Immune",]

### plots
Gene_correlation_TFs_high$is_immuneRelated = factor(Gene_correlation_TFs_high$is_immuneRelated)
p1 = ggplot(Gene_correlation_TFs_high, aes(is_immuneRelated,r.adj_NI)) + #volcanoplot with log2Foldchange versus pvalue
  geom_violin(aes(),alpha=0.2)+
  stat_summary(fun=mean, color="#de2d26",geom="point",fill="#de2d26",shape=95,size = 2)+
  ylab("R square")+
  geom_jitter(aes(color = TE_cluster),width = 0.3,na.rm=TRUE)+
  scale_color_manual(values = colorSet)+
  geom_text(aes(label=Label))+
  #geom_dotplot(binaxis='y', stackdir='center',binwidth=50,
  #             stackratio=1.5, dotsize=1,aes(col=TE_cluster,fill=TE_cluster))+
  #  geom_jitter(aes(col=TE_cluster,fill=TE_cluster,group=TE_cluster),position=position_jitter(0.2),na.rm=TRUE)+
  stat_compare_means(method = "t.test",aes(group = is_immuneRelated,label = paste0("p = ", ..p.format..)))+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text=element_text(colour="black",size=rel(1),angle = 0),
    axis.text.x=element_text(colour="black",vjust=.5),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.background = element_blank(),
    legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

### plots
Gene_correlation_TFs$is_immuneRelated = factor(Gene_correlation_TFs$is_immuneRelated)
p1 = ggplot(Gene_correlation_TFs, aes(is_immuneRelated,r.adj_NI)) + #volcanoplot with log2Foldchange versus pvalue
  geom_violin(aes(),alpha=0.2)+
  stat_summary(fun=mean, color="#de2d26",geom="point",fill="#de2d26",shape=95,size = 2)+
  ylab("R square")+
  geom_jitter(aes(color = TE_cluster),width = 0.3,na.rm=TRUE)+
  scale_color_manual(values = colorSet)+
  geom_text(aes(label=Label))+
  #geom_dotplot(binaxis='y', stackdir='center',binwidth=50,
  #             stackratio=1.5, dotsize=1,aes(col=TE_cluster,fill=TE_cluster))+
  #  geom_jitter(aes(col=TE_cluster,fill=TE_cluster,group=TE_cluster),position=position_jitter(0.2),na.rm=TRUE)+
  stat_compare_means(method = "t.test",aes(group = is_immuneRelated,label = paste0("p = ", ..p.format..)))+
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.text=element_text(colour="black",size=rel(1),angle = 0),
      axis.text.x=element_text(colour="black",vjust=.5),
      axis.title=element_text(colour="black",size=rel(1)),
      legend.position="right",
      legend.background = element_blank(),
      legend.text=element_text(size=rel(1)))
glist[[Order]] = p1
Order = Order + 1

pdf("EMC-Figure6F_predictiveModel5_2022_3_8.pdf",     
    width = 7,        # 5 x 300 pixels
    height = 4,
    pointsize = 10)
do.call("grid.arrange",c(glist,ncol=2))
dev.off()