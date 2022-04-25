###     Author:     Xun Chen
###     Date:       2022/2/26
###     Contact:    xunchen85@gmail.com

##### load pacakges
library(gplots)
library(ggplot2)
library(grid)
library(ggrepel)
library(dplyr)
library(splitstackshape)
library(logspline)
library(ggbeeswarm)
library(ggpubr)
library(multtest)
library(ggpmisc)
library(gridExtra)
library(cowplot)
library(ggforce)
library("viridis") 
library(RColorBrewer)
library(tidyverse)
library(scales)

## specify the color
colorset_TE_cluster = c("a" = "#7C7C7C","Group a" = "#7C7C7C", "b" = "#bf812d","Group b" = "#bf812d","NI" = "#beaed4","Not_candidate" = "#beaed4")
colorset = c("FluPeaks" = "#99000d","NIPeaks" = "#4d4d4d","NoPeaks"="#1f78b4","SharedPeaks"="#6a3d9a","WithPeaks"="#33a02c")
colorset_fluStatus = c("Flu" = "#99000d","NI" = "#4d4d4d","Share"="#6a3d9a")
Cluster_color = c("Cluster_1"="#4d004b", "Cluster_2"="#810f7c", "Cluster_3"="#88419d", "Cluster_4"="#8c6bb1", "Cluster_5"="#8c96c6","Cluster_others"="#bfd3e6","No"="#f7fcfd")
Cluster_color = c("Cluster_1"="#4d004b", "Cluster_2"="#810f7c", "Cluster_3"="#88419d", "Cluster_4"="#8c6bb1", "Cluster_5"="#8c96c6","Cluster_others"="#bfd3e6","No"="white")
PieFill = c("Fluonly"="#0570b0","NIonly"="#0570b0","Share"="#ece7f2")
IndividualGroup_color = c("Group a" = "#7C7C7C","Group b" = "#bf812d","a" = "#7C7C7C","b" = "#bf812d","NI" = "#beaed4","Not_candidate" = "#beaed4")

order_Clusters_NI = rev(c("WithTopConsensus Share","WithTopConsensus NIonly","WithConsensus Share","WithConsensus NIonly","No Share","No NIonly"))
color_Clusters_NI = c("WithTopConsensus Share"="#08306b","WithTopConsensus NIonly"="#08306b","WithConsensus Share"="#9ecae1","WithConsensus NIonly"="#9ecae1","No Share"="#f7fbff","No NIonly"="#f7fbff")
order_Clusters_Flu = rev(c("WithTopConsensus Share","WithTopConsensus Fluonly","WithConsensus Share","WithConsensus Fluonly","No Share","No Fluonly"))
color_Clusters_Flu = c("WithTopConsensus Share"="#08306b","WithTopConsensus Fluonly"="#08306b","WithConsensus Share"="#9ecae1","WithConsensus Fluonly"="#9ecae1","No Share"="#f7fbff","No Fluonly"="#f7fbff")


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

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#### specify my theme for ggplot
my_theme_1_noLegend = theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
                            plot.margin=unit(c(10,10,10,10),"mm"),
                            axis.line = element_line(colour = "black"),
                            axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                            axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                            axis.text=element_text(colour="black",size=rel(1),angle = 0),
                            axis.title=element_text(colour="black",size=rel(1.2)),
                            legend.position="none",
                            legend.text=element_text(size=rel(0.8)),
                            axis.text.x=element_text(angle=90,hjust = 0.95,vjust = .5),
                            legend.key = element_rect(fill = NA))
my_theme_1_Legend = theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(), 
                          plot.margin=unit(c(10,10,10,10),"mm"),
                          axis.line = element_line(colour = "black"),
                          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                          axis.text=element_text(colour="black",size=rel(1),angle = 0),
                          axis.title=element_text(colour="black",size=rel(1.2)),
                          legend.position="right",
                          legend.text=element_text(size=rel(0.8)),
                          axis.text.x=element_text(angle=90,hjust = 0.95,vjust = .5),
                          legend.key = element_rect(fill = NA))

my_theme_2_noLegend = theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
                            plot.margin=unit(c(10,10,10,10),"mm"),
                            axis.line = element_line(colour = "black"),
                            axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                            axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                            axis.text=element_text(colour="black",size=rel(1),angle = 0),
                            axis.title=element_text(colour="black",size=rel(1.2)),
                            legend.position="none",
                            legend.text=element_text(size=rel(0.8)),
                            #axis.text.x=element_text(angle=90,size=rel(0.7),hjust = 0.95,vjust = .5),
                            legend.key = element_rect(fill = NA))
###########################
####### load files
## Input
# input 1: TE list
TElist = read.delim("inputs/TElist.sep.27",header=T,sep="")
# input 3: bed.out file
INTname = "2020_12_21_Aggregate_noEU37_100.bed.out"                    ############ motif analysis I excluded EU37
# input 4: TE consensus length
TEconsensusLen = read.delim("inputs/76subfamily_all_Consensus_len.out",header=F,sep="")

## TE consensus length
colnames(TEconsensusLen) = c("Consensus","Consensus_Len_max")
TEconsensusLen$Consensus = as.character(TEconsensusLen$Consensus)
TEconsensusLen[nrow(TEconsensusLen)+1,] = NA
TEconsensusLen[nrow(TEconsensusLen),1] = "LTR39"
TEconsensusLen[nrow(TEconsensusLen),2] = 783
TEconsensusLen[nrow(TEconsensusLen)+1,] = NA
TEconsensusLen[nrow(TEconsensusLen),1] = "MSTB1"
TEconsensusLen[nrow(TEconsensusLen),2] = 432

Cluster_Consensus_Flu = read.csv(paste("inputs/","EMC_Flu_2020_12_21_Aggregate_noEU37_100.bed.out_allConsensus",sep=""))
Cluster_Consensus_NI = read.csv(paste("inputs/","EMC_NI_2020_12_21_Aggregate_noEU37_100.bed.out_allConsensus",sep=""))

###### Summary File
### TE instances
TEinstance_summary = read.delim("inputs/Summary.expAndCentroid.enriched.table",sep="",header=F)
TEinstance_summary_header = read.delim("inputs/Summary.expAndCentroid.table.header",sep="",header=T)
colnames(TEinstance_summary) = colnames(TEinstance_summary_header)

###################### same as the first part in Final-3-TEmotifPlot_2020_12_23.R
Motif_Flu = read.delim(paste("inputs/","EMC_Flu_","2020_12_21_Aggregate_noEU37_100.bed.out",sep=""),header=F,sep="") 
colnames(Motif_Flu) = c("File","Chr","Start_instance","End_instance","TEinfo","Chr2","Start_summit","End_summit","TEinfo_summit","Score","Strand","Posi","TE_coordinate_1bp","Strand2","Consensus","Ref_posi","Consensus_posi")
Motif_Flu$TE_coordinate_0bp = paste(Motif_Flu$Chr,":",Motif_Flu$Start_instance,"-",Motif_Flu$End_instance,sep="")
Motif_NI = read.delim(paste("inputs/","EMC_NI_","2020_12_21_Aggregate_noEU37_100.bed.out",sep=""),header=F,sep="") 
colnames(Motif_NI) = c("File","Chr","Start_instance","End_instance","TEinfo","Chr2","Start_summit","End_summit","TEinfo_summit","Score","Strand","Posi","TE_coordinate_1bp","Strand2","Consensus","Ref_posi","Consensus_posi")
Motif_NI$TE_coordinate_0bp = paste(Motif_NI$Chr,":",Motif_NI$Start_instance,"-",Motif_NI$End_instance,sep="")

Motif_cluster_Flu = read.csv(paste("inputs/","EMC_Flu_","2020_12_21_Aggregate_noEU37_100.bed.out","_allCluster",sep=""))
Motif_cluster_Flu$TEcoordinate_0bp = paste(Motif_cluster_Flu$Chr,":",Motif_cluster_Flu$Start_instance,"-",Motif_cluster_Flu$End_instance,sep="")
Motif_cluster_Flu = Motif_cluster_Flu[,c("TEcoordinate_0bp","TEsubfamily","TEinstance","TEfamily","TEgroup","Cluster_50bp","Win_50bp","Cluster_100bp","Win_100bp","Cluster_200bp","Win_200bp")]
Motif_cluster_NI = read.csv(paste("inputs/","EMC_NI_","2020_12_21_Aggregate_noEU37_100.bed.out","_allCluster",sep=""))
Motif_cluster_NI$TEcoordinate_0bp = paste(Motif_cluster_NI$Chr,":",Motif_cluster_NI$Start_instance,"-",Motif_cluster_NI$End_instance,sep="")
Motif_cluster_NI = Motif_cluster_NI[,c("TEcoordinate_0bp","TEsubfamily","TEinstance","TEfamily","TEgroup","Cluster_50bp","Win_50bp","Cluster_100bp","Win_100bp","Cluster_200bp","Win_200bp")]

Motif_Flu_full = merge(Motif_Flu,Motif_cluster_Flu,by.x="TE_coordinate_0bp",by.y="TEcoordinate_0bp",all.x=T)
Motif_NI_full = merge(Motif_NI,Motif_cluster_NI,by.x="TE_coordinate_0bp",by.y="TEcoordinate_0bp",all.x=T)
Motif_Flu_full = Motif_Flu_full[!duplicated(Motif_Flu_full$TE_coordinate_0bp),]
Motif_NI_full = Motif_NI_full[!duplicated(Motif_NI_full$TE_coordinate_0bp),]
rm(Motif_Flu)
rm(Motif_NI)

###### preparing the dataset
Motif_Flu_full = data.frame(cSplit(Motif_Flu_full,"TEinfo_summit",sep=":",type.convert = as.character))
TEinstance_summary_Flu = TEinstance_summary[TEinstance_summary$TE_family %in% as.character(TElist[TElist$TE_cluster!="NI",]$TEfamily),]
TEinstance_summary_Flu = merge(TEinstance_summary_Flu,Motif_Flu_full,by.x="TE_instance",by.y="TEinfo_summit_1",all.x=T)

###### version one: include any TEs overlap with peaks (Flu and NI peaks)
### including EU37, but the peak regions are still based on the data without EU37
#TEinstance_summary_Flu_withATAC_withEU37 = TEinstance_summary_Flu[TEinstance_summary_Flu$Counts.centroidATAC>0,]

### without EU37
TEinstance_summary_Flu_withATAC = TEinstance_summary_Flu[!is.na(TEinstance_summary_Flu$TEinfo_summit_2),]
####################### 
TEinstance_summary_Flu_withATAC$H3K27ac_overlap = ifelse(TEinstance_summary_Flu_withATAC$Counts.centroidH3K27ac>0,"Yes","No")
TEinstance_summary_Flu_withATAC$H3K27me3_overlap = ifelse(TEinstance_summary_Flu_withATAC$Counts.centroidH3K27me3>0,"Yes","No")
TEinstance_summary_Flu_withATAC$H3K4me1_overlap = ifelse(TEinstance_summary_Flu_withATAC$Counts.centroidH3K4me1>0,"Yes","No")
TEinstance_summary_Flu_withATAC$H3K4me3_overlap = ifelse(TEinstance_summary_Flu_withATAC$Counts.centroidH3K4me3>0,"Yes","No")
TEinstance_summary_Flu_withATAC$H3K27ac_H3K4me1_overlap = ifelse(TEinstance_summary_Flu_withATAC$Counts.centroidH3K27ac>0 & TEinstance_summary_Flu_withATAC$Counts.centroidH3K4me1>0,"Yes","No")
TEinstance_summary_Flu_withATAC$TEexp_overlap = ifelse(TEinstance_summary_Flu_withATAC$log2FoldChange.RNA>=1 & TEinstance_summary_Flu_withATAC$padj.RNA<=0.05,"Yes","No")

####################### TE DNA methylation
TEinstance_methyl = read.delim("inputs/Summary.5mC.enriched.table",sep="",header=F)
TEinstance_methyl_header = read.delim("inputs/Summary.5mC.table.header",sep="",header=T)
colnames(TEinstance_methyl) = colnames(TEinstance_methyl_header)
TEinstance_methyl = TEinstance_methyl[,c("TE_instance","Count.5mC","Flu_mean.5mC")]
nrow(TEinstance_summary_Flu_withATAC[!is.na(TEinstance_summary_Flu_withATAC$Cluster_100bp)&TEinstance_summary_Flu_withATAC$TE_family=="MLT2B3"&TEinstance_summary_Flu_withATAC$Cluster_100bp=="Cluster_3",])
nrow(TEinstance_summary_Flu_withATAC[!is.na(TEinstance_summary_Flu_withATAC$Cluster_100bp)&TEinstance_summary_Flu_withATAC$TE_family=="MLT2B3"&TEinstance_summary_Flu_withATAC$Cluster_100bp=="Cluster_4",])

###############################
############################### EMC_final_figure 5A and EMC_final_supplementary_figure
Filter = "all"
### without EU37
TEinstance_summary_Flu_withATAC = TEinstance_summary_Flu[!is.na(TEinstance_summary_Flu$TEinfo_summit_2) & TEinstance_summary_Flu$TEinfo_summit_3>0,]
####################### 
TEinstance_summary_Flu_withATAC$H3K27ac_overlap = ifelse(TEinstance_summary_Flu_withATAC$Flu.centroidH3K27ac!="Flu:0:0:0" & TEinstance_summary_Flu_withATAC$Flu.centroidH3K27ac !=0,"Yes","No")
TEinstance_summary_Flu_withATAC$H3K27me3_overlap = ifelse(TEinstance_summary_Flu_withATAC$Flu.centroidH3K27me3!="Flu:0:0:0" & TEinstance_summary_Flu_withATAC$Flu.centroidH3K27me3!=0,"Yes","No")
TEinstance_summary_Flu_withATAC$H3K4me1_overlap = ifelse(TEinstance_summary_Flu_withATAC$Flu.centroidH3K4me1!="Flu:0:0:0" & TEinstance_summary_Flu_withATAC$Flu.centroidH3K4me1!=0,"Yes","No")
TEinstance_summary_Flu_withATAC$H3K4me3_overlap = ifelse(TEinstance_summary_Flu_withATAC$Flu.centroidH3K4me3!="Flu:0:0:0" & TEinstance_summary_Flu_withATAC$Flu.centroidH3K4me3!=0,"Yes","No")
TEinstance_summary_Flu_withATAC$H3K27ac_H3K4me1_overlap = ifelse(TEinstance_summary_Flu_withATAC$Flu.centroidH3K27ac!="Flu:0:0:0" & TEinstance_summary_Flu_withATAC$Flu.centroidH3K4me1!="Flu:0:0:0" & TEinstance_summary_Flu_withATAC$Flu.centroidH3K27ac!=0 & TEinstance_summary_Flu_withATAC$Flu.centroidH3K4me1!=0,"Yes","No")
TEinstance_summary_Flu_withATAC$TEexp_overlap = ifelse(TEinstance_summary_Flu_withATAC$log2FoldChange.RNA>=1 & TEinstance_summary_Flu_withATAC$padj.RNA<=0.05,"Yes","No")
TEinstance_summary_Flu_withATAC$Flu.centroidATAC_tmp = TEinstance_summary_Flu_withATAC$Flu.centroidATAC
TEinstance_summary_Flu_withATAC = data.frame(cSplit(TEinstance_summary_Flu_withATAC,"Flu.centroidATAC_tmp",sep=":",type.convert = as.character))

TEinstance_summary_Flu_withATAC = TEinstance_summary_Flu_withATAC[!is.na(TEinstance_summary_Flu_withATAC$TE_family),]
####################### Nearest genes
### nearest genes
NearestGene = read.delim("inputs/ATAC.inTE_76.closest.bed",sep="",header=F)
NearestGene = cSplit(NearestGene,"V10",sep=",")
NearestGene = cSplit(NearestGene,"V4",sep=":")

### keep nearest genes within 100k
NearestGene = NearestGene[abs(NearestGene$V13)<=100000,]

### if it is up-regulated gene or not (log2FC >=0.5)
NearestGene$Geneexp_overlap = ifelse(NearestGene$V10_13>=0.5&NearestGene$V10_16<=0.05,"Yes","No")

NearestGene_kept = NearestGene[,c("V4_02","V10_01","V10_06","V10_13","V10_16","Geneexp_overlap")]
colnames(NearestGene_kept) = c("TE_instance","GeneID","GeneName","log2FC","pValue_gene","Geneexp_overlap")

####################### Combined table
TEinstance_summary_Flu_withATAC = merge(TEinstance_summary_Flu_withATAC,NearestGene_kept,by="TE_instance",all.x=T)
TEinstance_summary_Flu_withATAC = merge(TEinstance_summary_Flu_withATAC,TEinstance_methyl,by="TE_instance",all.x=T)

####################### cluster corrected
TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected = ifelse(is.na(TEinstance_summary_Flu_withATAC$Cluster_100bp),"Without_consensus",as.character(TEinstance_summary_Flu_withATAC$Cluster_100bp))
TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected_min10 = ifelse(TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected=="Cluster_1"| 
                                                                         TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected=="Cluster_2" |
                                                                         TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected=="Cluster_3" |
                                                                         TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected=="Cluster_4" |
                                                                         TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected=="Cluster_5" ,TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected,"Others")
TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected_min10 = ifelse(is.na(TEinstance_summary_Flu_withATAC$Cluster_100bp),"Without_consensus",as.character(TEinstance_summary_Flu_withATAC$Cluster_100bp_corrected_min10))

####################### calculate the proportion
### Sum of each cluster (we group peak regions with a number of less than 10 instances)
TEinstance_summary_Flu_withATAC_summary = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp) %>% dplyr::summarise(n = n(),na.rm = TRUE))
TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary$Cluster_100bp),"Without_consensus",as.character(TEinstance_summary_Flu_withATAC_summary$Cluster_100bp))
TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected_min10 = ifelse(TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected=="Cluster_1"| 
                                                                                 TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected=="Cluster_2" |
                                                                                 TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected=="Cluster_3" |
                                                                                 TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected=="Cluster_4" |
                                                                                 TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected=="Cluster_5" ,TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected,"Others")
TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected_min10 = ifelse(TEinstance_summary_Flu_withATAC_summary$n<10|TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected =="No"|TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected =="Without_consensus","Others",TEinstance_summary_Flu_withATAC_summary$Cluster_100bp_corrected_min10)
TEinstance_summary_Flu_withATAC_summary$ID = paste(TEinstance_summary_Flu_withATAC_summary$TE_family,TEinstance_summary_Flu_withATAC_summary$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary[TEinstance_summary_Flu_withATAC_summary$TE_family=="MLT2B3",]

#### H3K27ac
TEinstance_summary_Flu_withATAC_summary_tmp = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp,H3K27ac_overlap) %>% dplyr::summarise(n = n(),na.rm=TRUE))
TEinstance_summary_Flu_withATAC_summary_tmp$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp_yes = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K27ac_overlap=="Yes",]
TEinstance_summary_Flu_withATAC_summary_tmp_no = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K27ac_overlap=="No",]
TEinstance_summary_Flu_withATAC_summary_tmp2 = merge(TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,by="ID",all=T)
TEinstance_summary_Flu_withATAC_summary_tmp3 = TEinstance_summary_Flu_withATAC_summary_tmp2[,c(1,5,10)]
colnames(TEinstance_summary_Flu_withATAC_summary_tmp3) = c("ID",colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[4],colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[9])
TEinstance_summary_Flu_withATAC_summary_tmp3[,2] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,2]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,2])
TEinstance_summary_Flu_withATAC_summary_tmp3[,3] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,3]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,3])
rm(TEinstance_summary_Flu_withATAC_summary_tmp,TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,TEinstance_summary_Flu_withATAC_summary_tmp2)
TEinstance_summary_Flu_withATAC_summary = merge(TEinstance_summary_Flu_withATAC_summary,TEinstance_summary_Flu_withATAC_summary_tmp3,all.x=T,by="ID")

#### H3K4me1
TEinstance_summary_Flu_withATAC_summary_tmp = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp,H3K4me1_overlap) %>% dplyr::summarise(n = n(),na.rm=TRUE))
TEinstance_summary_Flu_withATAC_summary_tmp$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp_yes = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K4me1_overlap=="Yes",]
TEinstance_summary_Flu_withATAC_summary_tmp_no = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K4me1_overlap=="No",]
TEinstance_summary_Flu_withATAC_summary_tmp2 = merge(TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,by="ID",all=T)
TEinstance_summary_Flu_withATAC_summary_tmp3 = TEinstance_summary_Flu_withATAC_summary_tmp2[,c(1,5,10)]
colnames(TEinstance_summary_Flu_withATAC_summary_tmp3) = c("ID",colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[4],colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[9])
TEinstance_summary_Flu_withATAC_summary_tmp3[,2] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,2]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,2])
TEinstance_summary_Flu_withATAC_summary_tmp3[,3] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,3]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,3])
rm(TEinstance_summary_Flu_withATAC_summary_tmp,TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,TEinstance_summary_Flu_withATAC_summary_tmp2)
TEinstance_summary_Flu_withATAC_summary = merge(TEinstance_summary_Flu_withATAC_summary,TEinstance_summary_Flu_withATAC_summary_tmp3,all.x=T,by="ID")

#### H3K4me3
TEinstance_summary_Flu_withATAC_summary_tmp = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp,H3K4me3_overlap) %>% dplyr::summarise(n = n(),na.rm=TRUE))
TEinstance_summary_Flu_withATAC_summary_tmp$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp_yes = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K4me3_overlap=="Yes",]
TEinstance_summary_Flu_withATAC_summary_tmp_no = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K4me3_overlap=="No",]
TEinstance_summary_Flu_withATAC_summary_tmp2 = merge(TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,by="ID",all=T)
TEinstance_summary_Flu_withATAC_summary_tmp3 = TEinstance_summary_Flu_withATAC_summary_tmp2[,c(1,5,10)]
colnames(TEinstance_summary_Flu_withATAC_summary_tmp3) = c("ID",colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[4],colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[9])
TEinstance_summary_Flu_withATAC_summary_tmp3[,2] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,2]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,2])
TEinstance_summary_Flu_withATAC_summary_tmp3[,3] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,3]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,3])
rm(TEinstance_summary_Flu_withATAC_summary_tmp,TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,TEinstance_summary_Flu_withATAC_summary_tmp2)
TEinstance_summary_Flu_withATAC_summary = merge(TEinstance_summary_Flu_withATAC_summary,TEinstance_summary_Flu_withATAC_summary_tmp3,all.x=T,by="ID")

#### H3K27me3
TEinstance_summary_Flu_withATAC_summary_tmp = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp,H3K27me3_overlap) %>% dplyr::summarise(n = n(),na.rm=TRUE))
TEinstance_summary_Flu_withATAC_summary_tmp$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp_yes = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K27me3_overlap=="Yes",]
TEinstance_summary_Flu_withATAC_summary_tmp_no = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K27me3_overlap=="No",]
TEinstance_summary_Flu_withATAC_summary_tmp2 = merge(TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,by="ID",all=T)
TEinstance_summary_Flu_withATAC_summary_tmp3 = TEinstance_summary_Flu_withATAC_summary_tmp2[,c(1,5,10)]
colnames(TEinstance_summary_Flu_withATAC_summary_tmp3) = c("ID",colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[4],colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[9])
TEinstance_summary_Flu_withATAC_summary_tmp3[,2] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,2]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,2])
TEinstance_summary_Flu_withATAC_summary_tmp3[,3] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,3]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,3])
rm(TEinstance_summary_Flu_withATAC_summary_tmp,TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,TEinstance_summary_Flu_withATAC_summary_tmp2)
TEinstance_summary_Flu_withATAC_summary = merge(TEinstance_summary_Flu_withATAC_summary,TEinstance_summary_Flu_withATAC_summary_tmp3,all.x=T,by="ID")

#### H3K27ac_H3K4me1
TEinstance_summary_Flu_withATAC_summary_tmp = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp,H3K27ac_H3K4me1_overlap) %>% dplyr::summarise(n = n(),na.rm=TRUE))
TEinstance_summary_Flu_withATAC_summary_tmp$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp_yes = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K27ac_H3K4me1_overlap=="Yes",]
TEinstance_summary_Flu_withATAC_summary_tmp_no = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$H3K27ac_H3K4me1_overlap=="No",]
TEinstance_summary_Flu_withATAC_summary_tmp2 = merge(TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,by="ID",all=T)
TEinstance_summary_Flu_withATAC_summary_tmp3 = TEinstance_summary_Flu_withATAC_summary_tmp2[,c(1,5,10)]
colnames(TEinstance_summary_Flu_withATAC_summary_tmp3) = c("ID",colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[4],colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[9])
TEinstance_summary_Flu_withATAC_summary_tmp3[,2] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,2]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,2])
TEinstance_summary_Flu_withATAC_summary_tmp3[,3] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,3]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,3])
rm(TEinstance_summary_Flu_withATAC_summary_tmp,TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,TEinstance_summary_Flu_withATAC_summary_tmp2)
TEinstance_summary_Flu_withATAC_summary = merge(TEinstance_summary_Flu_withATAC_summary,TEinstance_summary_Flu_withATAC_summary_tmp3,all.x=T,by="ID")

#### TE exp.
TEinstance_summary_Flu_withATAC_summary_tmp = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp,TEexp_overlap) %>% dplyr::summarise(n = n(),na.rm=TRUE))
TEinstance_summary_Flu_withATAC_summary_tmp$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp_yes = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$TEexp_overlap=="Yes",]
TEinstance_summary_Flu_withATAC_summary_tmp_no = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$TEexp_overlap=="No",]
TEinstance_summary_Flu_withATAC_summary_tmp2 = merge(TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,by="ID",all=T)
TEinstance_summary_Flu_withATAC_summary_tmp3 = TEinstance_summary_Flu_withATAC_summary_tmp2[,c(1,5,10)]
colnames(TEinstance_summary_Flu_withATAC_summary_tmp3) = c("ID",colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[4],colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[9])
TEinstance_summary_Flu_withATAC_summary_tmp3[,2] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,2]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,2])
TEinstance_summary_Flu_withATAC_summary_tmp3[,3] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,3]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,3])
rm(TEinstance_summary_Flu_withATAC_summary_tmp,TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,TEinstance_summary_Flu_withATAC_summary_tmp2)
TEinstance_summary_Flu_withATAC_summary = merge(TEinstance_summary_Flu_withATAC_summary,TEinstance_summary_Flu_withATAC_summary_tmp3,all.x=T,by="ID")

#### Gene exp.
TEinstance_summary_Flu_withATAC_summary_tmp = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp,Geneexp_overlap) %>% dplyr::summarise(n = n(),na.rm=TRUE))
TEinstance_summary_Flu_withATAC_summary_tmp$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp_yes = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$Geneexp_overlap=="Yes",]
TEinstance_summary_Flu_withATAC_summary_tmp_no = TEinstance_summary_Flu_withATAC_summary_tmp[TEinstance_summary_Flu_withATAC_summary_tmp$Geneexp_overlap=="No",]
TEinstance_summary_Flu_withATAC_summary_tmp2 = merge(TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,by="ID",all=T)
TEinstance_summary_Flu_withATAC_summary_tmp3 = TEinstance_summary_Flu_withATAC_summary_tmp2[,c(1,5,10)]
colnames(TEinstance_summary_Flu_withATAC_summary_tmp3) = c("ID",colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[4],colnames(TEinstance_summary_Flu_withATAC_summary_tmp2)[9])
TEinstance_summary_Flu_withATAC_summary_tmp3[,2] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,2]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,2])
TEinstance_summary_Flu_withATAC_summary_tmp3[,3] = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_tmp3[,3]),0,TEinstance_summary_Flu_withATAC_summary_tmp3[,3])
rm(TEinstance_summary_Flu_withATAC_summary_tmp,TEinstance_summary_Flu_withATAC_summary_tmp_yes,TEinstance_summary_Flu_withATAC_summary_tmp_no,TEinstance_summary_Flu_withATAC_summary_tmp2)
TEinstance_summary_Flu_withATAC_summary = merge(TEinstance_summary_Flu_withATAC_summary,TEinstance_summary_Flu_withATAC_summary_tmp3,all.x=T,by="ID")

#### DNA methylation 2021_10_26
TEinstance_summary_Flu_withATAC_summary_tmp1 = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp_corrected_min10) %>% dplyr::summarize(Mean = mean(Flu_mean.5mC, na.rm=TRUE)))
TEinstance_summary_Flu_withATAC_summary_tmp1$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp1$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp1$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp2 = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family,Cluster_100bp_corrected_min10) %>% dplyr::summarize(Mean = mean(Count.5mC, na.rm=TRUE)))
TEinstance_summary_Flu_withATAC_summary_tmp2$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp2$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp2$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp_all1 = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family) %>% dplyr::summarize(Mean = mean(Flu_mean.5mC, na.rm=TRUE)))
TEinstance_summary_Flu_withATAC_summary_tmp_all1$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp_all1$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp_all1$Cluster_100bp)
TEinstance_summary_Flu_withATAC_summary_tmp_all2 = data.frame(TEinstance_summary_Flu_withATAC %>% group_by(TE_family) %>% dplyr::summarize(Mean = mean(Count.5mC, na.rm=TRUE)))
TEinstance_summary_Flu_withATAC_summary_tmp_all2$ID = paste(TEinstance_summary_Flu_withATAC_summary_tmp_all2$TE_family,TEinstance_summary_Flu_withATAC_summary_tmp_all2$Cluster_100bp)
### subfamily level
TEinstance_summary_Flu_withATAC_summary_methy_subfamily = merge(TEinstance_summary_Flu_withATAC_summary_tmp_all1,TEinstance_summary_Flu_withATAC_summary_tmp_all2,by="ID",all=T)
colnames(TEinstance_summary_Flu_withATAC_summary_methy_subfamily) = c("ID","TE_family.x","DNA_methyl.Level","TE_family.y","DNA_methyl.Count")
### region level
TEinstance_summary_Flu_withATAC_summary_methy_region = merge(TEinstance_summary_Flu_withATAC_summary_tmp1,TEinstance_summary_Flu_withATAC_summary_tmp2,by="ID",all=T)
colnames(TEinstance_summary_Flu_withATAC_summary_methy_region) = c("ID","TE_family.x","Cluster_100bp.x","DNA_methyl.Level","TE_family.y","Cluster_100bp.y","DNA_methyl.Count")
rm(TEinstance_summary_Flu_withATAC_summary_tmp1,TEinstance_summary_Flu_withATAC_summary_tmp2,TEinstance_summary_Flu_withATAC_summary_tmp_all1,TEinstance_summary_Flu_withATAC_summary_tmp_all2)


################ Sum up across peak regions
# region level
TEinstance_summary_Flu_withATAC_summary$ID_final = (TEinstance_summary_Flu_withATAC_summary$TE_family) 

# subfamily level
TEinstance_summary_Flu_withATAC_summary_final = data.frame(TEinstance_summary_Flu_withATAC_summary %>% 
                                                             group_by(ID_final) %>% 
                                                             summarise(Counts_total= sum(n,na.rm=TRUE),H3K27ac = sum(H3K27ac_overlap.x,na.rm=TRUE),H3K27ac_no = sum(H3K27ac_overlap.y,na.rm=TRUE),
                                                                       H3K4me3 = sum(H3K4me3_overlap.x,na.rm=TRUE),H3K4me3_no = sum(H3K4me3_overlap.y,na.rm=TRUE),
                                                                       H3K4me1 = sum(H3K4me1_overlap.x,na.rm=TRUE),H3K4me1_no = sum(H3K4me1_overlap.y,na.rm=TRUE),
                                                                       H3K27me3 = sum(H3K27me3_overlap.x,na.rm=TRUE),H3K27me3_no = sum(H3K27me3_overlap.y,na.rm=TRUE),
                                                                       H3K27ac_H3K4me1 = sum(H3K27ac_H3K4me1_overlap.x,na.rm=TRUE),H3K27ac_H3K4me1_no = sum(H3K27ac_H3K4me1_overlap.y,na.rm=TRUE),
                                                                       TEexp = sum(TEexp_overlap.x,na.rm=TRUE),TEexp_no = sum(TEexp_overlap.y,na.rm=TRUE),
                                                                       Geneexp = sum(Geneexp_overlap.x,na.rm=TRUE),Geneexp_no = sum(Geneexp_overlap.y,na.rm=TRUE)))

TEinstance_summary_Flu_withATAC_summary_final_plot = data.frame("ID_final"=NA,"Counts_total"=NA,"Counts_overlap"=NA,"Counts_notoverlap"=NA,"Mark"=NA,"perC" = NA,"Min"=NA,"Max"=NA,"perC_Zscaled"=NA,"perC_scaled_1_100"=NA)
i = 3
for(i in seq(3,ncol(TEinstance_summary_Flu_withATAC_summary_final),2)){
  TEinstance_summary_Flu_withATAC_summary_final_tmp = TEinstance_summary_Flu_withATAC_summary_final[,c(1,2,i,(i+1))]
  TEinstance_summary_Flu_withATAC_summary_final_tmp$Mark = colnames(TEinstance_summary_Flu_withATAC_summary_final_tmp)[3]
  TEinstance_summary_Flu_withATAC_summary_final_tmp$perC = TEinstance_summary_Flu_withATAC_summary_final_tmp[,3]/(TEinstance_summary_Flu_withATAC_summary_final_tmp[,3]+TEinstance_summary_Flu_withATAC_summary_final_tmp[,4])
  TEinstance_summary_Flu_withATAC_summary_final_tmp$Min = min(TEinstance_summary_Flu_withATAC_summary_final_tmp$perC)
  TEinstance_summary_Flu_withATAC_summary_final_tmp$Max = max(TEinstance_summary_Flu_withATAC_summary_final_tmp$perC)
  TEinstance_summary_Flu_withATAC_summary_final_tmp$perC_Zscaled = scale(TEinstance_summary_Flu_withATAC_summary_final_tmp$perC)
  TEinstance_summary_Flu_withATAC_summary_final_tmp$perC_scaled_1_100 = rescale(TEinstance_summary_Flu_withATAC_summary_final_tmp$perC, to = c(0.01, 1)) 
  
  colnames(TEinstance_summary_Flu_withATAC_summary_final_tmp)[3] = "Counts_overlap"
  colnames(TEinstance_summary_Flu_withATAC_summary_final_tmp)[4] = "Counts_notoverlap"
  
  if (nrow(TEinstance_summary_Flu_withATAC_summary_final_plot) == 1){
    TEinstance_summary_Flu_withATAC_summary_final_plot = TEinstance_summary_Flu_withATAC_summary_final_tmp
  } else {
    TEinstance_summary_Flu_withATAC_summary_final_plot = rbind(TEinstance_summary_Flu_withATAC_summary_final_plot,TEinstance_summary_Flu_withATAC_summary_final_tmp)
  }
}

## 
TEinstance_summary_Flu_withATAC_summary_final_plot$ID_final2 = TEinstance_summary_Flu_withATAC_summary_final_plot$ID_final
TEinstance_summary_Flu_withATAC_summary_final_plot = data.frame(cSplit(TEinstance_summary_Flu_withATAC_summary_final_plot,"ID_final2",sep=" ",type.convert = as.character))

Subfamily_FluOrder = TElist[TElist$TE_cluster!="NI",]
Subfamily_FluOrder = Subfamily_FluOrder[order(Subfamily_FluOrder$Reordered_by_motif),]
Subfamily_FluOrder = as.character(Subfamily_FluOrder$TEfamily)
TEinstance_summary_Flu_withATAC_summary_final_plot$ID_final = factor(TEinstance_summary_Flu_withATAC_summary_final_plot$ID_final,levels=rev(Subfamily_FluOrder))

Mark_names = c("H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K27ac_H3K4me1","TEexp","Geneexp")
TEinstance_summary_Flu_withATAC_summary_final_plot$Mark = factor(TEinstance_summary_Flu_withATAC_summary_final_plot$Mark,levels=Mark_names)
TEinstance_summary_Flu_withATAC_summary_final_plot$ID_final2_2 = TEinstance_summary_Flu_withATAC_summary_final_plot$ID_final
levels(TEinstance_summary_Flu_withATAC_summary_final_plot$ID_final2_2) = c((2:length(levels(TEinstance_summary_Flu_withATAC_summary_final_plot$ID_final2_2))-1),2,2) 
TEinstance_summary_Flu_withATAC_summary_final_plot$Mark2 = TEinstance_summary_Flu_withATAC_summary_final_plot$Mark
levels(TEinstance_summary_Flu_withATAC_summary_final_plot$Mark2) = c(1:(length(levels(TEinstance_summary_Flu_withATAC_summary_final_plot$Mark2))-1),1,1)
TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2 = ifelse(TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap>0 & TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap <=10,1,NA)
TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2 = ifelse(TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap>10 & TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap <=50,2,TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2)
TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2 = ifelse(TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap>50 & TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap <=100,3,TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2)
TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2 = ifelse(TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap>100 & TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap <=250,4,TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2)
TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2 = ifelse(TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap>250 & TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap <=500,5,TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2)
TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2 = ifelse(TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap>500 & TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap <=1000,6,TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2)
TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2 = ifelse(TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap>1000,7,TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2)
TEinstance_summary_Flu_withATAC_summary_final_plot$Shape = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2),"absent","present")
TEinstance_summary_Flu_withATAC_summary_final_plot$Shape = factor(TEinstance_summary_Flu_withATAC_summary_final_plot$Shape)
TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap3 = ifelse(is.na(TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2),5,TEinstance_summary_Flu_withATAC_summary_final_plot$Counts_overlap2)
TEinstance_summary_Flu_withATAC_summary_final_plot = merge(TEinstance_summary_Flu_withATAC_summary_final_plot,TElist,by.x="ID_final",by.y="TEfamily",all.x=T)
TEinstance_summary_Flu_withATAC_summary_methy_subfamily = merge(TEinstance_summary_Flu_withATAC_summary_methy_subfamily,TElist,by.x="TE_family.x",by.y="TEfamily",all.x=T)
# range

TEinstance_summary_Flu_withATAC_summary_final_range = TEinstance_summary_Flu_withATAC_summary_final_plot[!duplicated(TEinstance_summary_Flu_withATAC_summary_final_plot$Mark,TEinstance_summary_Flu_withATAC_summary_final_plot$Min,TEinstance_summary_Flu_withATAC_summary_final_plot$Max),c("Mark","Min","Max")]
TEinstance_summary_Flu_withATAC_summary_final_range$Filter = Filter
TEinstance_summary_Flu_withATAC_summary_final_range$Mark = as.character(TEinstance_summary_Flu_withATAC_summary_final_range$Mark)
TEinstance_summary_Flu_withATAC_summary_final_range[nrow(TEinstance_summary_Flu_withATAC_summary_final_range)+1,] = NA
TEinstance_summary_Flu_withATAC_summary_final_range[nrow(TEinstance_summary_Flu_withATAC_summary_final_range),]$Mark = "DNA_methyl"
TEinstance_summary_Flu_withATAC_summary_final_range[nrow(TEinstance_summary_Flu_withATAC_summary_final_range),]$Min = min(TEinstance_summary_Flu_withATAC_summary_methy_subfamily$DNA_methyl.Level)
TEinstance_summary_Flu_withATAC_summary_final_range[nrow(TEinstance_summary_Flu_withATAC_summary_final_range),]$Max = max(TEinstance_summary_Flu_withATAC_summary_methy_subfamily$DNA_methyl.Level)
TEinstance_summary_Flu_withATAC_summary_final_range[nrow(TEinstance_summary_Flu_withATAC_summary_final_range),]$Filter = Filter
TEinstance_summary_Flu_withATAC_summary_final_rangeAll = TEinstance_summary_Flu_withATAC_summary_final_range

Shapes = c("present"=19,"absent"=4)

TEinstance_summary_Flu_withATAC_summary_final_plot2 = TEinstance_summary_Flu_withATAC_summary_final_plot[,c("ID_final","Mark","Mark2","perC_scaled_1_100","Counts_overlap3","Shape")]
TEinstance_summary_Flu_withATAC_summary_final_plot2$perC_F = TEinstance_summary_Flu_withATAC_summary_final_plot2$perC_scaled_1_100*100
TEinstance_summary_Flu_withATAC_summary_final_plot2$CountsF = TEinstance_summary_Flu_withATAC_summary_final_plot2$Counts_overlap3

### Add the methylation
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2 = TEinstance_summary_Flu_withATAC_summary_methy_subfamily[,1:5]
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$perC_F = rescale(TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Level, to = c(0.01, 1))
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$perC_F = TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$perC_F*100
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF = ifelse(TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count<=1,1,NA)
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF = ifelse(TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count>1 & TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count<=5,2,TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF)
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF = ifelse(TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count>5 & TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count<=10,3,TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF)
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF = ifelse(TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count>10 & TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count<=25,4,TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF)
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF = ifelse(TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count>25 & TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count<=50,5,TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF)
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF = ifelse(TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$DNA_methyl.Count>50,6,TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$CountsF)
colnames(TEinstance_summary_Flu_withATAC_summary_methy_subfamily2) = c("ID_final_tmp","ID_final","perC_scaled_1_100",
                                                                       "TE_family2","Counts_overlap3","perC_F","CountsF")
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$ID_final = gsub(" $","",TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$ID_final)
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$Mark = "DNA_methyl"
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$Mark2 = 7
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2$Shape = "present"
TEinstance_summary_Flu_withATAC_summary_methy_subfamily2 = TEinstance_summary_Flu_withATAC_summary_methy_subfamily2[,c("ID_final","Mark","Mark2","perC_scaled_1_100","Counts_overlap3","Shape","perC_F","CountsF")]

### combined the table for plot
TEinstance_summary_Flu_withATAC_summary_final_plot2$ID_final = as.character(TEinstance_summary_Flu_withATAC_summary_final_plot2$ID_final)
TEinstance_summary_Flu_withATAC_summary_final_plot2$Mark2 = as.numeric(as.character(TEinstance_summary_Flu_withATAC_summary_final_plot2$Mark2))
TEinstance_summary_Flu_withATAC_summary_final_plot2 = rbind(TEinstance_summary_Flu_withATAC_summary_final_plot2,TEinstance_summary_Flu_withATAC_summary_methy_subfamily2)

### make factors
TEinstance_summary_Flu_withATAC_summary_final_plot2$ID_final = factor(TEinstance_summary_Flu_withATAC_summary_final_plot2$ID_final,levels=rev(Subfamily_FluOrder))
### order different marks
Mark_names = c("Geneexp","TEexp","DNA_methyl","H3K27ac","H3K4me1","H3K4me3","H3K27me3","H3K27ac_H3K4me1")
Mark_names2 = c("Geneexp","TEexp","DNA_methyl","H3K27ac","H3K4me1","H3K4me3","H3K27ac_H3K4me1","H3K27me3")

TEinstance_summary_Flu_withATAC_summary_final_plot2$Mark = factor(TEinstance_summary_Flu_withATAC_summary_final_plot2$Mark,levels=Mark_names)

TEinstance_summary_Flu_withATAC_summary_final_plot2$ID_final2_2 = TEinstance_summary_Flu_withATAC_summary_final_plot2$ID_final
levels(TEinstance_summary_Flu_withATAC_summary_final_plot2$ID_final2_2) = c((2:length(levels(TEinstance_summary_Flu_withATAC_summary_final_plot2$ID_final2_2))-1),2,2) 

Mark_orders = c("repressed","array")
for (Mark_order in Mark_orders){
  if (Mark_order == "repressed") {
    TEinstance_summary_Flu_withATAC_summary_final_plot2$Mark = factor(TEinstance_summary_Flu_withATAC_summary_final_plot2$Mark,levels=Mark_names2)
  } else {
    TEinstance_summary_Flu_withATAC_summary_final_plot2$Mark = factor(TEinstance_summary_Flu_withATAC_summary_final_plot2$Mark,levels=Mark_names)
  }
  gplist = list()
  Order = 1
  p1 = ggplot(TEinstance_summary_Flu_withATAC_summary_final_plot2,aes(x = Mark, y = ID_final)) +        ## global aes
    geom_point(aes(color=perC_F,size = CountsF,shape=Shape))  +    ## geom_point for circle illusion
    scale_size(range = c(1,5.5)) + #geom_text(aes(label = pValue.adj_label),col = "grey44") +
    scale_color_gradient2(high="#67000d",mid = "white", low = "#08306b")+
    scale_shape_manual(values = Shapes)+
    xlab("")+
    ylab("")+
    geom_hline(aes(yintercept = as.numeric(as.character(ID_final2_2))+0.5),color="grey44")+
    geom_hline(aes(yintercept = 15+0.5),color="black",size=0.5)+
    geom_vline(aes(xintercept = as.numeric(as.character(Mark2))+0.5),color="grey44")+
    theme(panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black",size = 0.8,fill=NA),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text=element_text(colour="black",size=rel(.8),angle = 0),
          axis.text.x=element_text(colour="black",angle = 90,hjust = 1,vjust = .5),
          axis.title=element_text(colour="black",size=rel(1.5)),
          legend.position="right",
          legend.background = element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.text=element_text(size=rel(1.2)))
  gplist[[Order]] = p1
  Order = Order + 1 
  p1 = ggplot(TEinstance_summary_Flu_withATAC_summary_final_plot2,aes(x = Mark, y = ID_final)) +        ## global aes
    geom_point(aes(color=perC_F,size = CountsF,shape=Shape))  +    ## geom_point for circle illusion
    scale_size(range = c(1,5.5)) + #geom_text(aes(label = pValue.adj_label),col = "grey44") +
    scale_color_gradient2(high="#08306b",mid = "white", low = "#08306b")+
    scale_shape_manual(values = Shapes)+
    xlab("")+
    ylab("")+
    geom_hline(aes(yintercept = as.numeric(as.character(ID_final2_2))+0.5),color="grey44")+
    geom_hline(aes(yintercept = 15+0.5),color="black",size=0.5)+
    geom_vline(aes(xintercept = as.numeric(as.character(Mark2))+0.5),color="grey44")+
    theme(panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black",size = 0.8,fill=NA),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text=element_text(colour="black",size=rel(.8),angle = 0),
          axis.text.x=element_text(colour="black",angle = 90,hjust = 1,vjust = .5),
          axis.title=element_text(colour="black",size=rel(1.5)),
          legend.position="right",
          legend.background = element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.text=element_text(size=rel(1.2)))
  gplist[[Order]] = p1
  Order = Order + 1 
  p1 = ggplot(TEinstance_summary_Flu_withATAC_summary_final_plot2,aes(x = Mark, y = ID_final)) +        ## global aes
    geom_point(aes(color=perC_F,size = CountsF,shape=Shape))  +    ## geom_point for circle illusion
    scale_size(range = c(1,5.5)) + #geom_text(aes(label = pValue.adj_label),col = "grey44") +
    scale_color_gradient2(high="#00441b",mid = "white", low = "#08306b")+
    scale_shape_manual(values = Shapes)+
    xlab("")+
    ylab("")+
    geom_hline(aes(yintercept = as.numeric(as.character(ID_final2_2))+0.5),color="grey44")+
    geom_hline(aes(yintercept = 15+0.5),color="black",size=0.5)+
    geom_vline(aes(xintercept = as.numeric(as.character(Mark2))+0.5),color="grey44")+
    theme(panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black",size = 0.8,fill=NA),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text=element_text(colour="black",size=rel(.8),angle = 0),
          axis.text.x=element_text(colour="black",angle = 90,hjust = 1,vjust = .5),
          axis.title=element_text(colour="black",size=rel(1.5)),
          legend.position="right",
          legend.background = element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.text=element_text(size=rel(1.2)))
  gplist[[Order]] = p1
  Order = Order + 1 
  pdf(paste("EMC-Figure4C-2021_10_26_",Filter,"-",Mark_order,".pdf",sep=""),    # create PNG for the heat map        
      width = 4.1*3,        # 5 x 300 pixels
      height = 8,
      pointsize = 10,useDingbats = F )        # smaller font size
  do.call("grid.arrange",c(gplist,ncol=3))
  dev.off()
}
