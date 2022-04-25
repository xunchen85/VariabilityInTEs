###     Author:     Xun Chen
###     Date:       2022/3/7
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
library(ggpmisc)
library(gridExtra)
library(cowplot)
library(ggforce)
library("viridis") 
library(RColorBrewer)
library(tidyverse)

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
# TE list
TElist = read.delim("inputs/TElist.sep.27",header=T,sep="")
# bed.out file
INTname = "2020_12_21_Aggregate_noEU37_100.bed.out"                    ############ motif analysis I excluded EU37
# TE consensus length
TEconsensusLen = read.delim("inputs/76subfamily_all_Consensus_len.out",header=F,sep="")

# global TF activity
TF_activity = read.delim("inputs/TF_activity_results.txt",sep="",header=T)

## TE consensus length
colnames(TEconsensusLen) = c("Consensus","Consensus_Len_max")
TEconsensusLen$Consensus = as.character(TEconsensusLen$Consensus)
# add LTR39
TEconsensusLen[nrow(TEconsensusLen)+1,] = NA
TEconsensusLen[nrow(TEconsensusLen),1] = "LTR39"
TEconsensusLen[nrow(TEconsensusLen),2] = 783
# add MSTB1
TEconsensusLen[nrow(TEconsensusLen)+1,] = NA
TEconsensusLen[nrow(TEconsensusLen),1] = "MSTB1"
TEconsensusLen[nrow(TEconsensusLen),2] = 432

Cluster_Consensus_Flu = read.csv(paste("inputs/","EMC_Flu_",INTname,"_allConsensus",sep=""))
Cluster_Consensus_NI = read.csv(paste("inputs/","EMC_NI_",INTname,"_allConsensus",sep=""))

Motif_Flu = read.delim(paste("inputs/","EMC_Flu_",INTname,sep=""),header=F,sep="") 

### counts in enhanced subfamilies
Motif_Flu_enhanced = Motif_Flu
Motif_Flu_enhanced = cSplit(Motif_Flu_enhanced,"V5",sep=":",type.convert = as.character)
Motif_Flu_enhanced = Motif_Flu_enhanced[Motif_Flu_enhanced$V5_1 %in% TElist[TElist$Group == "Flu_enriched",]$TEfamily,]

colnames(Motif_Flu) = c("File","Chr","Start_instance","End_instance","TEinfo","Chr2","Start_summit","End_summit","TEinfo_summit","Score","Strand","Posi","TE_coordinate_1bp","Strand2","Consensus","Ref_posi","Consensus_posi")
Motif_Flu$TE_coordinate_0bp = paste(Motif_Flu$Chr,":",Motif_Flu$Start_instance,"-",Motif_Flu$End_instance,sep="")
Motif_NI = read.delim(paste("inputs/","EMC_NI_",INTname,sep=""),header=F,sep="") 
colnames(Motif_NI) = c("File","Chr","Start_instance","End_instance","TEinfo","Chr2","Start_summit","End_summit","TEinfo_summit","Score","Strand","Posi","TE_coordinate_1bp","Strand2","Consensus","Ref_posi","Consensus_posi")
Motif_NI$TE_coordinate_0bp = paste(Motif_NI$Chr,":",Motif_NI$Start_instance,"-",Motif_NI$End_instance,sep="")

Motif_cluster_Flu = read.csv(paste("inputs/","EMC_Flu_",INTname,"_allCluster",sep=""))
Motif_cluster_Flu$TEcoordinate_0bp = paste(Motif_cluster_Flu$Chr,":",Motif_cluster_Flu$Start_instance,"-",Motif_cluster_Flu$End_instance,sep="")
Motif_cluster_Flu = Motif_cluster_Flu[,c("TEcoordinate_0bp","TEsubfamily","TEinstance","TEfamily","TEgroup","Cluster_50bp","Win_50bp","Cluster_100bp","Win_100bp","Cluster_200bp","Win_200bp")]
Motif_cluster_NI = read.csv(paste("inputs/","EMC_NI_",INTname,"_allCluster",sep=""))
Motif_cluster_NI$TEcoordinate_0bp = paste(Motif_cluster_NI$Chr,":",Motif_cluster_NI$Start_instance,"-",Motif_cluster_NI$End_instance,sep="")
Motif_cluster_NI = Motif_cluster_NI[,c("TEcoordinate_0bp","TEsubfamily","TEinstance","TEfamily","TEgroup","Cluster_50bp","Win_50bp","Cluster_100bp","Win_100bp","Cluster_200bp","Win_200bp")]

Motif_Flu_full = merge(Motif_Flu,Motif_cluster_Flu,by.x="TE_coordinate_0bp",by.y="TEcoordinate_0bp",all.x=T)
Motif_NI_full = merge(Motif_NI,Motif_cluster_NI,by.x="TE_coordinate_0bp",by.y="TEcoordinate_0bp",all.x=T)
Motif_Flu_full = Motif_Flu_full[!duplicated(Motif_Flu_full$TE_coordinate_0bp),]
Motif_NI_full = Motif_NI_full[!duplicated(Motif_NI_full$TE_coordinate_0bp),]

### counts in reduced subfamilies
Motif_NI_reduced = Motif_NI
Motif_NI_reduced = cSplit(Motif_NI_reduced,"TEinfo",sep=":",type.convert = as.character)
Motif_NI_reduced = Motif_NI_reduced[Motif_NI_reduced$TEinfo_1 %in% TElist[TElist$Group == "NI_enriched",]$TEfamily,]
head(Motif_NI_reduced)

### delete
rm(Motif_Flu)
rm(Motif_NI)

####################################
################################# Section 1
##### 2.1 Flu or NI group peak clustering
peakTypes = c("Flupeak","NIpeak")
peakType = "NIpeak"
TEgroup = "Fluenriched"


TopCount = 5
peakType = "Flupeak"
WinSize = "100bp"
Proportion_t = 0.2
Proportion_sum = 0.5

Motif_data_Flu = read.delim(paste("inputs/","EMC_Flu_",INTname,".fimo.uniqueCount_100bpCluster.summary",sep=""),sep="",header=F)
Motif_data_NI = read.delim(paste("inputs/","EMC_NI_",INTname,".fimo.uniqueCount_100bpCluster.summary",sep=""),sep="",header=F)

for (peakType in peakTypes){
  if (peakType == "Flupeak") {
    TEgroup = "Fluenriched"
    Motif_data = Motif_data_Flu
    colnames(Motif_data) = c("motifCount","motifName","TEsubfamily","Consensus","Cluster","Win_position")
    Motif_data$Cluster_all = "No"
    Motif_data$Cluster_all = ifelse(Motif_data$Cluster == "Cluster_1" |Motif_data$Cluster == "Cluster_2" |Motif_data$Cluster == "Cluster_3" |Motif_data$Cluster == "Cluster_4" |Motif_data$Cluster == "Cluster_5",as.character(Motif_data$Cluster),Motif_data$Cluster_all)
    Motif_data$Cluster_all = ifelse(grepl("Cluster_",Motif_data$Cluster) & Motif_data$Cluster_all == "No","Cluster_others",Motif_data$Cluster_all)
    Motif_data$UniqueID = paste(Motif_data$TEsubfamily,Motif_data$Cluster_all)
    ## cluster info
    Motif_cluster_tmp = Motif_cluster_Flu
    if (WinSize == "100bp"){
      Motif_cluster_tmp$Cluster_tmp = Motif_cluster_tmp$Cluster_100bp
      Motif_cluster_tmp$Win_all = Motif_cluster_tmp$Win_100bp
      
    } else if (WinSize == "200bp"){
      Motif_cluster_tmp$Cluster_tmp = Motif_cluster_tmp$Cluster_200bp
      Motif_cluster_tmp$Win_all = Motif_cluster_tmp$Win_200bp
      
    } else if (WinSize == "50bp"){
      Motif_cluster_tmp$Cluster_tmp = Motif_cluster_tmp$Cluster_50bp
      Motif_cluster_tmp$Win_all = Motif_cluster_tmp$Win_50bp
    }
    Motif_cluster_tmp$Cluster_all = "No"
    Motif_cluster_tmp$Cluster_all = ifelse(Motif_cluster_tmp$Cluster_tmp == "Cluster_1" |Motif_cluster_tmp$Cluster_tmp == "Cluster_2" |Motif_cluster_tmp$Cluster_tmp == "Cluster_3" |Motif_cluster_tmp$Cluster_tmp == "Cluster_4" |Motif_cluster_tmp$Cluster_tmp == "Cluster_5",as.character(Motif_cluster_tmp$Cluster_tmp),Motif_cluster_tmp$Cluster_all)
    Motif_cluster_tmp$Cluster_all = ifelse(grepl("Cluster_",Motif_cluster_tmp$Cluster_tmp) & Motif_cluster_tmp$Cluster_all == "No","Cluster_others",Motif_cluster_tmp$Cluster_all)
    
    ## count
    Motif_totalCount = data.frame(Motif_cluster_tmp %>% group_by(TEsubfamily,Cluster_all) %>% dplyr::count())
    Motif_totalCount$UniqueID = paste(Motif_totalCount$TEsubfamily,Motif_totalCount$Cluster_all)
    Motif_data = merge(Motif_data,Motif_totalCount,by="UniqueID",all.x=T)
    Motif_data$Proportion = as.numeric(as.character(Motif_data$motifCount))/as.numeric(as.character(Motif_data$n))
    unique(Motif_data$UniqueID)
    ## Filtering record without clusters
    Motif_data_filtered = Motif_data[!is.na(Motif_data$Cluster_all.x),]
    
    ## keep Flu-enriched TEs
    Motif_data_filtered = Motif_data_filtered[as.character(Motif_data_filtered$TEsubfamily.x) %in% as.character(TElist[TElist$TE_cluster!="NI",]$TEfamily),]
  } else if (peakType == "NIpeak"){
    TEgroup = "NIenriched"
    Motif_data = Motif_data_NI
    colnames(Motif_data) = c("motifCount","motifName","TEsubfamily","Consensus","Cluster","Win_position")
    Motif_data$Cluster_all = "No"
    Motif_data$Cluster_all = ifelse(Motif_data$Cluster == "Cluster_1" |Motif_data$Cluster == "Cluster_2" |Motif_data$Cluster == "Cluster_3" |Motif_data$Cluster == "Cluster_4" |Motif_data$Cluster == "Cluster_5",as.character(Motif_data$Cluster),Motif_data$Cluster_all)
    Motif_data$Cluster_all = ifelse(grepl("Cluster_",Motif_data$Cluster) & Motif_data$Cluster_all == "No","Cluster_others",Motif_data$Cluster_all)
    Motif_data$UniqueID = paste(Motif_data$TEsubfamily,Motif_data$Cluster_all)
    
    ## cluster info
    Motif_cluster_tmp = Motif_cluster_NI
    if (WinSize == "100bp"){
      Motif_cluster_tmp$Cluster_tmp = Motif_cluster_tmp$Cluster_100bp
      Motif_cluster_tmp$Win_all = Motif_cluster_tmp$Win_100bp
      
    } else if (WinSize == "200bp"){
      Motif_cluster_tmp$Cluster_tmp = Motif_cluster_tmp$Cluster_200bp
      Motif_cluster_tmp$Win_all = Motif_cluster_tmp$Win_200bp
      
    } else if (WinSize == "50bp"){
      Motif_cluster_tmp$Cluster_tmp = Motif_cluster_tmp$Cluster_50bp
      Motif_cluster_tmp$Win_all = Motif_cluster_tmp$Win_50bp
    }
    Motif_cluster_tmp$Cluster_all = "No"
    Motif_cluster_tmp$Cluster_all = ifelse(Motif_cluster_tmp$Cluster_tmp == "Cluster_1" |Motif_cluster_tmp$Cluster_tmp == "Cluster_2" |Motif_cluster_tmp$Cluster_tmp == "Cluster_3" |Motif_cluster_tmp$Cluster_tmp == "Cluster_4" |Motif_cluster_tmp$Cluster_tmp == "Cluster_5",as.character(Motif_cluster_tmp$Cluster_tmp),Motif_cluster_tmp$Cluster_all)
    Motif_cluster_tmp$Cluster_all = ifelse(grepl("Cluster_",Motif_cluster_tmp$Cluster_tmp) & Motif_cluster_tmp$Cluster_all == "No","Cluster_others",Motif_cluster_tmp$Cluster_all)
    
    ## count
    Motif_totalCount = data.frame(Motif_cluster_tmp %>% group_by(TEsubfamily,Cluster_all) %>% dplyr::count())
    Motif_totalCount$UniqueID = paste(Motif_totalCount$TEsubfamily,Motif_totalCount$Cluster_all)
    
    Motif_data = merge(Motif_data,Motif_totalCount,by="UniqueID",all.x=T)
    Motif_data$Proportion = as.numeric(as.character(Motif_data$motifCount))/as.numeric(as.character(Motif_data$n))
    
    ## Filtering record without clusters
    Motif_data_filtered = Motif_data[!is.na(Motif_data$Cluster_all.x),]
    
    ## keep NI-enriched TEs
    Motif_data_filtered = Motif_data_filtered[as.character(Motif_data_filtered$TEsubfamily.x) %in% as.character(TElist[TElist$TE_cluster=="NI",]$TEfamily),]
  } 

  Motif_data_filtered = merge(Motif_data_filtered,TElist,by.x="TEsubfamily.x",by.y="TEfamily",all.x=T)
  Motif_data_filtered = Motif_data_filtered[order(Motif_data_filtered$Reordered_by_motif,Motif_data_filtered$Cluster_all.x,-Motif_data_filtered$Proportion),]
  
  ## Kept top 5 per cluster
  Motif_data_filtered_top10 = data.frame(Motif_data_filtered %>% group_by(UniqueID) %>% slice_max(order_by = Proportion, n = TopCount))
  Motif_data_filtered_top10$ColumnOrder = c(1:nrow(Motif_data_filtered_top10))
  Motif_data_filtered_top10 = data.frame(Motif_data_filtered_top10 %>% group_by(UniqueID) %>% slice_max(order_by = ColumnOrder, n = TopCount))
  Motif_data_filtered_top10$UniqueID = factor(Motif_data_filtered_top10$UniqueID,levels=(unique(as.character(Motif_data_filtered$UniqueID))))
  Motif_data_filtered_top10$Order = Motif_data_filtered_top10$UniqueID
  Motif_data_filtered_top10$Line_posi = NA
  
  ## Filtering record with motif proportion < 0.2
  Motif_data_filtered_top10 = Motif_data_filtered_top10[Motif_data_filtered_top10$Proportion>=Proportion_t,]
  
  ## Filtering by the sum proportion
  Motif_data_filtered_top10_tmp = aggregate(Proportion ~ motifName, Motif_data_filtered_top10, sum)
  Motif_data_filtered_top10_tmp = Motif_data_filtered_top10_tmp[Motif_data_filtered_top10_tmp$Proportion>=Proportion_sum,]
  Motif_data_filtered_top10 = Motif_data_filtered_top10[Motif_data_filtered_top10$motifName %in% Motif_data_filtered_top10_tmp$motifName,]
  levels(Motif_data_filtered_top10$Order) = c(1:length(levels(Motif_data_filtered_top10$UniqueID)))
  rm(Motif_data_filtered_top10_tmp)
  
  ## y axis label
  Motif_data_filtered_top10[1,]$Line_posi = as.numeric(as.character(Motif_data_filtered_top10[1,]$Order))
  for (i in 2:nrow(Motif_data_filtered_top10)){
    if (as.character(Motif_data_filtered_top10[i,]$TEsubfamily.x) != as.character(Motif_data_filtered_top10[i-1,]$TEsubfamily.x)){
      Motif_data_filtered_top10[i,]$Line_posi = as.numeric(as.character(Motif_data_filtered_top10[i,]$Order))
    }
  }
  Motif_data_filtered_top10$Label_posi = NA
  Motif_data_filtered_top10$Label_posi = ifelse(!is.na(Motif_data_filtered_top10$Line_posi),as.character(Motif_data_filtered_top10$TEsubfamily.x),Motif_data_filtered_top10$Label_posi)
  Motif_data_filtered_top10_y = Motif_data_filtered_top10[!is.na(Motif_data_filtered_top10$Line_posi),]
  
  ## plot
  TEsubfamily_order = unique(as.character(Motif_data_filtered_top10[order(Motif_data_filtered_top10$Reordered_by_motif),]$TEsubfamily.x))
  Motif_data_filtered_top10$TEsubfamily.x = factor(Motif_data_filtered_top10$TEsubfamily.x,levels=TEsubfamily_order)
  Motif_data_filtered_top10$Cluster_all.x = factor(Motif_data_filtered_top10$Cluster_all.x,levels = rev(unique(as.character(Motif_data_filtered_top10$Cluster_all.x))))

  ## Add up
  Motif_data_filtered_top10_sum = data.frame(Motif_data_filtered_top10 %>% group_by(TEsubfamily.x,motifName,Reordered_by_motif,TE_cluster) %>% 
                                               dplyr::summarise(motifCount = sum(motifCount), totalInstance = sum(n)))
  Motif_data_filtered_top10_sum$ID = paste(Motif_data_filtered_top10_sum$TEsubfamily.x,Motif_data_filtered_top10_sum$motifName)
  
  ## achieve unique and top cluster
  Motif_data_filtered_top10_sum_tmp = Motif_data_filtered_top10[!duplicated(paste(Motif_data_filtered_top10$TEsubfamily.x,Motif_data_filtered_top10$motifName,sep="")),c("TEsubfamily.x","motifName","Cluster_all.y")]
  Motif_data_filtered_top10_sum_tmp$ID = paste(Motif_data_filtered_top10_sum_tmp$TEsubfamily.x,Motif_data_filtered_top10_sum_tmp$motifName)
  Motif_data_filtered_top10_sum_tmp = Motif_data_filtered_top10_sum_tmp[,c("ID","Cluster_all.y")]
  
  ## combine 
  Motif_data_filtered_top10_sum = merge(Motif_data_filtered_top10_sum,Motif_data_filtered_top10_sum_tmp,by="ID",all.x=T)
  ## re-calculate the proportion
  Motif_data_filtered_top10_sum$Proportion = Motif_data_filtered_top10_sum$motifCount/Motif_data_filtered_top10_sum$totalInstance
  rm(Motif_data_filtered_top10_sum_tmp)
  ## order by TE family
  Motif_data_filtered_top10_sum$TEsubfamily.x = factor(Motif_data_filtered_top10_sum$TEsubfamily.x,levels=rev(TEsubfamily_order))
  
  ### filter out by the total instances (>= 50)
  Motif_data_filtered_top10_sum_tmp = data.frame(Motif_data_filtered_top10_sum %>% group_by(motifName) %>% 
                                                   dplyr::summarise(motifTotalCount = sum(motifCount), totalInstanceCount = sum(totalInstance)))
  Motif_data_filtered_top10_sum_tmp = Motif_data_filtered_top10_sum_tmp[Motif_data_filtered_top10_sum_tmp$totalInstanceCount>=50,]
  
  ### filter by the total count of instances (>= 5)
  Motif_data_filtered_top10_sum_filterByCounts = Motif_data_filtered_top10_sum[Motif_data_filtered_top10_sum$totalInstance>=5,]
  
  ### further filter by total count of motifs
  Motif_data_filtered_top10_sum_filterByCounts = Motif_data_filtered_top10_sum_filterByCounts[Motif_data_filtered_top10_sum_filterByCounts$motifName %in% Motif_data_filtered_top10_sum_tmp$motifName,]
  
  ###
  Motif_data_filtered_top10_sum_filterByCounts$motifName = factor(Motif_data_filtered_top10_sum_filterByCounts$motifName)
  p1 = ggplot(Motif_data_filtered_top10_sum_filterByCounts, aes(y = TEsubfamily.x, x = motifName)) +        ## global aes
    geom_point(aes(color=Cluster_all.y,size =Proportion),shape=19)  +    ## geom_point for circle illusion (previously i use shape 16)
    scale_size(range = c(1,4))+             ## to tune the size of circles
    ylab("TE family")+
    xlab("TF")+
    scale_color_manual(values=Cluster_color)+
    theme(panel.grid.minor = element_line(colour = "blue", size = 2),
      plot.margin=unit(c(10,10,10,10),"mm"),
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.text=element_text(colour="black",size=rel(0.8),angle = 0),
      axis.title=element_text(colour="black",size=rel(1)),
      legend.position="right",
      legend.text=element_text(size=rel(0.8)),
      axis.text.x=element_text(angle=90,hjust = 0.95,size=rel(0.75),vjust = .5),
      legend.key = element_rect(fill = NA))
  
  ## TF activity heatmap
  TF_activity_sub = TF_activity[TF_activity$TF_name %in% Motif_data_filtered_top10_sum_filterByCounts$motifName,]
  TF_activity_sub$Col = 1
  TF_activity_sub$TF_name = factor(TF_activity_sub$TF_name,levels=levels(Motif_data_filtered_top10_sum_filterByCounts$motifName))
  TF_activity_final = data.frame("label"=levels(Motif_data_filtered_top10_sum_filterByCounts$motifName),order = 1)
  TF_activity_final$order = 1:nrow(TF_activity_final)
  TF_activity_final = merge(TF_activity_final,TF_activity_sub,by.x="label",by.y="TF_name",all.x=T)
  p2 = ggplot(TF_activity_final, aes(label,Col)) +
    geom_tile(aes(fill = Mean_activity),color= "white",size=0.1) +
    scale_fill_gradient2(low="#5e4fa2",high = "#9e0142")+
    scale_x_discrete(drop=FALSE)+
    xlab("")+
    ylab("")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text=element_text(colour="black",size=rel(1),angle = 0),
          axis.text.y=element_blank(),
          axis.text.x=element_text(colour="black",size=rel(1),angle = 90,vjust=.5, hjust=1),
          axis.title=element_text(colour="black",size=rel(1)),
          legend.position="right",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1.2)))
  
  ## For each region separately
  p3 = ggplot(Motif_data_filtered_top10[Motif_data_filtered_top10$n>=5 & Motif_data_filtered_top10$motifName %in% Motif_data_filtered_top10_sum_tmp$motifName,], aes(x = motifName, y = (Order))) +        ## global aes
    geom_hline(data=Motif_data_filtered_top10_y,aes(yintercept = as.character(Line_posi)),color = "white")+
    geom_point(aes(color=Cluster_all.x,fill=Cluster_all.x,size =Proportion),shape=19)  +    ## geom_point for circle illusion
    scale_size(range = c(1,4))+             ## to tune the size of circles
    ylab("TE family")+
    xlab("TF")+
    scale_color_manual(values=Cluster_color)+
    scale_fill_manual(values=Cluster_color)+
    scale_y_discrete(limits=rev,labels= Motif_data_filtered_top10_y$Label_posi,breaks = c(as.character(Motif_data_filtered_top10_y$Line_posi)),drop=TRUE)+
    theme(panel.background = element_rect(fill = "#e0e0e0",
                                          colour = "#e0e0e0"),
    plot.margin=unit(c(10,10,10,10),"mm"),
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.title=element_text(colour="black",size=rel(1)),
    legend.position="right",
    legend.text=element_text(size=rel(0.8)),
    axis.text=element_text(colour="black",size=rel(0.8),angle = 0),
    axis.text.x=element_text(angle=90,hjust = 0.95,vjust = .5),
    legend.key = element_rect(fill = NA))
  
  pdf(paste("EMC-Figure5D_",peakType,"_2022_4_20_motif.pdf",sep=""),    # create PNG for the heat map        
      width = 18,        # 5 x 300 pixels
      height = 7,
      pointsize = 10 )        # smaller font size
  grid.draw(p1)
  dev.off()
  pdf(paste("EMC-Figure5D_",peakType,"_2022_4_20_TFactivity.pdf",sep=""),    # create PNG for the heat map        
      width = 18,        # 5 x 300 pixels
      height = 2,
      pointsize = 10 )        # smaller font size
  grid.draw(p2)
  dev.off()
  pdf(paste("EMC-Figure5D_",peakType,"_2022_4_20_motifPerRegion.pdf",sep=""),    # create PNG for the heat map        
      width = 18,        # 5 x 300 pixels
      height = 15,
      pointsize = 10 )        # smaller font size
  grid.draw(p3)
  dev.off()
}
