###     Author:     Xun Chen
###     Date:       2021/4/8
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
library("ggpmisc")
library(gtools)
library(FSA)

####
setwd(dir="/Users/xchen/OneDrive - Kyoto University/Projects_Bourque/EMC_project_2020_12_1/Final_Rscripts_2022_4_20/")

### TE info
TElist = read.delim("inputs/TElist.sep.27",header=T,sep="")

### individual info
IndividualInfo = read.csv("inputs/H1N1_Viral_Load_5_27.csv")
IndividualInfo = IndividualInfo[!is.na(IndividualInfo$Group_individual),]
IndividualInfo$Group_individual_color = factor(IndividualInfo$Group_individual)
IndividualInfo$Sample_ID2 = gsub("_Flu|_NI","",IndividualInfo$Sample_ID)
levels(IndividualInfo$Group_individual_color)<-c("#377eb8","#4daf4a","#737373","#ff7f00")

############################# section 1 (create the csv file, only need to run it once)
types = c("ATAC","H3K27ac","H3K27me3","H3K4me1","H3K4me3")
colorset = c("Flu" = "#99000d","NI" = "#4d4d4d")
my_palette <- colorRampPalette(c("#377eb8","white", "#e41a1c"))(n = 100)

############################# section 2
type = "ATAC"
TE_type_list_mark_peaks = read.csv(paste("inputs/","EMC_flu-",type,"-","2021_4_2","-1.csv",sep=""))
TE_type_list_mark_peaks_beforeFiltering = TE_type_list_mark_peaks
TE_type_list_mark_peaks = TE_type_list_mark_peaks[,-1]
TE_type_list_mark_peaks$TEinfo = TE_type_list_mark_peaks$TE_family
TE_type_list_mark_peaks = data.frame(cSplit(TE_type_list_mark_peaks,"TEinfo",sep=":"))
TE_type_list_mark_peaks_Flu = TE_type_list_mark_peaks[TE_type_list_mark_peaks$TEinfo_1 %in% as.character(TElist[TElist$TE_cluster != "NI",]$TEfamily),]
TE_type_list_mark_peaks_NI = TE_type_list_mark_peaks[TE_type_list_mark_peaks$TEinfo_1 %in% as.character(TElist[TElist$TE_cluster == "NI",]$TEfamily),]

#### 6 fold change significantly diff. (paired t test)
enrichedStatuses = c("Fluenriched","NIenriched")
enrichedStatus = "Fluenriched"
for (enrichedStatus in enrichedStatuses){
  if (enrichedStatus == "Fluenriched") {
    TE_type_list_mark_peaks_tmp = data.frame(TE_type_list_mark_peaks_Flu)
    mycolors = c("#b15928","#1f78b4")
  } else {
    TE_type_list_mark_peaks_tmp = data.frame(TE_type_list_mark_peaks_NI)
    mycolors = c("#33a02c")
  }
  
  ### data for the heatmap, data: fold enrichment
  if (nrow(TE_type_list_mark_peaks_tmp)<2) {
    next
  }
  
  TE_type_list_mark_peaks_tmp = merge(TE_type_list_mark_peaks_tmp,TElist,by.x="TEinfo_1",by.y="TEfamily",all.x=T)
  TE_type_list_mark_peaks_tmp = TE_type_list_mark_peaks_tmp[order(TE_type_list_mark_peaks_tmp$Reordered_by_motif),]
  rownames(TE_type_list_mark_peaks_tmp) = TE_type_list_mark_peaks_tmp$TEinfo_1
  
  # extract FC
  input3 = TE_type_list_mark_peaks_tmp[,grepl("FC$",colnames(TE_type_list_mark_peaks_tmp))]
  colnames(input3) = gsub("_peaks.narrowPeak_summit.TE_FC","",colnames(input3))
  colnames(input3) = gsub("_peaks.broadPeak_summit.TE_FC","",colnames(input3))
  # column color
  c_color = as.character(IndividualInfo[match(colnames(input3),IndividualInfo$Sample_ID),]$Group_individual_color)
  
  # heatmap color
  my_palette <- colorRampPalette(c("#377eb8","white", "#e41a1c"))(n = 100)
  # convert to matrix format with log2 conversion
  input3 = as.matrix(log2(input3))
  # replace the -inf value using the mean
  if ("LTR10E" %in% rownames(input3)){
    input3["LTR10E","EU37_Flu"] = mean(input3["LTR10E",grepl("_Flu",colnames(input3)) & !grepl("EU37_Flu",colnames(input3))])
  }
  
  ###### clustering of Flu samples
  input3_Flu = input3[,grepl("_Flu",colnames(input3))]
  c_color = as.character(IndividualInfo[match(colnames(input3_Flu),IndividualInfo$Sample_ID),]$Group_individual_color)
  
  pdf(paste("EMC-Figure3C-",enrichedStatus,"-Flu",".pdf",sep=""),    # create PNG for the heat map        
      width = 6,        # 5 x 300 pixels
      height = 8,
      pointsize = 10
  )        # smaller font size
  heatmap.2(input3_Flu,
            main = "", # heat map title
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(12,9),     # widens margins around plot
            key.xlab=NULL,
            hclustfun = function(x) hclust(x,method = 'complete'),
            distfun = function(x) dist(x,method = 'euclidean'),
            ColSideColors=c_color,
            col=my_palette,       # use on color palette defined earlier 
            Rowv="NA",
  )     
  dev.off()
  
  ###### clustering of NI samples
  input3_NI = input3[,grepl("_NI",colnames(input3))]
  c_color = as.character(IndividualInfo[match(colnames(input3_NI),gsub("_Flu","_NI",as.character(IndividualInfo$Sample_ID))),]$Group_individual_color)
  pdf(paste("EMC_Figure3C-",enrichedStatus,"-NI",".pdf",sep=""),    # create PNG for the heat map        
      width = 6,        # 5 x 300 pixels
      height = 8,
      pointsize = 10
  )        # smaller font size
  heatmap.2(input3_NI,
            main = "", # heat map title
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(12,9),     # widens margins around plot
            key.xlab=NULL,
            ColSideColors=c_color,
            col=my_palette,       # use on color palette defined earlier 
            hclustfun = function(x) hclust(x,method = 'complete'),
            distfun = function(x) dist(x,method = 'euclidean'),
            Rowv="NA",
  )     
  dev.off()
}
  

