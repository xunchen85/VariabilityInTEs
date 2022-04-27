###     Author:     Xun Chen
###     Date:       2021/4/5
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
Date_input = "2021_4_2"
Date = "2022_4_20"

sampleInfo = read.delim("inputs/sampleList_chipseq",header=F,sep="")
colnames(sampleInfo) = c("sampleID","totalStrength","totalPeakCount")
types = c("ATAC","H3K27ac","H3K27me3","H3K4me1","H3K4me3")
colorset = c("Flu" = "#99000d","NI" = "#4d4d4d")

############################# section 1 (create the csv file, only need to run it once)
TE_type_list = list.files(path="inputs/EMC_Flu_counts/")
TE_type_list = TE_type_list[grepl("AF|EU",TE_type_list)]
TE_type_list = TE_type_list[grepl("counts$",TE_type_list)]

type = "ATAC"

# preparing table 1
for (type in types){
  if (type == "ATAC"){
    TE_type_list_mark = TE_type_list[grepl("NI_peaks|Flu_peaks",TE_type_list)]
  } else if (type == "H3K27ac") {
    TE_type_list_mark = TE_type_list[grepl("H3K27ac",TE_type_list)]
  } else if (type == "H3K27me3") {
  TE_type_list_mark = TE_type_list[grepl("H3K27me3",TE_type_list)]
  } else if (type == "H3K4me1") {
  TE_type_list_mark = TE_type_list[grepl("H3K4me1",TE_type_list)]
  } else if (type == "H3K4me3") {
  TE_type_list_mark = TE_type_list[grepl("H3K4me3",TE_type_list)]
  }
  for (j in 1:length(TE_type_list_mark)){
    peaks = read.delim(TE_type_list_mark[j],sep="",header=T)
    #### Extract non-shuffle peaks
    peaks_nonshuffle = peaks[,!grepl("shuffle",colnames(peaks))]
    peaks_samples = gsub("TE","",colnames(peaks_nonshuffle))
    peaks_nonshuffle[,4:6] = NA
    colnames(peaks_nonshuffle)[4] = paste(colnames(peaks_nonshuffle)[3],"_shuffle.mean",sep="")
    colnames(peaks_nonshuffle)[5] = paste(colnames(peaks_nonshuffle)[3],"_FC",sep="")
    colnames(peaks_nonshuffle)[6] = paste(colnames(peaks_nonshuffle)[3],"_pvalue",sep="")
    peaks_shuffle_tmp1 = peaks[,grepl(paste(peaks_samples[3],"shuffle",sep=""),colnames(peaks))]
    
    ##### p value were further adjusted
    for (i in 1:nrow(peaks_nonshuffle)){
      peaks_nonshuffle[i,4] = mean(as.numeric(as.character(peaks_shuffle_tmp1[i,])),na.rm=T)
      peaks_nonshuffle[i,5] = as.numeric(as.character(peaks_nonshuffle[i,3]))/mean(as.numeric(as.character(peaks_shuffle_tmp1[i,])),na.rm=T)
      peaks_nonshuffle[i,6] = 2*mean(as.numeric(as.character(peaks_shuffle_tmp1[i,])) >= as.numeric(as.character(peaks_nonshuffle[i,3])))
    }
    if (j == 1){
      TE_type_list_mark_peaks = peaks_nonshuffle
    } else {
      TE_type_list_mark_peaks = cbind(TE_type_list_mark_peaks,peaks_nonshuffle[,3:6])
    }
  }
  write.csv(TE_type_list_mark_peaks,file = paste("EMC_flu-",type,"-",Date,"-1.csv",sep=""))
}

############################# section 2
type = "ATAC"
for (type in types){
  TE_type_list_mark_peaks = read.csv(paste("inputs/","EMC_flu-",type,"-",Date_input,"-1.csv",sep=""))
  TE_type_list_mark_peaks_beforeFiltering = TE_type_list_mark_peaks
  TE_type_list_mark_peaks = TE_type_list_mark_peaks[,-1]
  
  ### 1 Normalized the counts by the number of peaks
  # 1.1
  CountTypes = c("Maximum","Mean")
  CountType = "Mean"
  
  # 1.2 obtain mean peak counts 
  Samples = colnames(TE_type_list_mark_peaks)[grepl("summit.TE$",colnames(TE_type_list_mark_peaks))]
  if (type == "ATAC"){
    FluMean = mean(sampleInfo[!grepl("_H3K",sampleInfo$sampleID) & grepl("_Flu",sampleInfo$sampleID),3])
    NIMean = mean(sampleInfo[!grepl("_H3K",sampleInfo$sampleID) & grepl("_NI",sampleInfo$sampleID),3])
    BothMean = mean(sampleInfo[!grepl("_H3K",sampleInfo$sampleID),3])
  } else if (type == "H3K27ac") {
    FluMean = mean(sampleInfo[grepl("H3K27ac",sampleInfo$sampleID) & grepl("_Flu",sampleInfo$sampleID),3])
    NIMean = mean(sampleInfo[grepl("H3K27ac",sampleInfo$sampleID) & grepl("_NI",sampleInfo$sampleID),3])
    BothMean = mean(sampleInfo[grepl("H3K27ac",sampleInfo$sampleID),3])
  } else if (type == "H3K27me3") {
    FluMean = mean(sampleInfo[grepl("H3K27me3",sampleInfo$sampleID) & grepl("_Flu",sampleInfo$sampleID),3])
    NIMean = mean(sampleInfo[grepl("H3K27me3",sampleInfo$sampleID) & grepl("_NI",sampleInfo$sampleID),3])
    BothMean = mean(sampleInfo[grepl("H3K27me3",sampleInfo$sampleID),3])
  } else if (type == "H3K4me1") {
    FluMean = mean(sampleInfo[grepl("H3K4me1",sampleInfo$sampleID) & grepl("_Flu",sampleInfo$sampleID),3])
    NIMean = mean(sampleInfo[grepl("H3K4me1",sampleInfo$sampleID) & grepl("_NI",sampleInfo$sampleID),3])
    BothMean = mean(sampleInfo[grepl("H3K4me1",sampleInfo$sampleID),3])
  } else if (type == "H3K4me3") {
    FluMean = mean(sampleInfo[grepl("H3K4me3",sampleInfo$sampleID) & grepl("_Flu",sampleInfo$sampleID),3])
    NIMean = mean(sampleInfo[grepl("H3K4me3",sampleInfo$sampleID) & grepl("_NI",sampleInfo$sampleID),3])
    BothMean = mean(sampleInfo[grepl("H3K4me3",sampleInfo$sampleID),3])
  }
  
  # 1.3 normalization
  Start1 = (ncol(TE_type_list_mark_peaks)+1)
  TE_type_list_mark_peaks[,Start1:(Start1-1+length(Samples))] = NA
  colnames(TE_type_list_mark_peaks)[Start1:(Start1-1+length(Samples))] = gsub("summit.TE$","summit.TE.normalizedBoth",Samples)
  Start2 = (ncol(TE_type_list_mark_peaks)+1)
  TE_type_list_mark_peaks[,Start2:(Start2-1+length(Samples))] = NA
  colnames(TE_type_list_mark_peaks)[Start2:(Start2-1+length(Samples))] = gsub("summit.TE$","summit.TE.normalizedEach",Samples)
  i=283
  for (i in Start1:ncol(TE_type_list_mark_peaks)){
    if (grepl("normalizedBoth",colnames(TE_type_list_mark_peaks)[i])){
      colTmp = which(colnames(TE_type_list_mark_peaks)==gsub(".normalizedBoth","",colnames(TE_type_list_mark_peaks)[i]))
      sampleTmp = gsub(".narrowPeak_summit.TE|.broadPeak_summit.TE","",colnames(TE_type_list_mark_peaks)[colTmp])
      sampleInfo[sampleInfo$sampleID == sampleTmp,3]
      TE_type_list_mark_peaks[,i] = (TE_type_list_mark_peaks[,colTmp]/sampleInfo[sampleInfo$sampleID == sampleTmp,3])* BothMean 
    } else if (grepl("normalizedEach",colnames(TE_type_list_mark_peaks)[i]) & grepl("_Flu",colnames(TE_type_list_mark_peaks)[i])){
      colTmp = which(colnames(TE_type_list_mark_peaks)==gsub(".normalizedEach","",colnames(TE_type_list_mark_peaks)[i]))
      sampleTmp = gsub(".narrowPeak_summit.TE|.broadPeak_summit.TE","",colnames(TE_type_list_mark_peaks)[colTmp])
      sampleInfo[sampleInfo$sampleID == sampleTmp,3]
      TE_type_list_mark_peaks[,i] = (TE_type_list_mark_peaks[,colTmp]/sampleInfo[sampleInfo$sampleID == sampleTmp,3])* FluMean 
    } else if (grepl("normalizedEach",colnames(TE_type_list_mark_peaks)[i]) & grepl("_NI",colnames(TE_type_list_mark_peaks)[i])){
      colTmp = which(colnames(TE_type_list_mark_peaks)==gsub(".normalizedEach","",colnames(TE_type_list_mark_peaks)[i]))
      sampleTmp = gsub(".narrowPeak_summit.TE|.broadPeak_summit.TE","",colnames(TE_type_list_mark_peaks)[colTmp])
      sampleInfo[sampleInfo$sampleID == sampleTmp,3]
      TE_type_list_mark_peaks[,i] = (TE_type_list_mark_peaks[,colTmp]/sampleInfo[sampleInfo$sampleID == sampleTmp,3])* NIMean
    }
  }
  
  ### t test
  FluSamples = colnames(TE_type_list_mark_peaks)[grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks)) & grepl("_Flu",colnames(TE_type_list_mark_peaks))]
  Start = (ncol(TE_type_list_mark_peaks)+1)
  TE_type_list_mark_peaks[Start:(Start-1+length(FluSamples))] = NA
  
  colnames(TE_type_list_mark_peaks)[Start:(Start-1+length(FluSamples))] = gsub("_Flu","",FluSamples)
  colnames(TE_type_list_mark_peaks)[Start:(Start-1+length(FluSamples))] = gsub("_peaks.narrowPeak_summit.TE.normalizedEach|_peaks.broadPeak_summit.TE.normalizedEach","",colnames(TE_type_list_mark_peaks)[Start:(Start-1+length(FluSamples))])
  
  for (FluSample in FluSamples){
    NISample = gsub("_Flu","_NI",FluSample)
    Sample = gsub("_Flu","",FluSample)
    Sample = gsub("_peaks.narrowPeak_summit.TE.normalizedEach|_peaks.broadPeak_summit.TE.normalizedEach","",Sample)
    TE_type_list_mark_peaks[,Sample] = ifelse(TE_type_list_mark_peaks[,NISample]>0,TE_type_list_mark_peaks[,FluSample]/TE_type_list_mark_peaks[,NISample],NA)
    TE_type_list_mark_peaks[,Sample]
  }
  TE_type_list_mark_peaks$pvalue = NA
  TE_type_list_mark_peaks$padj = NA
  
  for (j in 1:nrow(TE_type_list_mark_peaks)){
    TE_type_list_mark_peaks[j,]$pvalue = t.test(as.matrix(TE_type_list_mark_peaks[j,grepl("_Flu_",colnames(TE_type_list_mark_peaks)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks))]),as.matrix(TE_type_list_mark_peaks[j,grepl("_NI_",colnames(TE_type_list_mark_peaks)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks))]),paired = TRUE)$p.value
    TE_type_list_mark_peaks[j,]$padj = p.adjust(TE_type_list_mark_peaks[j,]$pvalue,n=nrow(TE_type_list_mark_peaks))
  }
  write.csv(TE_type_list_mark_peaks,file = paste("EMC_flu-",type,"-",Date,"-2-withPadj.csv",sep=""))
  ### 2 remove TE families have less than 20 instances genome-wide
  TE_type_list_mark_peaks = TE_type_list_mark_peaks[TE_type_list_mark_peaks$totalCounts>=20,]
  ### 3.1 remove satellite group
  TE_type_list_mark_peaks = TE_type_list_mark_peaks[!grepl("Satellite",TE_type_list_mark_peaks$TE_family),]
  ### 3.2 filter by the padj (pvalue 0.01)
  TE_type_list_mark_peaks = TE_type_list_mark_peaks[TE_type_list_mark_peaks$padj<=0.01 & !is.na(TE_type_list_mark_peaks$padj),]               
  
  ### 4 filter by the Counts (based on the normalization of Flu and NI separately)
  if (CountType == "Maximum") {
    # at least one sample have >=20 peaks
    TE_type_list_mark_peaks = TE_type_list_mark_peaks[apply(TE_type_list_mark_peaks[,grepl("_summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks))],1,max) >=20,]
  } else if (CountType == "Mean") {
    # every sample have >=20 peaks on average
    TE_type_list_mark_peaks = TE_type_list_mark_peaks[
      (apply(TE_type_list_mark_peaks[,grepl("_Flu_",colnames(TE_type_list_mark_peaks)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks))],1,mean)>=20 |
         apply(TE_type_list_mark_peaks[,grepl("_NI_",colnames(TE_type_list_mark_peaks)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks))],1,mean)>=20),]
  }
  ### 5 kept only the fold change higher than 1.5
  # 5.1 keep candidate TEs have more peaks in Flu
  TE_type_list_mark_peaks_Flu = TE_type_list_mark_peaks[
    apply((TE_type_list_mark_peaks[,grepl("_Flu_",colnames(TE_type_list_mark_peaks)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks))]),1,mean)/apply(TE_type_list_mark_peaks[,grepl("_NI_",colnames(TE_type_list_mark_peaks))& grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks))],1,mean)>=1.5,]
  # 5.2 keep candidate TEs have more peaks in NI
  TE_type_list_mark_peaks_NI = TE_type_list_mark_peaks[
    apply((TE_type_list_mark_peaks[,grepl("_NI_",colnames(TE_type_list_mark_peaks))& grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks))]),1,mean)/apply(TE_type_list_mark_peaks[,grepl("_Flu_",colnames(TE_type_list_mark_peaks))& grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks))],1,mean)>=1.5,]
  
  ### 6 fold change significantly diff. (paired t test)
  enrichedStatuses = c("Flu","NI")
  enrichedStatus = "Flu"
  for (enrichedStatus in enrichedStatuses){
    if (enrichedStatus == "Flu") {
      TE_type_list_mark_peaks_tmp = data.frame(TE_type_list_mark_peaks_Flu)
    } else {
      TE_type_list_mark_peaks_tmp = data.frame(TE_type_list_mark_peaks_NI)
    }
    if (nrow(TE_type_list_mark_peaks_tmp) == 0) {
      next
    }
    
    # log2 fold enrichment
    dataplot_violine = data.frame("TEname"=NA,"totalCount"=NA,"fold" = NA,"sampleID"=NA)
    TE_type_list_mark_peaks_tmp_FC = TE_type_list_mark_peaks_tmp[,grepl("FC$",colnames(TE_type_list_mark_peaks_tmp))]
    for (i in 1:ncol(TE_type_list_mark_peaks_tmp_FC)){
      data_tmp = TE_type_list_mark_peaks_tmp_FC[,c(1,2)]
      data_tmp$TEname = TE_type_list_mark_peaks_tmp$TE_family
      data_tmp$totalCount = TE_type_list_mark_peaks_tmp$totalCounts
      data_tmp = data_tmp[,c(-1,-2)]
      data_tmp$fold = TE_type_list_mark_peaks_tmp_FC[,i]
      data_tmp$sampleID = gsub("_peaks.narrowPeak_summit.TE_FC","",colnames(TE_type_list_mark_peaks_tmp_FC)[i])
      dataplot_violine = rbind(dataplot_violine,data_tmp)
    }
    
    dataplot_violine = dataplot_violine[-1,]
    dataplot_violine = cSplit(dataplot_violine,"sampleID",sep="_",type.convert = as.character)
    dataplot_violine = cSplit(dataplot_violine,"TEname",sep=":",type.convert = as.character)
    dataplot_violine$sampleID_2 = as.factor(dataplot_violine$sampleID_2)
    dataplot_violine$sampleID_1 = as.factor(dataplot_violine$sampleID_1)
    dataplot_violine$log2fold = log2(dataplot_violine$fold)
    dataplot_violine$log2fold = ifelse(dataplot_violine$log2fold=="-Inf",NA,dataplot_violine$log2fold)
    
    Sum = Summarize(log2fold ~ TEname_1 + sampleID_2, data=dataplot_violine,na.rm=TRUE)
    Sum$se = Sum$sd / sqrt(Sum$n)
    Sum = Sum[order(Sum$mean),]
    if (enrichedStatus == "Flu"){
      TE_order.list2 = Sum[Sum$sampleID_2=="NI",]$TEname_1
    } else {
      TE_order.list2 = Sum[Sum$sampleID_2=="Flu",]$TEname_1
    }
    
    Sum$TEname_1 = factor(Sum$TEname_1,levels = as.character(TE_order.list2))
    if (type == "ATAC"){
      pd = position_dodge(0)
      p2 = ggplot(Sum, aes(x  = TEname_1, y = mean, color = sampleID_2)) +
        geom_point(shape = 19,size  = 0.5,position = pd) +
        geom_errorbar(aes(ymin  = Q1,ymax  = Q3),width = 0.5,size  = 0.5,position = pd) +
        scale_color_manual(values=colorset)+
        coord_flip()+
        geom_hline(yintercept = 0, linetype="dotted", color = "grey44", size=1)+
        ylab("log2(fold enrichment)")+
        xlab("TE family")+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.text=element_text(colour="black",size=rel(.8),angle = 0),
              axis.title=element_text(colour="black",size=rel(1.5)),
              #        legend.position=c(0.8,0.8),
              #        legend.position="bottom",
              #legend.position="right",
              legend.position="none",
              legend.background = element_blank(),
              legend.text=element_text(size=rel(1.2)))
      pdf(paste("EMC-Figure2DE-",type,"-",enrichedStatus,"-",Date,"-enrichment",".pdf",sep=""),    # create PNG for the heat map        
          width = 4,
          height = 5,
          pointsize = 10)        # smaller font size
      grid.draw(p2)
      dev.off()
    }
    
    # normalized peak counts
    TE_type_list_mark_peaks_tmp$TE_family2 = TE_type_list_mark_peaks_tmp$TE_family
    TE_type_list_mark_peaks_tmp = cSplit(TE_type_list_mark_peaks_tmp,"TE_family2",sep=":")
    TE_type_list_mark_peaks_tmp$TE_family2 = paste(TE_type_list_mark_peaks_tmp$TE_family2_1,TE_type_list_mark_peaks_tmp$TE_family2_2,sep=":")
    TE_type_list_mark_peaks_tmp = data.frame(TE_type_list_mark_peaks_tmp)
    rownames(TE_type_list_mark_peaks_tmp) = TE_type_list_mark_peaks_tmp$TE_family2
    TE_type_list_mark_peaks_tmp2 = TE_type_list_mark_peaks_tmp[,grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaks_tmp))]
    colnames(TE_type_list_mark_peaks_tmp2) = gsub("_peaks.narrowPeak_summit.TE.normalizedEach|_peaks.broadPeak_summit.TE.normalizedEach","",colnames(TE_type_list_mark_peaks_tmp2))
    
    i=1
    Individual_abundance = data.frame("TE"=NA,"totalCounts"=NA,"Normalized.counts"=NA,"sampleID"=NA,"fluStatus"=NA)
    for (i in 1:ncol(TE_type_list_mark_peaks_tmp2)){
      Individual_abundance_tmp = TE_type_list_mark_peaks_tmp2[,c(1,i)]
      colnames(Individual_abundance_tmp)[2] = "Normalized.counts"
      Individual_abundance_tmp$TE = rownames(TE_type_list_mark_peaks_tmp2)
      Individual_abundance_tmp$totalCounts = TE_type_list_mark_peaks_tmp[i,2]
      Individual_abundance_tmp$sampleID = colnames(TE_type_list_mark_peaks_tmp2)[i]
      Individual_abundance_tmp$fluStatus = ifelse(grepl("Flu",colnames(TE_type_list_mark_peaks_tmp2)[i]),"Flu","NI")
      Individual_abundance_tmp = Individual_abundance_tmp[,c(3,4,2,5,6)]
      Individual_abundance = rbind(Individual_abundance,Individual_abundance_tmp)
    }
    Individual_abundance = Individual_abundance[-1,]

    ### plots
    # log abundance plot
    colorset = c("Flu" = "#99000d","NI" = "#4d4d4d")
    
    Individual_abundance$log10_count = log10(Individual_abundance$Normalized.counts)
    if (nrow(Individual_abundance[is.infinite(Individual_abundance$log10_count),])>0){
      Individual_abundance[is.infinite(Individual_abundance$log10_count),]$log10_count <-NA
    }
    plotData = cSplit(Individual_abundance,"TE",sep=":",type.convert = as.character)
    iris.summary <- aggregate(plotData$log10_count, by =list(plotData$TE_1,plotData$fluStatus), mean,na.rm=TRUE)
    colnames(iris.summary) = c("TE","fluStatus","mean")
    plotData$TE_1 = factor(plotData$TE_1,levels=TE_order.list2)
    if (type == "ATAC"){
      p3<-ggplot(plotData, aes(x=TE_1, y=log10_count, colour = fluStatus, fill=fluStatus)) + #volcanoplot with log2Foldchange versus pvalue
        geom_jitter(aes(colour = fluStatus),shape=19, size=1,  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))+
        geom_point(data=iris.summary, aes(x = TE, y = mean,color=fluStatus),size=4,shape=3,position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8))+
        scale_fill_manual(values=colorset)+
        scale_color_manual(values=colorset)+
        ylab("log10(# normalized accessible instances)") + 
        xlab("TE family") + 
        coord_flip()+
        theme(#panel.grid.major = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          
          axis.text=element_text(colour="black",size=rel(.8),angle = 0),
          axis.title=element_text(colour="black",size=rel(1.5)),
          legend.position="none",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1.2)))
      pdf(paste("EMC-Figure2DF-",type,"-",enrichedStatus,"-",Date,"-counts",".pdf",sep=""),    # create PNG for the heat map        
          width = 4,
          height = 5,
          pointsize = 10)        # smaller font size
      grid.draw(p3)
      dev.off()
    }

    ### obtain the TE order
    if (nrow(TE_type_list_mark_peaks_tmp)<2) {
      next
    }
    input3 = TE_type_list_mark_peaks_tmp[,grepl("FC$",colnames(TE_type_list_mark_peaks_tmp))]
    peaks_TEname = TE_type_list_mark_peaks_tmp[,c(1,2)]
    peaks_TEname = cSplit(peaks_TEname,"TE_family",sep=":",type.convert = as.character)
    peaks_TEname$TE_family = peaks_TEname$TE_family_1
    rownames(input3) = peaks_TEname$TE_family
    input3 = as.matrix(log2(input3))
    colnames(input3) = gsub("_peaks.narrowPeak_summit.TE_FC","",colnames(input3))
    colnames(input3) = gsub("_peaks.broadPeak_summit.TE_FC","",colnames(input3))
    my_palette <- colorRampPalette(c("#377eb8","white", "#e41a1c"))(n = 100)
    if ("LTR10E" %in% rownames(input3)){
      input3["LTR10E","EU37_Flu"] = mean(input3["LTR10E",grepl("_Flu",colnames(input3)) & !grepl("EU37_Flu",colnames(input3))])
    }
    TE_order = heatmap.2(input3,
                         main = "Correlation", # heat map title
                         #  cellnote = mat_data,  # same data set for cell labels
                         #  notecol="black",      # change font color of cell labels to black
                         density.info="none",  # turns off density plot inside color legend
                         trace="none",         # turns off trace lines inside the heat map
                         margins =c(12,9),     # widens margins around plot
                         #          labRow=NULL,
                         #          labCol=NULL,
                         #          xlab = NULL,
                         #          ylab = NULL,
                         key.xlab=NULL,
                         col=my_palette       # use on color palette defined earlier 
                         #          breaks=col_breaks,    # enable color transition at specified limits
                         #          RowSideColors=c_color,
                         #  col=redblue(16),
                         #  dendrogram="col",     # only draw a row dendrogram
                         #  Rowv="NA",
                         #  Colv="NA"            # turn off column clustering
    )     
    TE_order.list = rownames(input3)[TE_order$rowInd]
    write.table(TE_order.list,file = paste("EMC_flu-",type,"-",enrichedStatus,"-",Date,"-TE.order",sep=""))
  }
}

############################# section 3: combine marks
enrichedStatuses = c("Flu","NI")
types = c("ATAC","H3K27ac","H3K27me3","H3K4me1","H3K4me3")

type2 ="ATAC"
enrichedStatus = "Flu"
Sum_peakPadjTable = data.frame("mark"=NA,"TE"=NA,"normalizedCount_Flu_Mean","normalizedCount_NI_Mean","normalizedCount_Mean","meanFC"=NA,"log2FC"=NA,"padj"=NA,"padj_label"=NA)
for (enrichedStatus in enrichedStatuses){
  TE_order_peaksPadj = read.delim(paste("EMC_flu-",type2,"-",enrichedStatus,"-",Date,"-TE.order",sep=""),header=T,sep="")
  peakPadj = data.frame("mark"=NA,"TE" = NA,"meanFC"=NA,"log2FC"=NA,"padj" = NA,"padj_label"=NA)
  for (type in types){
    TE_type_list_mark_peaksPadj = read.csv(paste("EMC_flu-",type,"-",Date,"-2-withPadj.csv",sep=""))
    TE_type_list_mark_peaksPadj = data.frame(cSplit(TE_type_list_mark_peaksPadj,"TE_family",sep=":",type.convert = as.character))
    TE_type_list_mark_peaksPadj$mark = type
    TE_type_list_mark_peaksPadj$TE = TE_type_list_mark_peaksPadj$TE_family_1
    TE_type_list_mark_peaksPadj$meanFC = apply(TE_type_list_mark_peaksPadj[,!grepl("peaks",colnames(TE_type_list_mark_peaksPadj)) & grepl("AF|EU",colnames(TE_type_list_mark_peaksPadj))],1,mean,na.rm=T)
    TE_type_list_mark_peaksPadj$log2FC = log2(TE_type_list_mark_peaksPadj$meanFC)
    TE_type_list_mark_peaksPadj$normalizedCount_Flu_Mean = apply((TE_type_list_mark_peaksPadj[,grepl("_Flu_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaksPadj))]),1,mean)
    TE_type_list_mark_peaksPadj$normalizedCount_NI_Mean = apply((TE_type_list_mark_peaksPadj[,grepl("_NI_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaksPadj))]),1,mean)
    TE_type_list_mark_peaksPadj$normalizedCount_Mean = apply((TE_type_list_mark_peaksPadj[,grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaksPadj))]),1,mean)
    TE_type_list_mark_peaksPadj$padj_label = ifelse(TE_type_list_mark_peaksPadj$padj<=0.001,"***",NA)
    TE_type_list_mark_peaksPadj$padj_label = ifelse(TE_type_list_mark_peaksPadj$padj>0.001 & TE_type_list_mark_peaksPadj$padj<=0.01,"**",TE_type_list_mark_peaksPadj$padj_label)
    TE_type_list_mark_peaksPadj$padj_label = ifelse(TE_type_list_mark_peaksPadj$padj>0.01 & TE_type_list_mark_peaksPadj$padj<=0.05,"*",TE_type_list_mark_peaksPadj$padj_label) ## update on 2022/1/12
    if (enrichedStatus == "Flu"){
      Sum_peakPadj_tmp = TE_type_list_mark_peaksPadj[,c("mark","TE","normalizedCount_Flu_Mean","normalizedCount_NI_Mean","normalizedCount_Mean","meanFC","log2FC","padj","padj_label")]
      colnames(Sum_peakPadj_tmp)[c(1,3:9)] = paste(colnames(Sum_peakPadj_tmp[c(1,3:9)]),type,sep="_")
      if (nrow(Sum_peakPadjTable) == 1){
        Sum_peakPadjTable = Sum_peakPadj_tmp
      } else {
        Sum_peakPadjTable = merge(Sum_peakPadjTable,Sum_peakPadj_tmp,by="TE",all=T)
      }
    }
    TE_type_list_mark_peaksPadj = data.frame(TE_type_list_mark_peaksPadj[TE_type_list_mark_peaksPadj$TE_family_1 %in% TE_order_peaksPadj$x,])
    peakPadj_tmp = TE_type_list_mark_peaksPadj[,c("mark","TE","meanFC","log2FC","padj","padj_label")]
    peakPadj = rbind(peakPadj,peakPadj_tmp)
  }
  peakPadj = peakPadj[-1,]
  peakPadj$mark = factor(peakPadj$mark,levels=c("ATAC","H3K27ac","H3K4me1","H3K4me3","H3K27me3"))
  peakPadj$TE = factor(peakPadj$TE,levels=TE_order_peaksPadj$x) 
  ### another order
  tmp = peakPadj[peakPadj$mark=="H3K27ac",]
  tmp = tmp[order(tmp$padj),]
  peakPadj$TE = factor(peakPadj$TE,levels=tmp$TE)
  
  p<-ggplot(peakPadj, aes(mark,TE)) +
    geom_tile(aes(fill = log2FC),col = "white",size=0.5) + 
    geom_text(aes(label = padj_label),col = "grey44") +
    scale_fill_gradient2(low="#2166ac",mid = "white", high = "#b2182b",midpoint = 0)+
    xlab("")+
    ylab("Subfamily")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.text=element_text(colour="black",size=rel(.8),angle = 0),
          axis.text.x=element_text(colour="black",angle = 90,hjust = 0.95,vjust = .5),
          axis.title=element_text(colour="black",size=rel(1.5)),
          legend.position="right",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1.2)))
  
  pdf(paste("EMC-Figure2EG",type,"-",enrichedStatus,"-",Date,"-combined.pdf",sep=""),    # create PNG for the heat map        
      width = 3.5,        # 5 x 300 pixels
      height = 6,
      pointsize = 10 )        # smaller font size
  grid.draw(p)
  dev.off()
}
write.csv(Sum_peakPadjTable,file = paste("EMC_flu-",type,"-",enrichedStatus,"-",Date,"-AllMark-withPadj.csv",sep=""))

############################# section 4: organize the output as the supplementary
enrichedStatuses = c("Flu","NI")
types = c("ATAC","H3K27ac","H3K27me3","H3K4me1","H3K4me3")
enrichedStatus = "Flu"
Sum_peakPadjTable = data.frame("mark"=NA,"TE_family_3"=NA,"TE"=NA,"totalCounts"=NA,"actualCount_Flu_mean"=NA,"actualCount_Flu_sd"=NA,"normalizedCount_Flu_Mean"=NA,"normalizedCount_Flu_sd"=NA,"enrichmentLevel_Flu_Mean"=NA,"enrichmentLevel_Flu_sd"=NA,"actualCount_NI_mean"=NA,"actualCount_NI_sd"=NA,"normalizedCount_NI_Mean"=NA,"normalizedCount_NI_sd"=NA,"enrichmentLevel_NI_Mean"=NA,"enrichmentLevel_NI_sd"=NA,"meanFC"=NA,"log2FC"=NA,"padj"=NA,"padj_label"=NA)
for (enrichedStatus in enrichedStatuses){
  TE_order_peaksPadj = read.delim(paste("EMC_flu-ATAC-",enrichedStatus,"-",Date,"-TE.order",sep=""),header=T,sep="")
  peakPadj = data.frame("mark"=NA,"TE" = NA,"meanFC"=NA,"log2FC"=NA,"padj" = NA,"padj_label"=NA)
  for (type in types){
    TE_type_list_mark_peaksPadj = read.csv(paste("EMC_flu-",type,"-",Date,"-2-withPadj.csv",sep=""))
    TE_type_list_mark_peaksPadj = data.frame(cSplit(TE_type_list_mark_peaksPadj,"TE_family",sep=":",type.convert = as.character))
    TE_type_list_mark_peaksPadj$mark = type
    TE_type_list_mark_peaksPadj$TE = TE_type_list_mark_peaksPadj$TE_family_1
    TE_type_list_mark_peaksPadj$actualCount_Flu_mean = apply((TE_type_list_mark_peaksPadj[,grepl("_Flu_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE$",colnames(TE_type_list_mark_peaksPadj))]),1,mean)
    TE_type_list_mark_peaksPadj$actualCount_Flu_sd = apply((TE_type_list_mark_peaksPadj[,grepl("_Flu_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE$",colnames(TE_type_list_mark_peaksPadj))]),1,sd)
    TE_type_list_mark_peaksPadj$normalizedCount_Flu_Mean = apply((TE_type_list_mark_peaksPadj[,grepl("_Flu_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaksPadj))]),1,mean)
    TE_type_list_mark_peaksPadj$normalizedCount_Flu_sd = apply((TE_type_list_mark_peaksPadj[,grepl("_Flu_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaksPadj))]),1,sd)
    TE_type_list_mark_peaksPadj$enrichmentLevel_Flu_Mean = apply((TE_type_list_mark_peaksPadj[,grepl("_Flu_",colnames(TE_type_list_mark_peaksPadj)) & grepl("_summit.TE_FC",colnames(TE_type_list_mark_peaksPadj))]),1,mean)
    TE_type_list_mark_peaksPadj$enrichmentLevel_Flu_sd = apply((TE_type_list_mark_peaksPadj[,grepl("_Flu_",colnames(TE_type_list_mark_peaksPadj)) & grepl("_summit.TE_FC",colnames(TE_type_list_mark_peaksPadj))]),1,sd)
    TE_type_list_mark_peaksPadj$actualCount_NI_mean = apply((TE_type_list_mark_peaksPadj[,grepl("_NI_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE$",colnames(TE_type_list_mark_peaksPadj))]),1,mean)
    TE_type_list_mark_peaksPadj$actualCount_NI_sd = apply((TE_type_list_mark_peaksPadj[,grepl("_NI_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE$",colnames(TE_type_list_mark_peaksPadj))]),1,sd)
    TE_type_list_mark_peaksPadj$normalizedCount_NI_Mean = apply((TE_type_list_mark_peaksPadj[,grepl("_NI_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaksPadj))]),1,mean)
    TE_type_list_mark_peaksPadj$normalizedCount_NI_sd = apply((TE_type_list_mark_peaksPadj[,grepl("_NI_",colnames(TE_type_list_mark_peaksPadj)) & grepl("summit.TE.normalizedEach",colnames(TE_type_list_mark_peaksPadj))]),1,sd)
    TE_type_list_mark_peaksPadj$enrichmentLevel_NI_Mean = apply((TE_type_list_mark_peaksPadj[,grepl("_NI_",colnames(TE_type_list_mark_peaksPadj)) & grepl("_summit.TE_FC",colnames(TE_type_list_mark_peaksPadj))]),1,mean)
    TE_type_list_mark_peaksPadj$enrichmentLevel_NI_sd = apply((TE_type_list_mark_peaksPadj[,grepl("_NI_",colnames(TE_type_list_mark_peaksPadj)) & grepl("_summit.TE_FC",colnames(TE_type_list_mark_peaksPadj))]),1,sd)
    TE_type_list_mark_peaksPadj$meanFC = apply(TE_type_list_mark_peaksPadj[,!grepl("peaks",colnames(TE_type_list_mark_peaksPadj)) & grepl("AF|EU",colnames(TE_type_list_mark_peaksPadj))],1,mean,na.rm=T)
    TE_type_list_mark_peaksPadj$log2FC = log2(TE_type_list_mark_peaksPadj$meanFC)
    TE_type_list_mark_peaksPadj$padj_label = ifelse(TE_type_list_mark_peaksPadj$padj<=0.001,"***",NA)
    TE_type_list_mark_peaksPadj$padj_label = ifelse(TE_type_list_mark_peaksPadj$padj>0.001 & TE_type_list_mark_peaksPadj$padj<=0.01,"**",TE_type_list_mark_peaksPadj$padj_label)
    TE_type_list_mark_peaksPadj$padj_label = ifelse(TE_type_list_mark_peaksPadj$padj>0.01 & TE_type_list_mark_peaksPadj$padj<=0.05,"*",TE_type_list_mark_peaksPadj$padj_label)      ## update on 2022/1/12
    if (enrichedStatus == "Flu"){
      Sum_peakPadj_tmp = TE_type_list_mark_peaksPadj[,c("mark","TE_family_3","TE","totalCounts","actualCount_Flu_mean","actualCount_Flu_sd","normalizedCount_Flu_Mean","normalizedCount_Flu_sd","enrichmentLevel_Flu_Mean","enrichmentLevel_Flu_sd","actualCount_NI_mean","actualCount_NI_sd","normalizedCount_NI_Mean","normalizedCount_NI_sd","enrichmentLevel_NI_Mean","enrichmentLevel_NI_sd","meanFC","log2FC","padj","padj_label")]
      #colnames(Sum_peakPadj_tmp)[c(1,5:20)] = paste(colnames(Sum_peakPadj_tmp[c(1,5:20)]),type,sep="_")
      if (nrow(Sum_peakPadjTable) == 1){
        Sum_peakPadjTable = Sum_peakPadj_tmp
      } else {
        #Sum_peakPadjTable = merge(Sum_peakPadjTable,Sum_peakPadj_tmp,by="TE",all=T)
        Sum_peakPadjTable = rbind(Sum_peakPadjTable,Sum_peakPadj_tmp)
      }
    }
    head(Sum_peakPadj_tmp)
    TE_type_list_mark_peaksPadj = data.frame(TE_type_list_mark_peaksPadj[TE_type_list_mark_peaksPadj$TE_family_1 %in% TE_order_peaksPadj$x,])
    peakPadj_tmp = TE_type_list_mark_peaksPadj[,c("mark","TE","meanFC","log2FC","padj","padj_label")]
    peakPadj = rbind(peakPadj,peakPadj_tmp)
  }
  peakPadj = peakPadj[-1,]
  peakPadj$mark = factor(peakPadj$mark,levels=c("ATAC","H3K27ac","H3K4me1","H3K4me3","H3K27me3"))
  peakPadj$TE = factor(peakPadj$TE,levels=TE_order_peaksPadj$x) 
  ### another order
  tmp = peakPadj[peakPadj$mark=="H3K27ac",]
  tmp = tmp[order(tmp$padj),]
  peakPadj$TE = factor(peakPadj$TE,levels=tmp$TE)
  
}
write.csv(Sum_peakPadjTable,file = paste("EMC_supplementary_table5_2022_3_17.csv"))

