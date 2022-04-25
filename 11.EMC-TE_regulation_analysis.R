###     Author:     Xun Chen
###     Date:       2021/12/12
###     Contact:    xunchen85@gmail.com

library(ggplot2)
library(grid)
library(dplyr)
library(ggpmisc)
library(splitstackshape)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

######
data_summary_95quantile <- function(x) {
  m <- mean(x)
  ymin <- unname(quantile(x,probs=0.025))
  ymax <- unname(quantile(x,probs=0.975))
  return(c(y=m,ymin=ymin,ymax=ymax))
}
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

####### Set up enviroment;

setwd(dir="/Users/xchen/OneDrive - Kyoto University/Projects_Bourque/EMC_project_2020_12_1/Final_Rscripts_2022_4_20/")
TPM = read.csv("inputs/EMC_RNAseq_hg19.rename.CountM-TPM.csv")
#TPM = read.csv("EMC_RNAseq_hg19.rename.CountM-TPM-gene.csv")
Dif = read.delim("inputs/EMC_RNAseq_hg19.rename.CountM_pairedEnd_test_gene_TE_analysis.txt",header=T,sep="")
Gene_coordinates = read.delim("inputs/hg19.ensGene.Genename",header=F,sep="")
TElist = read.delim("inputs/TElist.sep.27",header=T,sep="")
# instances
TEinstance_summary = read.delim("inputs/Summary.expAndCentroid.enriched.table",sep="",header=F)
TEinstance_summary_header = read.delim("inputs/Summary.expAndCentroid.table.header",sep="",header=T)

#############################
############################# Section 1
####### TPM
TPM$Flu_mean = apply(TPM[,grepl("_Flu",colnames(TPM))],1,mean)
TPM$NI_mean = apply(TPM[,grepl("_NI",colnames(TPM))],1,mean)
TPM$Mean_max = ifelse(TPM$Flu_mean>=TPM$NI_mean,TPM$Flu_mean,TPM$NI_mean)
TPM = TPM[,c("name","len","Flu_mean","NI_mean","Mean_max")]
TPM_tpm1_500bp = TPM[TPM$len>=500 & (TPM$Flu_mean>=1 | TPM$NI_mean>=1),]
TPM_tpm2_500bp = TPM[TPM$len>=500 & (TPM$Flu_mean>=2 | TPM$NI_mean>=2),]
TPM_tpm0_500bp = TPM[TPM$len>=500 & (TPM$Flu_mean>0 | TPM$NI_mean>0),]

####### Gene differential expression
Dif$GeneName = rownames(Dif)
Dif = Dif[!grepl(":",Dif$GeneName),]

####### Coordinates
head(Gene_coordinates)
colnames(Gene_coordinates) = c("transcriptName","geneName","chr","strand","start","end","refGeneName")
Gene_coordinates$TSS = ifelse(Gene_coordinates$strand=="+",Gene_coordinates$start,Gene_coordinates$end)
Gene_coordinates_1 = Gene_coordinates[Gene_coordinates$strand == "+",]
Gene_coordinates_1 = Gene_coordinates_1[order(Gene_coordinates_1$start),]
Gene_coordinates_1 = Gene_coordinates_1[!duplicated(Gene_coordinates_1$geneName),]
Gene_coordinates_2 = Gene_coordinates[Gene_coordinates$strand == "-",]
Gene_coordinates_2 = Gene_coordinates_2[order(-Gene_coordinates_2$start),]
Gene_coordinates_2 = Gene_coordinates_2[!duplicated(Gene_coordinates_2$geneName),]
Gene_coordinates_unique = rbind(Gene_coordinates_1,Gene_coordinates_2)

####### Merge
Genes = merge(TPM_tpm1_500bp,Dif,all.x=T,by.x="name",by.y="GeneName")                #### filter by tmp1
Genes = merge(Gene_coordinates_unique,Genes,all.y=T,by.y="name",by.x="geneName")
####### format and write
Genes_bed = Genes
Genes_bed$anno = paste(Genes_bed$geneName,Genes_bed$chr,Genes_bed$strand,Genes_bed$start,Genes_bed$end,Genes_bed$refGeneName,Genes_bed$TSS,Genes_bed$len,Genes_bed$Flu_mean,Genes_bed$NI_mean,Genes_bed$Mean_max,Genes_bed$baseMean,Genes_bed$log2FoldChange,Genes_bed$lfcSE,Genes_bed$stat,Genes_bed$pvalue,Genes_bed$padj,sep=",")
Genes_bed$log2FC_padj = paste(Genes_bed$log2FoldChange,Genes_bed$pvalue,sep=",")
Genes_bed = Genes_bed[,c("chr","TSS","TSS","anno","log2FC_padj","strand")]
Genes_bed$TSS = ifelse(Genes_bed$strand=="+",Genes_bed$TSS-1,Genes_bed$TSS)
Genes_bed$TSS.1 = ifelse(Genes_bed$strand=="-",Genes_bed$TSS.1+1,Genes_bed$TSS.1)
### output the gene coordinates in BED
write.table(Genes_bed,file="EMC_highExp_genes_hg19_tss.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names = FALSE)

### All accessible instances
### TE instances
colnames(TEinstance_summary) = colnames(TEinstance_summary_header)
TEinstance_summary$is_ATAC = ifelse(TEinstance_summary$Counts.centroidATAC>0,"ATAC","No")
TEinstance_summary$is_K27ac = ifelse(TEinstance_summary$Counts.centroidH3K27ac>0,"K27ac","No")
TEinstance_summary$is_K4me1 = ifelse(TEinstance_summary$Counts.centroidH3K4me1>0,"K4me1","No")
TEinstance_summary$is_K4me3 = ifelse(TEinstance_summary$Counts.centroidH3K4me3>0,"K4me3","No")
TEinstance_summary$is_K27me3 = ifelse(TEinstance_summary$Counts.centroidH3K27me3>0,"K27me3","No")
TEinstance_summary = TEinstance_summary[,c("TE_instance","is_ATAC","is_K27ac","is_K4me1","is_K4me3","is_K27me3")]
TEinstance_summary = TEinstance_summary[TEinstance_summary$is_ATAC=="ATAC",]

### output instances with each mark
write.table(TEinstance_summary[TEinstance_summary$is_K27ac=="K27ac",]$TE_instance,file="TE_regulation.summary_2021_4_22_TEinstance_K27ac",col.names=FALSE,row.names=FALSE,sep="\t", quote = FALSE)
write.table(TEinstance_summary[TEinstance_summary$is_K4me1=="K4me1",]$TE_instance,file="TE_regulation.summary_2021_4_22_TEinstance_K4me1",col.names=FALSE,row.names=FALSE,sep="\t", quote = FALSE)
write.table(TEinstance_summary[TEinstance_summary$is_K4me3=="K4me3",]$TE_instance,file="TE_regulation.summary_2021_4_22_TEinstance_K4me3",col.names=FALSE,row.names=FALSE,sep="\t", quote = FALSE)
write.table(TEinstance_summary[TEinstance_summary$is_K27me3!="K27me3",]$TE_instance,file="TE_regulation.summary_2021_4_22_TEinstance_noK27me3",col.names=FALSE,row.names=FALSE,sep="\t", quote = FALSE)
write.table(TEinstance_summary[TEinstance_summary$is_K4me1=="K4me1" & TEinstance_summary$is_K27ac=="K27ac",]$TE_instance,file="TE_regulation.summary_2021_4_22_TEinstance_K4me1_K27ac",col.names=FALSE,row.names=FALSE,sep="\t", quote = FALSE)


################ section 2
### Read summary info
NearbyGenes = read.delim("inputs/TE_regulation.summary_2021_12_12_broad",sep="",header=F)

### read Window info
WindowSize = read.delim("inputs/TE_regulation.summary_2021_12_12_broad.windowSize",header=F,sep="")

### format
NearbyGenes$Group = ifelse(grepl("shuffle",NearbyGenes$V1),"Shuffle","NonShuffling")
NearbyGenes$Approach = ifelse(grepl("closest",NearbyGenes$V1),"Closest","Window")
colnames(NearbyGenes) = c("dataset","TE_cluster","Type","Window","up_sig_count","down_sig_count","not_sig_count","upregulated_perC","downregulated_perC","Group","Approach")

# filter out pairs longer than 1Mb
NearbyGenes = NearbyGenes[NearbyGenes$Window!= 0,]
NearbyGenes[NearbyGenes$Group!="Shuffle",]

# Only study specific peaks
NearbyGenes = NearbyGenes[NearbyGenes$Type!="Shared" & NearbyGenes$Type!="All",]
NearbyGenes = NearbyGenes[(NearbyGenes$TE_cluster == "NI" & NearbyGenes$Type != "Flu_specific") | (NearbyGenes$TE_cluster != "NI" & NearbyGenes$Type != "NI_specific"),]

# Analyze by using Window approach
NearbyGenes = NearbyGenes[NearbyGenes$Approach == "Window",]
# Windows
NearbyGenes$Window = factor(NearbyGenes$Window,levels=c(as.character(WindowSize$V1)))
# Corrected the Type of peaks
NearbyGenes$Type_modified = factor(ifelse(as.character(NearbyGenes$Type) == "All","All","Specific"))
NearbyGenes$total_count = NearbyGenes$up_sig_count+NearbyGenes$down_sig_count+NearbyGenes$not_sig_count

###################################################### plot
###### ATAC only
NearbyGenes_plot = NearbyGenes
NearbyGenes_plot = NearbyGenes_plot[!grepl("sort_K|sort_no|sort_MEF",NearbyGenes_plot$dataset),]
NearbyGenes_plot$Type_modified

### plots
lineType = c("All"="dotted","Specific"="solid")
colorSet = c("a" = "#7C7C7C", "b" = "#bf812d","NI" = "#beaed4")

p1 = ggplot(NearbyGenes_plot[grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type_modified=="Specific" & NearbyGenes_plot$Approach == "Window",], aes(x=Window,y=upregulated_perC*100,group=TE_cluster)) + 
  #stat_summary(fun = 'mean', colour = 'black', geom = 'line')+
  stat_summary(fun.data= data_summary_95quantile, geom = 'ribbon', alpha = 0.2, aes(fill=TE_cluster))+
  #stat_summary(fun.data= 'mean_cl_normal', geom = 'ribbon', alpha = 0.5, fun.args=(conf.int=0.95), fill="lightblue")+
  #geom_point(data =NearbyGenes[!grepl("shuffle",NearbyGenes$dataset),], aes(x=Window,y=upregulated_perC*100,color=TE_cluster,alpha=up_sig_count))+
  geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$TE_cluster == "a",],
            aes(x=Window,y=upregulated_perC*100,group=Type_modified,linetype=Type_modified),color="#7C7C7C")+
  geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$TE_cluster == "b",],
            aes(x=Window,y=upregulated_perC*100,group=Type_modified,linetype=Type_modified),color="#bf812d")+
  geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$TE_cluster == "NI",],
            aes(x=Window,y=upregulated_perC*100,group=Type_modified,linetype=Type_modified),color="#beaed4")+
  ylim(15,30)+
  ylab("% up-regulated genes")+
  scale_linetype_manual(values=lineType) +
  scale_color_manual(values=colorSet)+
  scale_fill_manual(values=colorSet)+
  my_theme_1_Legend


p2 = ggplot(NearbyGenes_plot[grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type_modified=="Specific" & NearbyGenes_plot$Approach == "Window",], aes(x=Window,y=downregulated_perC*100,group=TE_cluster)) + 
  stat_summary(fun.data= data_summary_95quantile, geom = 'ribbon', alpha = 0.2, aes(fill=TE_cluster))+
  #geom_point(data =NearbyGenes[!grepl("shuffle",NearbyGenes$dataset),], aes(x=Window,y=downregulated_perC*100,color=TE_cluster,alpha=down_sig_count))+
  geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$TE_cluster == "a",],
            aes(x=Window,y=downregulated_perC*100,group=Type_modified,linetype=Type_modified),color="#7C7C7C")+
  geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$TE_cluster == "b",],
            aes(x=Window,y=downregulated_perC*100,group=Type_modified,linetype=Type_modified),color="#bf812d")+
  geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$TE_cluster == "NI",],
            aes(x=Window,y=downregulated_perC*100,group=Type_modified,linetype=Type_modified),color="#beaed4")+
  ylim(15,30)+
  ylab("% down-regulated genes")+
  scale_linetype_manual(values=lineType) +
  scale_color_manual(values=colorSet)+
  scale_fill_manual(values=colorSet)+
  my_theme_1_Legend

pdf(paste("EMC-Figure5A_1_2021_12_12",".pdf",sep=""),                      
    width = 6,
    height = 4.5,
    pointsize = 10)        # smaller font size
grid.draw(p1)
dev.off()
pdf(paste("EMC-Figure5B_2_2021_12_12",".pdf",sep=""),                      
    width = 6,
    height = 4.5,
    pointsize = 10)        # smaller font size
grid.draw(p2)
dev.off()


############ separate by different peaks
NearbyGenes = read.delim("inputs/TE_regulation.summary_2021_12_12_broad",sep="",header=F)

### read Window info
WindowSize = read.delim("inputs/TE_regulation.summary_2021_12_12_broad.windowSize",header=F,sep="")

### format
NearbyGenes$Group = ifelse(grepl("shuffle",NearbyGenes$V1),"Shuffle","NonShuffling")
NearbyGenes$Approach = ifelse(grepl("closest",NearbyGenes$V1),"Closest","Window")
colnames(NearbyGenes) = c("dataset","TE_cluster","Type","Window","up_sig_count","down_sig_count","not_sig_count","upregulated_perC","downregulated_perC","Group","Approach")
# filter out pairs longer than 1Mb
NearbyGenes = NearbyGenes[NearbyGenes$Window!= 0,]

### different peak type
#NearbyGenes = NearbyGenes[(NearbyGenes$TE_cluster == "NI" & NearbyGenes$Type != "Flu_specific") | (NearbyGenes$TE_cluster != "NI" & NearbyGenes$Type != "NI_specific"),]

# Analyze by using Window approach
NearbyGenes = NearbyGenes[NearbyGenes$Approach == "Window",]
# Windows
NearbyGenes$Window = factor(NearbyGenes$Window,levels=c(as.character(WindowSize$V1)))
# Corrected the Type of peaks
NearbyGenes$Type_modified = factor(ifelse(as.character(NearbyGenes$Type) == "All","All","Specific"))
NearbyGenes$total_count = NearbyGenes$up_sig_count+NearbyGenes$down_sig_count+NearbyGenes$not_sig_count

NearbyGenes_plot = NearbyGenes

######## Different TE clusters
TEcluster = "NI"
for (TEcluster in c("a","b","NI")){
  NearbyGenes_plot = NearbyGenes
  NearbyGenes_plot = NearbyGenes_plot[NearbyGenes_plot$TE_cluster==TEcluster,]
  NearbyGenes_plot = NearbyGenes_plot[!grepl("sort_K|sort_noK|MEF2",NearbyGenes_plot$dataset),]
  
  NearbyGenes_plot$PointAlpha = ifelse(NearbyGenes_plot$total_count<=100,NearbyGenes_plot$total_count/100,1)
  NearbyGenes_plot$PointSize = ifelse(NearbyGenes_plot$total_count>=100,100,NearbyGenes_plot$total_count)
  
  ### plots
  lineType = c("All"="dotted","Specific"="solid")
  shapeType = c("All" = 19, "noMEF2" = 17,"MEF2" = 15)
  colorset2 = c("NI_specific" = "#636363","Flu_specific" = "#99000d","Shared"="#a6761d","lowCount"="green")
  colorSet = c("a" = "#7C7C7C", "b" = "#bf812d","NI" = "#beaed4","MEF2" = "#081d58","lowCount"="green")
  p1 = ggplot(NearbyGenes_plot[grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Approach == "Window",], aes(x=Window,y=upregulated_perC*100,group=Type)) + 
    stat_summary(fun.data= data_summary_95quantile, geom = 'ribbon', alpha = 0.2, aes(fill=Type))+
    geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "Flu_specific",],
              aes(x=Window,y=upregulated_perC*100,color=Type))+
    geom_point(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "Flu_specific",],
               aes(x=Window,y=upregulated_perC*100,color=Type))+
    geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "NI_specific",],
              aes(x=Window,y=upregulated_perC*100,color=Type))+
    geom_point(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "NI_specific",],
               aes(x=Window,y=upregulated_perC*100,color=Type))+
    geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "Shared",],
              aes(x=Window,y=upregulated_perC*100,color=Type))+
    geom_point(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "Shared",],
               aes(x=Window,y=upregulated_perC*100,color=Type))+
    ylab("% up-regulated genes")+
    ylim(0,50)+
    scale_color_manual(values=colorset2)+
    scale_shape_manual(values=shapeType)+
    scale_fill_manual(values=colorset2)+
    scale_size(range = c(1,4))+
    my_theme_1_Legend
  
  
  p2 = ggplot(NearbyGenes_plot[grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Approach == "Window",], aes(x=Window,y=downregulated_perC*100,group=Type)) + 
    stat_summary(fun.data= data_summary_95quantile, geom = 'ribbon', alpha = 0.2, aes(fill=Type))+
    geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "Flu_specific",],
              aes(x=Window,y=downregulated_perC*100,color=Type))+
    geom_point(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "Flu_specific",],
               aes(x=Window,y=downregulated_perC*100,color=Type))+
    geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "NI_specific",],
              aes(x=Window,y=downregulated_perC*100,color=Type))+
    geom_point(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "NI_specific",],
               aes(x=Window,y=downregulated_perC*100,color=Type))+
    geom_line(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "Shared",],
              aes(x=Window,y=downregulated_perC*100,color=Type))+
    geom_point(data = NearbyGenes_plot[!grepl("shuffle",NearbyGenes_plot$dataset) & NearbyGenes_plot$Type == "Shared",],
               aes(x=Window,y=downregulated_perC*100,color=Type))+
    ylab("% down-regulated genes")+
    ylim(0,50)+
    scale_color_manual(values=colorset2)+
    scale_shape_manual(values=shapeType)+
    scale_fill_manual(values=colorset2)+
    my_theme_1_Legend
  
  
  pdf(paste("EMC-Figure_S4A_",TEcluster,"_up_2021_12_12",".pdf",sep=""),                      
      width = 6,
      height = 4.5,
      pointsize = 10)        # smaller font size
  grid.draw(p1)
  dev.off()
  pdf(paste("EMC-Figure_S4A_",TEcluster,"_down_2021_12_12",".pdf",sep=""),                      
      width = 6,
      height = 4.5,
      pointsize = 10)        # smaller font size
  grid.draw(p2)
  dev.off()
} 
