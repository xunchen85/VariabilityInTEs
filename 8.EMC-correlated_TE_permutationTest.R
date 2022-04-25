###     Author:     Xun Chen
###     Date:       2022/4/20
###     Contact:    xunchen85@gmail.com


##### load pacakges
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(ggrepel)
library(dplyr)
library(splitstackshape)
library(gridExtra)


memory.limit(size=10240000)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
###############################
####### Set up enviroment;
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

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

####################################
TEgroups = c("DNA","LTR","LINE","SINE")
colorsetTE = c("DNA"="#1f78b4","LTR"="#e31a1c","LINE"="#33a02c","SINE"="#ff7f00")

##### read the table
Gene_exp_tpm_global_converted2 = read.csv("inputs/EMC_tpm_global_regression_2021_9_13.csv")
Gene_exp_tpm_TEs = Gene_exp_tpm_global_converted2[Gene_exp_tpm_global_converted2$TEinfo_3 %in% TEgroups,]
Gene_exp_tpm_TEs = Gene_exp_tpm_TEs[Gene_exp_tpm_TEs$TEinfo_3 != "Other",]
Gene_exp_tpm_TEs = Gene_exp_tpm_TEs[,!grepl("^EU",colnames(Gene_exp_tpm_TEs)) & !grepl("^AF",colnames(Gene_exp_tpm_TEs))]
Gene_exp_tpm_TEs$Group = ifelse(as.numeric(as.character(Gene_exp_tpm_TEs$log2_1_r.squared))>=0.3 & as.numeric(as.character(Gene_exp_tpm_TEs$log2_1_direction))>0,"Positive",NA)
Gene_exp_tpm_TEs$Group = ifelse(as.numeric(as.character(Gene_exp_tpm_TEs$log2_1_r.squared))>=0.3 & as.numeric(as.character(Gene_exp_tpm_TEs$log2_1_direction))<0,"Negative",Gene_exp_tpm_TEs$Group)

######################################## random sampling
# positive association
countUp_p = nrow(Gene_exp_tpm_TEs[Gene_exp_tpm_TEs$Group=="Positive" & !is.na(Gene_exp_tpm_TEs$Group),])
countDown_p = nrow(Gene_exp_tpm_TEs[Gene_exp_tpm_TEs$Group=="Negative" & !is.na(Gene_exp_tpm_TEs$Group),])

############################### making plots
Combination = "Positive"
Combinations = c("Negative","Positive")
Permutation_plot_points_combined = data.frame("Order"=NA,"Count"=NA,"TE"=NA,"pvalue"=NA,"pvalue_final"=NA,"Group"=NA)
SumTable = Gene_exp_tpm_TEs
j = 1
qlist = list()
for (Combination in Combinations){
  if (Combination == "Positive"){
    CountType = countUp_p
  } else if (Combination == "Negative"){
    CountType = countDown_p
  }
  
  ### Simulation
  PermuSum = data.frame("Order"=NA,"LTR"=NA,"DNA"=NA,"LINE"=NA,"SINE"=NA,"ERV1"=NA,"ERVK"=NA,"ERVL"=NA,"ERVL-MaLR"=NA)
  colorset = c("DNA"="#575A6C","LTR"="#B4C540","LINE"="#E0E2D2","SINE"="#3686C9","ERV1"="#78c679","ERVK"="#addd8e","ERVL"="#d9f0a3","ERVL.MaLR"="#f7fcb9")
  PermuSum = PermuSum[-1,]
  for (i in 1:10000){
    tmp1 = SumTable[sample(nrow(SumTable),CountType),]
    PermuSum1 = data.frame("Order"=NA,"LTR"=NA,"DNA"=NA,"LINE"=NA,"SINE"=NA,"ERV1"=NA,"ERVK"=NA,"ERVL"=NA,"ERVL-MaLR"=NA)
    PermuSum1[1,1] = i
    PermuSum1[1,2] = nrow(tmp1[tmp1$TEinfo_3=="LTR",])/nrow(SumTable[SumTable$TEinfo_3=="LTR",])
    PermuSum1[1,3] = nrow(tmp1[tmp1$TEinfo_3=="DNA",])/nrow(SumTable[SumTable$TEinfo_3=="DNA",])
    PermuSum1[1,4] = nrow(tmp1[tmp1$TEinfo_3=="LINE",])/nrow(SumTable[SumTable$TEinfo_3=="LINE",])
    PermuSum1[1,5] = nrow(tmp1[tmp1$TEinfo_3=="SINE",])/nrow(SumTable[SumTable$TEinfo_3=="SINE",])
    PermuSum1[1,6] = nrow(tmp1[tmp1$TEinfo_2=="ERV1",])/nrow(SumTable[SumTable$TEinfo_2=="ERV1",])
    PermuSum1[1,7] = nrow(tmp1[tmp1$TEinfo_2=="ERVK",])/nrow(SumTable[SumTable$TEinfo_2=="ERVK",])
    PermuSum1[1,8] = nrow(tmp1[tmp1$TEinfo_2=="ERVL",])/nrow(SumTable[SumTable$TEinfo_2=="ERVL",])
    PermuSum1[1,9] = nrow(tmp1[tmp1$TEinfo_2=="ERVL-MaLR",])/nrow(SumTable[SumTable$TEinfo_2=="ERVL-MaLR",])
    PermuSum = rbind(PermuSum,PermuSum1)
  }
  
  ### select for the simulated data
  Permutation_plot = data.frame("Order"=NA,"Count"=NA,"TE"=NA)
  Permutation_plot = Permutation_plot[-1,]
  for (i in 2:ncol(PermuSum)){
    Permutation_plot1 = PermuSum[,c(1,i)]
    Permutation_plot1$TE = colnames(PermuSum[i])
    colnames(Permutation_plot1) = c("Order","Count","TE")
    Permutation_plot = rbind(Permutation_plot,Permutation_plot1)
  }
  
  ### re-format and compute the proportion, actual data
  Permutation_plot_points = data.frame("Order"=NA,"Count"=NA,"TE"=NA)
  Permutation_plot_points = Permutation_plot_points[-1,]
  te = "LTR"
  for (te in levels(factor(Permutation_plot$TE))){
    Permutation_plot_points1 = data.frame("Order"=NA,"Count"=NA,"TE"=NA)
    Permutation_plot_points1[1,1] = 0
    te = gsub("ERVL.MaLR","ERVL-MaLR",te)
    if (Combination == "Positive"){
      if (te %in% SumTable$TEinfo_3){          #### SumTable[SumTable$Group=="Up_sig" is for p-value, should be changed to Group_padj for padjust
        Permutation_plot_points1[1,2] = nrow(SumTable[SumTable$Group=="Positive" & !is.na(SumTable$Group) & SumTable$TEinfo_3 == te,])/nrow(SumTable[SumTable$TEinfo_3 == te,])
      } else{
        Permutation_plot_points1[1,2] = nrow(SumTable[SumTable$Group=="Positive" & !is.na(SumTable$Group) & SumTable$TEinfo_2 == te,])/nrow(SumTable[SumTable$TEinfo_2 == te,])
      }
    } else if (Combination == "Negative"){
      if (te %in% SumTable$TEinfo_3){          #### SumTable[SumTable$Group=="Up_sig" is for p-value, should be changed to Group_padj for padjust
        Permutation_plot_points1[1,2] = nrow(SumTable[SumTable$Group=="Negative" & !is.na(SumTable$Group) & SumTable$TEinfo_3 == te,])/nrow(SumTable[SumTable$TEinfo_3 == te,])
      } else{
        Permutation_plot_points1[1,2] = nrow(SumTable[SumTable$Group=="Negative" & !is.na(SumTable$Group) & SumTable$TEinfo_2 == te,])/nrow(SumTable[SumTable$TEinfo_2 == te,])
      }
    } 
    Permutation_plot_points1[1,3] = te
    Permutation_plot_points = rbind(Permutation_plot_points,Permutation_plot_points1)
  }
  Permutation_plot_points$TE =  gsub("ERVL-MaLR","ERVL.MaLR",Permutation_plot_points$TE)
  ### z test
  te = "ERV1"
  Permutation_plot_points$pvalue = NA
  Permutation_plot_points$pvalue_final = NA
  for (te in levels(factor(Permutation_plot$TE))){
    Permutation_plot_points[Permutation_plot_points$TE == te,]$pvalue_final = 2*mean(Permutation_plot[Permutation_plot$TE == te,]$Count>=Permutation_plot_points[Permutation_plot_points$TE == te,]$Count)
    ttest = t.test(Permutation_plot[Permutation_plot$TE == te,]$Count, mu=Permutation_plot_points[Permutation_plot_points$TE == te,]$Count, alternative="less", conf.level=0.99)
    Permutation_plot_points[Permutation_plot_points$TE == te,]$pvalue = ttest$p.value
  }
  Permutation_plot_points$Group = Combination
  Permutation_plot_points_combined = rbind(Permutation_plot_points_combined,Permutation_plot_points)
  
  ### plot
  Permutation_plot$TE = factor(Permutation_plot$TE,levels=c("DNA","LINE","SINE","LTR","ERV1","ERVK","ERVL","ERVL.MaLR"))
  
  p<-ggplot(Permutation_plot, aes(TE, Count*100,group=TE))+
    xlab("TE group and superfamily")+
    ylab("% Permuted subfamilies")+
    geom_violin(aes(fill = TE,color = TE))+
    scale_fill_manual(values=colorset)+
    scale_color_manual(values=colorset)+
    stat_summary(fun.data=data_summary,geom="pointrange",color="black",shape=20)+
    geom_point(data = Permutation_plot_points, aes(TE, Count*100),shape=17,size = 2,color="#3f007d")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          
          axis.text=element_text(colour="black",size=rel(1),angle = 0),
          axis.text.x=element_text(colour="black",size=rel(1),angle = 45,vjust = 0.9),
          
          axis.title=element_text(colour="black",size=rel(1)),
          legend.position="none",
          legend.background = element_blank(),
          legend.text=element_text(size=rel(1.2)))
  qlist[[j]] <- ggplotGrob(p)
  j = j+1
}
pdf("EMC-Figure_S1C_permutation_test_2022_4_20.pdf",
    width = 4,
    height = 5,
    pointsize = 10)        # smaller font size
do.call("grid.arrange",c(qlist,ncol=1))
dev.off()
