###     Author:     Xun Chen
###     Date:       2021/9/12
###     Contact:    xunchen85@gmail.com

library(ggplot2)
library(grid)
library(dplyr)
library(ggpmisc)
library(splitstackshape)
library(logspline)
library(gridExtra)
library(ggrepel)

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
setwd(dir="/Users/xchen/OneDrive - Kyoto University/Projects_Bourque/EMC_project_2020_12_1/Final_Rscripts_2022_4_20/")
TEgroups = c("DNA","LTR","LINE","SINE")

colorsetTE = c("DNA"="#1f78b4","LTR"="#e31a1c","LINE"="#33a02c","SINE"="#ff7f00")

##################################### read log2FC
pValue = 0.001
pvaluetable<-read.delim("inputs/EMC_RNAseq_hg19.rename.CountM_pairedEnd_test_gene_TE_analysis.txt",sep="",header=T)
Viralload = read.csv("inputs/H1N1_Viral_Load_5_27.csv")

pvaluetable_TE<-pvaluetable[order(pvaluetable$padj),]
pvaluetable_TE2<-pvaluetable_TE[grepl(":",rownames(pvaluetable_TE)),]
pvaluetable_TE2$TEfamily = rownames(pvaluetable_TE2)
pvaluetable_TE2$TEfamily2 = rownames(pvaluetable_TE2)
pvaluetable_TE2 = data.frame(cSplit(pvaluetable_TE2,"TEfamily",sep=":"))
pvaluetable_TE2 = pvaluetable_TE2[!is.na(pvaluetable_TE2$pvalue),]

rownames(pvaluetable_TE2) = pvaluetable_TE2$TEfamily2
TEGroups = c("DNA","LINE","SINE","LTR","Other")                       ##### for hg19, it also include the SVA group
pvaluetable_TE2 = pvaluetable_TE2[pvaluetable_TE2$TEfamily_3 %in% TEGroups,]
mutateddf <- mutate(pvaluetable_TE2, sig=ifelse(pvaluetable_TE2$padj<=pValue & abs(pvaluetable_TE2$log2FoldChange)>=1, paste("padj<=",pValue,sep=""), "Not Sig")) #Will have different colors depending on significance
pvaluetable_TE3 = cbind(pvaluetable_TE2,mutateddf$sig)
colnames(pvaluetable_TE3)[5] = "sig"
pvaluetable_TE3[pvaluetable_TE3$TEfamily_3 == "Other",]

inputData = pvaluetable_TE3
inputData$Color = ifelse(inputData$padj<=pValue & inputData$log2FoldChange >= 1, paste("up_padj<",pValue,sep=""), "Not Sig")
inputData$Color = ifelse(inputData$padj<=pValue & inputData$log2FoldChange <= -1, paste("down_padj<",pValue,sep=""), inputData$Color)
inputData$Color = factor(inputData$Color)

#################################### read the table of log2FC
Gene_exp_tpm_global_converted2 = read.csv("inputs/EMC_tpm_global_regression_2021_9_13.csv")
Gene_exp_tpm_global_converted2$log2_1_r.squared_withDirection = ifelse(as.numeric(as.character(Gene_exp_tpm_global_converted2$log2_1_direction))>0,as.numeric(as.character(Gene_exp_tpm_global_converted2$log2_1_r.squared)),(-1 * as.numeric(as.character(Gene_exp_tpm_global_converted2$log2_1_r.squared))))
head(Gene_exp_tpm_global_converted2)
#################################### combine datasets and plot
inputData_correlation = merge(inputData,Gene_exp_tpm_global_converted2[,c("TEinfo2","log2_1_r.squared_withDirection")],by.x="TEfamily2",by.y="TEinfo2",all.x=T)
inputData_correlation = inputData_correlation[order(match(inputData_correlation$TEfamily2, inputData$TEfamily2)), ]
head(inputData_correlation)

p<-ggplot(inputData_correlation, aes(log2FoldChange, log2_1_r.squared_withDirection)) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Color),na.rm=TRUE,shape=19) + #add points colored by significance
  scale_color_manual(values=c("#2166ac","#4d4d4d", "#b2182b")) + 
  geom_hline(yintercept=c(-0.3,0,0.3), linetype='dashed')+
  geom_text_repel(data=head(inputData_correlation, 20), aes(label=head(as.character(inputData_correlation$TEfamily_1),20))) + #adding text for the top 20 genes
  xlab("R square")+
  ylab("log2FC")+
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
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.text=element_text(size=rel(1.2)))
pdf("EMC-Figure1D_log2FC_2021_10_23.pdf",                      
    width = 9.5,
    height = 6,
    pointsize = 10)        # smaller font size
grid.draw(p)
dev.off()

#################################### examples
Viralload$individualID = gsub("_Flu|_NI","",Viralload$Sample_ID)
Correlation_examples = data.frame(t(Gene_exp_tpm_global_converted2[Gene_exp_tpm_global_converted2$TEinfo_1 %in% c("LTR76","MER84-int","MER61F","PABL_A-int"),grepl("_log2",colnames(Gene_exp_tpm_global_converted2)) | grepl("TEinfo_1",colnames(Gene_exp_tpm_global_converted2))]))
colnames(Correlation_examples) = as.character(Correlation_examples["TEinfo_1",])
Correlation_examples = Correlation_examples[-1,]
Correlation_examples$individualID = gsub("_log2","",rownames(Correlation_examples))
Correlation_examples = merge(Correlation_examples,Viralload,by = "individualID",all.x = T)
my.formula = y~x

###
glist = list()
Order = 1
i = 2
for (i in seq(2,5,1)){
  p1 = ggplot(Correlation_examples,aes(x=ViralLoad_...,y=as.numeric(as.character(Correlation_examples[,i])))) + 
    geom_point(aes(),size=1)+
    xlab("Viral load (%)") + 
    ylab(paste("log2FC",colnames(Correlation_examples)[i])) +
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 p.digits = 3,
                 aes(label = paste(..rr.label.., sep = "*plain(\",\")~")), 
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
  i= i + 1
}

pdf("EMC-Figure1EF-log2FC_2021_10_23.pdf",                      
    width = 6,        # 5 x 300 pixels
    height = 6,
    pointsize = 10)
do.call("grid.arrange",c(glist,ncol=2))
dev.off()

################################### read the table of correlation
Gene_exp_tpm_global_converted2 = read.csv("inputs/EMC_tpm_global_regression_2021_9_13.csv")
Gene_exp_tpm_TEs = Gene_exp_tpm_global_converted2[Gene_exp_tpm_global_converted2$TEinfo_3 %in% TEgroups,]

############################# step 1
########### summary
Gene_exp_tpm_TEs_Flu =Gene_exp_tpm_TEs[Gene_exp_tpm_TEs$TEinfo_3!="Other" & as.numeric(as.character(Gene_exp_tpm_TEs$Flu_r.squared))>=0.3,]
Gene_exp_tpm_TEs_Flu$Direction =ifelse(as.numeric(as.character(Gene_exp_tpm_TEs_Flu$Flu_direction))>0,"positive","negative")
Gene_exp_tpm_TEs_Flu_sum = data.frame(Gene_exp_tpm_TEs_Flu %>% group_by(Direction) %>% dplyr::summarise(n = n()) %>% mutate(total = sum(n)) %>% mutate(freq = n/ sum(n)))
Gene_exp_tpm_TEs_Flu_sum$Group = "Flu"

Gene_exp_tpm_TEs_NI =Gene_exp_tpm_TEs[Gene_exp_tpm_TEs$TEinfo_3!="Other" & as.numeric(as.character(Gene_exp_tpm_TEs$NI_r.squared))>=0.3,]
Gene_exp_tpm_TEs_NI$Direction =ifelse(as.numeric(as.character(Gene_exp_tpm_TEs_NI$NI_direction))>0,"positive","negative")
Gene_exp_tpm_TEs_NI_sum = data.frame(Gene_exp_tpm_TEs_NI %>% group_by(Direction) %>% dplyr::summarise(n = n()) %>% mutate(total = sum(n)) %>% mutate(freq = n/ sum(n)))
Gene_exp_tpm_TEs_NI_sum$Group = "NI"

Gene_exp_tpm_TEs_log2FC =Gene_exp_tpm_TEs[Gene_exp_tpm_TEs$TEinfo_3!="Other" & as.numeric(as.character(Gene_exp_tpm_TEs$log2_1_r.squared))>=0.3,]
Gene_exp_tpm_TEs_log2FC$Direction =ifelse(as.numeric(as.character(Gene_exp_tpm_TEs_log2FC$log2_1_direction))>0,"positive","negative")
Gene_exp_tpm_TEs_log2FC_sum = data.frame(Gene_exp_tpm_TEs_log2FC %>% group_by(Direction) %>% dplyr::summarise(n = n()) %>% mutate(total = sum(n)) %>% mutate(freq = n/ sum(n)))
Gene_exp_tpm_TEs_log2FC_sum$Group = "log2FC"

### 
Summary = rbind(Gene_exp_tpm_TEs_Flu_sum,Gene_exp_tpm_TEs_NI_sum,Gene_exp_tpm_TEs_log2FC_sum)
Summary$Direction = ifelse(is.na(Summary$Direction),"Unknown",Summary$Direction)
Summary$Direction = factor(Summary$Direction,levels = c("Unknown","positive","negative"))
Summary$Group = factor(Summary$Group,levels=c("log2FC","NI","Flu"))

p = ggplot(Summary,aes(x=Group,y=freq*100,fill=Direction))+
  geom_col()+
  scale_fill_manual(values=c("white","#bdbdbd","#252525")) +
  geom_vline(xintercept = 0.3, size = 1,linetype = "dotted",color="black")+
  ylab("% correlated subfamilies")+
  xlab("")+
  coord_flip()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        axis.title=element_text(colour="black",size=rel(1.2)),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=rel(1.2)))
pdf(paste("EMC-Figure_S7A",".pdf",sep=""),     
    width = 6,        # 5 x 300 pixels
    height = 1.5,
    pointsize = 10)
grid.draw(p)
dev.off()

###################################
Combinations = c("Flu_r.squared","NI_r.squared","log2_1_r.squared")
### density
glist = list()
Order = 1
p<-ggplot(Gene_exp_tpm_TEs[Gene_exp_tpm_TEs$TEinfo_3!="Other",],aes(as.numeric(as.character(Flu_r.squared)),fill=TEinfo_3))+
  geom_density(alpha=0.4,adjust=0.4)+
  scale_fill_manual(values=colorsetTE) +
  scale_color_manual(values=colorsetTE) +
  xlim(0,0.8)+
  geom_vline(xintercept = 0.3, size = 1,linetype = "dotted",color="black")+
  ylab("Density")+
  xlab("R square")+
  ggtitle(gsub("_r.squared","",Combinations[Order])) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        axis.title=element_text(colour="black",size=rel(1.2)),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=rel(1.2)))
glist[[Order]] = p
Order = Order + 1
p<-ggplot(Gene_exp_tpm_TEs[Gene_exp_tpm_TEs$TEinfo_3!="Other",],aes(as.numeric(as.character(NI_r.squared)),fill=TEinfo_3))+
  geom_density(alpha=0.4,adjust=0.4)+
  scale_fill_manual(values=colorsetTE) +
  scale_color_manual(values=colorsetTE) +
  xlim(0,0.8)+
  geom_vline(xintercept = 0.3, size = 1,linetype = "dotted",color="black")+
  ylab("Density")+
  xlab("R square")+
  ggtitle(gsub("_r.squared","",Combinations[Order])) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        axis.title=element_text(colour="black",size=rel(1.2)),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=rel(1.2)))
glist[[Order]] = p
Order = Order + 1
p<-ggplot(Gene_exp_tpm_TEs[Gene_exp_tpm_TEs$TEinfo_3!="Other",],aes(as.numeric(as.character(log2_1_r.squared)),fill=TEinfo_3))+
  geom_density(alpha=0.4,adjust=0.4)+
  scale_fill_manual(values=colorsetTE) +
  scale_color_manual(values=colorsetTE) +
  xlim(0,0.8)+
  geom_vline(xintercept = 0.3, size = 1,linetype = "dotted",color="black")+
  ylab("Density")+
  xlab("R square")+
  ggtitle(gsub("_r.squared","",Combinations[Order])) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        axis.title=element_text(colour="black",size=rel(1.2)),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=rel(1.2)))
glist[[Order]] = p
pdf(paste("EMC-Figure6A_2021_9_12",".pdf",sep=""),     
    width = 6,        # 5 x 300 pixels
    height = 2* 3,
    pointsize = 10)
do.call("grid.arrange",c(glist,ncol=1))
dev.off()
