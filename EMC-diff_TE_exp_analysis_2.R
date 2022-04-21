library(gplots)
library(RColorBrewer)
library(devtools)
library(MASS)
library(ggplot2)
library(grid)
library(ggrepel)
library(dplyr)
library(ggbeeswarm)
library(splitstackshape)
library(gridExtra)

setwd(dir="/Users/xchen/OneDrive - Kyoto University/Projects_Bourque/EMC_project_2020_12_1/Final_Rscripts_2022_4_20/")

pValue = 0.001
pvaluetable<-read.delim("inputs/EMC_RNAseq_hg19.rename.CountM_pairedEnd_test_gene_TE_analysis.txt",sep="",header=T)
pvaluetable_TE<-pvaluetable[order(pvaluetable$padj),]
pvaluetable_TE2<-pvaluetable_TE[grepl(":",rownames(pvaluetable_TE)),]
pvaluetable_TE2$TEfamily = rownames(pvaluetable_TE2)
pvaluetable_TE2$TEfamily2 = rownames(pvaluetable_TE2)
pvaluetable_TE2 = data.frame(cSplit(pvaluetable_TE2,"TEfamily",sep=":",type.convert = as.character))
pvaluetable_TE2 = pvaluetable_TE2[!is.na(pvaluetable_TE2$pvalue),]

rownames(pvaluetable_TE2) = pvaluetable_TE2$TEfamily2
TEGroups = c("DNA","LINE","SINE","LTR","Other")                       ##### for hg19, it also include the SVA group
pvaluetable_TE2 = pvaluetable_TE2[pvaluetable_TE2$TEfamily_3 %in% TEGroups,]
mutateddf <- mutate(pvaluetable_TE2, sig=ifelse(pvaluetable_TE2$padj<=pValue & abs(pvaluetable_TE2$log2FoldChange)>=1, paste("padj<=",pValue,sep=""), "Not Sig")) #Will have different colors depending on significance
pvaluetable_TE3 = cbind(pvaluetable_TE2,mutateddf$sig)
colnames(pvaluetable_TE3)[5] = "sig"
pvaluetable_TE3[pvaluetable_TE3$TEfamily_3 == "Other",]

write.csv(pvaluetable_TE3,file="EMC_RNAseq_hg19.rename.CountM_supplementary_table_2021_6_24.csv")

inputData = pvaluetable_TE3
inputData$Color = ifelse(inputData$padj<=pValue & inputData$log2FoldChange >= 1, paste("up_padj<",pValue,sep=""), "Not Sig")
inputData$Color = ifelse(inputData$padj<=pValue & inputData$log2FoldChange <= -1, paste("down_padj<",pValue,sep=""), inputData$Color)
inputData$Color = factor(inputData$Color)

p<-ggplot(inputData, aes(log2FoldChange, -log10(padj))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Color),na.rm=TRUE,shape=19) + #add points colored by significance
  scale_color_manual(values=c("#2166ac","#4d4d4d", "#b2182b")) + 
  #ggtitle(outputFile) + #e.g. 'Volcanoplot DESeq2'
  geom_text_repel(data=head(inputData, 20), aes(label=head(as.character(inputData$TEfamily_1),20))) + #adding text for the top 20 genes
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
pdf("EMC_Figure1C_2021_6_24.pdf",                      
    width = 6,
    height = 4,
    pointsize = 10)        # smaller font size
grid.draw(p)
dev.off()


################ barplot EMC_Figure 1E.pdf
pvaluetable_TE<-pvaluetable[order(pvaluetable$padj),]
pvaluetable_TE2<-pvaluetable_TE[grepl(":",rownames(pvaluetable_TE)),]
pvaluetable_TE2$TEfamily = rownames(pvaluetable_TE2)
pvaluetable_TE2$TEfamily2 = rownames(pvaluetable_TE2)

pvaluetable_TE2 = data.frame(cSplit(pvaluetable_TE2,"TEfamily",sep=":"))

rownames(pvaluetable_TE2) = pvaluetable_TE2$TEfamily2

### not include SVA
TEGroups = c("DNA","LINE","SINE","LTR")
pvaluetable_TE2 = pvaluetable_TE2[pvaluetable_TE2$TEfamily_3 %in% TEGroups,]


pvaluetable_TE2$Color = ifelse(pvaluetable_TE2$padj<=pValue & pvaluetable_TE2$log2FoldChange >=1, "up 1log2FC", "Not Sig")
pvaluetable_TE2$Color = ifelse(is.na(pvaluetable_TE2$Color),"up not-sig",pvaluetable_TE2$Color)
pvaluetable_TE2$Color = ifelse(pvaluetable_TE2$padj<=pValue & (pvaluetable_TE2$log2FoldChange >=0.5 & pvaluetable_TE2$log2FoldChange < 1), "up 0.5log2FC", pvaluetable_TE2$Color)
pvaluetable_TE2$Color = ifelse(pvaluetable_TE2$log2FoldChange >0 & pvaluetable_TE2$Color == "Not Sig","up not-sig", pvaluetable_TE2$Color)
pvaluetable_TE2$Color = ifelse(pvaluetable_TE2$padj<=pValue & pvaluetable_TE2$log2FoldChange <= -1, "down 1log2FC", pvaluetable_TE2$Color)
pvaluetable_TE2$Color = ifelse(pvaluetable_TE2$padj<=pValue & (pvaluetable_TE2$log2FoldChange <=-0.5 & pvaluetable_TE2$log2FoldChange >-1), "down 0.5log2FC", pvaluetable_TE2$Color)
pvaluetable_TE2$Color = ifelse(pvaluetable_TE2$log2FoldChange <0 & pvaluetable_TE2$Color == "Not Sig","down not-sig", pvaluetable_TE2$Color)
pvaluetable_TE2$Color = factor(pvaluetable_TE2$Color)

pvaluetable_TE_summary = data.frame(pvaluetable_TE2 %>% group_by(TEfamily_3,Color) %>% dplyr::summarise(n = n(),na.rm=TRUE))
pvaluetable_TE_summary2 = data.frame(pvaluetable_TE2 %>% group_by(TEfamily_3) %>% dplyr::summarise(n = n(),na.rm=TRUE))
pvaluetable_TE_summary = merge(pvaluetable_TE_summary,pvaluetable_TE_summary2,by="TEfamily_3",all.x=T)
pvaluetable_TE_summary$perC = ifelse(grepl("up ",pvaluetable_TE_summary$Color),(pvaluetable_TE_summary$n.x/pvaluetable_TE_summary$n.y)*100,-(pvaluetable_TE_summary$n.x/pvaluetable_TE_summary$n.y)*100)
pvaluetable_TE_summary$Color = factor(pvaluetable_TE_summary$Color,levels=c("down 1log2FC", "down 0.5log2FC", "down not-sig", "up 1log2FC","up 0.5log2FC", "up not-sig"))
Colors = c("up 1log2FC"="#b2182b","up 0.5log2FC"="#d6604d","up not-sig"="#4d4d4d","down not-sig"="#4d4d4d","down 0.5log2FC"="#4393c3","down 1log2FC"="#2166ac")
pvaluetable_TE_summary$TEfamily_3 = factor(pvaluetable_TE_summary$TEfamily_3,levels=c("LTR","DNA","LINE","SINE"))
levels(pvaluetable_TE_summary$TEfamily_3) = c(paste("LTR(",as.character(pvaluetable_TE_summary[pvaluetable_TE_summary$TEfamily_3=="LTR",]$n.y)[1],")",sep=""),
                                              paste("DNA(",as.character(pvaluetable_TE_summary[pvaluetable_TE_summary$TEfamily_3=="DNA",]$n.y)[1],")",sep=""),
                                              paste("LINE(",as.character(pvaluetable_TE_summary[pvaluetable_TE_summary$TEfamily_3=="LINE",]$n.y)[1],")",sep=""),
                                              paste("SINE(",as.character(pvaluetable_TE_summary[pvaluetable_TE_summary$TEfamily_3=="SINE",]$n.y)[1],")",sep=""))

p<-ggplot(pvaluetable_TE_summary,aes(x=TEfamily_3,y=perC, fill=Color))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=Colors) + 
  geom_hline(yintercept = 0, size = 0.8)+
  coord_flip()+
  ylab("%")+
  xlab("TE group")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        
        axis.text=element_text(colour="black",size=rel(1.2),angle = 0),
        axis.title=element_text(colour="black",size=rel(1.5)),
        #        legend.position=c(0.8,0.8),
        panel.border = element_rect(colour = "black",size = 0.8,fill=NA),
        legend.position="right",
        legend.text=element_text(size=rel(1.2)))
pdf("EMC-Figure1C_2021_6_24_2.pdf",                      
    width = 6,
    height = 1.6,
    pointsize = 10)        # smaller font size
grid.draw(p)
dev.off()