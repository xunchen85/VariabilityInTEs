###     Author:     Xun Chen
###     Date:       2021/3/4
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
####### variables
Exclude_columns = 3
Plus_1 = "no"
my.formula = y~x
gene_foldChange_threadhold = 1
pValue = 0.001

####### load inputs
## sample information
Viralload = read.csv("inputs/H1N1_Viral_Load_5_27.csv")
Viralload$individualID = gsub("_Flu|_NI","",Viralload$Sample_ID)

## TE list
TElist = read.delim("inputs/TElist.sep.27",header=T,sep="")

## gene Symbol
hg19_geneName = read.delim("inputs/hg19.ensGene.Genename",header=F,sep="")
colnames(hg19_geneName) = c("transcriptID","geneID","chr","direction","start","end","geneName")
hg19_geneName = hg19_geneName[!duplicated(hg19_geneName[,c("geneID","geneName")]),]

## Subfamily level
FoldChange = read.delim("EMC_RNAseq_hg19.rename.CountM_pairedEnd_test_gene_TE_analysis.txt",sep="",header=T)
FoldChange$is_sig = ifelse(FoldChange$log2FoldChange>=1 & FoldChange$padj<=pValue,"sig_up","not_sig")
FoldChange$is_sig = ifelse(FoldChange$log2FoldChange>=1 & FoldChange$padj<=pValue,"sig_down",FoldChange$is_sig)

Gene_foldChange = FoldChange[!grepl(":",rownames(FoldChange)),]
TE_foldChange = FoldChange[grepl(":",rownames(FoldChange)),]

## individual level expression change
TE_exp = read.delim("EMC_RNAseq_hg19.rename.CountM_Normalized_PvalueOutput.txt",sep="",header=T)
Gene_exp = TE_exp[!grepl(":",rownames(TE_exp)),c(1:78)]
TE_exp = TE_exp[grepl(":",rownames(TE_exp)),c(1:78)]

## individual level expression change (TPM)
Both_exp_tpm = read.csv("inputs/EMC_RNAseq_hg19.rename.CountM-TPM.csv")
rownames(Both_exp_tpm) = Both_exp_tpm$name
Both_exp_tpm = Both_exp_tpm[,c(-1,-2)]
Gene_exp_tpm = Both_exp_tpm[!grepl(":",rownames(Both_exp_tpm)),c(1:78)]
TE_exp_tpm = Both_exp_tpm[grepl(":",rownames(Both_exp_tpm)),c(1:78)]

################################## section 1: global-wide Gene expression (genome-wide gene correlation analysis 2021_3_14)
Gene_exp_tpm_global = Both_exp_tpm
Gene_exp_tpm_global = Gene_exp_tpm_global[,colnames(Gene_exp_tpm_global) %in% c(paste(Viralload$individualID,"_NI",sep=""),paste(Viralload$individualID,"_Flu",sep=""))]
Gene_exp_tpm_global$meanFlu = apply(Gene_exp_tpm_global[grepl("_Flu",colnames(Gene_exp_tpm_global))],1,mean,na.rm=T)
Gene_exp_tpm_global$meanNI = apply(Gene_exp_tpm_global[grepl("_NI",colnames(Gene_exp_tpm_global))],1,mean,na.rm=T)

## highly expressed genes and TEs
Gene_exp_tpm_global = Gene_exp_tpm_global[Gene_exp_tpm_global$meanFlu>=1 | Gene_exp_tpm_global$meanNI >=1 | grepl(":",rownames(Gene_exp_tpm_global)),]
Gene_exp_tpm_global[nrow(Gene_exp_tpm_global)+1,] = NA

rownames(Gene_exp_tpm_global)[nrow(Gene_exp_tpm_global)]= "viralLoad"
for (i in 1:(ncol(Gene_exp_tpm_global)-2)){
  Gene_exp_tpm_global[nrow(Gene_exp_tpm_global),i] = as.numeric(as.character(Viralload[gsub("_Flu|_NI","",Viralload$Sample_ID) == gsub("_Flu|_NI","",colnames(Gene_exp_tpm_global)[i]),]$ViralLoad_...))
}  

## add the log2 fold change
Gene_exp_tpm_global[,((ncol(Gene_exp_tpm_global)+1):(ncol(Gene_exp_tpm_global)+length(Viralload$individualID)*2))] = NA
colnames(Gene_exp_tpm_global)[(ncol(Gene_exp_tpm_global)-length(Viralload$individualID)*2+1):ncol(Gene_exp_tpm_global)] = c(Viralload$individualID,paste(Viralload$individualID,"_log2",sep=""))

individualID = "AF04"
for(individualID in Viralload$individualID) {
  Gene_exp_tpm_global[,which(colnames(Gene_exp_tpm_global)== individualID)] = log2(Gene_exp_tpm_global[,which(colnames(Gene_exp_tpm_global)== paste(individualID,"_Flu",sep=""))]) - log2(Gene_exp_tpm_global[,which(colnames(Gene_exp_tpm_global)== paste(individualID,"_NI",sep=""))])
  Gene_exp_tpm_global[,paste(individualID,"_log2",sep="")] = log2(Gene_exp_tpm_global[,which(colnames(Gene_exp_tpm_global)== paste(individualID,"_Flu",sep=""))]+0.01) - log2(Gene_exp_tpm_global[,which(colnames(Gene_exp_tpm_global)== paste(individualID,"_NI",sep=""))]+0.01)
}
Gene_exp_tpm_global$TEinfo = rownames(Gene_exp_tpm_global)

### association between genes and viral loads
Gene_exp_tpm_global_converted = data.frame(t(Gene_exp_tpm_global))
Gene_exp_tpm_global_converted$individualID = rownames(Gene_exp_tpm_global_converted)

### add viral load
Viralload2 = Viralload
Viralload2$individualID = paste(Viralload2$individualID,"_Flu",sep="")
Viralload3 = Viralload
Viralload3$individualID = paste(Viralload3$individualID,"_NI",sep="")
Viralload4 = Viralload
Viralload4$individualID = paste(Viralload4$individualID,"_log2",sep="")
Viralload_converted = rbind(Viralload,Viralload2,Viralload3,Viralload4)
Viralload_converted = Viralload_converted[,c("individualID","ViralLoad_...","Age","Group_individual")]
rm(Viralload2,Viralload3,Viralload4)
Gene_exp_tpm_global_converted = merge(Gene_exp_tpm_global_converted,Viralload_converted,by="individualID",all.x=T)
Gene_exp_tpm_global_converted[(nrow(Gene_exp_tpm_global_converted)+1):(nrow(Gene_exp_tpm_global_converted)+20),] = NA
Gene_exp_tpm_global_converted[(nrow(Gene_exp_tpm_global_converted)-19):nrow(Gene_exp_tpm_global_converted),]$individualID = c("Flu_r.squared","Flu_r.adj","Flu_pValue","Flu_direction","Flu_corr",
                                                                                                                              "NI_r.squared","NI_r.adj","NI_pValue","NI_direction","NI_corr",
                                                                                                                              "log2_r.squared","log2_r.adj","log2_pValue","log2_direction","log2_corr",
                                                                                                                              "log2_1_r.squared","log2_1_r.adj","log2_1_pValue","log2_1_direction","log2_1_corr")
rownames(Gene_exp_tpm_global_converted) = Gene_exp_tpm_global_converted$individualID
Gene_exp_tpm_global_converted[,"ViralLoad_..."] = as.character(Gene_exp_tpm_global_converted[,"ViralLoad_..."])

i = 2
for (i in 2:(ncol(Gene_exp_tpm_global_converted)-4)){
  Gene_exp_tpm_global_converted[,i] = as.character(Gene_exp_tpm_global_converted[,i])
  if (!any(is.na(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% Viralload$individualID,i]))) | is.infinite(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% Viralload$individualID,i]))))){
    lm_tmp = summary(lm(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% Viralload$individualID,i]))
                        ~ as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% Viralload$individualID,"ViralLoad_..."]))))
    corr.test.out = cor.test(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% Viralload$individualID,i]))
                        , as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% Viralload$individualID,"ViralLoad_..."])))
    Gene_exp_tpm_global_converted["log2_r.squared",i] = lm_tmp$r.squared
    Gene_exp_tpm_global_converted["log2_r.adj",i] = lm_tmp$adj.r.squared
    Gene_exp_tpm_global_converted["log2_pValue",i] = lm_tmp$coefficients[,"Pr(>|t|)"][2]
    Gene_exp_tpm_global_converted["log2_direction",i] = lm_tmp$coefficients[2,1]
    Gene_exp_tpm_global_converted["log2_corr",i] = corr.test.out$estimate[[1]]
    
  }
  if (!any(is.na(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_log2",sep=""),i]))) | is.infinite(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_log2",sep=""),i]))))){
    lm_tmp = summary(lm(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_log2",sep=""),i]))
                        ~ as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_log2",sep=""),"ViralLoad_..."]))))
    corr.test.out = cor.test(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_log2",sep=""),i]))
                        ,as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_log2",sep=""),"ViralLoad_..."])))
    Gene_exp_tpm_global_converted["log2_1_r.squared",i] = lm_tmp$r.squared
    Gene_exp_tpm_global_converted["log2_1_r.adj",i] = lm_tmp$adj.r.squared
    Gene_exp_tpm_global_converted["log2_1_pValue",i] = lm_tmp$coefficients[,"Pr(>|t|)"][2]
    Gene_exp_tpm_global_converted["log2_1_direction",i] = lm_tmp$coefficients[2,1]
    Gene_exp_tpm_global_converted["log2_1_corr",i] = corr.test.out$estimate[[1]]
    
  }
  if (!any(is.na(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_NI",sep=""),i]))) | is.infinite(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_NI",sep=""),i]))))){
    lm_tmp = summary(lm(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_NI",sep=""),i]))
                        ~ as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_NI",sep=""),"ViralLoad_..."]))))
    corr.test.out = cor.test(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_NI",sep=""),i]))
                       ,as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_NI",sep=""),"ViralLoad_..."])))
    Gene_exp_tpm_global_converted["NI_r.squared",i] = lm_tmp$r.squared
    Gene_exp_tpm_global_converted["NI_r.adj",i] = lm_tmp$adj.r.squared
    Gene_exp_tpm_global_converted["NI_pValue",i] = lm_tmp$coefficients[,"Pr(>|t|)"][2]
    Gene_exp_tpm_global_converted["NI_direction",i] = lm_tmp$coefficients[2,1]
    Gene_exp_tpm_global_converted["NI_corr",i] = corr.test.out$estimate[[1]]
    
  }
  if (!any(is.na(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_Flu",sep=""),i]))) | is.infinite(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_Flu",sep=""),i]))))){
    lm_tmp = summary(lm(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_Flu",sep=""),i]))
                        ~ as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_Flu",sep=""),"ViralLoad_..."]))))
    corr.test.out = cor.test(as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_Flu",sep=""),i]))
                        , as.numeric(as.character(Gene_exp_tpm_global_converted[Gene_exp_tpm_global_converted$individualID %in% paste(Viralload$individualID,"_Flu",sep=""),"ViralLoad_..."])))
    Gene_exp_tpm_global_converted["Flu_r.squared",i] = lm_tmp$r.squared
    Gene_exp_tpm_global_converted["Flu_r.adj",i] = lm_tmp$adj.r.squared
    Gene_exp_tpm_global_converted["Flu_pValue",i] = lm_tmp$coefficients[,"Pr(>|t|)"][2]
    Gene_exp_tpm_global_converted["Flu_direction",i] = lm_tmp$coefficients[2,1]
    Gene_exp_tpm_global_converted["Flu_corr",i] = corr.test.out$estimate[[1]]
  }
}

Gene_exp_tpm_global_converted2 = data.frame(t(Gene_exp_tpm_global_converted))
Gene_exp_tpm_global_converted2$TEinfo2 = Gene_exp_tpm_global_converted2$TEinfo
Gene_exp_tpm_global_converted2 = data.frame(cSplit(Gene_exp_tpm_global_converted2,"TEinfo",sep=":"))
Gene_exp_tpm_global_converted2 = merge(Gene_exp_tpm_global_converted2,hg19_geneName,by.x="TEinfo_1",by.y="geneID",all.x=T)
Gene_exp_tpm_global_converted2$is_TEorGene = ifelse(grepl(":",Gene_exp_tpm_global_converted2$TEinfo2),"TE",NA)
Gene_exp_tpm_global_converted2$geneName = as.character(Gene_exp_tpm_global_converted2$geneName)
Gene_exp_tpm_global_converted2$is_TEorGene = ifelse(!is.na(Gene_exp_tpm_global_converted2$geneName),"Gene",Gene_exp_tpm_global_converted2$is_TEorGene)
write.csv(Gene_exp_tpm_global_converted2,file="EMC_tpm_global_regression_2021_9_13.csv")

################################## section 2: detection of correlated genes
Gene_exp_tpm_global_converted2 = read.csv("EMC_tpm_global_regression_2021_9_13.csv")

Gene_exp_tpm_global_converted2[Gene_exp_tpm_global_converted2$geneName=="TRIM25",]
Gene_exp_tpm_global_converted2 = Gene_exp_tpm_global_converted2[Gene_exp_tpm_global_converted2$is_TEorGene=="Gene" & !is.na(Gene_exp_tpm_global_converted2$is_TEorGene),]

### NI and log2FC
Gene_exp_tpm_global_converted2_NI = Gene_exp_tpm_global_converted2[as.numeric(as.character(Gene_exp_tpm_global_converted2$NI_r.squared))>=0.3 & Gene_exp_tpm_global_converted2$NI_r.squared != "NaN",]
Gene_exp_tpm_global_converted2_log2FC = Gene_exp_tpm_global_converted2[Gene_exp_tpm_global_converted2$TEinfo_1 %in% rownames(Gene_foldChange[Gene_foldChange$is_sig!="not_sig" & !is.na(Gene_foldChange$is_sig),]),]
Gene_exp_tpm_global_converted2_log2FC = Gene_exp_tpm_global_converted2_log2FC[as.numeric(as.character(Gene_exp_tpm_global_converted2_log2FC$log2_r.squared))>=0.3 & !is.na(Gene_exp_tpm_global_converted2_log2FC$log2_r.squared),]
write.csv(Gene_exp_tpm_global_converted2_NI,file="EMC_tpm_global_regression_2021_9_21_NIcorrelated.csv")
write.csv(Gene_exp_tpm_global_converted2_log2FC,file="EMC_tpm_global_regression_2021_9_21_log2FCcorrelated.csv")

### Submit to g:profiler
### Download the csv files
### rename the input csv files

### plot 
Type = "log2FC"
Source = "KEGG"
for (Type in c("NI","log2FC")) {
  Pathways = read.csv(paste("EMC_tpm_global_regression_2021_9_21_",Type,"correlated_gprofiler.csv",sep=""))
  Pathways$source = factor(gsub("GO:","GO_",as.character(Pathways$source)))
  for(Source in levels(Pathways$source)){
    Pathways_sub = Pathways[Pathways$source==Source,]
    Pathways_sub[,1:7]
    Pathways_sub = Pathways_sub[order(Pathways_sub$negative_log10_of_adjusted_p_value),]
    Pathways_sub$term_name = factor(Pathways_sub$term_name,levels=unique(Pathways_sub$term_name))
    Pathways_sub = Pathways_sub[tail(1:nrow(Pathways_sub),30),]
    p1 = ggplot(Pathways_sub, aes(x=negative_log10_of_adjusted_p_value,y=term_name)) +
      geom_col(aes(),fill="grey44")  +    ## geom_point for circle illusion
      geom_text(aes(x=0,label = term_name), hjust = 0)+
      xlab("-log10(padj)")+
      ylab("Pathway term")+
      #scale_x_discrete(position = "top") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
            axis.line.y = element_blank(),
            
            axis.text.x=element_text(colour="black",size=rel(1),vjust = 0.5,hjust = 0.4,angle = 0),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title=element_text(colour="black",size=rel(1)),
            #        legend.position=c(0.8,0.8),
            #        legend.position="bottom",
            legend.position="right",
            legend.background = element_blank(),
            legend.text=element_text(size=rel(1.2)))
    
    pdf(paste("EMC-Figure_S1B_",Type,"_correlated_gprofiler-",Source,".pdf",sep=""),    # create PNG for the heat map        
        width = 4,        # 5 x 300 pixels
        height = 7,
        pointsize = 10 )        # smaller font size
    grid.draw(p1)
    dev.off()
  }
}

################################## section 3 plot
Gene_exp_tpm_global_converted2 = read.csv("EMC_tpm_global_regression_2021_9_13.csv")
length(unique(Gene_exp_tpm_global_converted2$TEinfo_1))
#### immune TFs and TRIM28 and SETDB1
Gene_symbols = c("ENSG00000164342","ENSG00000127666","ENSG00000196664","ENSG00000101916",
                 "ENSG00000172936","ENSG00000107201","ENSG00000115267","ENSG00000088888",
                 "ENSG00000126456","ENSG00000185507","ENSG00000185436","ENSG00000243646",
                 "ENSG00000162434","ENSG00000105397","ENSG00000142166","ENSG00000159110",
                 "ENSG00000115415","ENSG00000170581","ENSG00000213928","ENSG00000130726","ENSG00000143379")

#### For IFN 1 pathways
Gene_exp_tpm_global_converted2_geneSymbol = Gene_exp_tpm_global_converted2[Gene_exp_tpm_global_converted2$type1_IFN != "-",]

#### For key immune regulators
#Gene_exp_tpm_global_converted2_geneSymbol = Gene_exp_tpm_global_converted2[Gene_exp_tpm_global_converted2$TEinfo2 %in% Gene_symbols,]
Gene_exp_tpm_global_converted2_geneSymbol = Gene_exp_tpm_global_converted2[(Gene_exp_tpm_global_converted2$TEinfo2 %in% c(Gene_symbols,"viralLoad","Age")) & !is.na(Gene_exp_tpm_global_converted2$TEinfo_1),]

####
Gene_exp_tpm_global_converted2_geneSymbol$geneName = as.character(Gene_exp_tpm_global_converted2_geneSymbol$geneName)
Gene_exp_tpm_global_converted2_geneSymbol[is.na(Gene_exp_tpm_global_converted2_geneSymbol$geneName),]$geneName = "viralLoad"
rownames(Gene_exp_tpm_global_converted2_geneSymbol) = Gene_exp_tpm_global_converted2_geneSymbol$geneName

Gene_exp_tpm_global_converted2_geneSymbol = data.frame(t(Gene_exp_tpm_global_converted2_geneSymbol[,-1]))
Gene_exp_tpm_global_converted2_geneSymbol$individualID = gsub("_Flu|_log2|_NI","",rownames(Gene_exp_tpm_global_converted2_geneSymbol))
Gene_exp_tpm_global_converted2_geneSymbol$individualID2 = rownames(Gene_exp_tpm_global_converted2_geneSymbol)
Gene_exp_tpm_global_converted2_geneSymbol = merge(Gene_exp_tpm_global_converted2_geneSymbol,Viralload,by="individualID",all.x=T)

Gene_exp_tpm_global_converted2_geneSymbol_kept = Gene_exp_tpm_global_converted2_geneSymbol[grepl("_NI",Gene_exp_tpm_global_converted2_geneSymbol$individualID2),]


### plots
qlist = list()
i = 1
geneName = "IFNAR1"
geneNames = c("TLR3","TICAM1","TLR7","TLR8","MYD88","DDX58","IFIH1","MAVS","IRF3","IRF7", "IFNLR1","IL10RB","IFNAR1","IFNAR2","JAK1","TYK2","STAT1","STAT2","IRF9","SETDB1","TRIM28")
for(geneName in geneNames){
  my.formula = y~x
  p1 = ggplot(Gene_exp_tpm_global_converted2_geneSymbol[grepl("_NI",Gene_exp_tpm_global_converted2_geneSymbol$individualID2),], aes(as.numeric(as.character(ViralLoad_...)),as.numeric(as.character(Gene_exp_tpm_global_converted2_geneSymbol[grepl("_NI",Gene_exp_tpm_global_converted2_geneSymbol$individualID2),which(colnames(Gene_exp_tpm_global_converted2_geneSymbol)==geneName)])))) +
    geom_point(aes(color = Age),size=2,shape=19) + 
    #geom_text(aes(label = sampleID_label)) +
    ylab(paste("Basal ",geneName," exp. (TPM)",sep=""))+
    xlab("Viral loads (%)")+
    #scale_color_gradient2(low="#f7fbff",high = "#08306b",na.value="black")+
    scale_color_gradient2(low="#f7fbff",high = "#08306b")+
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 p.digits = 3,
                 aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
                 parse = TRUE,size = 4.5,) +
    stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
    guides(color=guide_legend(title="Age (yr)"))+
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
      axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
      axis.title=element_text(colour="black",size=rel(1.4)),
      legend.key = element_rect(colour = "transparent", fill = "white"),
      #legend.position="right",
      legend.position = "none",
      legend.background = element_blank(),
      legend.text=element_text(size=rel(1)))
  qlist[[i]] <- ggplotGrob(p1)
  i = i +1
}
pdf(paste("EMC-Figure_S1A_","Immune_regulators_basal-TFs_2021_9_21",".pdf",sep=""),     
    width = 2.5 * 5,        # 5 x 300 pixels
    height = 3*5,
    pointsize = 10)
do.call("grid.arrange",c(qlist,ncol=5))
dev.off()

### log2 FC
### plots
qlist = list()
i = 1
geneName = "IFNAR1"
for(geneName in geneNames){
  my.formula = y~x
  p1 = ggplot(Gene_exp_tpm_global_converted2_geneSymbol[grepl("_log2",Gene_exp_tpm_global_converted2_geneSymbol$individualID2),], aes(as.numeric(as.character(ViralLoad_...)),as.numeric(as.character(Gene_exp_tpm_global_converted2_geneSymbol[grepl("_log2",Gene_exp_tpm_global_converted2_geneSymbol$individualID2),which(colnames(Gene_exp_tpm_global_converted2_geneSymbol)==geneName)])))) +
    geom_point(aes(color = Age),size=2,shape=19) + 
    #geom_text(aes(label = sampleID_label)) +
    ylab(paste("log2FC ",geneName,sep=""))+
    xlab("Viral loads (%)")+
    #scale_color_gradient2(low="#f7fbff",high = "#08306b",na.value="black")+
    scale_color_gradient2(low="#f7fbff",high = "#08306b")+
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black",formula = my.formula)+
    stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 p.digits = 3,
                 aes(label = paste(..rr.label..,stat(p.value.label), sep = "*plain(\",\")~")), 
                 parse = TRUE,size = 4.5,) +
    stat_quantile(quantiles = c(0.05,0.95),col="grey44")+
    guides(color=guide_legend(title="Age (yr)"))+
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.text=element_text(colour="black",size=rel(1.4),angle = 0),
      axis.text.x=element_text(colour="black",vjust=.5,angle = 0),
      axis.title=element_text(colour="black",size=rel(1.4)),
      legend.key = element_rect(colour = "transparent", fill = "white"),
      #legend.position="right",
      legend.position = "none",
      legend.background = element_blank(),
      legend.text=element_text(size=rel(1)))
  qlist[[i]] <- ggplotGrob(p1)
  i = i +1
}
pdf(paste("EMC-Figure_S1A_","Immune_regulators_log2FC-TFs_2021_9_21",".pdf",sep=""),     
    width = 2.5 * 5,        # 5 x 300 pixels
    height = 3* 5,
    pointsize = 10)
do.call("grid.arrange",c(qlist,ncol=5))
dev.off()
