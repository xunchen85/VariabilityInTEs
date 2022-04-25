###     Author:     Xun Chen
###     Date:       2020/5/10
###     Contact:    xunchen85@gmail.com

##### load pacakges
library(MASS)
library(dplyr)
library(splitstackshape)

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

###########################
####### load files
FileName = "EMC_Flu_2020_12_21_Aggregate_noEU37_100.bed.out"
TElist = read.delim2("inputs/TElist.sep.27",sep="",header=T)

## peak summit file
Peaks = read.delim2(paste("inputs/",FileName,sep=""),sep="",header=F)
head(Peaks)
Peaks = cSplit(Peaks,"V5",sep=":")
colnames(Peaks) = c("File","Chr","Start_instance","End_instance","Chr2","Start_summit","End_summit","Summit_info","Score","Strand","Posi","TEcoordinate_1bp","Strand2","Consensus","Ref_position","Consensus_posi","TEsubfamily","TEinstance","TEfamily","TEgroup")
WinSize = 100

## keep summit within Consensus
Peaks_filtered_kept = data.frame(Peaks[!is.na(Peaks$Consensus_posi),])
Peaks_filtered_kept$TEsubfamily = as.character(Peaks_filtered_kept$TEsubfamily)
Peaks_filtered_kept$Consensus = as.character(Peaks_filtered_kept$Consensus)
## 
Peaks_filtered_consensus = data.frame(Peaks_filtered_kept %>% group_by(File,TEsubfamily) %>% dplyr::count(Consensus,sort = TRUE))
Peaks_filtered_consensus$ID = paste(Peaks_filtered_consensus$File,Peaks_filtered_consensus$TEsubfamily)
Peaks_filtered_consensus = Peaks_filtered_consensus[order(-Peaks_filtered_consensus$n),]
Peaks_filtered_consensus = Peaks_filtered_consensus[!duplicated(Peaks_filtered_consensus$ID),]

##### all instances
Peaks_filtered = data.frame(Peaks_filtered_kept)
Peaks_filtered_cluster = data.frame("tmp1"=NA,"tmp2"=NA)
Row = 19
for (Row in 1:nrow(Peaks_filtered_consensus)){
  TEsubfamily = as.character(Peaks_filtered_consensus[Row,]$TEsubfamily)
  TEconsensus = as.character(Peaks_filtered_consensus[Row,]$Consensus)
  Peaks_filtered_tmp = data.frame(Peaks_filtered[as.character(Peaks_filtered$TEsubfamily)==TEsubfamily & as.character(Peaks_filtered$Consensus) == TEconsensus,])
  if (nrow(Peaks_filtered_tmp) == 0) {
    next
  }
  Peaks_filtered_tmp$Consensus_posi = as.numeric(as.character(gsub("_NA","",Peaks_filtered_tmp$Consensus_posi)))
  Peaks_filtered_tmp$Cluster = "No"
  Peaks_filtered_tmp$Win = 0

  k = 1
  GroupID = 1
  Max_posi = max(Peaks_filtered_tmp$Consensus_posi)
  
  ### winSize = 50 bp
  WinSize = 50
  for (k in 1:nrow(Peaks_filtered_tmp)) {
    tmpTable = data.frame(1:(max(Peaks_filtered_tmp$Consensus_posi)))
    tmpTable$SummitCountperWin = 0
    colnames(tmpTable)[1] = "Consensus_posi"
    for (i in (WinSize/2):(Max_posi-(WinSize/2))){
      tmpTable[i,]$SummitCountperWin = nrow(Peaks_filtered_tmp[Peaks_filtered_tmp$Cluster == "No" & Peaks_filtered_tmp$Consensus_posi >= (i -(WinSize/2)) & Peaks_filtered_tmp$Consensus_posi <= (i +(WinSize/2)),])
    }
    if (max(tmpTable$SummitCountperWin)>=5) {
      Max_consensus = tmpTable[tmpTable$SummitCountperWin == max(tmpTable$SummitCountperWin),]$Consensus_posi
      if (length(Max_consensus) > 1) {
        Max_consensus = Max_consensus[as.integer(length(Max_consensus)/2)]
      }
      Peaks_filtered_tmp$Cluster = ifelse(Peaks_filtered_tmp$Cluster == "No" & Peaks_filtered_tmp$Consensus_posi>=(Max_consensus - (WinSize/2)) & Peaks_filtered_tmp$Consensus_posi<=(Max_consensus + (WinSize/2)),paste("Cluster",GroupID,sep="_"),Peaks_filtered_tmp$Cluster)
      Peaks_filtered_tmp$Win = ifelse (Peaks_filtered_tmp$Cluster == paste("Cluster",GroupID,sep="_"),round(median(Peaks_filtered_tmp[Peaks_filtered_tmp$Cluster == paste("Cluster",GroupID,sep="_"),]$Consensus_posi),0),Peaks_filtered_tmp$Win)
      GroupID = GroupID + 1
    } else {
	    break
    }
  }
  Peaks_filtered_tmp$Cluster_50bp = Peaks_filtered_tmp$Cluster
  Peaks_filtered_tmp$Win_50bp = Peaks_filtered_tmp$Win
  
  ### winSize = 100 bp
  WinSize = 100
  Peaks_filtered_tmp$Cluster = "No"
  Peaks_filtered_tmp$Win = 0
  k = 1
  GroupID = 1
  Max_posi = max(Peaks_filtered_tmp$Consensus_posi)
  for (k in 1:nrow(Peaks_filtered_tmp)) {
    tmpTable = data.frame(1:(max(Peaks_filtered_tmp$Consensus_posi)))
    tmpTable$SummitCountperWin = 0
    colnames(tmpTable)[1] = "Consensus_posi"
    for (i in (WinSize/2):(Max_posi-(WinSize/2))){
      tmpTable[i,]$SummitCountperWin = nrow(Peaks_filtered_tmp[Peaks_filtered_tmp$Cluster == "No" & Peaks_filtered_tmp$Consensus_posi >= (i -(WinSize/2)) & Peaks_filtered_tmp$Consensus_posi <= (i +(WinSize/2)),])
    }
    if (max(tmpTable$SummitCountperWin)>=5) {
      Max_consensus = tmpTable[tmpTable$SummitCountperWin == max(tmpTable$SummitCountperWin),]$Consensus_posi
      if (length(Max_consensus) > 1) {
        Max_consensus = Max_consensus[as.integer(length(Max_consensus)/2)]
      }
      Peaks_filtered_tmp$Cluster = ifelse(Peaks_filtered_tmp$Cluster == "No" & Peaks_filtered_tmp$Consensus_posi>=(Max_consensus - (WinSize/2)) & Peaks_filtered_tmp$Consensus_posi<=(Max_consensus + (WinSize/2)),paste("Cluster",GroupID,sep="_"),Peaks_filtered_tmp$Cluster)
      Peaks_filtered_tmp$Win = ifelse (Peaks_filtered_tmp$Cluster == paste("Cluster",GroupID,sep="_"),round(median(Peaks_filtered_tmp[Peaks_filtered_tmp$Cluster == paste("Cluster",GroupID,sep="_"),]$Consensus_posi),0),Peaks_filtered_tmp$Win)
      GroupID = GroupID + 1
    } else {
	    break
    }
  }
  Peaks_filtered_tmp$Cluster_100bp = Peaks_filtered_tmp$Cluster
  Peaks_filtered_tmp$Win_100bp = Peaks_filtered_tmp$Win
  
  ### winSize = 200 bp
  WinSize = 200
  Peaks_filtered_tmp$Cluster = "No"
  Peaks_filtered_tmp$Win = 0
  k = 0
  GroupID = 1
  Max_posi = max(Peaks_filtered_tmp$Consensus_posi)
  for (k in 1: nrow(Peaks_filtered_tmp)) {
    tmpTable = data.frame(1:(max(Peaks_filtered_tmp$Consensus_posi)))
    tmpTable$SummitCountperWin = 0
    colnames(tmpTable)[1] = "Consensus_posi"
    for (i in (WinSize/2):(Max_posi-(WinSize/2))){
      tmpTable[i,]$SummitCountperWin = nrow(Peaks_filtered_tmp[Peaks_filtered_tmp$Cluster == "No" & Peaks_filtered_tmp$Consensus_posi >= (i -(WinSize/2)) & Peaks_filtered_tmp$Consensus_posi <= (i +(WinSize/2)),])
    }
    if (max(tmpTable$SummitCountperWin)>=5) {
      Max_consensus = tmpTable[tmpTable$SummitCountperWin == max(tmpTable$SummitCountperWin),]$Consensus_posi
      if (length(Max_consensus) > 1) {
        Max_consensus = Max_consensus[as.integer(length(Max_consensus)/2)]
      }
      Peaks_filtered_tmp$Cluster = ifelse(Peaks_filtered_tmp$Cluster == "No" & Peaks_filtered_tmp$Consensus_posi>=(Max_consensus - (WinSize/2)) & Peaks_filtered_tmp$Consensus_posi<=(Max_consensus + (WinSize/2)),paste("Cluster",GroupID,sep="_"),Peaks_filtered_tmp$Cluster)
      Peaks_filtered_tmp$Win = ifelse (Peaks_filtered_tmp$Cluster == paste("Cluster",GroupID,sep="_"),round(median(Peaks_filtered_tmp[Peaks_filtered_tmp$Cluster == paste("Cluster",GroupID,sep="_"),]$Consensus_posi),0),Peaks_filtered_tmp$Win)
      GroupID = GroupID + 1
    } else {
      break
    }
  }
  Peaks_filtered_tmp$Cluster_200bp = Peaks_filtered_tmp$Cluster
  Peaks_filtered_tmp$Win_200bp = Peaks_filtered_tmp$Win
  
  if (nrow(Peaks_filtered_cluster) == 1) {
    Peaks_filtered_cluster = Peaks_filtered_tmp
  } else {
    Peaks_filtered_cluster = rbind(Peaks_filtered_cluster,Peaks_filtered_tmp)
  }
  cat("all",Row,nrow(Peaks_filtered_consensus),"\n")
}
write.csv(Peaks_filtered_consensus,file=paste(FileName,"_allConsensus",sep=""))
write.csv(Peaks_filtered_cluster,file=paste(FileName,"_allCluster",sep=""))
