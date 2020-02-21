library(tidyverse)
library(BiocManager)
library(clusterProfiler)
library(ggthemes)
library(VennDiagram)
library(data.table)
library(pheatmap)
library(viridis)
tar<-fread("../ZJH_miRNA/data/all_miRNA_Targets.csv")
colnames(tar)[2]<-"miRNA"
expr<-fread("../ZJH_miRNA/data/miRNA_expression.csv",header = T)%>%as.data.frame()
anno<-data.frame(row.names=colnames(expr)[-c(1,2)],
                 Group = c(rep("ckd",3),rep("ctr",5)))
head(anno)
mat<-expr[,3:10]
rownames(mat)<-expr$miRNA
pheatmap(log(mat+1),
         annotation_col = anno,
         fontsize = 12,
         annotation_colors = list(Group=c(ckd="blue",ctr="#D95F02")),
         color = viridis(50),
         border_color = NA,scale = "row")
tar_data<-merge(tar,expr[,1:2],by="miRNA")
head(tar_data)
tar_list<-split(tar_data,f = tar_data$dataset,drop = F)
names(tar_list)
miRDB<-tar_list[[1]]
miRDB_list<-split(miRDB,f=miRDB$miRNA,drop = F)

miRWalk<-tar_list[[2]]
miRWalk_list<-split(miRWalk,f=miRWalk$miRNA,drop = F)
names(miRWalk_list)
TargetScan<-tar_list[[3]]
TargetScan_list<-split(TargetScan,f=TargetScan$miRNA,drop = F)
names(TargetScan_list)
setdiff(names(miRWalk_list),names(TargetScan_list))
miRWalk_list1<-miRWalk_list[-which(names(miRWalk_list)%in%setdiff(names(miRWalk_list),names(TargetScan_list)))]
miRDB_list1<-miRDB_list[-which(names(miRDB_list)%in%setdiff(names(miRWalk_list),names(TargetScan_list)))]
names(miRWalk_list1)
names(TargetScan_list)
names(miRDB_list1)

diff<-list()
veen_list<-list()
venn.plot<-list()
for (i in seq(length(miRWalk_list1))) {
 veen_list[[i]]<-list(miRWalk=miRWalk_list1[[i]]$target,
                      TargetScan=TargetScan_list[[i]]$target,
                      miRDB=miRDB_list1[[i]]$target)
 venn.plot[[i]] <- venn.diagram( veen_list[[i]] , NULL, 
                           fill=c("darkmagenta", "darkblue","yellow"), 
                           alpha=c(0.5,0.5,0.5), cex = 2, 
                           cat.fontface=4,  
                           main=paste0("overlap of target ( ",names(miRWalk_list1)[i]," )"))
 
 grid.draw(venn.plot[[i]]) 
 }

