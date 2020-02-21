library(tidyverse)
library(clusterProfiler)
library(ggthemes)
library(VennDiagram)
library(data.table)
library(pheatmap)
library(viridis)
options(stringsAsFactors = F)
tar<-fread("../ZJH_miRNA/data/all_miRNA_Targets.csv")
colnames(tar)[2]<-"miRNA"
expr<-fread("../ZJH_miRNA/data/miRNA_expression.csv",header = T)%>%as.data.frame()
anno<-data.frame(row.names=colnames(expr)[-c(1,2)],
                 Group = c(rep("ckd",3),rep("ctr",5)))
head(anno)
#heatmap
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
# veen plot
veen_list<-list()
venn.plot<-list()
overlap<-list()
overlap_genes<-list()
for (i in seq(length(miRWalk_list1))) {
        overlap[[i]]<-intersect(intersect(miRWalk_list1[[i]]$target,TargetScan_list[[i]]$target),miRDB_list1[[i]]$target)
        names(overlap)[i]<-names(miRWalk_list1)[i]
 veen_list[[i]]<-list(miRWalk=miRWalk_list1[[i]]$target,
                      TargetScan=TargetScan_list[[i]]$target,
                      miRDB=miRDB_list1[[i]]$target)
 venn.plot[[i]] <- venn.diagram( veen_list[[i]] , NULL, 
                           fill=c("darkmagenta", "darkblue","yellow"), 
                           alpha=c(0.5,0.5,0.5), cex = 2, 
                           cat.fontface=4,  
                           main=paste0("overlap of targets ( ",names(miRWalk_list1)[i]," )"))
 overlap_genes[[i]]<-data.frame(miRNA=rep(names(overlap)[i],length(overlap[[i]])),
                                target=overlap[[i]])
 grid.newpage()
 pdf(paste0("../ZJH_miRNA/results/veen_",names(miRWalk_list1)[i],".pdf"))
 grid.draw(venn.plot[[i]])
 dev.off()
 png(paste0("../ZJH_miRNA/results/veen_",names(miRWalk_list1)[i],".png"))          
 grid.draw(venn.plot[[i]])
 dev.off()
}
overlap_data<-bind_rows(overlap_genes)
overlap_res<-data.frame(miRNA=names(overlap),overlap_num=sapply(overlap, length))
write.csv(overlap_res,"../ZJH_miRNA/results/overlap_3database_summary.csv",row.names = F)
write.csv(overlap_data,"../ZJH_miRNA/results/overlap_3database_genes.csv",row.names = F)
#enrichment
down_miRNA<-tar_data%>%subset(.,log2FC<0)
up_miRNA<-tar_data%>%subset(.,log2FC>0)
down_miRNA_list<-levels(factor(down_miRNA$miRNA))
up_miRNA_list<-levels(factor(up_miRNA$miRNA))
up_target<-overlap[which(names(overlap)%in%up_miRNA_list)]%>%unlist()
down_target<-overlap[-which(names(overlap)%in%up_miRNA_list)]%>%unlist()
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
list_up=select(org.Hs.eg.db,keys = up_target,column = "ENTREZID", keytype = "SYMBOL",multiVals = "first")
list_down=select(org.Hs.eg.db,keys = down_target,column = "ENTREZID", keytype = "SYMBOL",multiVals = "first")
geneUp<-list_up[,2]
geneDown<-list_down[,2]
#organism:'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
ekkUp <- enrichKEGG(gene = geneUp,organism = 'hsa',pvalueCutoff = 0.05)
ekkDown <- enrichKEGG(gene = geneDown,organism = 'hsa',pvalueCutoff = 0.05)
write.csv(as.matrix(ekkUp@result), file="../ZJH_miRNA/results/KEGG_Up.csv",row.names = F)
write.table(as.matrix(ekkDown@result), file="../ZJH_miRNA/results/KEGG_Down.csv",row.names = F,sep = "\t",quote = F)

egoUp1 <- enrichGO(gene=geneUp,org.Hs.eg.db,ont="BP",pvalueCutoff=0.05,readable=TRUE)
egoDown1 <- enrichGO(gene=geneDown,org.Hs.eg.db,ont="BP",pvalueCutoff=0.05,readable=TRUE)
rUp1<-rep("BP",length(row_number(egoUp1@result$Description)))
rUp1<-as.data.frame(rUp1)
colnames(rUp1)<-("Class")
egoUp1<-bind_cols(egoUp1@result,rUp1)

rDown1<-rep("BP",length(row_number(egoDown1@result$Description)))
rDown1<-as.data.frame(rDown1)
colnames(rDown1)<-("Class")
egoDown1<-bind_cols(egoDown1@result,rDown1)

egoUp2 <- enrichGO(gene=geneUp,org.Hs.eg.db,ont="MF",pvalueCutoff=0.05,readable=TRUE)
egoDown2<- enrichGO(gene=geneDown,org.Hs.eg.db,ont="MF",pvalueCutoff=0.05,readable=TRUE)

rUp2<-rep("MF",length(row_number(egoUp2@result$Description)))
rUp2<-as.data.frame(rUp2)
colnames(rUp2)<-("Class")
egoUp2<-bind_cols(egoUp2@result,rUp2)

rDown2<-rep("MF",length(row_number(egoDown2@result$Description)))
rDown2<-as.data.frame(rDown2)
colnames(rDown2)<-("Class")
egoDown2<-bind_cols(egoDown2@result,rDown2)

egoUp3 <- enrichGO(gene=geneUp,org.Hs.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)
egoDown3<- enrichGO(gene=geneDown,org.Hs.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)

rUp3<-rep("CC",length(row_number(egoUp3@result$Description)))
rUp3<-as.data.frame(rUp3)
colnames(rUp3)<-("Class")
egoUp3<-bind_cols(egoUp3@result,rUp3)

rDown3<-rep("CC",length(row_number(egoDown3@result$Description)))
rDown3<-as.data.frame(rDown3)
colnames(rDown3)<-("Class")
egoDown3<-bind_cols(egoDown3@result,rDown3)

egoUp<-bind_rows(egoUp1,egoUp2,egoUp3)
egoDown<-bind_rows(egoDown1,egoDown2,egoDown3)
write.csv(as.matrix(egoUp), file="../ZJH_miRNA/results/GO_up.csv",row.names = F)
write.csv(as.matrix(egoDown), file="../ZJH_miRNA/results/GO_down.csv",row.names = F)
##plot GO
library(ggpubr)
library(dplyr)
library(ggthemes)
df<-read.csv("../ZJH_miRNA/results/plot_GO.csv")
BP<-subset(df,Class=="BP")
MF<-subset(df,Class=="MF")
CC<-subset(df,Class=="CC")
df2<-read.csv("../ZJH_miRNA/results/plot_GO.csv")
p<-ggbarplot(df2, x="Description", y="Count", fill = "Class", color = "white",
             palette = "ucscgb",
             sort.val = "asc",
             sort.by.groups=TRUE,
             x.text.angle=45, ylab = "GeneNumber",legend.title="GO_classes")
p+theme_hc()+theme(axis.text = element_text(size = 12),
                   axis.title  = element_text(size = 14))
pp<- ggbarplot(df2, x="Description", y="Count", 
               fill = "Class", color = "white", 
               width = 1,palette = c("blue","yellow","red"), 
               sort.val = "asc", sort.by.groups = T, 
               x.text.angle=90, ylab = "GeneNumber", 
               xlab = "Go_term", 
               legend.title="Class", 
               rotate=TRUE, ggtheme = theme_tufte())
#plotKEGG
df<-read.csv("../ZJH_miRNA/results/plot_KEGG.csv",header = T)
library(ggplot2)
library(ggthemes)
library(tidyverse)
library(Cairo)
library(ggpubr)
head(df)
df$Group<-factor(df$Group,levels = c("Up","Down"))
df$Term<-reorder(df$Term,df$GeneRatio)
p<-ggplot(df,aes(GeneRatio,Term))
pp<-p+geom_point(aes(size=GeneNumber,color=negLOG2.qvalue.,shape=Group))
p1<-pp+ scale_color_gradient(low = "blue",high = "red")+
        theme(axis.text = element_text(colour = "black",size = 12),
              legend.text = element_text(colour = "black",size = 12),
              panel.grid.major = element_line(colour = "gray88",linetype = "dashed"),
              panel.background=element_blank(),
              axis.ticks = element_line(size = 1, colour = "black"),
              axis.line = element_line(size = 1, colour = "black"),
              legend.position = "right",
              panel.border = element_rect(linetype = "solid", fill = NA,size = 2),
              axis.title = element_text(size = 12)
              ,title = element_text(size = 12))
p1+labs(x="Ratio")+theme_base()+
        theme(axis.text = element_text(colour = "black",size = 12),
              legend.text = element_text(colour = "black",size = 12),
              panel.grid.major = element_line(colour = "gray88",linetype = "dashed"),
              panel.background=element_blank(),
              axis.ticks = element_line(size = 1, colour = "black"),
              axis.line = element_line(size = 1, colour = "black"),
              legend.position = "right",
              panel.border = element_rect(linetype = "solid", fill = NA,size = 2),
              axis.title = element_text(size = 12)
              ,title = element_text(size = 12))

