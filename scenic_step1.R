#SCENIC####
#数据准备
library(Seurat)
str_scRNA<-readRDS("sub_str.rds")
Idents(str_scRNA)<-"label"
#DimPlot(str_scRNA,group.by = "label")
#注意矩阵转置
count<-t(as.matrix(str_scRNA@assays$RNA@layers$counts))
rownames(count)<-colnames(str_scRNA)
colnames(count)<-row.names(str_scRNA)
write.csv(count,file = "scenic/str/str_all_count.csv")