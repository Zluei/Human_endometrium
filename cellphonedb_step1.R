setwd("/mnt/My_disk/Zhang_leilei/Human_endo/scRNA_Seq/cpdb/v2")
scRNA<-readRDS("../../Human_LH_06112024_renamed_v3.rds")
LH2<-subset(scRNA,subset = sub_pup == "LH2")
count<-GetAssayData(LH2,slot = "counts")
write.table(count,"LH2/LH2_count.txt",sep = "\t",row.names = T,col.names = T,quote = F)
meta<-data.frame(cell = rownames(LH2@meta.data),label = LH2@meta.data$new_label)
write.table(meta,"LH2/LH2_meta.txt",sep = "\t",row.names = F,col.names = T,quote = F)

LH7<-subset(scRNA,subset = sub_pup == "LH7")
count<-GetAssayData(LH7,slot = "counts")
write.table(count,"LH7/LH7_count.txt",sep = "\t",row.names = T,col.names = T,quote = F)
meta<-data.frame(cell = rownames(LH7@meta.data),label = LH7@meta.data$new_label)
write.table(meta,"LH7/LH7_meta.txt",sep = "\t",row.names = F,col.names = T,quote = F)

LH9<-subset(scRNA,subset = sub_pup == "LH9")
count<-GetAssayData(LH9,slot = "counts")
write.table(count,"LH9/LH9_count.txt",sep = "\t",row.names = T,col.names = T,quote = F)
meta<-data.frame(cell = rownames(LH9@meta.data),label = LH9@meta.data$new_label)
write.table(meta,"LH9/LH9_meta.txt",sep = "\t",row.names = F,col.names = T,quote = F)

LH12<-subset(scRNA,subset = sub_pup == "LH12")
count<-GetAssayData(LH12,slot = "counts")
write.table(count,"LH12/LH12_count.txt",sep = "\t",row.names = T,col.names = T,quote = F)
meta<-data.frame(cell = rownames(LH12@meta.data),label = LH12@meta.data$new_label)
write.table(meta,"LH12/LH12_meta.txt",sep = "\t",row.names = F,col.names = T,quote = F)

