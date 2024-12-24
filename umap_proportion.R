setwd("/mnt/My_disk/Zhang_leilei/Human_endo/scRNA_Seq/")
scRNA<-readRDS("Human_LH_06112024_renamed_v3.rds")
#umap####
col = c("#8DD469",#cili
        "#702189",#endo
        "#FEA800",#epi
        "#007E67",#lym
        "#FEA88D",#mac
        "#9E7257",#nk
        "#D19FF6",#per
        "#79DACE",#proS
        "#F8887D"#str
)
col =  c("#E14AEB",#cili
         "#EB1D34",#endo
         "#FFFC04",#epi
         "#70C33B",#lym
         "#A30D59",#mac
         "#1BE6FF",#nk
         "#917C85",#per
         "white",#proS
         "#001AA6"#str "#6970FF"
)
DimPlot(scRNA,group.by = "new_label",label = F,cols = col)+
  theme( panel.background = element_rect(fill = 'black'))
DimPlot(scRNA,group.by = "new_label",label = T,split.by = "sub_pup")
#细胞占比####
library(Seurat)
library(dplyr)
library(cowplot)
scRNA@active.ident<-factor(scRNA$new_label,levels = unique(scRNA$new_label))
meta<-scRNA@meta.data
LH2<-meta[meta$sub_pup == "LH2",]
LH7<-meta[meta$sub_pup == "LH7",]
LH9<-meta[meta$sub_pup == "LH9",]
LH12<-meta[meta$sub_pup == "LH12",]
data.plot<-data.frame(LH2 = as.data.frame(table(LH2[,"new_label"]))[,2] / sum(as.data.frame(table(LH2[,"new_label"]))[,2]),
                      LH7 = as.data.frame(table(LH7[,"new_label"]))[,2] / sum(as.data.frame(table(LH7[,"new_label"]))[,2]),
                      LH9 = as.data.frame(table(LH9[,"new_label"]))[,2] / sum(as.data.frame(table(LH9[,"new_label"]))[,2]),
                      LH12 = as.data.frame(table(LH12[,"new_label"]))[,2] / sum(as.data.frame(table(LH12[,"new_label"]))[,2]),
                      row.names = as.data.frame(table(LH7[,"new_label"]))[,1])
data.plot<-data.plot * 100
data.plot["group",] = colnames(data.plot)
data.plot<-as.data.frame(t(data.plot))
data.plot$group<-factor(data.plot$group,levels = data.plot$group)

p1<-ggplot(as.data.frame(data.plot),aes(group,as.numeric(data.plot$proStromal),group = 1))+
  geom_point(color = "#79DACE",size = 2.5)+
  geom_line(color = "#79DACE",cex = 1.2)+
  scale_y_continuous(limits = c(0,12),breaks = c(0,4,8))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "Percentage(%)",x = "Group")+
  annotate("text", x=1.2, y=12, label="proStromal")

p2<-ggplot(as.data.frame(data.plot),aes(group,as.numeric(data.plot$Stromal),group = 1))+
  geom_point(color = "#F8887D",size = 2.5)+
  geom_line(color = "#F8887D",cex = 1.2)+
  scale_y_continuous(limits = c(0,100),breaks = c(0,40,80))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "Percentage(%)",x = "Group")+
  annotate("text", x=1, y=100, label="Stromal")
p3<-ggplot(as.data.frame(data.plot),aes(group,as.numeric(data.plot$Lymphocyte),group = 1))+
  geom_point(color = "#007E67",size = 2.5)+
  geom_line(color = "#007E67",cex = 1.2)+
  scale_y_continuous(limits = c(0,12),breaks = c(0,4,8))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "Percentage(%)",x = "Group")+
  annotate("text", x=1.5, y=12, label="Lymphocyte")

p4<-ggplot(as.data.frame(data.plot),aes(group,as.numeric(data.plot$Epithelial),group = 1))+
  geom_point(color = "#FEA800",size = 2.5)+
  geom_line(color = "#FEA800",cex = 1.2)+
  scale_y_continuous(limits = c(0,12),breaks = c(0,4,8))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "Percentage(%)",x = "Group")+
  annotate("text", x=1, y=12, label="Epithelial")

p5<-ggplot(as.data.frame(data.plot),aes(group,as.numeric(data.plot$Endothelial),group = 1))+
  geom_point(color = "#702189",size = 2.5)+
  geom_line(color = "#702189",cex = 1.2)+
  scale_y_continuous(limits = c(0,12),breaks = c(0,4,8))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "Percentage(%)",x = "Group")+
  annotate("text", x=1, y=12, label="Endothelial")
p6<-ggplot(as.data.frame(data.plot),aes(group,as.numeric(data.plot$`NK cell`),group = 1))+
  geom_point(color = "#9E7257",size = 2.5)+
  geom_line(color = "#9E7257",cex = 1.2)+
  scale_y_continuous(limits = c(0,5),breaks = c(0,2,4))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "Percentage(%)",x = "Group")+
  annotate("text", x=1, y=5, label="NK cell")
p7<-ggplot(as.data.frame(data.plot),aes(group,as.numeric(data.plot$Macrophages),group = 1))+
  geom_point(color = "#FEA88D",size = 2.5)+
  geom_line(color = "#FEA88D",cex = 1.2)+
  scale_y_continuous(limits = c(0,5),breaks = c(0,2,4))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "Percentage(%)",x = "Group")+
  annotate("text", x=1, y=5, label="Macrophages")
p8<-ggplot(as.data.frame(data.plot),aes(group,as.numeric(data.plot$Pericytes),group = 1))+
  geom_point(color = "#D19FF6",size = 2.5)+
  geom_line(color = "#D19FF6",cex = 1.2)+
  scale_y_continuous(limits = c(0,5),breaks = c(0,2,4))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "Percentage(%)",x = "Group")+
  annotate("text", x=1, y=5, label="Pericytes")

p9<-ggplot(as.data.frame(data.plot),aes(group,as.numeric(data.plot$`Ciliated cell`),group = 1))+
  geom_point(color = "#8DD469",size = 2.5)+
  geom_line(color = "#8DD469",cex = 1.2)+
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.4,0.8))+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  labs(y = "Percentage(%)",x = "Group")+
  annotate("text", x=1, y=1, label="Ciliated cell")
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,nrow = 2) 

library(reshape2)
meta<-scRNA@meta.data
data.plot<-melt(table(meta[,c("sub_pup","new_label")]))
data.plot$sub_pup<-factor(data.plot$sub_pup,levels = c("LH2","LH7","LH9","LH12"))
data.plot$new_label<-factor(data.plot$new_label,levels = data.plot[data.plot$sub_pup == "LH2","new_label"][order(data.plot[data.plot$sub_pup == "LH2","value"],decreasing = T)])

#data.plot$new_label<-factor(data.plot$new_label,levels = unique(scRNA$new_label))
ggplot(data.plot,mapping = aes(sub_pup,value,fill=new_label))+
  geom_bar(stat='identity',
           position='fill',width = 0.85) +
  labs(x = 'Group',y = 'frequnency') +
  theme(axis.title =element_text(size = 16),
        axis.text =element_text(size = 14, color = 'black'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = c("#F8887D",#str
                               "#79DACE",#proS
                               "#007E67",#lym
                               "#702189",#endo
                               "#D19FF6",#per
                               "#FEA88D",#mac
                               "#FEA800",#epi
                               "#8DD469",#cili
                               "#9E7257"#nk
  ))
