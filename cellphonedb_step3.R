##heatmap####
library(pheatmap)
library(reshape)
library("RColorBrewer")
library("viridis")
#library("dpler")
library("ape") 

cpdb_data<-read.table("LH7/plot/count_network.txt",sep = "\t",header = T)
sum<-aggregate(cpdb_data$count, by=list(type=cpdb_data$SOURCE),sum)
sx<-sum$type[order(sum$x,decreasing = T)]
heatmap_data <- cast(cpdb_data,SOURCE~TARGET)
row.names(heatmap_data)<-heatmap_data$SOURCE
heatmap_data<-heatmap_data[,-1]
heatmap_data<-heatmap_data[sx,sx]
ys <- c(seq(1,250,by=1))
my_col = c(colorRampPalette(colors = c("#FFECE5","#FFB6B2"))(length(ys)/2),
           colorRampPalette(colors = c("#FFB6B2","#EC4C53"))(length(ys)/2))
pheatmap(heatmap_data,border_color = NA,cluster_cols = F,cluster_rows = F,
         color = my_col,show_rownames = T,show_colnames = T,
         legend_breaks=seq(0,250,50),
         breaks=c(seq(1,250,by=1))
)

sx<-row.names(heatmap_data)
cpdb_data<-read.table("LH12/plot/count_network.txt",sep = "\t",header = T)
#da<-cpdb_data[order(cpdb_data$count,decreasing = T),]
heatmap_data <- cast(cpdb_data,SOURCE~TARGET)
row.names(heatmap_data)<-heatmap_data$SOURCE
heatmap_data<-heatmap_data[,-1]
heatmap_data<-heatmap_data[sx,sx]
pheatmap(heatmap_data,border_color = NA,cluster_cols = F,cluster_rows = F,
         color = my_col,
         legend_breaks=seq(0,250,50),
         breaks=c(seq(1,250,by=1))
)

##circos plot####

#install.packages("circlize")
library(Seurat)
library(circlize)
library(dplyr)
library(ggplot2)
D4<-read.table("~/project/ZhaoHui/D4_D8/cellphonedb/D4/plot/count_network.txt",sep = "\t",header = T)
cir_data<-D4
sx<-c("Mac_0", "Mono_1", "NK_2", "B_3", "trMac_4", "ProNK_5", "ProMac_6", "mNK_7", "Neutro_8", "cDC1_9", "cDC2_10", "pro_B_11", "pDC_12", "T_13", "NKT_14", "Mast_15", "B_16")
#sum<-aggregate(cir_data$count, by=list(type=cir_data$SOURCE),sum)
#sx<-sum$type[order(sum$x,decreasing = T)]
data<-readRDS("~/project/ZhaoHui/D4_D8/dim15_resolution0.35_cluster18.rds")
p<-DimPlot(data)
a<-ggplot_build(p)$data[[1]] %>% .[,c(1,5)] %>% duplicated.data.frame(.)
col<-ggplot_build(p)$data[[1]] %>% .[,c(1,5)]
col<-unique(col)
col<-col[order(col$group),]
#col<-col[c(15,7,1,11,5,6,8,2,10,13,3,14,12,9,17,4,16),]
color<-factor(col$colour,levels = col$colour)
chordDiagram(cir_data,
             annotationTrack = c("grid"),
             grid.col = col$colour,
             transparency = 0,
             order=sx,
             link.lwd = 0.5,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cir_data))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

sx<-sum$type[order(sum$x,decreasing = T)]
D8<-read.table("~/project/ZhaoHui/D4_D8/cellphonedb/D8/plot/count_network.txt",sep = "\t",header = T)
cir_data<-D8
chordDiagram(cir_data,
             annotationTrack = c("grid"),
             grid.col = col$colour,
             transparency = 0,
             order=sx,
             link.lwd = 0.5,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cir_data))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)