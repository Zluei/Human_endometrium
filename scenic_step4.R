#SCENIC####

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
# If your bioconductor version is previous to 4.0, see the section bellow

## Required
# BiocManager::install(c("AUCell", "RcisTarget"))
# devtools::install_github("aertslab/RcisTarget")
# BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# # To score the network on cells (i.e. run AUCell):
# BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# # For various visualizations and perform t-SNEs:
# 
# # To support paralell execution (not available in Windows):
# BiocManager::install(c("doMC", "doRNG"))
# # To export/visualize in http://scope.aertslab.org
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE,force = T)
# 
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("aertslab/SCENIC",force = T) 
library(SCENIC)
packageVersion("SCENIC")
library(SCopeLoomR)
library(ggplot2)
library(dplyr)

setwd("/mnt/My_disk/Zhang_leilei/Human_endo/scRNA_Seq/scenic/str/0722/")

# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(Seurat)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(pheatmap)
.regularise_df <- function(df, drop_single_values = TRUE) {
  if (ncol(df) == 0) df[['name']] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0)
      warning(paste('Dropping single category variables:'),
              paste(colnames(df)[k_singular], collapse=', '))
    df <- df[, !k_singular, drop=F]
    if (ncol(df) == 0) df[['name']] <- rownames(df)
  }
  return(df)
}
#str
setwd("str/")
str<-str_scRNA
str<-readRDS("sub_str.rds")
loom <- open_loom("all/str_all_count_SCENIC.loom")
DefaultAssay(str) <- "RNA"

obs <- .regularise_df(str@meta.data, drop_single_values = T)
str[["RNA"]] <- as(object = str[["RNA"]], Class = "Assay")
var <- .regularise_df(GetAssay(str, assay = "RNA")@meta.features, drop_single_values = T)
obsm <- NULL
reductions <- names(str@reductions)
if (length(reductions) > 0) {
  obsm <- sapply(
    reductions,
    function(name) as.matrix(Embeddings(str, reduction=name)),
    simplify = FALSE
  )
  names(obsm) <- paste0('X_', tolower(names(str@reductions)))
}


### Initialize settings
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name="RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- obsm
#cellClusters <- data.frame(clusters = obs[,c("new_label")],row.names = rownames(obs))
#cellClusters$clusters<-factor(cellClusters$clusters,levels = c("str_3","str_2","str_0","str_1"))

cellClusters <- data.frame(clusters = obs[,c("label")],row.names = rownames(obs))

close_loom(loom)

##

# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,1]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) {rowMeans(getAUC(regulonAUC)[,cells])})
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
# plot:
# options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
# hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
#                                    #row_names_gp=grid::gpar(fontsize=6),
#                                    show_row_names = F)) # row font size
hm<-ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, 
                            name="Regulon activity",show_row_names = T,cluster_columns = F)
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 30))###每一组细胞内活跃的TF

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC), "clusters"])
write.csv(rss,"all/Rugulon_Specificity_Score.csv")
## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
for (i in unique(cellClusters$clusters)) {
  p<-plotRSS_oneSet(rss, setName = i) # cluster ID
  pdf(paste0("all/rss/",i,"_rss.pdf"),width = 5,height = 6)
  print(p)
  dev.off()
}
#plotRSS_oneSet(rss, setName = "Mac_0") # cluster ID

par(mfrow=c(1,1))
#regulonsToPlot <- "Atf4(+)"
for (i in names(regulons)) {
  pdf(paste0("all/regulon_umap/",i,"_umap.pdf"),width = 5,height = 5)
  AUCell_plotTSNE(embeddings[["X_umap"]], exprMat_log, regulonAUC[i,], plots=c("AUC"), cex = .3,sub="X_umap")
  dev.off()
}
#AUCell_plotTSNE(embeddings[["X_umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .3,sub="X_umap")

# aucellApp <- AUCell_createViewerApp(auc=regulonAUC,
#                                     thresholds=regulonAucThresholds,
#                                     tSNE=embeddings[["X_umap"]], 
#                                     exprMat=exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)

# This function will be included in the next version of AUCell
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}
a<-as.data.frame(regulonAucThresholds)
regulonAucThresholds<-rownames(a)
names(regulonAucThresholds)<-a$regulonAucThresholds
binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)


thresholds <- regulonAucThresholds[intersect(names(regulonAucThresholds), rownames(regulonAUC))]
regulonsCells <- setNames(lapply(names(thresholds), 
                                 function(x) {
                                   trh <- thresholds[x]
                                   names(which(getAUC(regulonAUC)[x,]>trh))
                                 }),names(thresholds))

regulonActivity <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
class(binaryRegulonActivity) <- "matrix"  
cellClusters<-data.frame(cluster = cellClusters[colnames(binaryRegulonActivity),],
                         row.names = colnames(binaryRegulonActivity))
binaryRegulonActivity<-binaryRegulonActivity[,rownames(cellClusters)[order(cellClusters$cluster)]]

#col = c("#fccccb","#bdb5e1","#b0d992","#f9d580")
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
names(col) <- unique(cellClusters$cluster)
column_ann <- HeatmapAnnotation(
  celltype = cellClusters[order(cellClusters$cluster),],
  col = list(celltype = col)
)

ComplexHeatmap::Heatmap(binaryRegulonActivity, name="Binarized activity", 
                        col = c("white", "black"),
                        cluster_rows = TRUE, cluster_columns=F,
                        show_column_names = FALSE,
                        row_names_gp=grid::gpar(fontsize=6),use_raster = F,
                        top_annotation = column_ann) # row font size

