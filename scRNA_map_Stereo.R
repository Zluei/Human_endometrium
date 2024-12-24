library(Seurat)
library(SeuratDisk)
Convert("Spatial_Endo/new/LH12_cellbin_r.h5ad", 
        assay = "RNA", dest = "h5seurat", overwrite = TRUE)

data =  LoadH5Seurat("Spatial_Endo/new/LH12_cellbin_r.h5seurat",meta.data = FALSE, misc = FALSE)
data <- SCTransform(data, assay = "RNA", ncells = 3000, verbose = FALSE)

scRNA<-readRDS("scRNA_Seq/Human_LH_06112024_renamed_v3.rds")
scRNA <- UpdateSeuratObject(scRNA)
scRNA <- SCTransform(scRNA, assay = "RNA", ncells = 3000, verbose = FALSE)

anchors <- FindTransferAnchors(reference = scRNA, query = data, normalization.method = "SCT",
                               npcs = 50,recompute.residuals = FALSE)
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = scRNA$new_label,
                                  prediction.assay = F)

data <- AddMetaData(object = data, metadata = predictions.assay)

data@active.ident = factor(data$predicted.id)
p<-DimPlot(data,group.by = "predicted.id",
           reduction = 'spatial',
           # cols =  c("#E14AEB",#cili
           #           "#EB1D34",#endo
           #           "#FFFC04",#epi
           #           "#70C33B",#lym
           #           "#A30D59",#mac
           #           "#1BE6FF",#nk
           #           "#917C85",#per
           #           "white",#proS
           #           "#001AA6"#str "#6970FF"
           # ),
           # > unique(data$predicted.id)
           # [1] "Stromal"       "NK cell"       "proStromal"   
           # [4] "Epithelial"    "Macrophages"   "Lymphocyte"   
           # [7] "Pericytes"     "Endothelial"   "Ciliated cell"
           #cols = "#606060",
           cols = "lightgrey",#str
           cols.highlight = "#6970FF",
           cells.highlight = WhichCells(data,idents = 'Stromal'),
           sizes.highlight = 0.1 #str
           #sizes.highlight = 0.5
)+ theme( panel.background = element_rect(fill = 'black'))
pdf("../plot/LH9_cellbin_r_str.pdf",width = 20,height = 15)
print(p)
dev.off()
saveRDS(data,"Spatial_Endo/new/LH12_cellbin_r_annoed.rds")
