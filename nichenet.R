#nichenet#####
setwd("/mnt/My_disk/Zhang_leilei/Human_endo/scRNA_Seq/")
#devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)

data<-readRDS("Human_LH_04232024_renamed_v2.rds")
data@active.assay
data@active.ident<-factor(data$new_label,levels = unique(data$new_label))
#先验配体-受体加权网络模型#
ligand_target_matrix = readRDS("nichenet/ligand_target_matrix.rds")
lr_network = readRDS("nichenet/lr_network.rds")
weighted_networks = readRDS("nichenet/weighted_networks.rds") # interactions and their weights in the ligand-receptor + signaling network
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))# interactions and their weights in the gene regulatory network

# LH2<-subset(data,subset = sub_pup == "LH2")
# LH7<-subset(data,subset = sub_pup == "LH7")
# LH9<-subset(data,subset = sub_pup == "LH9")
# LH12<-subset(data,subset = sub_pup == "LH12")
# data<-LH12
#定义受体细胞
receiver = "Stromal"
#筛选受体
expressed_genes_receiver <- get_expressed_genes(receiver, data, pct = 0.05)
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
#根据受体获取潜在配体
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

#定义配体细胞
sender_celltypes <- c("NK cell")
# Use lapply to get the expressed genes of every sender cell type separately here
#获取配体细胞所表达的配体
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, data, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# ave <- AverageExpression(data) %>% as.data.frame()
# geneset_oi<-names(rowSums(ave)[order(rowSums(ave),decreasing = T)])[1:2000]
# geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
condition_oi <-  c("LH12")
condition_reference <- c("LH7","LH9")

seurat_obj_receiver <- subset(data, idents = receiver)
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "sub_pup",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 1) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]


background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) 

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") 
ligand_aupr_matrix<-data.frame(aupr_corrected = ligand_aupr_matrix$aupr_corrected,row.names = rownames(ligand_aupr_matrix)) %>% arrange(aupr_corrected)

vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr
# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target
# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 
vis_ligand_receptor_network<-vis_ligand_receptor_network[,rev(best_upstream_ligands)]
p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor
best_upstream_ligands_all %in% rownames(data) %>% table()

# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(data, new_label %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot

celltype_order <- levels(Idents(data)) 

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype,
  seurat_obj = data,
  condition_colname = "sub_pup",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "new_label",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands
)

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>%
  column_to_rownames("gene")

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), ])
rownames(vis_ligand_lfc)<-rev(best_upstream_ligands)
colnames(vis_ligand_lfc)<-sender_celltypes
p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc),
                                        high_color = "red",
                                        legend_title = "LFC")

p_lfc

figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  # p_lfc + theme(legend.position = "none",
  #               axis.title.y = element_blank()),
  p_ligand_receptor + theme(legend.position = "none",
                            axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+20, ncol(vis_ligand_lfc)+20,ncol(vis_ligand_target), ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  #ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_receptor)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(2, 2, 2))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot