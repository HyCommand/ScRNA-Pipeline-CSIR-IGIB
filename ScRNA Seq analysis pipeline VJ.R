#Libraries ----
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(celldex)
library(SingleR)
library(readxl)
library(pheatmap)
library(here)
library(Matrix)
library(monocle3)
library(DESeq2)
library(scCustomize)
library(remotes)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(ggforce)
library(org.Hs.eg.db)
library(CytoTRACE)
library(tidyr)
library(writexl)
set.seed(1234)
#Data Loading ----
matrix_WT <- ReadMtx(
  mtx = '',
  cells = '',
  features = ''
)

matrix_KO <- ReadMtx(mtx = '',
                     cells = '',
                     features = '')

#Seurat object creation ----

Ctrl_unfiltered <- CreateSeuratObject(counts = matrix_WT, project = "Dist", min.cells = 3, min.features = 200)
KO_unfiltered <- CreateSeuratObject(counts = matrix_KO, project = "Dist", min.cells = 3, min.features = 200)

#QC (Mito genes) ----

Ctrl_unfiltered[["percentMT"]] <- PercentageFeatureSet(Ctrl_unfiltered , pattern = "^MT-")
KO_unfiltered[["percentMT"]] <- PercentageFeatureSet(KO_unfiltered, pattern = "^MT-")

VlnPlot(Ctrl_unfiltered, features = c("nFeature_RNA", "percentMT")) 
Ctrl <- subset(Ctrl_unfiltered, subset = c(nFeature_RNA >2000 & percentMT < 7))
VlnPlot(KO_unfiltered, features = c("nFeature_RNA", "percentMT"))
KO <- subset(KO_unfiltered, subset = c(nFeature_RNA >2000 & percentMT < 7))

#Sub Sampling ----
s
ncell <- min(ncol(Ctrl), ncol(KO))
Ctrl_subset <- subset(Ctrl, cells = sample(colnames(Ctrl), ncell))
KO_subset <- subset(KO, cells = sample(colnames(KO), ncell))

#Normalization ----

Ctrl_work <- NormalizeData(Ctrl_subset)
Ctrl_work <- FindVariableFeatures(Ctrl_work)
Ctrl_work <- ScaleData (Ctrl_work)

Ctrl_work$condition <- "Control"

KO_work <- NormalizeData(KO_subset)
KO_work <- FindVariableFeatures(KO_work)
KO_work <- ScaleData(KO_work)

KO_work$condition <- "Knockout"

#Data Integration ----

data_list <- list(Ctrl_work,KO_work)
features_int <- SelectIntegrationFeatures(object.list = data_list)
anchor <- FindIntegrationAnchors(object.list = data_list, anchor.features = features_int)
Int <- IntegrateData(anchorset = anchor)

DefaultAssay(Int) <- "integrated"

Int <- ScaleData(Int)
Int <- RunPCA(Int , npcs = 50)
ElbowPlot(Int, reduction = 'pca', ndims = 50)
Int <- FindNeighbors(Int, reduction = 'pca', dims = 1:30)
Int <- FindClusters(Int,resolution = c(0.5))
Int <- RunUMAP(Int, reduction = 'pca', dims = 1:30,seed.use = 1234)
DimPlot(Int, reduction = 'umap', label = T)
DimPlot(Int, reduction = 'umap', group.by = "condition", label = T)
DimPlot(Int, reduction = 'umap', split.by = "condition", label = T, label.size = 5, repel = T)
label_data <- p$data %>%
  group_by(ident) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

p + geom_label(data = label_data,
               aes(x = umap_1, y = umap_2, label = ident),
               size       = 3.5,
               fontface   = "bold",
               fill       = alpha("white", 0.7),
               color      = "black",
               label.size = 0.2) + theme(legend.position = "none")
#Marker identification ----
DefaultAssay(Int) <- "RNA"
Int <- JoinLayers(Int)
markers <- FindAllMarkers(Int, only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(markers,"cluster_markers1.csv", row.names = F)
#Marker Check ----
genes_to_check = c('SFRP5','CHEK1','HSD11B2','CDH1','CLDN4',
                   'CLDN7','DSP','KRT8','KRT18','HAPLN1','MEOX1','SMTNL2','CDK1',
                   'PLK1','PCLAF','DLGAP5','ABHD10','SOX9','TFAP2A','RET','FOXD3',
                   'GADD45A','NGFR','SNAI2',"SOX10",'MITF','PMEL', 'DCT', 'ERBB3','MCAM' ,
                   'SNAI1','HEYL','TFAP2B','CA11', 'NEUROG1', 'NEUROG2', 'NEUROD1', 
                   'NEUROD4', 'ASCL1', 'PRRX1', 'FOXC1','TWIST1', 'COL1A1',
                   'COL3A1', 'CDH2', 'VIM', 'CDON' , 'SOX2','OLIG3', 'PAX6', 
                   'MSX1', 'MSX2', 'WNT1', 'WNT3A','ZIC1',
                   'ZIC4', 'ZIC3', 'MAFB', 'LMX1A', 'GDF7', 'PAK3')
FeaturePlot(Int,
            features   = genes_to_check,
            reduction  = "umap",
            order      = TRUE,         
            pt.size    = 0.2,
            ncol       = 2) &         
  scale_color_viridis_c() &
  theme(legend.position = "right")

#Annotation VJ ----
cluster_annotations <- c(
  "0"  = "Migratory NCCs",
  "1"  = "Schwann cell precursors",
  "2"  = "Transitioning Epithelial Cells",
  "3"  = "Cycling cells",
  "4"  = "NCC Derived Mesenchymal/Stromal Cells",
  "5"  = "iPSCs",
  #"6"  = "Aberrant Pre-Migratory Intermediate",
  "6"  = "KO specific population",
  "7"  = "Melanoblasts",
  #"8"  = "Aberrant Partial EMT State",
  "8"  = "KO specific population",
  "9"  = "Delaminating NCCs",
  "10" = "Neuronal progenitors",
  "11" = "Cycling cells",
  "12" = "iPSCs"
)


Int$cell_type <- dplyr::recode(
  as.character(Int$seurat_clusters),
  !!!cluster_annotations
)

table(Int$cell_type, Int$seurat_clusters)

DimPlot(Int, reduction = "umap",
        group.by = "cell_type",
        label = TRUE,
        repel = TRUE) +
  labs(title = "Annotated Cell Types") +
  theme(legend.text = element_text(size = 8))

DimPlot(Int, reduction = "umap",
        group.by = "cell_type",
        split.by = "condition",
        label = TRUE,
        repel = TRUE) +
  labs(title = "Annotated Cell Types by Condition") +
  theme(legend.text = element_text(size = 8))
ggsave("umap_annotated_split.pdf", width = 16, height = 8)

#Marker Representation----
DotPlot(Int,
        features  =genes_to_check,
        group.by  = "cell_type") +
  RotatedAxis() +
  scale_color_viridis_c() +
  labs(title = "Gene Expression Across Cell Types")

print(
  DotPlot(Int,
          features  = genes_to_check,
          group.by  = "cell_type",
          dot.scale = 8) +
    RotatedAxis() +
    scale_color_viridis_c() +
    labs(title = "Gene Expression Across Cell Types") +
    theme(
      axis.text.x  = element_text(size = 13, angle = 90),
      
      axis.text.y  = element_text(size = 13),
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 12),
      plot.title   = element_text(size = 15, face = "bold")
    )
)

genes_ordered <- c(
  "SOX2", "WNT1", "WNT3A",
  "CDH1", "CLDN4", "CLDN7", "DSP", "KRT8", "KRT18",
  "ZIC1", "ZIC3", "ZIC4", "MSX1", "MSX2", "PAX6", "OLIG3",
  "SNAI1", "SNAI2", "TWIST1", "PRRX1", "FOXC1", "HEYL",
  "SOX9", "TFAP2A", "TFAP2B", "FOXD3", "NGFR", "ERBB3", "RET",
  "HAPLN1", "MEOX1", "MCAM", "CDH2", "VIM",
  "COL1A1", "COL3A1",
  "SOX10","MITF", "PMEL", "DCT", "CDON", "MAFB",
  "SMTNL2", "SFRP5", "CA11", "HSD11B2",
  "NEUROG1", "NEUROG2", "NEUROD1", "NEUROD4", "ASCL1",
  "LMX1A", "GDF7", "PAK3",
  "CDK1", "PLK1", "PCLAF", "DLGAP5", "ABHD10",
  "GADD45A", "CHEK1"
)
cell_type_order <- c(
  "KO specific population",
  "NCC Derived Mesenchymal/Stromal Cells",
  "Neuronal progenitors",
  "Schwann cell precursors",
  "Melanoblasts",
  "Migratory NCCs",
  "Delaminating NCCs",
  "Cycling cells",
  "Transitioning Epithelial Cells",
  "iPSCs"
)

Int <- RenameIdents(Int, cluster_annotations)
Idents(Int) <- factor(Idents(Int), levels = cell_type_order)


Int$cell_type <- Idents(Int)

levels(Idents(Int))

avg_expr <- AverageExpression(Int, 
                              features = genes_ordered, 
                              assay    = "RNA",
                              group.by = "cell_type")$RNA


print(colnames(avg_expr))

gene_peak <- data.frame(
  gene         = rownames(avg_expr),
  peak_cluster = colnames(avg_expr)[apply(avg_expr, 1, which.max)]
)

gene_peak$peak_cluster <- factor(gene_peak$peak_cluster, levels = cell_type_order)
gene_peak_sorted       <- gene_peak[order(gene_peak$peak_cluster), ]

genes_diagonal <- gene_peak_sorted$gene
print(genes_diagonal)

genes_diagonal <- genes_diagonal[genes_diagonal %in% rownames(Int)]

DotPlot(Int,
        features = genes_diagonal) +
  RotatedAxis() +
  scale_color_viridis_c(option = "viridis", name = "Average\nExpression") +
  labs(title = "Gene Expression Across Cell Types", x = "Features", y = "Identity") +
  theme(
    axis.text.x = element_text(size = 14, angle = 80),
    axis.text.y = element_text(size = 14),
    plot.title  = element_text(size = 12, face = "bold", hjust = 0.5)
  )


#Comparison vs condition 

genes_compare <- c("LMX1A",'WNT3A','GDF7','ZIC4','ZIC1','MSX1','MSX2','OLIG3','PAX6','CDON','WNT1','CDH2','MAFB','PAK3','TWIST1','SNAI1','SNAI2','ERBB3','FOXD3','TFAP2A','ABHD10', 'PCLAF','CDK1','PLK1','HAPLN1','SOX10') 
pops_compare <- c("Migratory NCCs", "Delaminating NCCs", "KO specific population")

Int_sub_comparison <- subset(Int, idents = pops_compare)

Int_sub_comparison$group <- paste(Int_sub_comparison$cell_type, Int_sub_comparison$condition, sep = " | ")

group_order <- c(
  "Migratory NCCs | Control",
  "Migratory NCCs | Knockout",
  "KO specific population | Control",
  "KO specific population | Knockout",
  "Delaminating NCCs | Control",
  "Delaminating NCCs | Knockout"
)

Int_sub_comparison$group <- factor(Int_sub_comparison$group, levels = group_order)

genes_compare <- genes_compare[genes_compare %in% rownames(Int_sub_comparison)]

DotPlot(Int_sub_comparison,
        features = genes_compare,
        group.by = "group") +
  RotatedAxis() +
  scale_color_viridis_c(option  = "viridis",
                        name    = "Average\nExpression") +
  geom_hline(yintercept = c(2.5, 4.5),
             linetype   = "dashed",
             color      = "grey50",
             linewidth  = 0.4) +
  labs(
    title = "Gene Expression: Control vs Knockout by Population",
    x     = "Features",
    y     = NULL
  ) +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  ) +
 
  theme(axis.text.y = element_text(
    color = rep(c("#3C8DC4", "#3C8DC4",  
                  "#E06C3A", "#E06C3A",   
                  "#8B5EA6", "#8B5EA6"),  
                1)
  ))


#Cluster wise cell proportions----
composition_table <- table(Int$cell_type, Int$condition)
print(composition_table)

composition_prop <- prop.table(composition_table, margin = 1) * 100
print(round(composition_prop, 2))

comp_df <- as.data.frame(composition_table)
colnames(comp_df) <- c("Cluster", "Condition", "Count")

comp_df <- comp_df %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate(Proportion = Count / sum(Count) * 100)

ggplot(comp_df, aes(x = Cluster, y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Control" = "steelblue", "Knockout" = "firebrick")) +
  theme_classic() +
  labs(title = "Cell Composition per Cell type",
       x = "Cluster", y = "Proportion (%)") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(comp_df, aes(x = Cluster, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Control" = "steelblue", "Knockout" = "firebrick")) +
  theme_classic() +
  labs(title = "Cell Counts per Cluster",
       x = "Cluster", y = "Cell Count") + theme(axis.text.x = element_text(angle = 45, hjust = 1))  

#Monocle3 Analysis----
run_monocle_trajectory <- function(seurat_obj, condition_name) {
  
  gene_annotation <- as.data.frame(
    rownames(seurat_obj@assays$RNA),
    row.names = rownames(seurat_obj@assays$RNA)
  )
  colnames(gene_annotation) <- "gene_short_name"
  
  cell_metadata     <- seurat_obj@meta.data
  expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  
  reducedDim(cds, type = "UMAP") <- seurat_obj@reductions$umap@cell.embeddings
  
  cds <- cluster_cells(cds, reduction_method = "UMAP")

  cds <- learn_graph(cds,
                     use_partition = T,
                     learn_graph_control = list(
                       minimal_branch_len = 5
                     ),close_loop = F)
  
  get_earliest_principal_node <- function(cds, root_cell_type = "Transitioning Epithelial Cells") {
    
    cell_ids <- which(colData(cds)$cell_type == root_cell_type)
    
    closest_vertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    
    root_pr_nodes <- igraph::V(principal_graph(cds)$UMAP)$name[
      as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
    ]
    
    root_pr_nodes
  }
  
  cds <- order_cells(cds,
                     root_pr_nodes = get_earliest_principal_node(cds))
  
  
  p1 <- plot_cells(cds,
                   color_cells_by         = "cell_type",
                   label_cell_groups      = TRUE,
                   label_leaves           = TRUE,
                   label_branch_points    = TRUE,
                   graph_label_size       = 3,
                   cell_size              = 0.7,
                   trajectory_graph_color = "black",
                   trajectory_graph_segment_size = 0.8) +
    labs(title = paste("Trajectory -", condition_name)) +
    theme(legend.text     = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"))
  
  p2 <- plot_cells(cds,
                   color_cells_by         = "pseudotime",
                   label_cell_groups      = FALSE,
                   label_leaves           = TRUE,
                   label_branch_points    = FALSE,
                   graph_label_size       = 3,
                   cell_size              = 0.7,
                   trajectory_graph_color = "black",
                   trajectory_graph_segment_size = 0.8) +
    scale_color_viridis_c() +
    labs(title = paste("Pseudotime -", condition_name)) +
    theme(legend.text     = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"))
  
  return(list(cds = cds, plot_celltype = p1, plot_pseudotime = p2))
}


WT_int <- subset(Int, subset = condition == "Control")
KO_int <- subset(Int, subset = condition == "Knockout")

WT_traj <- run_monocle_trajectory(WT_int, "Control")
KO_traj <- run_monocle_trajectory(KO_int, "Knockout")

WT_traj$plot_celltype + KO_traj$plot_celltype +
  plot_annotation(
    title = "Trajectory Comparison: Control vs Knockout",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
  )

ggsave("trajectory_celltype_comparison.pdf", width = 18, height = 8)

WT_traj$plot_pseudotime + KO_traj$plot_pseudotime +
  plot_annotation(
    title = "Pseudotime Comparison: Control vs Knockout",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
  )

ggsave("trajectory_pseudotime_comparison.pdf", width = 18, height = 8)
#Pseudotime export to Excel ----
library(writexl)

export_pseudotime <- function(traj, condition_name) {
  cds  <- traj$cds
  umap <- reducedDims(cds)$UMAP
  pt   <- pseudotime(cds)
  
  df <- data.frame(
    cell_barcode   = colnames(cds),
    condition      = condition_name,
    seurat_cluster = as.character(colData(cds)$seurat_clusters),
    cell_type      = as.character(colData(cds)$cell_type),
    partition      = as.character(partitions(cds)),
    pseudotime     = pt,
    UMAP_1         = umap[, 1],
    UMAP_2         = umap[, 2],
    row.names      = NULL
  )
  
  # Cells in partitions disconnected from the root get Inf -> flag them, set NA
  df$reachable_from_root <- is.finite(df$pseudotime)
  df$pseudotime[!is.finite(df$pseudotime)] <- NA
  
  # 0-1 scaled within condition (useful for side-by-side comparison)
  rng <- range(df$pseudotime, na.rm = TRUE)
  df$pseudotime_scaled <- (df$pseudotime - rng[1]) / (rng[2] - rng[1])
  df
}

pt_WT <- export_pseudotime(WT_traj, "Control")
pt_KO <- export_pseudotime(KO_traj, "Knockout")
pt_all <- rbind(pt_WT, pt_KO)

safe_stat <- function(x, f) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)

summarise_pt <- function(df) {
  df %>%
    dplyr::group_by(condition, cell_type) %>%
    dplyr::summarise(
      n_cells   = dplyr::n(),
      n_with_pt = sum(!is.na(pseudotime)),
      pct_with_pt = round(100 * n_with_pt / n_cells, 1),
      mean_pt   = safe_stat(pseudotime, mean),
      median_pt = safe_stat(pseudotime, median),
      sd_pt     = safe_stat(pseudotime, sd),
      min_pt    = safe_stat(pseudotime, min),
      max_pt    = safe_stat(pseudotime, max),
      .groups = "drop"
    ) %>%
    dplyr::arrange(condition, median_pt)
}

write_xlsx(
  list(
    "All_cells"           = pt_all,
    "Control_cells"       = pt_WT,
    "Knockout_cells"      = pt_KO,
    "Summary_by_celltype" = summarise_pt(pt_all)
  ),
  path = "Pseudotime_Control_vs_Knockout.xlsx"
)
#Cluster wise DE (Bootstrapping iteration) ----
celltype_to_analyze <- "Cycling cells"
n_iterations <- 25
output_prefix <- paste0(celltype_to_analyze, "_iter", n_iterations)

cluster_subset <- subset(Int, subset = cell_type == celltype_to_analyze)
DefaultAssay(cluster_subset) <- "RNA"
cluster_subset <- JoinLayers(cluster_subset)
cat("\n--- Population Check: Cluster", celltype_to_analyze, "---\n")
pop_table <- table(cluster_subset$condition)
print(pop_table)

ctrl_cells <- colnames(cluster_subset)[cluster_subset$condition == "Control"]
ko_cells <- colnames(cluster_subset)[cluster_subset$condition == "Knockout"]
min_cells <- min(length(ctrl_cells), length(ko_cells))
cat("Subsampling to", min_cells, "cells per condition.\n")

all_iters_list <- list()

for (i in 1:n_iterations) {
  cat("Running Iteration:", i, "/", n_iterations, "...\n")
  set.seed(42 + i)
  balanced_cells <- c(
    sample(ctrl_cells, min_cells),
    sample(ko_cells, min_cells)
  )
  
  tmp_obj <- subset(cluster_subset, cells = balanced_cells)
  Idents(tmp_obj) <- "condition"
  
  de_res <- FindMarkers(
    tmp_obj,
    ident.1 = "Knockout",
    ident.2 = "Control",
    test.use = "wilcox",
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  
  de_res$gene <- rownames(de_res)
  de_res$iteration <- i
  all_iters_list[[i]] <- de_res
}

combined_results <- do.call(rbind, all_iters_list)
  consensus_summary <- combined_results %>%
  group_by(gene) %>%
  summarise(
    times_significant = sum(p_val_adj < 0.05),
    avg_log2FC = mean(avg_log2FC),
    avg_p_val_adj = mean(p_val_adj),
    consistency_score = (times_significant / n_iterations) * 100
  ) %>%
  arrange(desc(times_significant), desc(abs(avg_log2FC)))

write.csv(combined_results, paste0(output_prefix, "_all_iterations_raw.csv"), row.names = FALSE)
write.csv(consensus_summary, paste0(output_prefix, "_consensus_summary.csv"), row.names = FALSE)

top_robust_genes <- consensus_summary %>% 
  head(20) %>% 
  pull(gene)

if(length(top_robust_genes) > 0) {
  p <- DotPlot(cluster_subset, features = top_robust_genes, group.by = "condition") +
    coord_flip() +
    theme_bw() +
    labs(title = paste("Consensus Top Genes - ", celltype_to_analyze),
         subtitle = paste("Ranked by frequency across", n_iterations, "iterations"))
  
  
  ggsave(paste0(output_prefix, "_consensus_dotplot.pdf"), plot = p, width = 7, height = 8)
  
  plot_data <- consensus_summary
  plot_data$avg_p_val_adj[plot_data$avg_p_val_adj == 0] <- 1e-300
  plot_data$Significance <- "Not Significant"
  plot_data$Significance[plot_data$avg_log2FC > 0.25 & plot_data$avg_p_val_adj < 0.05] <- "Upregulated"
  plot_data$Significance[plot_data$avg_log2FC < -0.25 & plot_data$avg_p_val_adj < 0.05] <- "Downregulated"
  
  plot_data$label <- ifelse(plot_data$gene %in% top_robust_genes, plot_data$gene, NA)
  
  
  v_plot <- ggplot(plot_data, aes(x = avg_log2FC, y = -log10(avg_p_val_adj), color = Significance)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Upregulated" = "firebrick3", "Downregulated" = "navy", "Not Significant" = "grey80")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_text_repel(aes(label = label), size = 3, color = "black", max.overlaps = 50) +
    theme_bw() +
    labs(
      title = paste("Consensus Volcano Plot -", celltype_to_analyze),
      subtitle = paste("Averaged across", n_iterations, "iterations"),
      x = expression(Average~Log[2]~Fold~Change),
      y = expression(-Log[10]~(Average~Adjusted~P-value))
    )
  
  ggsave(paste0(output_prefix, "_consensus_volcano.pdf"), plot = v_plot, width = 7, height = 7)
  
}
p
v_plot

#GO Analysis ----
go_genes_KO <- consensus_summary %>%
  filter(times_significant > 0 & avg_log2FC > 0) %>%
  arrange(desc(consistency_score), desc(avg_log2FC)) %>%
  pull(gene)

go_genes_Ctrl <- consensus_summary %>%
  filter(times_significant > 0 & avg_log2FC < 0) %>%
  arrange(desc(consistency_score), desc(avg_log2FC)) %>%
  pull(gene)

cat("\nGenes upregulated in KO:", length(go_genes_KO), "\n")
cat("Genes upregulated in Control:", length(go_genes_Ctrl), "\n")


entrez_KO <- bitr(go_genes_KO,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Hs.eg.db)

entrez_Ctrl <- bitr(go_genes_Ctrl,
                    fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Hs.eg.db)

cat("KO genes mapped:", nrow(entrez_KO), "out of", length(go_genes_KO), "\n")
cat("Control genes mapped:", nrow(entrez_Ctrl), "out of", length(go_genes_Ctrl), "\n")

if(nrow(entrez_KO) >= 5) {
  
  GO_BP_KO <- enrichGO(gene          = entrez_KO$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.3,
                       readable      = TRUE)
  
  GO_MF_KO <- enrichGO(gene          = entrez_KO$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2,
                       readable      = TRUE)
  
  GO_CC_KO <- enrichGO(gene          = entrez_KO$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2,
                       readable      = TRUE)
  
  write.csv(as.data.frame(GO_BP_KO),
            paste0(output_prefix, "_GO_BP_KO.csv"), row.names = FALSE)
  write.csv(as.data.frame(GO_MF_KO),
            paste0(output_prefix, "_GO_MF_KO.csv"), row.names = FALSE)
  write.csv(as.data.frame(GO_CC_KO),
            paste0(output_prefix, "_GO_CC_KO.csv"), row.names = FALSE)
  
  cat("\nKO GO terms - BP:", nrow(as.data.frame(GO_BP_KO) %>% filter(p.adjust < 0.05)),
      "MF:", nrow(as.data.frame(GO_MF_KO) %>% filter(p.adjust < 0.05)),
      "CC:", nrow(as.data.frame(GO_CC_KO) %>% filter(p.adjust < 0.05)), "\n")
  
  if(nrow(GO_BP_KO) > 0) {
    dotplot(GO_BP_KO, showCategory = 15) +
      labs(title = paste(celltype_to_analyze, "- GO BP (KO Upregulated)"))
    ggsave(paste0(output_prefix, "_GO_BP_KO_dotplot.pdf"), width = 10, height = 8)
    
    emapplot(pairwise_termsim(GO_BP_KO), showCategory = 20) +
      labs(title = paste(celltype_to_analyze, "- GO BP Network (KO Upregulated)"))
    ggsave(paste0(output_prefix, "_GO_BP_KO_network.pdf"), width = 12, height = 10)
  }
  
  if(nrow(GO_MF_KO) > 0) {
    dotplot(GO_MF_KO, showCategory = 15) +
      labs(title = paste(celltype_to_analyze, "- GO MF (KO Upregulated)"))
    ggsave(paste0(output_prefix, "_GO_MF_KO_dotplot.pdf"), width = 10, height = 8)
  }
  
  if(nrow(GO_CC_KO) > 0) {
    dotplot(GO_CC_KO, showCategory = 15) +
      labs(title = paste(celltype_to_analyze, "- GO CC (KO Upregulated)"))
    ggsave(paste0(output_prefix, "_GO_CC_KO_dotplot.pdf"), width = 10, height = 8)
  }
  
} else {
  cat("Not enough KO genes mapped for GO analysis.\n")
}


if(nrow(entrez_Ctrl) >= 5) {
  
  GO_BP_Ctrl <- enrichGO(gene          = entrez_Ctrl$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.2,
                         readable      = TRUE)
  
  GO_MF_Ctrl <- enrichGO(gene          = entrez_Ctrl$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.2,
                         readable      = TRUE)
  
  GO_CC_Ctrl <- enrichGO(gene          = entrez_Ctrl$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.2,
                         readable      = TRUE)
  
  write.csv(as.data.frame(GO_BP_Ctrl),
            paste0(output_prefix, "_GO_BP_Ctrl.csv"), row.names = FALSE)
  write.csv(as.data.frame(GO_MF_Ctrl),
            paste0(output_prefix, "_GO_MF_Ctrl.csv"), row.names = FALSE)
  write.csv(as.data.frame(GO_CC_Ctrl),
            paste0(output_prefix, "_GO_CC_Ctrl.csv"), row.names = FALSE)
  
  cat("\nControl GO terms - BP:", nrow(as.data.frame(GO_BP_Ctrl) %>% filter(p.adjust < 0.05)),
      "MF:", nrow(as.data.frame(GO_MF_Ctrl) %>% filter(p.adjust < 0.05)),
      "CC:", nrow(as.data.frame(GO_CC_Ctrl) %>% filter(p.adjust < 0.05)), "\n")
  
  if(nrow(GO_BP_Ctrl) > 0) {
    dotplot(GO_BP_Ctrl, showCategory = 15) +
      labs(title = paste(celltype_to_analyze, "- GO BP (Control Upregulated)"))
    ggsave(paste0(output_prefix, "_GO_BP_Ctrl_dotplot.pdf"), width = 10, height = 8)
    
    emapplot(pairwise_termsim(GO_BP_Ctrl), showCategory = 20) +
      labs(title = paste(celltype_to_analyze, "- GO BP Network (Control Upregulated)"))
    ggsave(paste0(output_prefix, "_GO_BP_Ctrl_network.pdf"), width = 12, height = 10)
  }
  
  if(nrow(GO_MF_Ctrl) > 0) {
    dotplot(GO_MF_Ctrl, showCategory = 15) +
      labs(title = paste(celltype_to_analyze, "- GO MF (Control Upregulated)"))
    ggsave(paste0(output_prefix, "_GO_MF_Ctrl_dotplot.pdf"), width = 10, height = 8)
  }
  
  if(nrow(GO_CC_Ctrl) > 0) {
    dotplot(GO_CC_Ctrl, showCategory = 15) +
      labs(title = paste(celltype_to_analyze, "- GO CC (Control Upregulated)"))
    ggsave(paste0(output_prefix, "_GO_CC_Ctrl_dotplot.pdf"), width = 10, height = 8)
  }
  
} else {
  cat("Not enough Control genes mapped for GO analysis.\n")
}


if(exists("GO_BP_KO") & exists("GO_BP_Ctrl")) {
  if(nrow(GO_BP_KO) > 0 & nrow(GO_BP_Ctrl) > 0) {
    
    p_KO <- dotplot(GO_BP_KO, showCategory = 15) +
      labs(title = "KO Upregulated") +
      theme(plot.title = element_text(hjust = 0.5, color = "red"))
    
    p_Ctrl <- dotplot(GO_BP_Ctrl, showCategory = 15) +
      labs(title = "Control Upregulated") +
      theme(plot.title = element_text(hjust = 0.5, color = "blue"))
    
    p_KO + p_Ctrl +
      plot_annotation(title = paste(celltype_to_analyze, "- GO BP Comparison"))
    
    ggsave(paste0(output_prefix, "_GO_BP_comparison.pdf"), width = 18, height = 8)
  }
}


#DE + GO within conditions ----

celltype_focus       <- "KO specific population"
celltype_predecessor <- "Transitioning Epithelial Cells"
celltype_successor   <- "Migratory NCCs"

for(cond in c("Control", "Knockout")) {
  
  obj <- subset(Int, subset = condition == cond)
  
  focus_n <- sum(obj$cell_type == celltype_focus)
  pred_n  <- sum(obj$cell_type == celltype_predecessor)
  succ_n  <- sum(obj$cell_type == celltype_successor)
  
  max_focus    <- min(focus_n, pred_n * 2, succ_n * 2)
  min_per_half <- floor(max_focus / 2)
  
  cat("Condition:", cond, "\n")
  cat(sprintf("  %-45s %d cells\n", celltype_focus,       focus_n))
  cat(sprintf("  %-45s %d cells\n", celltype_predecessor, pred_n))
  cat(sprintf("  %-45s %d cells\n", celltype_successor,   succ_n))
  cat("  → Max usable N (focus cells):", max_focus, "\n")
  cat("  → N/2 sampled per neighbor:  ", min_per_half, "\n\n")
}


n_iterations <- 25   
min_cells    <- 62  

run_trajectory_DE <- function(seurat_obj,
                              condition_name,
                              focus,
                              predecessor,
                              successor,
                              n_iter,
                              n_cells) {
  
  output_prefix <- paste0(
    gsub(" ", "_", focus), "_vs_neighbors_",
    condition_name, "_iter", n_iter
  )
  
  cat("Running:", condition_name, "|", focus, "vs neighbors\n")
  cat("N per group:", n_cells, "| Iterations:", n_iter, "\n")

  
  obj <- subset(seurat_obj, subset = condition == condition_name)
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  
  focus_cells <- colnames(obj)[obj$cell_type == focus]
  pred_cells  <- colnames(obj)[obj$cell_type == predecessor]
  succ_cells  <- colnames(obj)[obj$cell_type == successor]
  
  half_n <- floor(n_cells / 2)
  
  if(length(focus_cells) < n_cells) {
    cat("WARNING: Focus has only", length(focus_cells), 
        "cells but n_cells =", n_cells, "— reduce min_cells\n")
    return(NULL)
  }
  if(length(pred_cells) < half_n) {
    cat("WARNING: Predecessor has only", length(pred_cells),
        "cells but need", half_n, "— reduce min_cells\n")
    return(NULL)
  }
  if(length(succ_cells) < half_n) {
    cat("WARNING: Successor has only", length(succ_cells),
        "cells but need", half_n, "— reduce min_cells\n")
    return(NULL)
  }
  
  
  all_iters_list <- list()
  
  for(i in 1:n_iter) {
    cat("Iteration:", i, "/", n_iter, "...\n")
    
    set.seed(42 + i)
    
    sampled_focus <- sample(focus_cells, n_cells)
    sampled_pred  <- sample(pred_cells,  half_n)
    sampled_succ  <- sample(succ_cells,  half_n)
    
    balanced_cells <- c(sampled_focus, sampled_pred, sampled_succ)
    
    tmp_obj <- subset(obj, cells = balanced_cells)
    tmp_obj$traj_group <- ifelse(tmp_obj$cell_type == focus,
                                 "Focus", "Neighbor")
    Idents(tmp_obj) <- "traj_group"
    
    de_res <- FindMarkers(
      tmp_obj,
      ident.1         = "Focus",
      ident.2         = "Neighbor",
      test.use        = "wilcox",
      min.pct         = 0.25,
      logfc.threshold = 0.25,
      verbose         = FALSE
    )
    
    de_res$gene      <- rownames(de_res)
    de_res$iteration <- i
    all_iters_list[[i]] <- de_res
  }
  
  combined_results <- do.call(rbind, all_iters_list)
  
  consensus_summary <- combined_results %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      times_significant = sum(p_val_adj < 0.05),
      avg_log2FC        = mean(avg_log2FC),
      avg_p_val_adj     = mean(p_val_adj),
      consistency_score = (times_significant / n_iter) * 100
    ) %>%
    dplyr::arrange(desc(times_significant), desc(abs(avg_log2FC)))
  
  write.csv(combined_results,
            paste0(output_prefix, "_all_iterations_raw.csv"), row.names = FALSE)
  write.csv(consensus_summary,
            paste0(output_prefix, "_consensus_summary.csv"), row.names = FALSE)
  
  
  top_genes <- consensus_summary %>%
    dplyr::filter(avg_log2FC > 0) %>%
    head(20) %>%
    dplyr::pull(gene)
  
  if(length(top_genes) > 0) {
    
    display_obj <- subset(obj,
                          subset = cell_type %in% c(focus, predecessor, successor))
    display_obj$cell_type <- factor(display_obj$cell_type,
                                    levels = c(predecessor, focus, successor))
    
    p <- DotPlot(display_obj,
                 features = top_genes,
                 group.by = "cell_type") +
      coord_flip() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title    = paste("Unique Genes —", focus),
           subtitle = paste(condition_name, "| N =", n_cells,
                            "| Iterations =", n_iter))
    
    ggsave(paste0(output_prefix, "_dotplot.pdf"),
           plot = p, width = 8, height = 9)
  }
  
  
  go_genes <- consensus_summary %>%
    dplyr::filter(times_significant > 0 & avg_log2FC > 0) %>%
    dplyr::pull(gene)
  
  cat("\nGenes upregulated in focus:", length(go_genes), "\n")
  
  if(length(go_genes) >= 5) {
    
    entrez <- bitr(go_genes,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)
    
    cat("Mapped to Entrez:", nrow(entrez), "\n")
    
    GO_BP <- enrichGO(gene          = entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)
    
    GO_MF <- enrichGO(gene          = entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)
    
    GO_CC <- enrichGO(gene          = entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)
    
    write.csv(as.data.frame(GO_BP),
              paste0(output_prefix, "_GO_BP.csv"), row.names = FALSE)
    write.csv(as.data.frame(GO_MF),
              paste0(output_prefix, "_GO_MF.csv"), row.names = FALSE)
    write.csv(as.data.frame(GO_CC),
              paste0(output_prefix, "_GO_CC.csv"), row.names = FALSE)
    
    cat("GO terms — BP:", nrow(GO_BP),
        "| MF:", nrow(GO_MF),
        "| CC:", nrow(GO_CC), "\n")
    
    if(nrow(GO_BP) > 0) {
      p_bp <- dotplot(GO_BP, showCategory = 15) +
        labs(title = paste(focus, "—", condition_name, "| GO BP"))
      ggsave(paste0(output_prefix, "_GO_BP_dotplot.pdf"),
             plot = p_bp, width = 10, height = 8)
      
      if(nrow(GO_BP) > 1) {
        p_net <- emapplot(pairwise_termsim(GO_BP), showCategory = 20) +
          labs(title = paste(focus, "—", condition_name, "| GO BP Network"))
        ggsave(paste0(output_prefix, "_GO_BP_network.pdf"),
               plot = p_net, width = 12, height = 10)
      }
    }
    
    if(nrow(GO_MF) > 0) {
      p_mf <- dotplot(GO_MF, showCategory = 15) +
        labs(title = paste(focus, "—", condition_name, "| GO MF"))
      ggsave(paste0(output_prefix, "_GO_MF_dotplot.pdf"),
             plot = p_mf, width = 10, height = 8)
    }
    
    if(nrow(GO_CC) > 0) {
      p_cc <- dotplot(GO_CC, showCategory = 15) +
        labs(title = paste(focus, "—", condition_name, "| GO CC"))
      ggsave(paste0(output_prefix, "_GO_CC_dotplot.pdf"),
             plot = p_cc, width = 10, height = 8)
    }
    
  } else {
    cat("Not enough genes for GO analysis.\n")
  }
  
  cat("\nComplete:", output_prefix, "\n")
  
  return(list(consensus = consensus_summary,
              combined  = combined_results))
}

results_Control <- run_trajectory_DE(
  seurat_obj     = Int,
  condition_name = "Control",
  focus          = celltype_focus,
  predecessor    = celltype_predecessor,
  successor      = celltype_successor,
  n_iter         = n_iterations,
  n_cells        = min_cells
)

results_KO <- run_trajectory_DE(
  seurat_obj     = Int,
  condition_name = "Knockout",
  focus          = celltype_focus,
  predecessor    = celltype_predecessor,
  successor      = celltype_successor,
  n_iter         = n_iterations,
  n_cells        = min_cells
)



if(!is.null(results_Control) & !is.null(results_KO)) {
  
  wt_unique <- results_Control$consensus %>%
    dplyr::filter(times_significant > 0 & avg_log2FC > 0) %>%
    dplyr::pull(gene)
  
  ko_unique <- results_KO$consensus %>%
    dplyr::filter(times_significant > 0 & avg_log2FC > 0) %>%
    dplyr::pull(gene)
  
  cat("Genes unique to WT:", length(setdiff(wt_unique, ko_unique)), "\n")
  cat("Genes unique to KO:", length(setdiff(ko_unique, wt_unique)), "\n")
  cat("Shared genes:",       length(intersect(wt_unique, ko_unique)), "\n")
  
  write.csv(data.frame(gene = setdiff(wt_unique, ko_unique)),
            paste0(gsub(" ", "_", celltype_focus), "_genes_lost_in_KO.csv"),
            row.names = FALSE)
  write.csv(data.frame(gene = setdiff(ko_unique, wt_unique)),
            paste0(gsub(" ", "_", celltype_focus), "_genes_gained_in_KO.csv"),
            row.names = FALSE)
  write.csv(data.frame(gene = intersect(wt_unique, ko_unique)),
            paste0(gsub(" ", "_", celltype_focus), "_genes_shared.csv"),
            row.names = FALSE)
}
#Against EMT genes----
gene_data <- read.csv('/Users/vijeyeshwar/CSIR IGIB/ScRNAseq/Seurat Archive/Gene List/Human_EMTome_Binned.csv')
clean_genes <- function(x) {
  unique(na.omit(x[x != ""]))
}
marker_lists <- list(
  epithelial    = clean_genes(gene_data$EPITHELIAL),
  mesenchymal   = clean_genes(gene_data$MESENCHYMAL),
  regulatory    = clean_genes(gene_data$REGULATORY),
  translational = clean_genes(gene_data$TRANSLATIONAL)
)
marker_lists_Fil <- list(
  epithelial    = intersect(marker_lists$epithelial, rownames(Int)),
  mesenchymal   = intersect(marker_lists$mesenchymal, rownames(Int)),
  regulatory    = intersect(marker_lists$regulatory, rownames(Int)),
  translational = intersect(marker_lists$translational, rownames(Int))
)
Emt_final <- c(
  marker_lists_Fil$epithelial,
  marker_lists_Fil$mesenchymal,
  marker_lists_Fil$regulatory,
  marker_lists_Fil$translational
)
Int$cell_condition <- paste0(Int$cell_type, "_", Int$condition)
unique_cell_types <- c(
  "Migratory NCCs",
  "Schwann cell precursors",
  "Transitioning Epithelial Cells",
  "Cycling cells",
  "NCC Derived Mesenchymal/Stromal Cells",
  "iPSCs",
  "KO specific population",
  "Delaminating NCCs",
  "Neuronal progenitors" )

conditions <- c("Control", "Knockout")
new_order <- paste0(rep(unique_cell_types, each = length(conditions)), "_", conditions)
Int$cell_condition <- factor(Int$cell_condition, levels = new_order)
Idents(Int) <- "cell_condition"
Int_downsampled <- subset(Int, downsample = 100) 
Int_clean <- subset(Int, cells = colnames(Int)[!is.na(Int$cell_condition)])
Idents(Int_clean) <- "cell_condition"
Int_downsampled <- subset(Int_clean, downsample = 100)

DoHeatmap(Int_downsampled, features = Emt_final, size = 3) + 
  ggtitle("Overall EMTome Transcriptional Profile by Cell Type") +
  theme(text = element_text(size = 8))
DoHeatmap(Int_downsampled, features = marker_lists_Fil$epithelial, size = 3) + 
  ggtitle("Epithelial Markers") +
  theme(text = element_text(size = 8))

DoHeatmap(Int_downsampled, features = marker_lists_Fil$mesenchymal, size = 3) + 
  ggtitle("Mesenchymal Markers") +
  theme(text = element_text(size = 8))

DoHeatmap(Int_downsampled, features = marker_lists_Fil$regulatory, size = 3) + 
  ggtitle("Regulatory Markers") +
  theme(text = element_text(size = 8))

DoHeatmap(Int_downsampled, features = marker_lists_Fil$translational, size = 3) + 
  ggtitle("Translational Markers") +
  theme(text = element_text(size = 8))

#EMT Scoring----
Int <- AddModuleScore(object = Int, 
                      features = list(marker_lists_Fil$epithelial), 
                      name = "Epithelial_Score",
                      ctrl = 10)
Int <- AddModuleScore(object = Int, 
                      features = list(marker_lists_Fil$mesenchymal), 
                      name = "Mesenchymal_Score",
                      ctrl = 10)
offset <- abs(min(c(Int$Epithelial_Score1, Int$Mesenchymal_Score1))) + 0.1

Int$Epi_shifted    <- Int$Epithelial_Score1   + offset
Int$Mes_shifted    <- Int$Mesenchymal_Score1  + offset

Int$EMT_ratio <- Int$Epi_shifted / Int$Mes_shifted
Int$EMT_log_ratio <- log2(Int$Epi_shifted / Int$Mes_shifted)


cluster_order <- Int@meta.data %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(mean_ratio = mean(EMT_ratio, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(mean_ratio) %>%
  dplyr::pull(cell_type)

Int$cell_type <- factor(Int$cell_type, levels = cluster_order)

plot_data_vln <- Int@meta.data %>%
  dplyr::filter(!(cell_type == "KO specific population" & condition == "Control"))
Int_vln <- subset(Int, 
                  cells = rownames(Int@meta.data %>%
                                     dplyr::filter(!(cell_type == "KO specific population" & condition == "Control"))))

p <- VlnPlot(Int,
             features  = "EMT_ratio",
             group.by  = "cell_type",
             split.by  = "condition",
             pt.size   = 0)

p +
  geom_hline(yintercept = 1,
             linetype = "dashed", color = "red", linewidth = 0.8) +
  stat_summary(aes(group = interaction(ident, split)),
               fun = mean, geom = "point",
               size = 2, color = "black",
               position = position_dodge(width = 0.9)) +
  labs(title = "Epithelial / Mesenchymal Score Ratio",
       y     = "E/M Ratio (>1 = Epithelial, <1 = Mesenchymal)")


paired_data <- Int@meta.data %>%
  dplyr::filter(!(cell_type == "KO specific population" & condition == "Control")) %>%
  dplyr::group_by(cell_type, condition) %>%
  dplyr::summarise(
    mean_ratio = mean(EMT_ratio, na.rm = TRUE),
    .groups    = "drop"
  )

delta_arrow <- paired_data %>%
  dplyr::mutate(cell_type = as.character(cell_type)) %>%
  dplyr::group_by(cell_type, condition) %>%
  dplyr::summarise(mean_ratio = mean(mean_ratio, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, 
                     values_from = mean_ratio) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      is.na(Control)     ~ "KO only",
      Knockout > Control ~ "Epithelial shift in KO",
      Knockout < Control ~ "Mesenchymal shift in KO",
      TRUE               ~ "No change"
    )
  )

delta_arrow$cell_type <- factor(
  delta_arrow$cell_type,
  levels = delta_arrow$cell_type[order(delta_arrow$Control, na.last = TRUE)]
)

ggplot(delta_arrow, aes(y = cell_type)) +
  geom_segment(aes(x = Control, xend = Knockout,
                   y = cell_type, yend = cell_type,
                   color = direction),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth = 0.8) +
  geom_point(aes(x = Control), shape = 21, 
             fill = "white", size = 2.5, stroke = 1) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "red", linewidth = 0.6) +
  scale_color_manual(values = c(
    "Epithelial shift in KO"  = "#E87D72",
    "Mesenchymal shift in KO" = "#53B4E0",
    "KO only"                 = "#9B59B6",
    "No change"               = "grey60"
  )) +
  theme_classic(base_size = 14) +        # controls all base text
  theme(
    legend.position    = "bottom",
    legend.text        = element_text(size = 13),
    legend.title       = element_text(size = 13),
    axis.text          = element_text(size = 13),
    axis.title         = element_text(size = 14),
    plot.title         = element_text(size = 15, face = "bold"),
    plot.subtitle      = element_text(size = 12),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4)
  ) +
  labs(
    title    = "EMT  Shift: Control → Knockout",
    subtitle = "Arrow tip = KO value  ·  Open circle = Control value  ·  Red line = equipoise",
    x        = "Mean E/M Ratio",
    y        = NULL,
    color    = NULL
  )
#Shannon Entropy----
DefaultAssay(Int) <- "RNA"
Int <- JoinLayers(Int, assay = "RNA")
set.seed(1234)
Int_sub <- subset(Int, cells = sample(Cells(Int), 3000))
counts_mat_sub <- as.matrix(Int_sub[["RNA"]]$counts)
counts_mat<-as.matrix(Int[["RNA"]]$counts)

Int_sub <- JoinLayers(Int_sub)
Int_sub <- NormalizeData(Int_sub)
Int_sub <- ScaleData(Int_sub)
Int_sub <- RunPCA(Int_sub)
ElbowPlot(Int_sub, reduction = "pca")
Int_sub <- FindNeighbors(Int_sub)
Int_sub <- FindClusters(Int_sub)
Int_sub <- RunUMAP(Int_sub, dims = 1:20, reduction = "pca")

Sys.setenv("R_MAX_VSIZE" = 100e9)
cytotrace_results <- CytoTRACE(counts_mat_sub)
Int <- AddMetaData(object = Int, 
                   metadata = cytotrace_results$CytoTRACE, 
                   col.name = "CytoTRACE_Score")
Int_sub <- AddMetaData(object = Int_sub, 
                       metadata = cytotrace_results$CytoTRACE, 
                       col.name = "CytoTRACE_Score")

plot_data <- Int_sub@meta.data
plot_data <- Int@meta.data

plot_data_filtered <- plot_data %>%
  filter(cell_type %in% c("Migratory NCCs",'iPSCs',"Melanoblasts","Aberrant Partial EMT State","Aberrant Pre-Migratory Intermediate"))

plot_data_clean <- plot_data_filtered %>%
  filter(!is.na(CytoTRACE_Score))

plot_data_clean <- plot_data_clean %>%
  mutate(
    
    cell_type = factor(cell_type, levels = c("iPSCs", "Migratory NCCs", "Melanoblasts", "Aberrant Partial EMT State","Aberrant Pre-Migratory Intermediate")),
    
    condition = factor(condition, levels = c("Control", "Knockout")) 
  )

ggplot(plot_data_clean, aes(x = cell_type, y = CytoTRACE_Score, fill = condition)) +
  
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", 
              position = position_dodge(width = 0.8)) +
  
  geom_boxplot(aes(group = interaction(cell_type, condition)), 
               width = 0.1, fill = "white", color = "black", outlier.shape = NA, 
               position = position_dodge(width = 0.8)) +
  
  theme_classic() +
  
  labs(title = "Transcriptomic Diversity by Cell Type",
       x = "Cell Type",
       y = "CytoTRACE Score",
       fill = "Condition") + 
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_text(face = "bold", color = "black"),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "right") + 
  
  scale_fill_manual(values = c("Control" = "#377EB8", "Knockout" = "#E41A1C"))

cells_with_scores <- rownames(Int@meta.data[!is.na(Int@meta.data$CytoTRACE_Score), ])
Int_scored_only <- subset(Int, cells = cells_with_scores)

FeaturePlot(Int_scored_only,
            features = 'CytoTRACE_Score', 
            cols = c("lightgrey", "red" ),split.by = "condition", 
            pt.size = 1) 

#Cell cycle score ----
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Int <- CellCycleScoring(Int, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
DimPlot(Int, reduction = "umap", group.by = "Phase", split.by = "condition")
VlnPlot(Int, features = c("S.Score", "G2M.Score"), group.by = "cell_type")
table(Int$Phase, Idents(Int))
prop.table(table(Int$Phase, Idents(Int)), margin = 2)

metadata <- Int@meta.data
ggplot(metadata, aes(x = cell_type, fill = Phase)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(y = "Proportion", x = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + facet_wrap(~ condition)
