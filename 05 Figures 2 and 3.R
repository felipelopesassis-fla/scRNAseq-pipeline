---
Author: Felipe Assis, PhD
title: "Figures 2 and 3"
date: '12-14-2022'
---

## Load Libraries 

library(Seurat)
library(patchwork)
library(metap)
library(ggplot2)
library(dplyr)
library(cowplot)
library(Platypus)
library(scales)
library(dittoSeq)

### Figure 2B

Moderna.integrated <- RenameIdents(Moderna.integrated, `7` = "Naive-C1")
Moderna.integrated <- RenameIdents(Moderna.integrated, `6` = "Naive-C7") 
Moderna.integrated <- RenameIdents(Moderna.integrated, `5` = "Switched MBC-C5")
Moderna.integrated <- RenameIdents(Moderna.integrated, `4` = "PB-C6")
Moderna.integrated <- RenameIdents(Moderna.integrated, `3` = "Switched MBC-C4") 
Moderna.integrated <- RenameIdents(Moderna.integrated, `2` = "unswitched MBC-C3")
Moderna.integrated <- RenameIdents(Moderna.integrated, `1` = "Naive-C2")
Moderna.integrated <- RenameIdents(Moderna.integrated, `0` = "Naive-C1")

p1 <- DimPlot(Moderna.integrated, reduction = "umap") 

### Figure 2C

## Cluster annotation using RNA transcriptome (RNA) and protein surface (ADT) data

# Protein surface normalization
DefaultAssay(Moderna.integrated) <- "ADT"
Moderna.integrated <- NormalizeData(Moderna.integrated, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

## Set the lists of B cell markers 

rna.features <- c( "CD19", "MS4A1", "CD38", "CD27", "IGHG1", "IGHD", "IGHM", "IGHA1") 
prot.features <- c("CD19-PROT", "CD20-PROT", "CD21-PROT", "CD27-PROT", "CD38-PROT", "IgGFc-PROT", "IgD-PROT", "IgM-PROT") 

p2 <- VlnPlot(Moderna.integrated, assay = "RNA", feature = rna.features, stack = TRUE) + NoLegend()
p3 <- VlnPlot(Moderna.integrated, assay = "ADT", feature = prot.features, stack = TRUE) + NoLegend()
p2+p3 

### Figure 2D

## Find differentially expressed genes (DEGs) 
moderna.markers <- FindAllMarkers(Moderna.integrated, assay = "RNA", only.pos = TRUE)

## Remove Ig genes and find Top10 DEGs 
moderna.markers = moderna.markers[!grepl("IG",moderna.markers$gene),]
moderna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> Top10 

## Remove duplicate genes if needed and generate the Figure 
Top10 <- Top10[-c(46),]   

p4 <- DotPlot(Moderna.integrated, assay = "RNA", features = Top10$gene, cols = "RdBu", dot.scale = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  theme(axis.text.y = element_text(size = 10)) + coord_flip()


### Figure 2F

## Retrieve the Ids from the selected cells

Idents(Moderna.integrated) <- "orig.ident"

v1D14 <- WhichCells(Moderna.integrated, idents = "S2pos.v1D14")
v2D6 <- WhichCells(Moderna.integrated, idents = "S2pos.v2D6")
v2D9 <- WhichCells(Moderna.integrated, idents = "S2pos.v2D9")
v2D14 <- WhichCells(Moderna.integrated, idents = "S2pos.v2D14")
v2D28 <- WhichCells(Moderna.integrated, idents = "S2pos.v2D28")
M6 <- WhichCells(Moderna.integrated, idents = "S2pos.M6")

p5a <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v1D14)
p5b <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v2D6)
p5c <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v2D9)
p5d <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v2D14)
p5e <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v2D28)
p5f <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = M6)

########### Analysis of switched MBC-C5 (re-clustering) 

cluster5 <- subset(Moderna.integrated, idents = "switched MBC-C5")

DefaultAssay(cluster5) <- "RNA"

## Normalize, scale and find clusters 

cluster5 <- NormalizeData(cluster5, normalization.method = "LogNormalize", scale.factor = 10000)
cluster5 <- ScaleData(cluster5, verbose = FALSE)
cluster5 <- RunPCA(cluster5, npcs = 15, verbose = FALSE)
cluster5 <- RunUMAP(cluster5, reduction = "pca", dims = 1:15)
cluster5 <- FindNeighbors(cluster5, reduction = "pca", dims = 1:15)
cluster5 <- FindClusters(cluster5, graph.name = "RNA_snn", resolution = 0.5)

## Data visualization  

DimPlot(cluster5, reduction = "umap") 

## Rename clusters

cluster5 <- RenameIdents(cluster5, `4` = "MBC-SC5.5")
cluster5 <- RenameIdents(cluster5, `3` = "MBC-SC5.4")
cluster5 <- RenameIdents(cluster5, `2` = "MBC-SC5.3")
cluster5 <- RenameIdents(cluster5, `1` = "MBC-SC5.2")
cluster5 <- RenameIdents(cluster5, `0` = "MBC-SC5.1")

### Figure 3A

colors.cluster5 <- c("#999999", "#e41a1c", "#ff7f00", "#f781bf", "#377eb8")

p6 <- DimPlot(cluster5, reduction = "umap", cols = colors.cluster5) + 
  theme(legend.text = element_text(size = 16)) +
  theme(axis.text = element_text(size = 18)) +
  theme(text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes=list(shape = 15, size = 6))) 

### Figure 3B

C5.markers <- FindAllMarkers(cluster5, assay = "RNA", only.pos = FALSE)
C5.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> Top10

Top10 <- Top10[-c(8),] #remove duplicates if needed

p7 <- dittoHeatmap(cluster5, assay = "RNA", slot = "data", order.by = "seurat_clusters", annot.by = "seurat_clusters", 
             annot.colors = colors.cluster5, 
             heatmap.colors = colorRampPalette(c("blue", "white", "firebrick"))(100), 
             breaks = seq(-3, 3, 0.06), genes = Top10$gene, scaled.to.max = FALSE)

### Figure 3C

## Normalize and scale ADT

DefaultAssay(cluster5) <- "ADT"

cluster5 <- NormalizeData(cluster5, normalization.method = 'CLR', scale.factor = 10000) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

## List of B cell markers and plots 

feat.adt.c5 <- c("CD19-PROT", "CD20-PROT", "CD21-PROT", "CD27-PROT", "CD38-PROT", "CD95-PROT", "CD11c-PROT", "CD71-PROT", "IgGFc-PROT", "IgD-PROT", "IgM-PROT")
feat.rna.c5 <- c("CD19", "MS4A1", "CD27", "CD38",  "FAS", "ITGAX", "TFRC", "IGHG1","IGHD", "IGHM", "IGHA1") 

p8 <- VlnPlot(cluster5, assay = "ADT", feature = feat.adt.c5, stack = TRUE) + NoLegend() + 
  theme(text = element_text(size = 14)) + theme(axis.text.x = element_text(size = 12))

p9 <- VlnPlot(cluster5, assay = "RNA", feature = feat.rna.c5, stack = TRUE) + NoLegend() +
  theme(text = element_text(size = 14)) + theme(axis.text.x = element_text(size = 11))

### Figure 3D

## Retrieve the Ids from the selected cells/timepoints   

Idents(Moderna.integrated) <- "orig.ident"

p10a <- DimPlot(cluster5, reduction = "umap", cells.highlight = v1D14, sizes.highlight = 0.2) + 
  theme(text = element_text(size = 18)) + ggtitle("v1D14") + NoLegend() 
p10b <- DimPlot(cluster5, reduction = "umap", cells.highlight = v2D6, sizes.highlight = 0.2) + 
  theme(text = element_text(size = 18)) + ggtitle("v2D6") + NoLegend()
p10c <- DimPlot(cluster5, reduction = "umap", cells.highlight = v2D9, sizes.highlight = 0.2) + 
  theme(text = element_text(size = 18)) + ggtitle("v2D9") + NoLegend()
p10d <- DimPlot(cluster5, reduction = "umap", cells.highlight = v2D14, sizes.highlight = 0.2) + 
  theme(text = element_text(size = 18)) + ggtitle("v2D14") + NoLegend() 
p10e <- DimPlot(cluster5, reduction = "umap", cells.highlight = v2D28, sizes.highlight = 0.2) + 
  theme(text = element_text(size = 18)) + ggtitle("v2D28") + NoLegend() 
p10f <- DimPlot(cluster5, reduction = "umap", cells.highlight = M6, sizes.highlight = 0.2) + 
  theme(text = element_text(size = 18)) + ggtitle("M6") + NoLegend()

### Figure 3E

## Retrieve cell ids

D28.cells <- WhichCells(cluster5, idents = "v2D28")
M6.cells <- WhichCells(cluster5, idents = "M6")

## Get the counts for RBD and S1 for each timepoint  

S1.D28 <- FetchData(cluster5, "S1-PROT", slot = "count", cells = D28.cells)
S1.M6 <- FetchData(cluster5, "S1-PROT", slot = "count", cells = M6.cells)
RBD.D28 <- FetchData(cluster5, "RBD-PROT", slot = "count", cells = D28.cells)
RBD.M6 <- FetchData(cluster5, "RBD-PROT", slot = "count", cells = M6.cells)

## Get the Threshold for each protein using Threshold-seq (https://cm.jefferson.edu/threshold-seq/)   
## Filter S1/RBD negative cells out

s1.d28.pos <- filter(S1.D28, S1 >= 5)
s1.m6.pos <- filter(S1.M6, S1 >= 150)
rbd.m6.pos <- filter(RBD.M6, RBD >= 75 )
rbd.d28.pos <- filter(RBD.D28, RBD >= 31)

## Ids of S1/RBD positive cells to highlight and plots  

rbd.m6.highlight <- factor(c(rownames(rbd.m6.pos)))
rbd.d28.highlight <- factor(c(rownames(rbd.d28.pos)))

s1.m6.highlight <- factor(c(rownames(s1.m6.pos)))
s1.d28.highlight <- factor(c(rownames(s1.d28.pos)))

Idents(cluster5) <- "timepoint"

p11a <- DimPlot(cluster5, reduction = "umap", cells.highlight = list(v1D28 = s1.d28.highlight, M6 = s1.m6.highlight), 
                sizes.highlight = 0.5, cols.highlight = c("#b2182b", "#2166ac", "lightgrey"), ) +
  theme(text = element_text(size = 18)) + ggtitle("S1+ cells") 

P11b <- DimPlot(cluster5, reduction = "umap", cells.highlight = list(v1D28 = rbd.d28.highlight, M6 = rbd.m6.highlight), 
              sizes.highlight = 0.5, cols.highlight = c("#b2182b", "#2166ac", "lightgrey"), ) +
  theme(text = element_text(size = 18)) + ggtitle("RBD+ cells") 

### Figure 3F

## load packages 

library(monocle3)
library(SeuratWrappers)

## Convert Seurat object into monocle dataset

cds2 <- as.cell_data_set(cluster5)
cds2 <- cluster_cells(cds2)
p12 <- plot_cells(cds, show_trajectory_graph = FALSE)
p13 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p12, p13)

## Fit a principal graph within each partition

cds2 <- learn_graph(cds2)
plot_cells(cds2, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, color_cells_by = "partition")

## Specify the root nodes of the trajectory graph

cds2 <- order_cells(cds2, reduction_method = "UMAP")
plot_cells(cds2, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE, label_roots = TRUE, graph_label_size = 10, ) + 
  theme(text = element_text(size = 18)) + 
  theme(legend.text = element_text(size = 16)) +
  theme(axis.text = element_text(size = 18)) +
  theme(text = element_text(size = 18))

### END ###
