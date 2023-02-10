---
Author: 'Felipe Assis, PhD'
title: 'scRNAseq analysis - Tracking B Cell Responses to the SARS-CoV-2 mRNA-1273 Vaccine'
date: '12-14-2022'
---
##Load Libraries 
library(Seurat)
library(ggplot2)
library(dplyr)
library(dittoSeq)

Moderna.integrated <- readRDS("Seurat/Moderna.integrated.rds")

#### Figure 2B
Moderna.integrated <- RenameIdents(Moderna.integrated, `7` = "Naive-C1")
Moderna.integrated <- RenameIdents(Moderna.integrated, `6` = "Naive-C7") 
Moderna.integrated <- RenameIdents(Moderna.integrated, `5` = "Switched MBC-C5")
Moderna.integrated <- RenameIdents(Moderna.integrated, `4` = "PB-C6")
Moderna.integrated <- RenameIdents(Moderna.integrated, `3` = "Switched MBC-C4") 
Moderna.integrated <- RenameIdents(Moderna.integrated, `2` = "Unswitched MBC-C3")
Moderna.integrated <- RenameIdents(Moderna.integrated, `1` = "Naive-C2")
Moderna.integrated <- RenameIdents(Moderna.integrated, `0` = "Naive-C1")

p1 <- DimPlot(Moderna.integrated, pt.size = 0.5, reduction = "umap") 
pdf("Figures/Figure_2B.pdf",width=8,height=4,useDingbats=FALSE)
print(p1)
dev.off()

#### Figure 2C
### Cluster annotation using RNA transcriptome (RNA) and protein surface (ADT) data
## Protein surface normalization
DefaultAssay(Moderna.integrated) <- "ADT"

Moderna.integrated <- NormalizeData(Moderna.integrated, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

## Set the lists of B cell markers 

rna.features <- c( "CD19", "MS4A1", "CD38", "CD27", "IGHG1", "IGHD", "IGHM", "IGHA1") 
prot.features <- c("CD19-PROT", "CD20-PROT", "CD21-PROT", "CD27-PROT", "CD38-PROT", "IgGFc-PROT", "IgD-PROT", "IgM-PROT") 

p2 <- VlnPlot(Moderna.integrated,  assay = "RNA", feature = rna.features, stack = TRUE) + NoLegend()
p3 <- VlnPlot(Moderna.integrated, assay = "ADT", feature = prot.features, stack = TRUE) + NoLegend()

pdf("Figures/Figure2C.pdf",width=10,height=4,useDingbats=FALSE)
print(p2+p3)
dev.off() 

#### Figure 2D
## Find all DEGs for each cluster 
moderna.markers <- FindAllMarkers(Moderna.integrated, assay = "RNA", only.pos = TRUE)

# Remove Ig genes and filter Top 10 DEGs
moderna.markers = moderna.markers[!grepl("IG",moderna.markers$gene),]

Top10 <- moderna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 

# Remove duplicate genes if needed 
Top10 <- Top10[-c(48,41,12),] #Check the list and make sure the numbers match.      

p4 <- DotPlot(Moderna.integrated, assay = "RNA", features = Top10$gene, cols = "RdBu", dot.scale = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  theme(axis.text.y = element_text(size = 10)) + coord_flip()
pdf("Figures/Figure2D.pdf",width=5,height=12,useDingbats=FALSE)
print(p4)
dev.off() 

#### Figure 2F
### Retrieve the Ids from the selected cells

Idents(Moderna.integrated) <- "orig.ident"

v1D14.S2pos.ids <- WhichCells(Moderna.integrated, idents = "S2pos.D1D14")
v2D6.S2pos.ids <- WhichCells(Moderna.integrated, idents = "S2pos.D2D6")
v2D9.S2pos.ids <- WhichCells(Moderna.integrated, idents = "S2pos.D2D9")
v2D14.S2pos.ids <- WhichCells(Moderna.integrated, idents = "S2pos.D2D14")
v2D28.S2pos.ids <- WhichCells(Moderna.integrated, idents = "S2pos.D2D28")
M6.S2pos.ids <- WhichCells(Moderna.integrated, idents = "S2pos.D2M6")

p5a <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v1D14.S2pos.ids) + NoLegend() + ggtitle("V1D14")
p5b <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v2D6.S2pos.ids) + NoLegend() + ggtitle("V2D6")
p5c <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v2D9.S2pos.ids) + NoLegend() + ggtitle("V2D9")
p5d <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v2D14.S2pos.ids) + NoLegend() + ggtitle("V2D14")
p5e <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = v2D28.S2pos.ids) + NoLegend() + ggtitle("V2D28")
p5f <- DimPlot(Moderna.integrated, reduction = "umap", cells.highlight = list("S-2P+ B cells" = M6.S2pos.ids)) + ggtitle("M6")

pdf("Figures/Figure2F.pdf",width=12,height=6,useDingbats=FALSE)
print(p5a+p5b+p5c+p5d+p5e+p5f)
dev.off() 

########### Analysis of switched MBC-C5 (re-clustering) 

Idents(Moderna.integrated) <- "clusters"

cluster5 <- subset(Moderna.integrated, idents = "Switched MBC-C5")

DefaultAssay(cluster5) <- "RNA"

# Normalize, scale and find clusters 
cluster5 <- NormalizeData(cluster5, normalization.method = "LogNormalize")
cluster5 <- ScaleData(cluster5, verbose = FALSE)
cluster5 <- RunPCA(cluster5, npcs = 15, verbose = FALSE)
cluster5 <- RunUMAP(cluster5, reduction = "pca", dims = 1:15)
cluster5 <- FindNeighbors(cluster5, reduction = "pca", dims = 1:15)
cluster5 <- FindClusters(cluster5, graph.name = "RNA_snn", resolution = 0.5)

# Data visualization

DimPlot(cluster5, reduction = "umap")

# Rename Idents based on gene expression analysis 
cluster5 <- RenameIdents(cluster5, `4` = "Atypical MBC-SC5.5")
cluster5 <- RenameIdents(cluster5, `3` = "Atypical MBC-SC5.4")
cluster5 <- RenameIdents(cluster5, `2` = "Activated MBC-SC5.3")
cluster5 <- RenameIdents(cluster5, `1` = "Resting MBC-SC5.2")
cluster5 <- RenameIdents(cluster5, `0` = "Mixed MBC-SC5.1")
cluster5$clusters <- Idents(cluster5)

saveRDS(cluster5, file = "cluster5.rds")

#### Figure 3A
Idents(cluster5) <- "clusters"
colors.cluster5 <- c("#999999", "#e41a1c", "#ff7f00", "#f781bf", "#4575b4")

p6 <- DimPlot(cluster5, reduction = "umap", cols = colors.cluster5) + 
  theme(legend.text = element_text(size = 16)) +
  theme(axis.text = element_text(size = 18)) +
  theme(text = element_text(size = 18)) 
p6$labels$colour <- "MBC-C5 subclusters"
p6
pdf("Figures/Figure3A.pdf",width=10,height=6,useDingbats=FALSE)
print(p6)
dev.off() 

#### Figure 3B
## Find all DEGs for each cluster 

C5.markers <- FindAllMarkers(cluster5, assay = "RNA", only.pos = FALSE)
C5.markers = C5.markers[!grepl("IG",C5.markers$gene),]

# Top 10 DEGs and Remove duplicate genes if needed
Top10.C5 <- C5.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
Top10.C5 <- Top10.C5[!duplicated(Top10.C5$gene), ] #remove duplicates if needed
genes <- Top10.C5$gene

#Order timepoints for good plotting
p7 <- dittoHeatmap(cluster5, assay = "RNA", slot = "data",order.by = "clusters", 
                   annot.by = c("clusters"), annot.colors = colors.cluster5,
                   heatmap.colors = colorRampPalette(c("darkblue", "white", "firebrick"))(100), 
                   breaks = seq(-3, 3, 0.06), genes = genes, scaled.to.max = FALSE)

pdf("Figures/Figure3B.pdf",width=10,height=10,useDingbats=FALSE)
print(p7)
dev.off() 

#### Figure 3C
# Normalize and scale ADT

DefaultAssay(cluster5) <- "ADT"

cluster5 <- NormalizeData(cluster5, normalization.method = 'CLR', scale.factor = 10000) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')

feat.adt.c5 <- c("CD19-PROT", "CD20-PROT", "CD21-PROT", "CD27-PROT", "CD38-PROT", "CD95-PROT", "CD11c-PROT", "CD71-PROT", "IgGFc-PROT", "IgD-PROT", "IgM-PROT")
feat.rna.c5 <- c("CD19", "MS4A1", "CD27", "CD38",  "FAS", "ITGAX", "TFRC", "IGHG1","IGHD", "IGHM", "IGHA1") 

p8 <- VlnPlot(cluster5, assay = "ADT", feature = feat.adt.c5, stack = TRUE) + NoLegend() + 
  theme(text = element_text(size = 14)) + theme(axis.text.x = element_text(size = 12))

p9 <- VlnPlot(cluster5, assay = "RNA", feature = feat.rna.c5, stack = TRUE) + NoLegend() +
  theme(text = element_text(size = 14)) + theme(axis.text.x = element_text(size = 11))

pdf("Figures/Figure3C.pdf",width=20,height=5,useDingbats=FALSE)
print(p8+p9)
dev.off() 

#### Figure 3D
### Retrieve the Ids from the selected cells/timepoints   
Idents(cluster5) <- "orig.ident"

v1D14.S2pos.ids.c5 <- WhichCells(cluster5, idents = "S2pos.D1D14")
v2D6.S2pos.ids.c5 <- WhichCells(cluster5, idents = "S2pos.D2D6")
v2D9.S2pos.ids.c5 <- WhichCells(cluster5, idents = "S2pos.D2D9")
v2D14.S2pos.ids.c5 <- WhichCells(cluster5, idents = "S2pos.D2D14")
v2D28.S2pos.ids.c5 <- WhichCells(cluster5, idents = "S2pos.D2D28")
M6.S2pos.ids.c5 <- WhichCells(cluster5, idents = "S2pos.D2M6")

p10a <- DimPlot(cluster5, reduction = "umap", cells.highlight = v1D14.S2pos.ids.c5, sizes.highlight = 0.2) + ggtitle("v1D14") + NoLegend() 
p10b <- DimPlot(cluster5, reduction = "umap", cells.highlight = v2D6.S2pos.ids.c5, sizes.highlight = 0.2) + ggtitle("v2D6") + NoLegend()
p10c <- DimPlot(cluster5, reduction = "umap", cells.highlight = v2D9.S2pos.ids.c5, sizes.highlight = 0.2) + ggtitle("v2D9") + NoLegend()
p10d <- DimPlot(cluster5, reduction = "umap", cells.highlight = v2D14.S2pos.ids.c5, sizes.highlight = 0.2) + ggtitle("v2D14") + NoLegend() 
p10e <- DimPlot(cluster5, reduction = "umap", cells.highlight = v2D28.S2pos.ids.c5, sizes.highlight = 0.2) + ggtitle("v2D28") + NoLegend() 
p10f <- DimPlot(cluster5, reduction = "umap", cells.highlight = list("S-2P+ B cells" = M6.S2pos.ids.c5), sizes.highlight = 0.2) + ggtitle("M6") 
                
pdf("Figures/Figure3D.pdf",width=12,height=6,useDingbats=FALSE)
print(p10a+p10b+p10c+p10d+p10e+p10f)
dev.off() 

#### Figure 3E
# Retrieve features count data from selected cells 
DefaultAssay(cluster5) <- "ADT"
Idents(cluster5) <- "timepoint"

S1.D28.c5.data <- FetchData(cluster5, "S1-PROT", slot = "counts", cells = v2D28.S2pos.ids.c5)
S1.M6.c5.data <- FetchData(cluster5, "S1-PROT", slot = "counts", cells = M6.S2pos.ids.c5)
RBD.D28.c5.data <- FetchData(cluster5, "RBD-PROT", slot = "counts", cells = v2D28.S2pos.ids.c5)
RBD.M6.c5.data <- FetchData(cluster5, "RBD-PROT", slot = "counts", cells = M6.S2pos.ids.c5)

s1.pos.d28.c5 <- filter(S1.D28.c5.data, `S1-PROT` >= 10)
s1.pos.m6.c5 <- filter(S1.M6.c5.data, `S1-PROT` >= 10)
rbd.pos.d28.c5 <- filter(RBD.D28.c5.data, `RBD-PROT` >= 10)
rbd.pos.m6.pos <- filter(RBD.M6.c5.data, `RBD-PROT` >= 10)

# Set the Ids of S1/RBD positive cells to highlight  

s1.d28.highlight <- factor(rownames(s1.pos.d28.c5))
s1.m6.highlight <- factor(rownames(s1.pos.m6.c5))
rbd.d28.highlight <- factor(rownames(rbd.pos.d28.c5))
rbd.m6.highlight <- factor(rownames(rbd.pos.m6.pos))

p11a <- DimPlot(cluster5, reduction = "umap", cells.highlight = list(v2D28 = rownames(s1.pos.d28.c5), M6 = rownames(s1.pos.m6.c5)), 
                sizes.highlight = 0.5, cols.highlight = c("#b2182b", "#2166ac", "lightgrey")) + NoLegend() + 
  theme(text = element_text(size = 18)) + ggtitle("S1+ cells") 

p11b <- DimPlot(cluster5, reduction = "umap", cells.highlight = list(v2D28 = rownames(rbd.pos.d28.c5), M6 = rownames(rbd.pos.m6.pos)), 
              sizes.highlight = 0.5, cols.highlight = c("#b2182b", "#2166ac", "lightgrey")) +
  theme(text = element_text(size = 18)) + ggtitle("RBD+ cells") 

pdf("Figures/Figure3E.pdf",width=12,height=4,useDingbats=FALSE)
print(p11a+p11b)
dev.off() 

#### Figure 3F

# load packages 

library(monocle3)
library(SeuratWrappers)

# Convert Seurat object into monocle dataset

cds <- as.cell_data_set(cluster5, assay = "RNA", reductions = "umap")
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Fit a principal graph within each partition
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, color_cells_by = "partition")

# Specify the root nodes of the trajectory graph
cds <- order_cells(cds, reduction_method = "UMAP")
p12 <- plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE, label_roots = TRUE, 
           graph_label_size = 7) +
  theme(text = element_text(size = 18)) + 
  theme(legend.text = element_text(size = 16)) +
  theme(axis.text = element_text(size = 18)) +
  theme(text = element_text(size = 18))
p12
pdf("Figures/Figure3F.pdf",width=5,height=3,useDingbats=FALSE)
print(p12)
dev.off() 

##Figure 2E
colors.moderna.7 <- c("#ED4648", "#A6780A", "#408002", "#1EAB6F", "#7F52FF", "#1A95E0", "#F218C3") ## All cells
colors.moderna.4 <- c("#408002", "#1EAB6F", "#7F52FF", "#1A95E0") ## Non-naive cells

fig2e <- table(Idents(Moderna.integrated), Moderna.integrated$Cells, Moderna.integrated$Donor, Moderna.integrated$timepoint)
df1 <- as.data.frame(fig2e)

## Rename variables and order the levels 
df1 <-  df1 %>% rename("Clusters" = "Var1", "Cells" = "Var2", 
                           "Donor" = "Var3", "Timepoint" = "Var4", 
                           "Frequency" = "Freq")

df1$Timepoint <- factor(df1$Timepoint, levels = c("v1D14", "v2D6", "v2D9", "v2D14", "v2D28", "v2M6"))
df1$Cells <- factor(df1$Cells, levels = c("Plasmablast", "S2P-", "S2P+"))
df1$Clusters <- factor(df1$Clusters, levels = c("Naive-C1", "Naive-C2", "Unswitched MBC-C3",
"Switched MBC-C4", "Switched MBC-C5", "PB-C6", "Naive-C7"))

df1 <- df1[!(df1$Cells == "Plasmablast" & df1$Timepoint == "v2D14"),]
df1 <- df1[!(df1$Cells == "Plasmablast" & df1$Timepoint == "v2D28"),]
df1 <- df1[!(df1$Cells == "Plasmablast" & df1$Timepoint == "v2M6"),]

df2 <- subset(df1, Clusters == "Unswitched MBC-C3" | Clusters == "Switched MBC-C4" | 
                Clusters == "Switched MBC-C5" | Clusters == "PB-C6")

### Figure 2E
# all cells
g <- ggplot(df1, aes(x=Timepoint, y=Frequency, fill=Clusters)) +
  geom_bar(position="fill", stat="identity")  + 
  scale_fill_manual(values = colors.moderna.7) +
  facet_grid(~Cells, scales ="free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size = 20))
pdf("Figures/Figure2E.pdf",width=10,height=4,useDingbats=FALSE)
print(g)
dev.off()
# non-naive B cells 
g <- ggplot(df2, aes(x=Timepoint, y=Frequency, fill=Clusters)) +
  geom_bar(position="fill", stat="identity")  + 
  scale_fill_manual(values = colors.moderna.4) +
  facet_grid(~Cells, scales ="free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size = 20))
pdf("Figures/Figure2E_non-naive.pdf",width=10,height=4,useDingbats=FALSE)
print(g)
dev.off()
## Figure S1B 
#all cells
g <- ggplot(df1, aes(x=Timepoint, y=Frequency, fill=Clusters)) +
  geom_bar(position="fill", stat="identity")  + 
  scale_fill_manual(values = colors.moderna.7) +
  facet_grid(Donor ~ Cells, scales ="free_x", ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size = 20))
pdf("Figures/FigureS2B.pdf",width=10,height=5,useDingbats=FALSE)
print(g)
dev.off()
# non-naive B cells
g <- ggplot(df2, aes(x=Timepoint, y=Frequency, fill=Clusters)) +
  geom_bar(position="fill", stat="identity")  + 
  scale_fill_manual(values = colors.moderna.4) +
  facet_grid(Donor ~ Cells, scales ="free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size = 20)) +
pdf("Figures/FigureS2B_non-naive.pdf",width=10,height=5,useDingbats=FALSE)
print(g)
dev.off()
##### Figure 4A and S2
# Load library

library(VennDiagram)

# Load clone analysis dataset

across_clones <- readRDS("all_public_analyzed_filtered_cloned_data.rds")
across_clone_clean <- across_clones %>% select(58,17,95,121,122,126,128,134,5) %>% #Select columns
  group_by(clone_id) %>%
  mutate(count = n()) #add clone_id counts column 

##Figure 4
# Generate three datasets containing counts of clones within each sorted cell   
set1 <- filter(across_clone_clean, group == "S-2P+") %>% select(7) 
set2 <- filter(across_clone_clean, group == "S-2P-") %>% select(7)
set3 <- filter(across_clone_clean, group == "Plasmablast") %>% select(7)
# Make a list of clone Ids
set1 <- as.vector(set1$clone_id, mode = "list")
set2 <- as.vector(set2$clone_id, mode = "list")
set3 <- as.vector(set3$clone_id, mode = "list")
# Main Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c(
    expression('S-2P'^'+'),
    expression('S-2P'^'-'),
    ("PBs")),
  filename = '#venn_diagramm.png',
  output=TRUE,
  # Output features
  imagetype="png",
  height = 480, 
  width = 480, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 3,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontfamily = "sans", 
  # Set names
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.pos = c(-10, 10, 180), # best setup for 3 sets. Position in degrees. 
  cat.dist = c(0.08, 0.08, 0.08),
  cat.fontfamily = "sans",
  rotation = 1,
  na = "none"
)

##Figure S2
# Generate datasets for each individual containing counts of clones within each sorted cell 
set4 <- filter(across_clone_clean, group == "S2P+" & Donor == "VAC716") %>% select(7)
set5 <- filter(across_clone_clean, group == "S2P-" & Donor == "VAC716") %>% select(7)
set6 <- filter(across_clone_clean, group == "Plasmablast" & Donor == "VAC716") %>% select(7)

set7 <- filter(across_clone_clean, group == "S2P+" & Donor == "VAC611") %>% select(7)
set8 <- filter(across_clone_clean, group == "S2P-" & Donor == "VAC611") %>% select(7)
set9 <- filter(across_clone_clean, group == "Plasmablast" & Donor == "VAC611") %>% select(7)

set10 <- filter(across_clone_clean, group == "S2P+" & Donor == "VAC003") %>% select(7)
set11 <- filter(across_clone_clean, group == "S2P-" & Donor == "VAC003") %>% select(7)
set12 <- filter(across_clone_clean, group == "Plasmablast" & Donor == "VAC003") %>% select(7)

# Make a list of clone Ids
set4 <- as.vector(set4$clone_id, mode = "list")
set5 <- as.vector(set5$clone_id, mode = "list")
set6 <- as.vector(set6$clone_id, mode = "list")

set7 <- as.vector(set7$clone_id, mode = "list")
set8 <- as.vector(set8$clone_id, mode = "list")
set9 <- as.vector(set9$clone_id, mode = "list")

set10 <- as.vector(set10$clone_id, mode = "list")
set11 <- as.vector(set11$clone_id, mode = "list")
set12 <- as.vector(set12$clone_id, mode = "list")

# Chart VAC716
venn.diagram(
  x = list(set4, set5, set6),
  category.names = c(
    expression('S-2P'^'+'),
    expression('S-2P'^'-'),
    ("PBs")),
  filename = '#VAC716 venn_diagramm.png',
  output=TRUE,
  # Output features
  imagetype="png",
  height = 480, 
  width = 480, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 3,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontfamily = "sans", 
  # Set names
  main = "VAC716",
  main.pos = c(0.1,1.2),
  main.cex = 0.8,
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.pos = c(-10, 10, 180), # best setup for 3 sets. Position in degrees. 
  cat.dist = c(0.08, 0.08, 0.08),
  cat.fontfamily = "sans",
  rotation = 1,
  na = "none"
)

# Chart VAC611
venn.diagram(
  x = list(set7, set8, set9),
  category.names = c(
    expression('S-2P'^'+'),
    expression('S-2P'^'-'),
    ("PBs")),
  filename = '#VAC611 venn_diagramm.png',
  output=TRUE,
  # Output features
  imagetype="png",
  height = 480, 
  width = 480, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 3,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontfamily = "sans", 
  # Set names
  main = "VAC611",
  main.pos = c(0.1,1.2),
  main.cex = 0.8,
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.pos = c(-10, 10, 180), # best setup for 3 sets. Position in degrees. 
  cat.dist = c(0.08, 0.08, 0.08),
  cat.fontfamily = "sans",
  rotation = 1,
  na = "none"
)

# Chart VAC003
venn.diagram(
  x = list(set10, set11, set12),
  category.names = c(
    expression('S-2P'^'+'),
    expression('S-2P'^'-'),
    ("PBs")),
  filename = '#VAC003 venn_diagramm.png',
  output=TRUE,
  # Output features
  imagetype="png",
  height = 480, 
  width = 480, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 3,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontfamily = "sans", 
  # Set names
  main = "VAC003",
  main.pos = c(0.1,1.2),
  main.cex = 0.8,
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.pos = c(-10, 10, 180), # best setup for 3 sets. Position in degrees. 
  cat.dist = c(0.08, 0.08, 0.08),
  cat.fontfamily = "sans",
  rotation = 1,
  na = "none"
)

##### END #####