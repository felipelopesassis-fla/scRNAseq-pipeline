---
Author: 'Felipe Assis, PhD'
title: 'scRNAseq - Integration and clustering'
date: '12-14-2022'
---

## Load Libraries 

library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)

## Merge Seurat Objects (MSO) 

MSO.merged <- merge(S2pos.v1D14, y = c(S2neg.v1D14, Plasmablast.v1D14, S2pos.v2D6, S2neg.v2D6, Plasmablast.v2D6,
                                            S2pos.v2D9, S2neg.v2D9, Plasmablast.v2D9, S2pos.v2D14, S2neg.v2D14, 
                                            S2pos.v2D28, S2neg.v2D28, S2pos.M6, S2neg.M6))

## Split the dataset into a list of Seurat objects

MSO.merged <- SplitObject(MSO.merged, split.by = "orig.ident") 

## Normalize and identify variable features for each dataset independently

MSO.list <- lapply(X = MSO.merged, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE, )
  x <- FindVariableFeatures(x, verbose = FALSE)
})

## Integration of large data set.
# Select features for downstream integration, and run PCA on each object in the list. 

features <- SelectIntegrationFeatures(object.list = Moderna.list, )
MSO.list <- lapply(X = MSO.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find anchors between the Seurat objects and create the "integrated" data assay.

anchors <- FindIntegrationAnchors(object.list = MSO.list, reduction = "rpca", dims = 1:50)
 
Moderna.integrated <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize", dims = 1:50)

### It is recommended removing immunoglobulin genes from the list of variable genes since they could influence the clustering due their 
### high expression in specific cells. For that, run the code below. Otherwise, proceed to ScaleData.

load("QC_features_meta.RData") 

biotypes_excl = unique(features_meta[["gene_biotype"]])[grepl(pattern="^IG_|^TR_", x=unique(features_meta[["gene_biotype"]]))]
remove.genes = features_meta[["external_gene_name"]][features_meta[["gene_biotype"]] %in% biotypes_excl]

Moderna.integrated <- NormalizeData(Moderna.integrated, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
Moderna.integrated <- FindVariableFeatures(Moderna.integrated, assay = "RNA", selection.method = "vst", nfeatures = 2000)

# Remove IG genes from list of variable genes
bool.remove.genes = Moderna.integrated@assays$RNA@var.features %in% remove.genes
Moderna.integrated@assays$RNA@var.features = Moderna.integrated@assays$RNA@var.features[!bool.remove.genes]

## Scaling, Dimensional Reduction and Clustering. 

DefaultAssay(Moderna.integrated) <- "RNA"
Moderna.integrated <- ScaleData(Moderna.integrated, verbose = FALSE)
Moderna.integrated <- RunPCA(Moderna.integrated, npcs = 50, verbose = FALSE) 
Moderna.integrated <- RunUMAP(Moderna.integrated, reduction = "pca", dims = 1:50)
Moderna.integrated <- FindNeighbors(Moderna.integrated, reduction = "pca", dims = 1:50)
Moderna.integrated <- FindClusters(Moderna.integrated, graph.name = "RNA_snn", resolution = 0.3)

## General visualization 

DimPlot(Moderna.integrated, reduction = "umap") + 
  theme(axis.text = element_text(size = 18)) +
  theme(text = element_text(size = 18)) 

saveRDS(Moderna.integrated, file = "Moderna.integrated.rds")

### END ###
