---
Author: 'Felipe Assis, PhD'
title: 'Pool 3 data pre-processing'
date: '12-14-2022'
---

## Load Libraries
library("Seurat") #load Seurat 3.1
library("dplyr")
library("matrixStats")
library('tidyverse')

### Import the data

M1C1.data <- Read10X(data.dir = "Moderna6mo_multi1/outs/multi/count/raw_feature_bc_matrix")
M1C2.data <- Read10X(data.dir = "Moderna6mo_multi2/outs/multi/count/raw_feature_bc_matrix")
M1C3.data <- Read10X(data.dir = "Moderna6mo_multi3/outs/multi/count/raw_feature_bc_matrix")
M1C4.data <- Read10X(data.dir = "Moderna6mo_multi4/outs/multi/count/raw_feature_bc_matrix")

## Initialize the Seurat object with the raw (non-normalized data).

M1C1 <- CreateSeuratObject(counts = M1C1.data$`Gene Expression`, project = "M1C1")
M1C2 <- CreateSeuratObject(counts = M1C2.data$`Gene Expression`, project = "M1C2")
M1C3 <- CreateSeuratObject(counts = M1C3.data$`Gene Expression`, project = "M1C3")
M1C4 <- CreateSeuratObject(counts = M1C4.data$`Gene Expression`, project = "M1C4")

## Prepare ADT
M1C1[["ADT"]] <- CreateAssayObject(counts = M1C1.data$`Antibody Capture`)
M1C2[["ADT"]] <- CreateAssayObject(counts = M1C2.data$`Antibody Capture`)
M1C3[["ADT"]] <- CreateAssayObject(counts = M1C3.data$`Antibody Capture`)
M1C4[["ADT"]] <- CreateAssayObject(counts = M1C4.data$`Antibody Capture`)

## Add metadata to distinguish each object 
# Batch
M1C1$Batch <- rep("M1C1", length(colnames(M1C1)))
M1C2$Batch <- rep("M1C2", length(colnames(M1C2)))
M1C3$Batch <- rep("M1C3", length(colnames(M1C3)))
M1C4$Batch <- rep("M1C4", length(colnames(M1C4)))

# Orig.ident
M1C1$orig.ident <- "S2pos.v2M6"
M1C2$orig.ident <- "S2pos.v2M6"
M1C3$orig.ident <- "S2neg.v2M6"
M1C4$orig.ident <- "S2neg.v2M6"

## Rename Barcode for each object to match barcodes in SNP file

M1C1 <- RenameCells(M1C1, new.names = paste(substr(colnames(M1C1), start = 1, stop = 17),"1", sep = ""))
M1C2 <- RenameCells(M1C2, new.names = paste(substr(colnames(M1C2), start = 1, stop = 17),"2", sep = ""))
M1C3 <- RenameCells(M1C3, new.names = paste(substr(colnames(M1C3), start = 1, stop = 17),"3", sep = ""))
M1C4 <- RenameCells(M1C4, new.names = paste(substr(colnames(M1C4), start = 1, stop = 17),"4", sep = ""))

## Add new metadata (SNP <- DEMUXLET) 

M1C1_demuxbestList = list()
for(i in 1:length("1")){
  M1C1_demuxbestList[[i]] = read.table(paste("v2multi_config_1.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C1_demuxbestList[[i]]$NewBarcode = paste(substr(M1C1_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"1", sep = "")
}
M1C1_demuxbestdf <- plyr::ldply(M1C1_demuxbestList, data.frame)
length(which(colnames(M1C1) %in% M1C1_demuxbestdf$NewBarcode))
setdiff(colnames(M1C1), M1C1_demuxbestdf$NewBarcode)
rownames(M1C1_demuxbestdf) <- M1C1_demuxbestdf$NewBarcode
M1C1 <- AddMetaData(M1C1, metadata = M1C1_demuxbestdf[colnames(M1C1),])

M1C2_demuxbestList = list()
for(i in 1:length("1")){
  M1C2_demuxbestList[[i]] = read.table(paste("v2multi_config_2.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C2_demuxbestList[[i]]$NewBarcode = paste(substr(M1C2_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"2", sep = "")
}
M1C2_demuxbestdf <- plyr::ldply(M1C2_demuxbestList, data.frame)
length(which(colnames(M1C2) %in% M1C2_demuxbestdf$NewBarcode))
setdiff(colnames(M1C2), M1C2_demuxbestdf$NewBarcode)
rownames(M1C2_demuxbestdf) <- M1C2_demuxbestdf$NewBarcode
M1C2 <- AddMetaData(M1C2, metadata = M1C2_demuxbestdf[colnames(M1C2),])

M1C3_demuxbestList = list()
for(i in 1:length("1")){
  M1C3_demuxbestList[[i]] = read.table(paste("v2multi_config_3.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C3_demuxbestList[[i]]$NewBarcode = paste(substr(M1C3_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"3", sep = "")
}
M1C3_demuxbestdf <- plyr::ldply(M1C3_demuxbestList, data.frame)
length(which(colnames(M1C3) %in% M1C3_demuxbestdf$NewBarcode))
setdiff(colnames(M1C3), M1C3_demuxbestdf$NewBarcode)
rownames(M1C3_demuxbestdf) <- M1C3_demuxbestdf$NewBarcode
M1C3 <- AddMetaData(M1C3, metadata = M1C3_demuxbestdf[colnames(M1C3),])

M1C4_demuxbestList = list()
for(i in 1:length("1")){
  M1C4_demuxbestList[[i]] = read.table(paste("v2multi_config_4.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C4_demuxbestList[[i]]$NewBarcode = paste(substr(M1C4_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"4", sep = "")
}
M1C4_demuxbestdf <- plyr::ldply(M1C4_demuxbestList, data.frame)
length(which(colnames(M1C4) %in% M1C4_demuxbestdf$NewBarcode))
setdiff(colnames(M1C4), M1C4_demuxbestdf$NewBarcode)
rownames(M1C4_demuxbestdf) <- M1C4_demuxbestdf$NewBarcode
M1C4 <- AddMetaData(M1C4, metadata = M1C4_demuxbestdf[colnames(M1C4),])

## Rename the cells again (optional), to prepare cells to incorporate vdj 

M1C1 <- RenameCells(M1C1, new.names = paste(substr(colnames(M1C1), start = 1, stop = 17),"<sample1><pool3>", sep = ""))
M1C2 <- RenameCells(M1C2, new.names = paste(substr(colnames(M1C2), start = 1, stop = 17),"<sample2><pool3>", sep = ""))
M1C3 <- RenameCells(M1C3, new.names = paste(substr(colnames(M1C3), start = 1, stop = 17),"<sample3><pool3>", sep = ""))
M1C4 <- RenameCells(M1C4, new.names = paste(substr(colnames(M1C4), start = 1, stop = 17),"<sample4><pool3>", sep = ""))

## Merge objects

S2pos.M6 <- merge(M1C1, y = M1C2)
S2neg.M6 <- merge(M1C3, y = M1C4)

## Remove abmormally high count cells

S2pos.M6[["percent.mt"]] <- PercentageFeatureSet(S2pos.v2M6, pattern = "^MT-")
S2neg.M6[["percent.mt"]] <- PercentageFeatureSet(S2neg.v2M6, pattern = "^MT-")

S2pos.M6 <- subset(S2pos.M6, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                       nFeature_RNA < 4000 & percent.mt < 10)
S2neg.M6 <- subset(S2neg.M6, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                       nFeature_RNA < 4000 & percent.mt < 10)

## Adittional metadata

S2pos.M6$Donor = sapply(strsplit(as.character(S2pos.M6$BEST.GUESS),split = ","),'[',1)
S2pos.M6$Donor = sapply(strsplit(as.character(S2pos.M6$Donor),split = "_"),'[',1)
S2pos.M6$Sample = paste(S2pos.M6$Batch, S2pos.M6$Donor)
S2pos.M6$timepoint <- "M6"
S2pos.M6$Cells <- "S2P+"

S2neg.M6$Donor = sapply(strsplit(as.character(S2neg.M6$BEST.GUESS),split = ","),'[',1)
S2neg.M6$Donor = sapply(strsplit(as.character(S2neg.M6$Donor),split = "_"),'[',1)
S2neg.M6$Sample = paste(S2neg.M6$Batch, S2neg.M6$Donor)
S2neg.M6$timepoint <- "M6"
S2neg.M6$Cells <- "S2P-"

## Save final files 

saveRDS(S2pos.M6, file = "S2pos.M6.rds")
saveRDS(S2neg.M6, file = "S2neg.M6.rds")

### End ###
