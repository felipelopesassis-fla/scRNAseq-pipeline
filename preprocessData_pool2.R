#this is tested using R 3.6.1 on a high-performance comupting node with 16 cores and at least 160 gb or ram. 
library("Seurat") #load Seurat 3.1
library("dplyr")
library("matrixStats")
library('tidyverse')

### Import the data
M1C1.data <- Read10X(data.dir = "Moderna_P2_multi1/outs/multi/count/raw_feature_bc_matrix")
M1C2.data <- Read10X(data.dir = "Moderna_P2_multi2/outs/multi/count/raw_feature_bc_matrix")
M1C3.data <- Read10X(data.dir = "Moderna_P2_multi3/outs/multi/count/raw_feature_bc_matrix")
M1C4.data <- Read10X(data.dir = "Moderna_P2_multi4/outs/multi/count/raw_feature_bc_matrix")
M1C5.data <- Read10X(data.dir = "Moderna_P2_multi5/outs/multi/count/raw_feature_bc_matrix")
M1C6.data <- Read10X(data.dir = "Moderna_P2_multi6/outs/multi/count/raw_feature_bc_matrix")
M1C7.data <- Read10X(data.dir = "Moderna_P2_multi7/outs/multi/count/raw_feature_bc_matrix")
M1C8.data <- Read10X(data.dir = "Moderna_P2_multi8/outs/multi/count/raw_feature_bc_matrix")
M1C9.data <- Read10X(data.dir = "Moderna_P2_multi9/outs/multi/count/raw_feature_bc_matrix")
M1C10.data <- Read10X(data.dir = "Moderna_P2_multi10/outs/multi/count/raw_feature_bc_matrix")
M1C11.data <- Read10X(data.dir = "Moderna_P2_multi11/outs/multi/count/raw_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data)
M1C1 <- CreateSeuratObject(counts = M1C1.data$`Gene Expression`, project = "M1C1", assay = "RNA", min.feature = 5)
M1C2 <- CreateSeuratObject(counts = M1C2.data$`Gene Expression`, project = "M1C2", assay = "RNA", min.feature = 5)
M1C3 <- CreateSeuratObject(counts = M1C3.data$`Gene Expression`, project = "M1C3", assay = "RNA", min.feature = 5)
M1C4 <- CreateSeuratObject(counts = M1C4.data$`Gene Expression`, project = "M1C4", assay = "RNA", min.feature = 5)
M1C5 <- CreateSeuratObject(counts = M1C5.data$`Gene Expression`, project = "M1C5", assay = "RNA", min.feature = 5)
M1C6 <- CreateSeuratObject(counts = M1C6.data$`Gene Expression`, project = "M1C6", assay = "RNA", min.feature = 5)
M1C7 <- CreateSeuratObject(counts = M1C7.data$`Gene Expression`, project = "M1C7", assay = "RNA", min.feature = 5)
M1C8 <- CreateSeuratObject(counts = M1C8.data$`Gene Expression`, project = "M1C8", assay = "RNA", min.feature = 5)
M1C9 <- CreateSeuratObject(counts = M1C9.data$`Gene Expression`, project = "M1C9", assay = "RNA", min.feature = 5)
M1C10 <- CreateSeuratObject(counts = M1C10.data$`Gene Expression`, project = "M1C10", assay = "RNA", min.feature = 5)
M1C11 <- CreateSeuratObject(counts = M1C11.data$`Gene Expression`, project = "M1C11", assay = "RNA", min.feature = 5)

#########Preparing ADT
M1C1[["ADT"]] <- CreateAssayObject(counts = M1C1.data$`Antibody Capture`[1:19,colnames(M1C1)])
M1C2[["ADT"]] <- CreateAssayObject(counts = M1C2.data$`Antibody Capture`[1:19,colnames(M1C2)])
M1C3[["ADT"]] <- CreateAssayObject(counts = M1C3.data$`Antibody Capture`[1:19,colnames(M1C3)])
M1C4[["ADT"]] <- CreateAssayObject(counts = M1C4.data$`Antibody Capture`[1:19,colnames(M1C4)])
M1C5[["ADT"]] <- CreateAssayObject(counts = M1C5.data$`Antibody Capture`[1:19,colnames(M1C5)])
M1C6[["ADT"]] <- CreateAssayObject(counts = M1C6.data$`Antibody Capture`[1:19,colnames(M1C6)])
M1C7[["ADT"]] <- CreateAssayObject(counts = M1C7.data$`Antibody Capture`[1:19,colnames(M1C7)])
M1C8[["ADT"]] <- CreateAssayObject(counts = M1C8.data$`Antibody Capture`[1:19,colnames(M1C8)])
M1C9[["ADT"]] <- CreateAssayObject(counts = M1C9.data$`Antibody Capture`[1:19,colnames(M1C9)])
M1C10[["ADT"]] <- CreateAssayObject(counts = M1C10.data$`Antibody Capture`[1:19,colnames(M1C10)])
M1C11[["ADT"]] <- CreateAssayObject(counts = M1C11.data$`Antibody Capture`[1:19,colnames(M1C11)])

#########Preparing HTO
M1C1[["HTO"]] <- CreateAssayObject(counts = M1C1.data$`Antibody Capture`[20:21,colnames(M1C1)])
M1C2[["HTO"]] <- CreateAssayObject(counts = M1C2.data$`Antibody Capture`[20:21,colnames(M1C2)])
M1C3[["HTO"]] <- CreateAssayObject(counts = M1C3.data$`Antibody Capture`[20:21,colnames(M1C3)])
M1C4[["HTO"]] <- CreateAssayObject(counts = M1C4.data$`Antibody Capture`[20:21,colnames(M1C4)])
M1C5[["HTO"]] <- CreateAssayObject(counts = M1C5.data$`Antibody Capture`[20:21,colnames(M1C5)])
M1C6[["HTO"]] <- CreateAssayObject(counts = M1C6.data$`Antibody Capture`[20:21,colnames(M1C6)])
M1C7[["HTO"]] <- CreateAssayObject(counts = M1C7.data$`Antibody Capture`[20:22,colnames(M1C7)])
M1C8[["HTO"]] <- CreateAssayObject(counts = M1C8.data$`Antibody Capture`[20:22,colnames(M1C8)])
M1C9[["HTO"]] <- CreateAssayObject(counts = M1C9.data$`Antibody Capture`[20:22,colnames(M1C9)])
M1C10[["HTO"]] <- CreateAssayObject(counts = M1C10.data$`Antibody Capture`[20:22,colnames(M1C10)])
M1C11[["HTO"]] <- CreateAssayObject(counts = M1C11.data$`Antibody Capture`[20:22,colnames(M1C11)])

######## Add metadata to distinguish each object 
#batch
M1C1$Batch <- rep("M1C1", length(colnames(M1C1)))
M1C2$Batch <- rep("M1C2", length(colnames(M1C2)))
M1C3$Batch <- rep("M1C3", length(colnames(M1C3)))
M1C4$Batch <- rep("M1C4", length(colnames(M1C4)))
M1C5$Batch <- rep("M1C5", length(colnames(M1C5)))
M1C6$Batch <- rep("M1C6", length(colnames(M1C6)))
M1C7$Batch <- rep("M1C7", length(colnames(M1C7)))
M1C8$Batch <- rep("M1C8", length(colnames(M1C8)))
M1C9$Batch <- rep("M1C9", length(colnames(M1C9)))
M1C10$Batch <- rep("M1C10", length(colnames(M1C10)))
M1C11$Batch <- rep("M1C11", length(colnames(M1C11)))

#orig.ident
M1C1$orig.ident <- "S2pos.v2D6"
M1C2$orig.ident <- "S2pos.v2D6"
M1C3$orig.ident <- "S2neg.v2D6"
M1C4$orig.ident <- "S2neg.v2D6"
M1C5$orig.ident <- "Plasmablast.v2D6"
M1C6$orig.ident <- "Plasmablast.v2D6"
M1C7$orig.ident <- "S2pos.v1D14"
M1C8$orig.ident <- "S2neg.v1D14"
M1C9$orig.ident <- "S2neg.v1D14"
M1C10$orig.ident <- "Plasmablast.v1D14"
M1C11$orig.ident <- "Plasmablast.v1D14"

#########Rename Barcode for each object to match barcodes in SNP file

M1C1 <- RenameCells(M1C1, new.names = paste(substr(colnames(M1C1), start = 1, stop = 17),"1", sep = ""))
M1C2 <- RenameCells(M1C2, new.names = paste(substr(colnames(M1C2), start = 1, stop = 17),"2", sep = ""))
M1C3 <- RenameCells(M1C3, new.names = paste(substr(colnames(M1C3), start = 1, stop = 17),"3", sep = ""))
M1C4 <- RenameCells(M1C4, new.names = paste(substr(colnames(M1C4), start = 1, stop = 17),"4", sep = ""))
M1C5 <- RenameCells(M1C5, new.names = paste(substr(colnames(M1C5), start = 1, stop = 17),"5", sep = ""))
M1C6 <- RenameCells(M1C6, new.names = paste(substr(colnames(M1C6), start = 1, stop = 17),"6", sep = ""))
M1C7 <- RenameCells(M1C7, new.names = paste(substr(colnames(M1C7), start = 1, stop = 17),"7", sep = ""))
M1C8 <- RenameCells(M1C8, new.names = paste(substr(colnames(M1C8), start = 1, stop = 17),"8", sep = ""))
M1C9 <- RenameCells(M1C9, new.names = paste(substr(colnames(M1C9), start = 1, stop = 17),"9", sep = ""))
M1C10 <- RenameCells(M1C10, new.names = paste(substr(colnames(M1C10), start = 1, stop = 17),"10", sep = ""))
M1C11 <- RenameCells(M1C11, new.names = paste(substr(colnames(M1C11), start = 1, stop = 17),"11", sep = ""))

##### Add new metadata (SNP <- DEMUXLET) 

M1C1_demuxbestList = list()
for(i in 1:length("1")){
  M1C1_demuxbestList[[i]] = read.table(paste("v3multi_config_1.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C2_demuxbestList[[i]] = read.table(paste("v3multi_config_2.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C3_demuxbestList[[i]] = read.table(paste("v3multi_config_3.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C4_demuxbestList[[i]] = read.table(paste("v3multi_config_4.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C4_demuxbestList[[i]]$NewBarcode = paste(substr(M1C4_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"4", sep = "")
}
M1C4_demuxbestdf <- plyr::ldply(M1C4_demuxbestList, data.frame)
length(which(colnames(M1C4) %in% M1C4_demuxbestdf$NewBarcode))
setdiff(colnames(M1C4), M1C4_demuxbestdf$NewBarcode)
rownames(M1C4_demuxbestdf) <- M1C4_demuxbestdf$NewBarcode
M1C4 <- AddMetaData(M1C4, metadata = M1C4_demuxbestdf[colnames(M1C4),])

M1C5_demuxbestList = list()
for(i in 1:length("1")){
  M1C5_demuxbestList[[i]] = read.table(paste("v3multi_config_5.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C5_demuxbestList[[i]]$NewBarcode = paste(substr(M1C5_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"5", sep = "")
}
M1C5_demuxbestdf <- plyr::ldply(M1C5_demuxbestList, data.frame)
length(which(colnames(M1C5) %in% M1C5_demuxbestdf$NewBarcode))
setdiff(colnames(M1C5), M1C5_demuxbestdf$NewBarcode)
rownames(M1C5_demuxbestdf) <- M1C5_demuxbestdf$NewBarcode
M1C5 <- AddMetaData(M1C5, metadata = M1C5_demuxbestdf[colnames(M1C5),])

M1C6_demuxbestList = list()
for(i in 1:length("1")){
  M1C6_demuxbestList[[i]] = read.table(paste("v3multi_config_6.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C6_demuxbestList[[i]]$NewBarcode = paste(substr(M1C6_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"6", sep = "")
}
M1C6_demuxbestdf <- plyr::ldply(M1C6_demuxbestList, data.frame)
length(which(colnames(M1C6) %in% M1C6_demuxbestdf$NewBarcode))
setdiff(colnames(M1C6), M1C6_demuxbestdf$NewBarcode)
rownames(M1C6_demuxbestdf) <- M1C6_demuxbestdf$NewBarcode
M1C6 <- AddMetaData(M1C6, metadata = M1C6_demuxbestdf[colnames(M1C6),])

M1C7_demuxbestList = list()
for(i in 1:length("1")){
  M1C7_demuxbestList[[i]] = read.table(paste("v3multi_config_7.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C7_demuxbestList[[i]]$NewBarcode = paste(substr(M1C7_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"7", sep = "")
}
M1C7_demuxbestdf <- plyr::ldply(M1C7_demuxbestList, data.frame)
length(which(colnames(M1C7) %in% M1C7_demuxbestdf$NewBarcode))
setdiff(colnames(M1C7), M1C7_demuxbestdf$NewBarcode)
rownames(M1C7_demuxbestdf) <- M1C7_demuxbestdf$NewBarcode
M1C7 <- AddMetaData(M1C7, metadata = M1C7_demuxbestdf[colnames(M1C7),])

M1C8_demuxbestList = list()
for(i in 1:length("1")){
  M1C8_demuxbestList[[i]] = read.table(paste("v3multi_config_8.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C8_demuxbestList[[i]]$NewBarcode = paste(substr(M1C8_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"8", sep = "")
}
M1C8_demuxbestdf <- plyr::ldply(M1C8_demuxbestList, data.frame)
length(which(colnames(M1C8) %in% M1C8_demuxbestdf$NewBarcode))
setdiff(colnames(M1C8), M1C8_demuxbestdf$NewBarcode)
rownames(M1C8_demuxbestdf) <- M1C8_demuxbestdf$NewBarcode
M1C8 <- AddMetaData(M1C8, metadata = M1C8_demuxbestdf[colnames(M1C8),])

M1C9_demuxbestList = list()
for(i in 1:length("1")){
  M1C9_demuxbestList[[i]] = read.table(paste("v3multi_config_9.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C9_demuxbestList[[i]]$NewBarcode = paste(substr(M1C9_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"9", sep = "")
}
M1C9_demuxbestdf <- plyr::ldply(M1C9_demuxbestList, data.frame)
length(which(colnames(M1C9) %in% M1C9_demuxbestdf$NewBarcode))
setdiff(colnames(M1C9), M1C9_demuxbestdf$NewBarcode)
rownames(M1C9_demuxbestdf) <- M1C9_demuxbestdf$NewBarcode
M1C9 <- AddMetaData(M1C9, metadata = M1C9_demuxbestdf[colnames(M1C9),])

M1C10_demuxbestList = list()
for(i in 1:length("1")){
  M1C10_demuxbestList[[i]] = read.table(paste("v3multi_config_10.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C10_demuxbestList[[i]]$NewBarcode = paste(substr(M1C10_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"10", sep = "")
}
M1C10_demuxbestdf <- plyr::ldply(M1C10_demuxbestList, data.frame)
length(which(colnames(M1C10) %in% M1C10_demuxbestdf$NewBarcode))
setdiff(colnames(M1C10), M1C10_demuxbestdf$NewBarcode)
rownames(M1C10_demuxbestdf) <- M1C10_demuxbestdf$NewBarcode
M1C10 <- AddMetaData(M1C10, metadata = M1C10_demuxbestdf[colnames(M1C10),])

M1C11_demuxbestList = list()
for(i in 1:length("1")){
  M1C11_demuxbestList[[i]] = read.table(paste("v3multi_config_11.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C11_demuxbestList[[i]]$NewBarcode = paste(substr(M1C11_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"11", sep = "")
}
M1C11_demuxbestdf <- plyr::ldply(M1C11_demuxbestList, data.frame)
length(which(colnames(M1C11) %in% M1C11_demuxbestdf$NewBarcode))
setdiff(colnames(M1C11), M1C11_demuxbestdf$NewBarcode)
rownames(M1C11_demuxbestdf) <- M1C11_demuxbestdf$NewBarcode
M1C11 <- AddMetaData(M1C11, metadata = M1C11_demuxbestdf[colnames(M1C11),])

#rename the cells again (optional), to prepare cells to incorporate vdj 

M1C1 <- RenameCells(M1C1, new.names = paste(substr(colnames(M1C1), start = 1, stop = 17),"<sample1><pool2>", sep = ""))
M1C2 <- RenameCells(M1C2, new.names = paste(substr(colnames(M1C2), start = 1, stop = 17),"<sample2><pool2>", sep = ""))
M1C3 <- RenameCells(M1C3, new.names = paste(substr(colnames(M1C3), start = 1, stop = 17),"<sample3><pool2>", sep = ""))
M1C4 <- RenameCells(M1C4, new.names = paste(substr(colnames(M1C4), start = 1, stop = 17),"<sample4><pool2>", sep = ""))
M1C5 <- RenameCells(M1C5, new.names = paste(substr(colnames(M1C5), start = 1, stop = 17),"<sample5><pool2>", sep = ""))
M1C6 <- RenameCells(M1C6, new.names = paste(substr(colnames(M1C6), start = 1, stop = 17),"<sample6><pool2>", sep = ""))
M1C7 <- RenameCells(M1C7, new.names = paste(substr(colnames(M1C7), start = 1, stop = 17),"<sample7><pool2>", sep = ""))
M1C8 <- RenameCells(M1C8, new.names = paste(substr(colnames(M1C8), start = 1, stop = 17),"<sample8><pool2>", sep = ""))
M1C9 <- RenameCells(M1C9, new.names = paste(substr(colnames(M1C9), start = 1, stop = 17),"<sample9><pool2>", sep = ""))
M1C10 <- RenameCells(M1C10, new.names = paste(substr(colnames(M1C10), start = 1, stop = 17),"<sample10><pool2>", sep = ""))
M1C11 <- RenameCells(M1C11, new.names = paste(substr(colnames(M1C11), start = 1, stop = 17),"<sample11><pool2>", sep = ""))

## use HTOdemux to identify negatives droplets to use for DSB norm
# merge objects
S2pos.v2D6 <- merge(M1C1, y = M1C2)
S2neg.v2D6 <- merge(M1C3, y = M1C4)
Plasmablast.v2D6 <- merge(M1C5, y = M1C6)
S2pos.v1D14 <- M1C7
S2neg.v1D14<- merge(M1C8, y = M1C9)
Plasmablast.v1D14 <- merge(M1C10, y = M1C11)

# remove abnormally high count cells
S2pos.v2D6 <- subset(S2pos.v2D6, subset = nCount_ADT < 30000 & nCount_HTO < 15001)
S2neg.v2D6 <- subset(S2neg.v2D6, subset = nCount_ADT < 30000 & nCount_HTO < 15001)
Plasmablast.v2D6 <- subset(Plasmablast.v2D6, subset = nCount_ADT < 30000 & nCount_HTO < 15001)
S2pos.v1D14 <- subset(S2pos.v1D14, subset = nCount_ADT < 30000 & nCount_HTO < 15001)
S2neg.v1D14 <- subset(S2neg.v1D14, subset = nCount_ADT < 30000 & nCount_HTO < 15001)
Plasmablast.v1D14 <- subset(Plasmablast.v1D14, subset = nCount_ADT < 30000 & nCount_HTO < 15001)

S2pos.v2D6[["percent.mt"]] <- PercentageFeatureSet(S2pos.v2D6, pattern = "^MT-")
S2neg.v2D6[["percent.mt"]] <- PercentageFeatureSet(S2neg.v2D6, pattern = "^MT-")
Plasmablast.v2D6[["percent.mt"]] <- PercentageFeatureSet(Plasmablast.v2D6, pattern = "^MT-")
S2pos.v1D14[["percent.mt"]] <- PercentageFeatureSet(S2pos.v1D14, pattern = "^MT-")
S2neg.v1D14[["percent.mt"]] <- PercentageFeatureSet(S2neg.v1D14, pattern = "^MT-")
Plasmablast.v1D14[["percent.mt"]] <- PercentageFeatureSet(Plasmablast.v1D14, pattern = "^MT-")

### HTO demultiplex

S2pos.v2D6 <- NormalizeData(S2pos.v2D6, assay = "HTO", normalization.method = "CLR", margin = 2)
S2pos.v2D6 <- ScaleData(S2pos.v2D6, assay = "HTO", model.use = "linear")
S2pos.v2D6 = HTODemux(S2pos.v2D6)
Idents(S2pos.v2D6) <- "hash.ID"
S2pos.v2D6 <- subset(S2pos.v2D6, idents = "HTO9-PROT")
S2pos.v2D6 <- subset(S2pos.v2D6, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                       nFeature_RNA < 4000 & percent.mt < 10)

S2neg.v2D6 <- NormalizeData(S2neg.v2D6, assay = "HTO", normalization.method = "CLR", margin = 2)
S2neg.v2D6 <- ScaleData(S2neg.v2D6, assay = "HTO", model.use = "linear")
S2neg.v2D6 = HTODemux(S2neg.v2D6)
Idents(S2neg.v2D6) <- "hash.ID"
S2neg.v2D6 <- subset(S2neg.v2D6, idents = c("HTO9-PROT"))
S2neg.v2D6 <- subset(S2neg.v2D6, subset = DROPLET.TYPE == "SNG" & 
                       nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

Plasmablast.v2D6 <- NormalizeData(Plasmablast.v2D6, assay = "HTO", normalization.method = "CLR", margin = 2)
Plasmablast.v2D6 <- ScaleData(Plasmablast.v2D6, assay = "HTO", model.use = "linear")
Plasmablast.v2D6 = HTODemux(Plasmablast.v2D6)
Idents(Plasmablast.v2D6) <- "hash.ID"
Plasmablast.v2D6 <- subset(Plasmablast.v2D6, idents = c("HTO9-PROT"))
Plasmablast.v2D6 <- subset(Plasmablast.v2D6, subset = DROPLET.TYPE == "SNG" &
                             nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

S2pos.v1D14 <- NormalizeData(S2pos.v1D14, assay = "HTO", normalization.method = "CLR", margin = 2)
S2pos.v1D14 <- ScaleData(S2pos.v1D14, assay = "HTO", model.use = "linear")
S2pos.v1D14 = HTODemux(S2pos.v1D14)
Idents(S2pos.v1D14) <- "hash.ID"
S2pos.v1D14 <- subset(S2pos.v1D14, idents = c("HTO3-PROT", "HTO4-PROT"))
S2pos.v1D14 <- subset(S2pos.v1D14, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                        nFeature_RNA < 4000 & percent.mt < 10)

S2neg.v1D14 <- NormalizeData(S2neg.v1D14, assay = "HTO", normalization.method = "CLR", margin = 2)
S2neg.v1D14 <- ScaleData(S2neg.v1D14, assay = "HTO", model.use = "linear")
S2neg.v1D14 = HTODemux(S2neg.v1D14)
Idents(S2neg.v1D14) <- "hash.ID"
S2neg.v1D14 <- subset(S2neg.v1D14, idents = c("HTO3-PROT", "HTO4-PROT"))
S2neg.v1D14 <- subset(S2neg.v1D14, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                        nFeature_RNA < 4000 & percent.mt < 10)

Plasmablast.v1D14 <- NormalizeData(Plasmablast.v1D14, assay = "HTO", normalization.method = "CLR", margin = 2)
Plasmablast.v1D14 <- ScaleData(Plasmablast.v1D14, assay = "HTO", model.use = "linear")
Plasmablast.v1D14 = MULTIseqDemux(Plasmablast.v1D14, assay = "HTO", quantile = 0.5, autoThresh = FALSE)
table(Plasmablast.v1D14$MULTI_ID)
Idents(Plasmablast.v1D14) <- "MULTI_ID"
Plasmablast.v1D14 <- subset(Plasmablast.v1D14, idents = c("HTO3-PROT", "HTO4-PROT"))
Plasmablast.v1D14 <- subset(Plasmablast.v1D14, subset = DROPLET.TYPE == "SNG" & 
                              nFeature_RNA > 100 & nFeature_RNA < 4000 & percent.mt < 10)

## add metadata
S2pos.v2D6$Donor = sapply(strsplit(as.character(S2pos.v2D6$BEST.GUESS),split = ","),'[',1)
S2pos.v2D6$Donor = sapply(strsplit(as.character(S2pos.v2D6$Donor),split = "_"),'[',1)
S2pos.v2D6$Sample = paste(S2pos.v2D6$Batch, S2pos.v2D6$Donor)
S2pos.v2D6$timepoint = "v2D6"
S2pos.v2D6$Cells <- "S2P+"

S2neg.v2D6$Donor = sapply(strsplit(as.character(S2neg.v2D6$BEST.GUESS),split = ","),'[',1)
S2neg.v2D6$Donor = sapply(strsplit(as.character(S2neg.v2D6$Donor),split = "_"),'[',1)
S2neg.v2D6$Sample = paste(S2neg.v2D6$Batch, S2neg.v2D6$Donor)
S2neg.v2D6$timepoint = "v2D6"
S2neg.v2D6$Cells <- "S2P-"

S2pos.v1D14$Donor = sapply(strsplit(as.character(S2pos.v1D14$BEST.GUESS),split = ","),'[',1)
S2pos.v1D14$Donor = sapply(strsplit(as.character(S2pos.v1D14$Donor),split = "_"),'[',1)
S2pos.v1D14$Sample = paste(S2pos.v1D14$Batch, S2pos.v1D14$Donor)
S2pos.v1D14$timepoint = "v1D14"
S2pos.v1D14$Cells <- "S2P+"

S2neg.v1D14$Donor = sapply(strsplit(as.character(S2neg.v1D14$BEST.GUESS),split = ","),'[',1)
S2neg.v1D14$Donor = sapply(strsplit(as.character(S2neg.v1D14$Donor),split = "_"),'[',1)
S2neg.v1D14$Sample = paste(S2neg.v1D14$Batch, S2neg.v1D14$Donor)
S2neg.v1D14$timepoint = "v1D14"
S2neg.v1D14$Cells <- "S2P-"

Plasmablast.v2D6$Donor = sapply(strsplit(as.character(Plasmablast.v2D6$BEST.GUESS),split = ","),'[',1)
Plasmablast.v2D6$Donor = sapply(strsplit(as.character(Plasmablast.v2D6$Donor),split = "_"),'[',1)
Plasmablast.v2D6$Sample = paste(Plasmablast.v2D6$Batch, Plasmablast.v2D6$Donor)
Plasmablast.v2D6$timepoint = "v2D6"
Plasmablast.v2D6$Cells <- "Plasmablast"

Plasmablast.v1D14$Donor = sapply(strsplit(as.character(Plasmablast.v1D14$BEST.GUESS),split = ","),'[',1)
Plasmablast.v1D14$Donor = sapply(strsplit(as.character(Plasmablast.v1D14$Donor),split = "_"),'[',1)
Plasmablast.v1D14$Sample = paste(Plasmablast.v1D14$Batch, Plasmablast.v1D14$Donor)
Plasmablast.v1D14$timepoint = "v1D14"
Plasmablast.v1D14$Cells <- "Plasmablast"

#save objects for clustering
saveRDS(S2pos.v2D6, file = "S2pos.v2D6.rds")
saveRDS(S2neg.v2D6, file = "S2neg.v2D6.rds")
saveRDS(S2pos.v1D14, file = "S2pos.v1D14.rds")
saveRDS(S2neg.v1D14, file = "S2neg.v1D14.rds")
saveRDS(Plasmablast.v2D6, file = "Plasmablast.v2D6.rds")
saveRDS(Plasmablast.v1D14, file = "Plasmablast.v1D14.rds")

### End