---
Author: 'Felipe Assis, PhD'
title: 'Pool 1 data pre-processing'
date: '12-14-2022'
---
## Load Libraries
library("Seurat") 
library("dplyr")
library("matrixStats")
library('tidyverse')

### Import the data

M1C1.data <- Read10X(data.dir = "Moderna_P1_multi1/outs/multi/count/raw_feature_bc_matrix")
M1C2.data <- Read10X(data.dir = "Moderna_P1_multi2/outs/multi/count/raw_feature_bc_matrix")
M1C3.data <- Read10X(data.dir = "Moderna_P1_multi3/outs/multi/count/raw_feature_bc_matrix")
M1C4.data <- Read10X(data.dir = "Moderna_P1_multi4/outs/multi/count/raw_feature_bc_matrix")
M1C5.data <- Read10X(data.dir = "Moderna_P1_multi5/outs/multi/count/raw_feature_bc_matrix")
M1C6.data <- Read10X(data.dir = "Moderna_P1_multi6/outs/multi/count/raw_feature_bc_matrix")
M1C7.data <- Read10X(data.dir = "Moderna_P1_multi7/outs/multi/count/raw_feature_bc_matrix")
M1C8.data <- Read10X(data.dir = "Moderna_P1_multi8/outs/multi/count/raw_feature_bc_matrix")
M1C9.data <- Read10X(data.dir = "Moderna_P1_multi9/outs/multi/count/raw_feature_bc_matrix")
M1C10.data <- Read10X(data.dir = "Moderna_P1_multi10/outs/multi/count/raw_feature_bc_matrix")
M1C11.data <- Read10X(data.dir = "Moderna_P1_multi11/outs/multi/count/raw_feature_bc_matrix")
M1C12.data <- Read10X(data.dir = "Moderna_P1_multi12/outs/multi/count/raw_feature_bc_matrix")

### Initialize the Seurat object with the raw (non-normalized data).
M1C1 <- CreateSeuratObject(counts = M1C1.data$`Gene Expression`, project = "M1C1", assay = "RNA", min.features = 5)
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
M1C12 <- CreateSeuratObject(counts = M1C12.data$`Gene Expression`, project = "M1C12", assay = "RNA")

  ### Create ADT Assay
M1C1[["ADT"]] <- CreateAssayObject(counts = M1C1.data$`Antibody Capture`[1:18,colnames(M1C1)]) # RBD and S1, but no CD11c, CD307e and CD95 ADTs
M1C2[["ADT"]] <- CreateAssayObject(counts = M1C2.data$`Antibody Capture`[1:18,colnames(M1C2)]) # RBD and S1, but no CD11c, CD307e and CD95 ADTs
M1C3[["ADT"]] <- CreateAssayObject(counts = M1C3.data$`Antibody Capture`[1:18,colnames(M1C3)]) # RBD and S1, but no CD11c, CD307e and CD95 ADTs
M1C4[["ADT"]] <- CreateAssayObject(counts = M1C4.data$`Antibody Capture`[1:19,colnames(M1C4)]) # no RBD and S1 ADTs
M1C5[["ADT"]] <- CreateAssayObject(counts = M1C5.data$`Antibody Capture`[1:19,colnames(M1C5)]) # no RBD and S1 ADTs
M1C6[["ADT"]] <- CreateAssayObject(counts = M1C6.data$`Antibody Capture`[1:19,colnames(M1C6)]) # no RBD and S1 ADTs
M1C7[["ADT"]] <- CreateAssayObject(counts = M1C7.data$`Antibody Capture`[1:19,colnames(M1C7)]) # no RBD and S1 ADTs
M1C8[["ADT"]] <- CreateAssayObject(counts = M1C8.data$`Antibody Capture`[1:19,colnames(M1C8)]) # no RBD and S1 ADTs
M1C9[["ADT"]] <- CreateAssayObject(counts = M1C9.data$`Antibody Capture`[1:19,colnames(M1C9)]) # no RBD and S1 ADTs
M1C10[["ADT"]] <- CreateAssayObject(counts = M1C10.data$`Antibody Capture`[1:19,colnames(M1C10)])# no RBD and S1 ADTs
M1C11[["ADT"]] <- CreateAssayObject(counts = M1C11.data$`Antibody Capture`[1:19,colnames(M1C11)])# no RBD and S1 ADTs
M1C12[["ADT"]] <- CreateAssayObject(counts = M1C12.data$`Antibody Capture`[1:19,colnames(M1C12)])# no RBD and S1 ADTs

  ### Create HTO Assay
M1C1[["HTO"]] <- CreateAssayObject(counts = M1C1.data$`Antibody Capture`[19:20,colnames(M1C1)])
M1C2[["HTO"]] <- CreateAssayObject(counts = M1C2.data$`Antibody Capture`[19:20,colnames(M1C2)])
M1C3[["HTO"]] <- CreateAssayObject(counts = M1C3.data$`Antibody Capture`[19:20,colnames(M1C3)])
M1C4[["HTO"]] <- CreateAssayObject(counts = M1C4.data$`Antibody Capture`[20:21,colnames(M1C4)])
M1C5[["HTO"]] <- CreateAssayObject(counts = M1C5.data$`Antibody Capture`[20:21,colnames(M1C5)])
M1C6[["HTO"]] <- CreateAssayObject(counts = M1C6.data$`Antibody Capture`[20:21,colnames(M1C6)])
M1C7[["HTO"]] <- CreateAssayObject(counts = M1C7.data$`Antibody Capture`[20:21,colnames(M1C7)])
M1C8[["HTO"]] <- CreateAssayObject(counts = M1C8.data$`Antibody Capture`[20:21,colnames(M1C8)])
M1C9[["HTO"]] <- CreateAssayObject(counts = M1C9.data$`Antibody Capture`[20:21,colnames(M1C9)])
M1C10[["HTO"]] <- CreateAssayObject(counts = M1C10.data$`Antibody Capture`[20:21,colnames(M1C10)])
M1C11[["HTO"]] <- CreateAssayObject(counts = M1C11.data$`Antibody Capture`[20:21,colnames(M1C11)])
M1C12[["HTO"]] <- CreateAssayObject(counts = M1C12.data$`Antibody Capture`[20:21,colnames(M1C12)])

#### Add metadata to distinguish each object 
### Batch
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
M1C12$Batch <- rep("M1C12", length(colnames(M1C12)))

### orig.ident
M1C1$orig.ident <- "S2pos.v2D28"
M1C2$orig.ident <- "S2pos.v2D28"
M1C3$orig.ident <- "S2neg.v2D28"
M1C4$orig.ident <- "S2pos.v2D14"
M1C5$orig.ident <- "S2pos.v2D14"
M1C6$orig.ident <- "S2neg.v2D14"
M1C7$orig.ident <- "S2neg.v2D14"
M1C8$orig.ident <- "S2pos.v2D9"
M1C9$orig.ident <- "S2pos.v2D9"
M1C10$orig.ident <- "S2neg.v2D9"
M1C11$orig.ident <- "S2neg.v2D9"
M1C12$orig.ident <- "Plasmablast.v2D9"

### Rename barcode for each object to match barcodes in SNP file

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
M1C12 <- RenameCells(M1C12, new.names = paste(substr(colnames(M1C12), start = 1, stop = 17),"12", sep = ""))

##### Add new metadata (SNP <- DEMUXLET) 

M1C1_demuxbestList = list()
for(i in 1:length("1")){
  M1C1_demuxbestList[[i]] = read.table(paste("multi_config_1.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C2_demuxbestList[[i]] = read.table(paste("multi_config_2.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C3_demuxbestList[[i]] = read.table(paste("multi_config_3.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C4_demuxbestList[[i]] = read.table(paste("multi_config_4.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C5_demuxbestList[[i]] = read.table(paste("multi_config_5.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C6_demuxbestList[[i]] = read.table(paste("multi_config_6.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C7_demuxbestList[[i]] = read.table(paste("multi_config_7.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C8_demuxbestList[[i]] = read.table(paste("multi_config_8.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C9_demuxbestList[[i]] = read.table(paste("multi_config_9.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C10_demuxbestList[[i]] = read.table(paste("multi_config_10.best", sep = ""), sep = "\t", header = TRUE)
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
  M1C11_demuxbestList[[i]] = read.table(paste("multi_config_11.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C11_demuxbestList[[i]]$NewBarcode = paste(substr(M1C11_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"11", sep = "")
}
M1C11_demuxbestdf <- plyr::ldply(M1C11_demuxbestList, data.frame)
length(which(colnames(M1C11) %in% M1C11_demuxbestdf$NewBarcode))
setdiff(colnames(M1C11), M1C11_demuxbestdf$NewBarcode)
rownames(M1C11_demuxbestdf) <- M1C11_demuxbestdf$NewBarcode
M1C11 <- AddMetaData(M1C11, metadata = M1C11_demuxbestdf[colnames(M1C11),])

M1C12_demuxbestList = list()
for(i in 1:length("1")){
  M1C12_demuxbestList[[i]] = read.table(paste("multi_config_12.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length("1")){
  M1C12_demuxbestList[[i]]$NewBarcode = paste(substr(M1C12_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),"12", sep = "")
}
M1C12_demuxbestdf <- plyr::ldply(M1C12_demuxbestList, data.frame)
length(which(colnames(M1C12) %in% M1C12_demuxbestdf$NewBarcode))
setdiff(colnames(M1C12), M1C12_demuxbestdf$NewBarcode)
rownames(M1C12_demuxbestdf) <- M1C12_demuxbestdf$NewBarcode
M1C12 <- AddMetaData(M1C12, metadata = M1C12_demuxbestdf[colnames(M1C12),])

### Rename the cells again to match barcodes from VDJ data.

M1C1 <- RenameCells(M1C1, new.names = paste(substr(colnames(M1C1), start = 1, stop = 17),"<sample1><pool1>", sep = ""))
M1C2 <- RenameCells(M1C2, new.names = paste(substr(colnames(M1C2), start = 1, stop = 17),"<sample2><pool1>", sep = ""))
M1C3 <- RenameCells(M1C3, new.names = paste(substr(colnames(M1C3), start = 1, stop = 17),"<sample3><pool1>", sep = ""))
M1C4 <- RenameCells(M1C4, new.names = paste(substr(colnames(M1C4), start = 1, stop = 17),"<sample4><pool1>", sep = ""))
M1C5 <- RenameCells(M1C5, new.names = paste(substr(colnames(M1C5), start = 1, stop = 17),"<sample5><pool1>", sep = ""))
M1C6 <- RenameCells(M1C6, new.names = paste(substr(colnames(M1C6), start = 1, stop = 17),"<sample6><pool1>", sep = ""))
M1C7 <- RenameCells(M1C7, new.names = paste(substr(colnames(M1C7), start = 1, stop = 17),"<sample7><pool1>", sep = ""))
M1C8 <- RenameCells(M1C8, new.names = paste(substr(colnames(M1C8), start = 1, stop = 17),"<sample8><pool1>", sep = ""))
M1C9 <- RenameCells(M1C9, new.names = paste(substr(colnames(M1C9), start = 1, stop = 17),"<sample9><pool1>", sep = ""))
M1C10 <- RenameCells(M1C10, new.names = paste(substr(colnames(M1C10), start = 1, stop = 17),"<sample10><pool1>", sep = ""))
M1C11 <- RenameCells(M1C11, new.names = paste(substr(colnames(M1C11), start = 1, stop = 17),"<sample11><pool1>", sep = ""))
M1C12 <- RenameCells(M1C12, new.names = paste(substr(colnames(M1C12), start = 1, stop = 17),"<sample12><pool1>", sep = ""))

### merge objects
S2pos.v2D28 <- merge(M1C1, y = M1C2)
S2neg.v2D28 <- M1C3
S2pos.v2D14 <- merge(M1C4, y = M1C5)
S2neg.v2D14 <- merge(M1C6, y = M1C7)
S2pos.v2D9 <- merge(M1C8, y = M1C9)
S2neg.v2D9 <- merge(M1C10, y = M1C11)
Plasmablast.v2D9 <- M1C12

# remove abmormally high count cells
S2pos.v2D28 <- subset(S2pos.v2D28, subset = nCount_ADT < 30000 & nCount_HTO < 15000)
S2neg.v2D28 <- subset(S2neg.v2D28, subset = nCount_ADT < 30000 & nCount_HTO < 15000)
S2pos.v2D14 <- subset(S2pos.v2D14, subset = nCount_ADT < 30000 & nCount_HTO < 15000)
S2neg.v2D14 <- subset(S2neg.v2D14, subset = nCount_ADT < 30000 & nCount_HTO < 15000)
S2pos.v2D9 <- subset(S2pos.v2D9, subset = nCount_ADT < 30000 & nCount_HTO < 15000)
S2neg.v2D9 <- subset(S2neg.v2D9, subset = nCount_ADT < 30000 & nCount_HTO < 15000)
Plasmablast.v2D9 <- subset(Plasmablast.v2D9, subset = nCount_ADT < 30000 & nCount_HTO < 15000)

S2pos.v2D28[["percent.mt"]] <- PercentageFeatureSet(S2pos.v2D28, pattern = "^MT-")
S2neg.v2D28[["percent.mt"]] <- PercentageFeatureSet(S2neg.v2D28, pattern = "^MT-")
S2pos.v2D14[["percent.mt"]] <- PercentageFeatureSet(S2pos.v2D14, pattern = "^MT-")
S2neg.v2D14[["percent.mt"]] <- PercentageFeatureSet(S2neg.v2D14, pattern = "^MT-")
S2pos.v2D9[["percent.mt"]] <- PercentageFeatureSet(S2pos.v2D9, pattern = "^MT-")
S2neg.v2D9[["percent.mt"]] <- PercentageFeatureSet(S2neg.v2D9, pattern = "^MT-")
Plasmablast.v2D9[["percent.mt"]] <- PercentageFeatureSet(Plasmablast.v2D9, pattern = "^MT-")

### HTO demultiplex

S2pos.v2D28 <- NormalizeData(S2pos.v2D28, assay = "HTO", normalization.method = "CLR", margin = 2)
S2pos.v2D28 <- ScaleData(S2pos.v2D28, assay = "HTO", model.use = "linear")
S2pos.v2D28 = MULTIseqDemux(S2pos.v2D28, autoThresh = FALSE, quantile = 0.5)
Idents(S2pos.v2D28) <- "MULTI_ID"
S2pos.v2D28 <- subset(S2pos.v2D28, idents = c("HTO3-PROT"))
S2pos.v2D28 <- subset(S2pos.v2D28, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                        nFeature_RNA < 4000 & percent.mt < 10)

S2neg.v2D28 <- NormalizeData(S2neg.v2D28, assay = "HTO", normalization.method = "CLR", margin = 2)
S2neg.v2D28 <- ScaleData(S2neg.v2D28, assay = "HTO", model.use = "linear")
S2neg.v2D28 = MULTIseqDemux(S2neg.v2D28, autoThresh = FALSE, quantile = 0.5)
Idents(S2neg.v2D28) <- "MULTI_ID"
S2neg.v2D28 <- subset(S2neg.v2D28, idents = c("HTO3-PROT"))
S2neg.v2D28 <- subset(S2neg.v2D28, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                        nFeature_RNA < 4000 & percent.mt < 10)

S2pos.v2D14 <- NormalizeData(S2pos.v2D14, assay = "HTO", normalization.method = "CLR", margin = 2)
S2pos.v2D14 <- ScaleData(S2pos.v2D14, assay = "HTO", model.use = "linear")
S2pos.v2D14 = MULTIseqDemux(S2pos.v2D14, autoThresh = FALSE, quantile = 0.5)
table(S2pos.v2D14$MULTI_ID)
S2pos.v2D14 <- subset(S2pos.v2D14, idents = c("HTO5-PROT"))
S2pos.v2D14 <- subset(S2pos.v2D14, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                        nFeature_RNA < 4000 & percent.mt < 10)

S2neg.v2D14 <- NormalizeData(S2neg.v2D14, assay = "HTO", normalization.method = "CLR")
S2neg.v2D14 <- ScaleData(S2neg.v2D14, assay = "HTO", model.use = "linear")
S2neg.v2D14 = MULTIseqDemux(S2neg.v2D14, autoThresh = FALSE, quantile = 0.5)
Idents(S2neg.v2D14) <- "MULTI_ID"
S2neg.v2D14 <- subset(S2neg.v2D14, idents = c("HTO5-PROT"))
S2neg.v2D14 <- subset(S2neg.v2D14, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                        nFeature_RNA < 4000 & percent.mt < 10)

S2pos.v2D9 <- NormalizeData(S2pos.v2D9, assay = "HTO", normalization.method = "CLR", margin = 2)
S2pos.v2D9 <- ScaleData(S2pos.v2D9, assay = "HTO", model.use = "linear")
S2pos.v2D9 = MULTIseqDemux(S2pos.v2D9, assay = "HTO", autoThresh = FALSE, quantile = 0.5)
Idents(S2pos.v2D9) <- "MULTI_ID"
S2pos.v2D9 <- subset(S2pos.v2D9, idents = c("HTO7-PROT"))
S2pos.v2D9 <- subset(S2pos.v2D9, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                       nFeature_RNA < 4000 & percent.mt < 10)

S2neg.v2D9 <- NormalizeData(S2neg.v2D9, assay = "HTO", normalization.method = "CLR")
S2neg.v2D9 <- ScaleData(S2neg.v2D9, assay = "HTO", model.use = "linear")
S2neg.v2D9 = MULTIseqDemux(S2neg.v2D9, assay = "HTO", autoThresh = FALSE, quantile = 0.5)
Idents(S2neg.v2D9) <- "MULTI_ID"
S2neg.v2D9 <- subset(S2neg.v2D9, idents = c("HTO7-PROT"))
S2neg.v2D9 <- subset(S2neg.v2D9, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 200 & 
                       nFeature_RNA < 4000 & percent.mt < 10)

Plasmablast.v2D9 <- NormalizeData(Plasmablast.v2D9, assay = "HTO", normalization.method = "CLR", margin = 2)
Plasmablast.v2D9 <- ScaleData(Plasmablast.v2D9, assay = "HTO", model.use = "linear")
Plasmablast.v2D9 = MULTIseqDemux(Plasmablast.v2D9, autoThresh = FALSE, quantile = 0.5)
Idents(Plasmablast.v2D9) <- "MULTI_ID"
Plasmablast.v2D9 <- subset(Plasmablast.v2D9, idents = c("HTO7-PROT"))
Plasmablast.v2D9.100 <- subset(Plasmablast.v2D9, subset = DROPLET.TYPE == "SNG" & nFeature_RNA > 100 & 
                             nFeature_RNA < 4000 & percent.mt < 10)

## add metadata
S2pos.v2D28$Donor = sapply(strsplit(as.character(S2pos.v2D28$BEST.GUESS),split = ","),'[',1)
S2pos.v2D28$Donor = sapply(strsplit(as.character(S2pos.v2D28$Donor),split = "_"),'[',1)
S2pos.v2D28$Sample = paste(S2pos.v2D28$Batch, S2pos.v2D28$Donor)
S2pos.v2D28$timepoint = "v2D28"
S2pos.v2D28$Cells = "S2P+"

S2neg.v2D28$Donor = sapply(strsplit(as.character(S2neg.v2D28$BEST.GUESS),split = ","),'[',1)
S2neg.v2D28$Donor = sapply(strsplit(as.character(S2neg.v2D28$Donor),split = "_"),'[',1)
S2neg.v2D28$Sample = paste(S2neg.v2D28$Batch, S2neg.v2D28$Donor)
S2neg.v2D28$timepoint = "v2D28"
S2neg.v2D28$Cells = "S2P-"

S2pos.v2D14$Donor = sapply(strsplit(as.character(S2pos.v2D14$BEST.GUESS),split = ","),'[',1)
S2pos.v2D14$Donor = sapply(strsplit(as.character(S2pos.v2D14$Donor),split = "_"),'[',1)
S2pos.v2D14$Sample = paste(S2pos.v2D14$Batch, S2pos.v2D14$Donor)
S2pos.v2D14$timepoint = "v2D14"
S2pos.v2D14$Cells = "S2P+"

S2neg.v2D14$Donor = sapply(strsplit(as.character(S2neg.v2D14$BEST.GUESS),split = ","),'[',1)
S2neg.v2D14$Donor = sapply(strsplit(as.character(S2neg.v2D14$Donor),split = "_"),'[',1)
S2neg.v2D14$Sample = paste(S2neg.v2D14$Batch, S2neg.v2D14$Donor)
S2neg.v2D14$timepoint = "v2D14"
S2neg.v2D14$Cells = "S2P-"

S2pos.v2D9$Donor = sapply(strsplit(as.character(S2pos.v2D9$BEST.GUESS),split = ","),'[',1)
S2pos.v2D9$Donor = sapply(strsplit(as.character(S2pos.v2D9$Donor),split = "_"),'[',1)
S2pos.v2D9$Sample = paste(S2pos.v2D9$Batch, S2pos.v2D9$Donor)
S2pos.v2D9$timepoint = "v2D9"
S2pos.v2D9$Cells = "S2P+"

S2neg.v2D9$Donor = sapply(strsplit(as.character(S2neg.v2D9$BEST.GUESS),split = ","),'[',1)
S2neg.v2D9$Donor = sapply(strsplit(as.character(S2neg.v2D9$Donor),split = "_"),'[',1)
S2neg.v2D9$Sample = paste(S2neg.v2D9$Batch, S2neg.v2D9$Donor)
S2neg.v2D9$timepoint = "v2D9"
S2neg.v2D9$Cells = "S2P-"

Plasmablast.v2D9$Donor = sapply(strsplit(as.character(Plasmablast.v2D9$BEST.GUESS),split = ","),'[',1)
Plasmablast.v2D9$Donor = sapply(strsplit(as.character(Plasmablast.v2D9$Donor),split = "_"),'[',1)
Plasmablast.v2D9$Sample = paste(Plasmablast.v2D9$Batch, Plasmablast.v2D9$Donor)
Plasmablast.v2D9$timepoint = "v2D9"
Plasmablast.v2D9Cells = "Plasmablast"

#save objects for clustering
saveRDS(S2pos.v2D28, file = "S2pos.v2D28.rds")
saveRDS(S2neg.v2D28, file = "S2neg.v2D28.rds")
saveRDS(S2pos.v2D14, file = "S2pos.v2D14.rds")
saveRDS(S2neg.v2D14, file = "S2neg.v2D14.rds")
saveRDS(S2pos.v2D9, file = "S2pos.v2D9.rds")
saveRDS(S2neg.v2D9, file = "S2neg.v2D9.rds")
saveRDS(Plasmablast.v2D9, file = "Plasmablast.v2D9.rds")

##### End
