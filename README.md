# Scripts for scRNAseq analysis performed at "Tracking B Cell Responses to the SARS-CoV-2 mRNA-1273 Vaccine" study

Felipe Lopes de Assis felipe.lopesdeassis@nih.gov 12/14/2022

# Dependencies

R 4.1.0 

R packages:

library(Seurat)

library(patchwork)

library(ggplot2)

library(dplyr)

library(scales)

library(RColorBrewer)

library(dittoSeq)

# Download raw data files into appropriate place and identify each file properly 

mkdir raw.data

data_source = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219098

accession number: GSE219098

# Import each raw data file using Seurat Read10X funnction       

raw.data <- Read10X(data.dir = "raw.data/raw_feature_bc_matrix") 

# Initialize the Seurat objects with the raw data (non-normalized data), and perform SNP and HTO demultiplexing:   

R scripts: 

01.Pool1.preprocessing.R

02.Pool2.preprocessing.R

03.Pool3.preprocessing.R

# Combine and integrate datasets into one file

Rscript 04.Integration&clustering.R

# Once all of the above are successfully run (congratulations!), make figures:

Rscript 05.Figures.R

Figure drafts will be in the analysis/figures/drafts folder. These are touched up with Adobe Illustrator to create the final figure panels.
