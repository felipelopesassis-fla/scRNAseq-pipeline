# Scripts for scRNAseq analysis performed at "Tracking B Cell Responses to the SARS-CoV-2 mRNA-1273 Vaccine" study
Felipe Lopes de Assis felipe.lopesdeassis@nih.gov 12/14/2022

Dependencies

R 4.1.0 

# R packages:
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)
library(dittoSeq)

Processing Assis et al. 2023 raw data

To process data from the B CELL TRAJECTORY study cohort, move to data/Trajectory and run

# download all raw data files into appropriate place

# Import and identify each file properly      

raw.data <- Read10X(data.dir = "raw_feature_bc_matrix") 

# Initialize the Seurat objects with the raw data (non-normalized data), and perform SNP 
# and HTO demultiplexing:   

Rscripts: 

01.Pool1.preprocessing.R
02.Pool2.preprocessing.R
03.Pool3.preprocessing.R

# Combine and integrate datasets into one file

Rscript 04.Integration.and.clustering.R

# Once all of the above are successfully run (congratulations!), make figures:

Rscript 05.Figures.R

Figure drafts will be in the analysis/figures/drafts folder. These are touched up with Adobe Illustrator to create the final figure panels.
