### Lung Adenocarcinoma PanCancer Atlas RNA Enrichment ###
## Should be adopted dependent on the analysis in question ##
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") } else { 
    (BiocManager::install("GSVA", version = "3.8")) }
library("GSVA")
library("GSEABase")
library("Biobase", "genefilter", "limma", "RColorBrewer")
source("./R_scripts/Functions/functions.R")

setwd("~/DataShare/TCGA_RNA_Analysis/Input/")
rna <- read.csv("luad_tcga_pancancer_rna.csv", header = T)

