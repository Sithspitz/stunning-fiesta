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
mut2 <- read.csv("KRAS_STK11_plus_dbl_tidy_vs_other.csv", header = T)


mergev1 <- merge(rna, mut2, by = "Patient.ID")
mergev2 <- mergev1
mergev3 <- mergev1
mergev4 <- mergev1
extracting_other_muts <- mergev3[mergev3 == "Other_Mut"]

Intermediate <- droplevels(subset(mergev4, Mutation_Status != "Other_Mut"))
### THIS DROPLEVELS SCRIPT REALLY IS GOD HERE, WRITE UP PROPERLY WITH THIS SCRIPT!! ###