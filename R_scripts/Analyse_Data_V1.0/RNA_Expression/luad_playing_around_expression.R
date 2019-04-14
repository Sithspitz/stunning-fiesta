### Lung Adenocarcinoma PanCancer Atlas RNA Expression Playing Around ###
source("./R_scripts/Functions/functions.R")
source("./R_scripts/Functions/GeneSets.R")
setwd("~/DataShare/TCGA_RNA_Analysis/Testing/")
rna <- read.csv("luad_tcga_pancancer_rna.csv", header = T)

rna$DataType <- rna$DataType
rna$DataType <- as.factor(rna$DataType)


# Label Gene Signatures
rna$th17 <- ifelse((rna$Hugo_Symbol %in% Th17_genes), T, F)

th17datatemp <- CalGeneMeanRNASeqz(rna, th17)
### sort out calgenemeanrnaseqz
### play around with the CRC microarray data and see how it works
### check how i can apply to my data from this