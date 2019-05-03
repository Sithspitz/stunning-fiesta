### Lung Adenocarcinoma PanCancer Atlas RNA Enrichment ###
## Should be adopted dependent on the analysis in question ##
mypackages <- c("GSEABase", "GSVA", "Biobase", "genefilter",
                "limma", "RColorBrewer", "GSVAdata")
lapply(mypackages, library, character.only = T)
source("./R_scripts/Functions/functions.R")

setwd("~/DataShare/TCGA_RNA_Analysis/Input/")
rna <- read.csv("luad_tcga_pancancer_rna.csv", header = T)
mut <- read.csv("Mut_tidy_vs_other.csv", header = T)


# Merge and clean up Entrez Ids
m1 <- merge(rna, mut, by = "Patient.ID")
Intermediate2 <- m1
Intermediate2$Entrez_Gene_Id <- NULL 

# Droplevels RNASeqZScore and remove the datatype column
Intermediate3 <- droplevels(subset(Intermediate2, DataType == "RNASeqEL"))
Intermediate4 <- Intermediate3
Intermediate4$DataType <- NULL


### I have managed to process the data to how I want it 
### But need to manipulate it futher to ensure that the input is correct for GSVA
### Is there anyway of doing this without losing information?
### Perhaps look at creating an ExpressionSet object to keep metadata?
### Need to transpose the data into the correct data.matrix in paper




