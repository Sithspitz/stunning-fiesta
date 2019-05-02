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
mut <- read.csv("KRAS_STK11_plus_dbl_tidy_vs_other.csv", header = T)


# Droplevels !Other_Mut #
## Gives just KRAS, STK11 and double MT ##
m1 <- merge(rna, mut, by = "Patient.ID")
Intermediate <- droplevels(subset(m1, Mutation_Status != "Other_Mut"))

# Remove Entrez_Gene_Id #
Intermediate2 <- Intermediate
Intermediate2$Entrez_Gene_Id <- NULL 

# Droplevels RNASeqZScore #
## Gives just the non-zscore RSEM values ##
Intermediate3 <- droplevels(subset(Intermediate2, DataType == "RNASeqEL"))



#### TEST GVSA WITH FAKE DATA ####
## This is the non-S2 gene set enriched example ##
not_s2_enriched_1 <- read.csv("test_gene_set.csv", header = T, row.names = 1)
not_s2_enriched_2 <- data.matrix(not_s2_enriched_1) ## has to be a matrix
gmt_gene_set <- getGmt("test_enrich_Set_gmt_idea.txt") ## Using a gmt gene set

not_s2_enriched_output <- gsva(not_s2_enriched_2, gmt_gene_set, method = "gsva")

## This is the S2 gene set enriched example ##
s2_enriched_1 <- read.csv("test_gene_set_s2_enriched.csv", header = T, row.names = 1)
s2_enriched_2 <- data.matrix(s2_enriched_1)

s2_enriched_output <- gsva(s2_enriched_2, gmt_gene_set, method = "gsva")
s2_enriched_output_ssgsea <- gsva(s2_enriched_2, gmt_gene_set, method = "ssgsea")

## This is the S2 high and S3 low gene set enriched example ##

s2_hi_s3_lo_1 <- read.csv("test_gene_set_s2_enriched_s3_low.csv", header = T, row.names = 1)
s2_hi_s3_lo_1 <- data.matrix(s2_hi_s3_lo_1)

s2_hi_s3_lo_output <- gsva(s2_hi_s3_lo_11, gmt_gene_set, method = "gsva")
s2_hi_s3_lo_output_ssgsea <- gsva(s2_hi_s3_lo_1, gmt_gene_set, method = "ssgsea")


#### FINISHED END OF LAST WEEK BUT WHAT DO THESE VALUES ACTUALLY MEAN? 
#### AND WHAT ON EARTH DO I DO NEXT...?

