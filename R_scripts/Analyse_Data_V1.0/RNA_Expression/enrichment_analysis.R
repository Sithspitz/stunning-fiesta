#### Lung Adenocarcinoma PanCancer Atlas RNA Enrichment ####
### Should be adopted dependent on the analysis in question ###

mypackages <- c("GSEABase", "GSVA", "Biobase", "genefilter",
                "limma", "RColorBrewer", "GSVAdata")
lapply(mypackages, library, character.only = T)
source("./R_scripts/Functions/functions.R")



### Data Processing Stages Pre-Enrichment Analysis ###
## To align RNA expression data and MTs and get into a format for GSVA ##
# Can skip this section and load the data later on in 'Enrichment Analysis' section #

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
Intermediate5 <- Intermediate4

# Subset just the Mut and the WT and drop the Mutation_Status column
my_Mut_subset <- droplevels(subset(Intermediate5, Mutation_Status == "Mut"))
my_WT_subset <- droplevels(subset(Intermediate5, Mutation_Status == "Other_Mut"))
my_Mut_subset2 <- my_Mut_subset
my_WT_subset2 <- my_WT_subset

my_Mut_subset2$Mutation_Status <- NULL
my_WT_subset2$Mutation_Status <- NULL

# Data wrangling to the correct format #

# Preparation and factoring
my_Mut_subset3 <- my_Mut_subset2
my_WT_subset3 <- my_WT_subset2

dofactor <- c("Patient.ID", "Hugo_Symbol")
Mut_factored <- factorthese(my_Mut_subset3, dofactor)
WT_factored <- factorthese(my_WT_subset3, dofactor)

# Removing duplicate values (including Hugo Symbol blanks) from the dataset 
Mut_extra <- Mut_factored
WT_extra <- WT_factored
Mut_extra$extra <- as.factor(paste(Mut_extra$Patient.ID, Mut_extra$Hugo_Symbol, sep = ","))
WT_extra$extra <- as.factor(paste(WT_extra$Patient.ID, WT_extra$Hugo_Symbol, sep = ","))

dropped_Mut <- Mut_extra[!duplicated(Mut_extra$extra), ]
dropped_WT <- WT_extra[!duplicated(WT_extra$extra), ]
dropped_Mut$extra <- NULL
dropped_WT$extra <- NULL

# Recast, remove the uneeded first row and the 'Hugo_Symbol' label 
# Export the total MT vs Other MT Data
setwd("~/DataShare/TCGA_RNA_Analysis/Input/gsva_correct_format/")

mutant_expression_data <- recast(dropped_Mut, Hugo_Symbol ~ Patient.ID, id.var = 1:2)
mutant_expression_data <- mutant_expression_data[-1, ]
removed_mutation_expression_data <- mutant_expression_data
colnames(removed_mutation_expression_data)[which(names(removed_mutation_expression_data)
                                                 == "Hugo_Symbol")] <- ""

no_mutation_expression_data <- recast(dropped_WT, Hugo_Symbol ~ Patient.ID, id.var = 1:2)
no_mutation_expression_data <- no_mutation_expression_data[-1, ]
removed_no_mutation_expression_data <- no_mutation_expression_data
colnames(removed_no_mutation_expression_data)[which(names(removed_no_mutation_expression_data)
                                                 == "Hugo_Symbol")] <- ""

write.csv(removed_mutation_expression_data, "total_STK11_and_KRAS_MT_rna_seq_EL_expression_data.csv", row.names = F)
write.csv(removed_no_mutation_expression_data, "total_no_STK11_or_KRAS_MT_rna_seq_EL_expression_data.csv", row.names = F)



### Enrichment Analysis ###

# Import Expression Data
## Make into 'Data Matrix' for GSVA
setwd("~/DataShare/TCGA_RNA_Analysis/Input/gsva_correct_format/")
no_MT_data <- read.csv("total_no_STK11_or_KRAS_MT_rna_seq_EL_expression_data.csv", header = T, row.names = 1)
MT_data <- read.csv("total_STK11_and_KRAS_MT_rna_seq_EL_expression_data.csv", header = T, row.names = 1)

no_MT_data_matrix <- as.matrix(no_MT_data)
MT_data_matrix <- as.matrix(MT_data)

# Import Gene Set (GMT Text Format, save an Excel doc as .txt)
setwd("H:/My Documents/R_WD/stunning-fiesta/R_scripts/Gene_Sets/th17_enhanced_geneset")
th17_geneset_gmt <- getGmt("th17_enhanced_geneset_gmt.txt")

# Run the test
## At the moment am doing GSVA and SSGSEA until Boris advise
setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_gsva_th17_enhanced_signature")

gsva_no_MT_enrichment_output <- gsva(no_MT_data_matrix, th17_geneset_gmt, method = "gsva")
gsva_MT_enrichment_output <- gsva(MT_data_matrix, th17_geneset_gmt, method = "gsva")

write.csv(gsva_no_MT_enrichment_output, "no_mut_gsva_output.csv", row.names = T)
write.csv(gsva_MT_enrichment_output, "mut_gsva_output.csv", row.names = T)

ssgsea_no_MT_enrichment_output <- gsva(no_MT_data_matrix, th17_geneset_gmt, method = "ssgsea")
ssgsea_MT_enrichment_output <- gsva(MT_data_matrix, th17_geneset_gmt, method = "ssgsea")

write.csv(ssgsea_no_MT_enrichment_output, "no_mut_ssgsea_output.csv", row.names = T)
write.csv(ssgsea_MT_enrichment_output, "mut_ssgsea_output.csv", row.names = T)

