#### Lung Adenocarcinoma PanCancer Atlas RNA Enrichment Pre-Processing ####
### Is the scripts from 'enrichment_analysis.R' but on their own ###
mypackages <- c("GSEABase", "GSVA", "Biobase", "genefilter",
                "limma", "RColorBrewer", "GSVAdata")
lapply(mypackages, library, character.only = T)
source("./R_scripts/Functions/functions.R")

## Preprocessing the total MT vs total other_MT datasets ##

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

# Recast, transpose and export the total MT vs Other MT Data
setwd("~/DataShare/TCGA_RNA_Analysis/Input/gsva_correct_format/")

mutant_expression_data <- recast(dropped_Mut, Patient.ID ~ Hugo_Symbol, id.var = 1:2)
mutant_expression_data$Var.2 <- NULL

no_mutation_expression_data <- recast(dropped_WT, Patient.ID ~ Hugo_Symbol, id.var = 1:2)
no_mutation_expression_data$Var.2 <- NULL

mutant_expression_data_transposed <- t(mutant_expression_data)
colnames(mutant_expression_data_transposed) <- mutant_expression_data_transposed[1, ]
mutant_expression_data_transposed <- mutant_expression_data_transposed[-1, ]

no_mutation_expression_data_transposed <- t(no_mutation_expression_data)
colnames(no_mutation_expression_data_transposed) <- no_mutation_expression_data_transposed[1, ]
no_mutation_expression_data_transposed <- no_mutation_expression_data_transposed[-1, ]

write.csv(mutant_expression_data_transposed, "total_STK11_and_KRAS_MT_rna_seq_EL_expression_data.csv", row.names = F)
write.csv(no_mutation_expression_data_transposed, "total_no_STK11_or_KRAS_MT_rna_seq_EL_expression_data.csv", row.names = F)
