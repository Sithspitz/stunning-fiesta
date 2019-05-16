#### Lung Adenocarcinoma PanCancer Atlas RNA Enrichment ####
### Should be adopted dependent on the analysis in question ###
## For a functional example of these scripts look at th17_enrichment_analysis.R ##

mypackages <- c("GSEABase", "GSVA", "Biobase", "genefilter",
                "limma", "RColorBrewer", "GSVAdata", "scales", "dplyr")
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
## Or just leave it and use the entire dataset 
my_Mut_subset <- droplevels(subset(Intermediate5, Mutation_Status == "Mut"))
my_WT_subset <- droplevels(subset(Intermediate5, Mutation_Status == "Other_Mut"))
my_Mut_subset2 <- my_Mut_subset
my_WT_subset2 <- my_WT_subset

my_Mut_subset2$Mutation_Status <- NULL
my_WT_subset2$Mutation_Status <- NULL
Intermediate5$Mutation_Status <- NULL

# Data wrangling to the correct format #

# Preparation and factoring
my_Mut_subset3 <- my_Mut_subset2
my_WT_subset3 <- my_WT_subset2

dofactor <- c("Patient.ID", "Hugo_Symbol")
Mut_factored <- factorthese(my_Mut_subset3, dofactor)
WT_factored <- factorthese(my_WT_subset3, dofactor)
Intermediate6 <- factorthese(Intermediate5, dofactor)

# Removing duplicate values (including Hugo Symbol blanks) from the dataset 
Mut_extra <- Mut_factored
WT_extra <- WT_factored
Inter_extra <- Intermediate6
Mut_extra$extra <- as.factor(paste(Mut_extra$Patient.ID, Mut_extra$Hugo_Symbol, sep = ","))
WT_extra$extra <- as.factor(paste(WT_extra$Patient.ID, WT_extra$Hugo_Symbol, sep = ","))
Inter_extra$extra <- as.factor(paste(Inter_extra$Patient.ID, Inter_extra$Hugo_Symbol, sep = ","))

dropped_Mut <- Mut_extra[!duplicated(Mut_extra$extra), ]
dropped_WT <- WT_extra[!duplicated(WT_extra$extra), ]
dropped_Inter <- Inter_extra[!duplicated(Inter_extra$extra), ]
dropped_Mut$extra <- NULL
dropped_WT$extra <- NULL
dropped_Inter$extra <- NULL

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

inter_expression_data <- recast(dropped_Inter, Hugo_Symbol ~ Patient.ID, id.var = 1:2)
inter_expression_data <- inter_expression_data[-1, ]
removed_inter_expression_data <- inter_expression_data
colnames(removed_inter_expression_data)[which(names(removed_inter_expression_data) == "Hugo_Symbol")] <- ""

write.csv(removed_mutation_expression_data, "total_STK11_and_KRAS_MT_rna_seq_EL_expression_data.csv", row.names = F)
write.csv(removed_no_mutation_expression_data, "total_no_STK11_or_KRAS_MT_rna_seq_EL_expression_data.csv", row.names = F)
write.csv(removed_inter_expression_data, "total_rna_seq_EL_expression_data.csv" , row.names = F)


### Enrichment Analysis ###

# Import Expression Data
## Make into 'Data Matrix' for GSVA
setwd("~/DataShare/TCGA_RNA_Analysis/Input/gsva_correct_format/")
no_MT_data <- read.csv("total_no_STK11_or_KRAS_MT_rna_seq_EL_expression_data.csv", header = T, row.names = 1)
MT_data <- read.csv("total_STK11_and_KRAS_MT_rna_seq_EL_expression_data.csv", header = T, row.names = 1)
total_data <- read.csv("total_rna_seq_EL_expression_data.csv", header = T, row.names = 1)

no_MT_data_matrix <- as.matrix(no_MT_data)
MT_data_matrix <- as.matrix(MT_data)
total_data_matrix <- as.matrix(total_data)

# Import Gene Set (GMT Text Format, save an Excel doc as .txt)
setwd("H:/My Documents/R_WD/stunning-fiesta/R_scripts/Gene_Sets/th17_enhanced_geneset")
th17_geneset_gmt <- getGmt("th17_enhanced_geneset_gmt.txt")

# Run the test
## At the moment am doing GSVA SSGSEA and PLAGE until Boris advise
setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_gsva_th17_enhanced_signature")

gsva_no_MT_enrichment_output <- gsva(no_MT_data_matrix, th17_geneset_gmt, method = "gsva")
gsva_MT_enrichment_output <- gsva(MT_data_matrix, th17_geneset_gmt, method = "gsva")
gsva_total_enrichment_output <- gsva(total_data_matrix, th17_geneset_gmt, method = "gsva")

write.csv(gsva_no_MT_enrichment_output, "no_mut_gsva_output.csv", row.names = T)
write.csv(gsva_MT_enrichment_output, "mut_gsva_output.csv", row.names = T)
write.csv(gsva_total_enrichment_output, "gsva_output.csv", row.names = T)

ssgsea_no_MT_enrichment_output <- gsva(no_MT_data_matrix, th17_geneset_gmt, method = "ssgsea")
ssgsea_MT_enrichment_output <- gsva(MT_data_matrix, th17_geneset_gmt, method = "ssgsea")
ssgsea_total_enrichment_output <- gsva(total_data_matrix, th17_geneset_gmt, method = "ssgsea")

write.csv(ssgsea_no_MT_enrichment_output, "no_mut_ssgsea_output.csv", row.names = T)
write.csv(ssgsea_MT_enrichment_output, "mut_ssgsea_output.csv", row.names = T)
write.csv(ssgsea_total_enrichment_output, "ssgsea_output.csv", row.names = T)

plage_no_MT_enrichment_output <- gsva(no_MT_data_matrix, th17_geneset_gmt, method = "plage")
plage_MT_enrichment_output <- gsva(MT_data_matrix, th17_geneset_gmt, method = "plage")
plage_total_enrichment_output <- gsva(total_data_matrix, th17_geneset_gmt, method = "plage")

write.csv(plage_no_MT_enrichment_output, "no_mut_plage_output.csv", row.names = T)
write.csv(plage_MT_enrichment_output, "mut_plage_output.csv", row.names = T)
write.csv(plage_total_enrichment_output, "plage_output.csv", row.names = T)

# Unionise Enrichment Output Data 
## For the WT vs MT data ONLY
gsva_no_MT_enrichment_transposed <- as.data.frame(t(gsva_no_MT_enrichment_output))
gsva_no_MT_enrichment_transposed["mutation"] <- paste("WT")

gsva_MT_enrichment_output_transposed <- as.data.frame(t(gsva_MT_enrichment_output))
gsva_MT_enrichment_output_transposed["mutation"] <- paste("MT")

gsva_union <- union(gsva_no_MT_enrichment_transposed, gsva_MT_enrichment_output_transposed)

ssgsea_no_MT_enrichment_transposed <- as.data.frame(t(ssgsea_no_MT_enrichment_output))
ssgsea_no_MT_enrichment_transposed["mutation"] <- paste("WT")

ssgsea_MT_enrichment_transposed <- as.data.frame(t(ssgsea_MT_enrichment_output))
ssgsea_MT_enrichment_transposed["mutation"] <- paste("MT")

ssgsea_union <- union(ssgsea_no_MT_enrichment_transposed, ssgsea_MT_enrichment_transposed)

plage_no_MT_enrichment_transposed <- as.data.frame(t(plage_no_MT_enrichment_output))
plage_no_MT_enrichment_transposed["mutation"] <- paste("WT")

plage_MT_enrichment_output_transposed <- as.data.frame(t(plage_MT_enrichment_output))
plage_MT_enrichment_output_transposed["mutation"] <- paste("MT")

plage_union <- union(plage_no_MT_enrichment_transposed, plage_MT_enrichment_output_transposed)

setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_gsva_th17_enhanced_signature")
write.csv(plage_union, "plage_union_output.csv", row.names = F)
write.csv(ssgsea_union, "ssgsea_union_output.csv", row.names = F)
write.csv(gsva_union, "gsva_union_output.csv", row.names = F)

# Align total data to mutational subtype
setwd("M:/Richard B/Analysis/2019/May 2019/TCGA_RNA_Analysis/Input")
mut <- read.csv("Mut_tidy_vs_other.csv", header = T)

gsva_total_output_transposed <- as.data.frame((t(gsva_total_enrichment_output)))
gsva_total_output_transposed <- rownames_to_column(gsva_total_output_transposed, var = "Patient.ID")
gsva_total_output_mutation_merge <- merge(gsva_total_output_transposed, mut)
gsva_total_output_mutation_merge$Patient.ID <- NULL

ssgsea_total_output_transposed <- as.data.frame((t(ssgsea_total_enrichment_output)))
ssgsea_total_output_transposed <- rownames_to_column(ssgsea_total_output_transposed, var = "Patient.ID")
ssgsea_total_output_mutation_merge <- merge(ssgsea_total_output_transposed, mut)
ssgsea_total_output_mutation_merge$Patient.ID <- NULL

plage_total_output_transposed <- as.data.frame((t(plage_total_enrichment_output)))
plage_total_output_transposed <- rownames_to_column(plage_total_output_transposed, var = "Patient.ID")
plage_total_output_mutation_merge <- merge(plage_total_output_transposed, mut)
plage_total_output_mutation_merge$Patient.ID <- NULL

setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_gsva_th17_enhanced_signature")
write.csv(gsva_total_output_mutation_merge, "gsva_total_merged.csv", row.names = F)
write.csv(ssgsea_total_output_mutation_merge, "ssgsea_total_merged.csv", row.names = F)
write.csv(plage_total_output_mutation_merge, "plage_total_merged.csv", row.names = F)


### Stats Comparison and Plotting ###

# Subsetting total data
### gsva_no_mt_subset <- subset(gsva_no_MT_enrichment_output, rownames(gsva_MT_enrichment_output) %in% "Th17_Genes")

gsva_total_output_wt <- droplevels(subset(gsva_total_output_mutation_merge, Mutation_Status == "Other_Mut"))
gsva_total_output_wt$Mutation_Status <- NULL
gsva_total_output_wt_transposed <- as.numeric(t(gsva_total_output_wt))
gsva_total_output_mt <- droplevels(subset(gsva_total_output_mutation_merge, Mutation_Status == "Mut"))
gsva_total_output_mt$Mutation_Status <- NULL
gsva_total_output_mt_transposed <- as.numeric(t(gsva_total_output_mt))

ssgsea_total_output_wt <- droplevels(subset(ssgsea_total_output_mutation_merge, Mutation_Status == "Other_Mut"))
ssgsea_total_output_wt$Mutation_Status <- NULL
ssgsea_total_output_wt_transposed <- as.numeric(t(ssgsea_total_output_wt))
ssgsea_total_output_mt <- droplevels(subset(ssgsea_total_output_mutation_merge, Mutation_Status == "Mut"))
ssgsea_total_output_wt$Mutation_Status <- NULL
ssgsea_total_output_mt_transposed <- as.numeric(t(ssgsea_total_output_mt))

plage_total_output_wt <- droplevels(subset(plage_total_output_mutation_merge, Mutation_Status == "Other_Mut"))
plage_total_output_wt$Mutation_Status <- NULL
plage_total_output_wt_transposed <- as.numeric(t(plage_total_output_wt))
plage_total_output_mt <- droplevels(subset(plage_total_output_mutation_merge, Mutation_Status == "Mut"))
plage_total_output_mt$Mutation_Status <- NULL
plage_total_output_mt_transposed <- as.numeric(t(plage_total_output_mt))

# Wilcox Test Old - on the pre-MT group sorted scores
stat_result_gsva_output <- wilcox.test(gsva_no_MT_enrichment_output, gsva_MT_enrichment_output)
stat_result_ssgsea_output <- wilcox.test(ssgsea_no_MT_enrichment_output, ssgsea_MT_enrichment_output)
stat_result_plage_output <- wilcox.test(plage_no_MT_enrichment_output, plage_MT_enrichment_output)
stat_result_gsva_output
stat_result_ssgsea_output
stat_result_plage_output

pVal_1_gsva <- wilcox.test(gsva_no_MT_enrichment_output, gsva_MT_enrichment_output)$p.value
pVal_2_gsva <- format(round(pVal_1_gsva, 4), nsmall = 4)
pVal_3_gsva <- paste("p = ", pVal_2_gsva, sep = "")

pVal_1_ssgsea <- wilcox.test(ssgsea_no_MT_enrichment_output, ssgsea_MT_enrichment_output)$p.value
pVal_2_ssgsea <- format(round(pVal_1_ssgsea, 4), nsmall = 4)
pVal_3_ssgsea <- paste("p = ", pVal_2_ssgsea, sep = "")

pVal_1_plage <- wilcox.test(plage_no_MT_enrichment_output, plage_MT_enrichment_output)$p.value
pVal_2_plage <- format(round(pVal_1_plage, 4), nsmall = 4)
pVal_3_plage <- paste("p = ", pVal_2_plage, sep = "")

# Wilcox Test New - on post-MT group sorted scored
stat_result_gsva_total_output <- wilcox.test(gsva_total_output_wt_transposed, gsva_total_output_mt_transposed)
stat_result_ssgsea_total_output <- wilcox.test(ssgsea_total_output_wt_transposed, ssgsea_total_output_mt_transposed)
stat_result_plage_total_output <- wilcox.test(plage_total_output_wt_transposed, plage_total_output_mt_transposed)
stat_result_gsva_total_output
stat_result_ssgsea_total_output
stat_result_plage_total_output

pVal_1_gsva_total <- wilcox.test(gsva_total_output_wt_transposed, gsva_total_output_mt_transposed)$p.value
pVal_2_gsva_total <- format(round(pVal_1_gsva_total, 4), nsmall = 4)
pVal_3_gsva_total <- paste("p = ", pVal_2_gsva_total, sep = "")

pVal_1_ssgsea_total <- wilcox.test(ssgsea_total_output_wt_transposed, ssgsea_total_output_mt_transposed)$p.value
pVal_2_ssgsea_total <- format(round(pVal_1_ssgsea_total, 4), nsmall = 4)
pVal_3_ssgsea_total <- paste("p = ", pVal_2_ssgsea_total, sep = "")

pVal_1_plage_total <- wilcox.test(plage_total_output_wt_transposed, plage_total_output_mt_transposed)$p.value
pVal_2_plage_total <- format(round(pVal_1_plage_total, 4), nsmall = 4)
pVal_3_plage_total <- paste("p = ", pVal_2_plage_total, sep = "")

# GSVA Violin Plot - old comparison
setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_gsva_th17_enhanced_signature")
cbcols <- c("WT" = "#0000FF",
            "MT" = "#FF0000")

gsva_union$mutation <-factor(gsva_union$mutation, levels = c("WT", "MT"))
my_comparisons <- list(c("WT", "MT"))

cairo_pdf("./gsva_violin.pdf")
violin_1 <- ggplot(gsva_union, aes(x = mutation, y = Th17_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(mutation, fill = mutation),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.5, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "GSVA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("GSVA Enrichment of Th17 Enhanced Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_gsva, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()
# Button to delete the violin file: unlink("gsva_violin.pdf")

# SSGSEA Violin Plot - old comparison
setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_gsva_th17_enhanced_signature")
cbcols <- c("WT" = "#0000FF",
            "MT" = "#FF0000")

ssgsea_union$mutation <-factor(ssgsea_union$mutation, levels = c("WT", "MT"))
my_comparisons <- list(c("WT", "MT"))

cairo_pdf("./ssgsea_violin.pdf")
violin_1 <- ggplot(ssgsea_union, aes(x = mutation, y = Th17_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(mutation, fill = mutation),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.5, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "SSGSEA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("SSGSEA Enrichment of Th17 Enhanced Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_ssgsea, x = 2.4, y = -0.2, size = 4)
violin_1
dev.off()

# Plage Violin Plot - old comparison
setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_gsva_th17_enhanced_signature")
cbcols <- c("WT" = "#0000FF",
            "MT" = "#FF0000")

plage_union$mutation <-factor(plage_union$mutation, levels = c("WT", "MT"))
my_comparisons <- list(c("WT", "MT"))

cairo_pdf("./plage_violin.pdf")
violin_1 <- ggplot(plage_union, aes(x = mutation, y = Th17_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(mutation, fill = mutation),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.5, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "plage Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("plage Enrichment of Th17 Enhanced Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_plage, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# GSVA Violin Plot - new comparison
setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_gsva_th17_enhanced_signature")
cbcols <- c("Other_Mut" = "#0000FF",
            "Mut" = "#FF0000")

gsva_total_output_mutation_merge$Mutation_Status <-factor(gsva_total_output_mutation_merge$Mutation_Status,
                                                          levels = c("Other_Mut", "Mut"))
my_comparisons <- list(c("Other_Mut", "Mut"))

cairo_pdf("./gsva_total_violin.pdf")
violin_1 <- ggplot(gsva_total_output_mutation_merge, aes(x = Mutation_Status, y = Th17_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.5, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "GSVA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("GSVA Enrichment of Th17 Enhanced Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_gsva_total, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# ssgsea Violin Plot - new comparison
setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_ssgsea_th17_enhanced_signature")
cbcols <- c("Other_Mut" = "#0000FF",
            "Mut" = "#FF0000")

ssgsea_total_output_mutation_merge$Mutation_Status <-factor(ssgsea_total_output_mutation_merge$Mutation_Status,
                                                          levels = c("Other_Mut", "Mut"))
my_comparisons <- list(c("Other_Mut", "Mut"))

cairo_pdf("./ssgsea_total_violin.pdf")
violin_1 <- ggplot(ssgsea_total_output_mutation_merge, aes(x = Mutation_Status, y = Th17_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.5, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "ssgsea Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("SSGSEA Enrichment of Th17 Enhanced Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_ssgsea_total, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# plage Violin Plot - new comparison
setwd("~/Datashare/TCGA_RNA_Analysis/Output/test_plage_th17_enhanced_signature")
cbcols <- c("Other_Mut" = "#0000FF",
            "Mut" = "#FF0000")

plage_total_output_mutation_merge$Mutation_Status <-factor(plage_total_output_mutation_merge$Mutation_Status,
                                                          levels = c("Other_Mut", "Mut"))
my_comparisons <- list(c("Other_Mut", "Mut"))

cairo_pdf("./plage_total_violin.pdf")
violin_1 <- ggplot(plage_total_output_mutation_merge, aes(x = Mutation_Status, y = Th17_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.5, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "plage Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("PLAGE Enrichment of Th17 Enhanced Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_plage_total, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()
