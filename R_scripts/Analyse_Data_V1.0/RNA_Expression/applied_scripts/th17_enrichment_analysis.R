#### ----------------------------------------------------------------------------------- ####
############################### Th17 Enrichment Analysis ####################################
#### ----------------------------------------------------------------------------------- ####

##### This script is in this repo as an example script ######################################

mypackages <- c("GSEABase", "GSVA", "Biobase", "genefilter",
                "limma", "RColorBrewer", "GSVAdata", "scales", "dplyr")
lapply(mypackages, library, character.only = T)
source("M:/Richard B/R_WD/jack_richard_functions/functions.R")



### Enrichment Analysis ###

# Import Expression Data
## Make into 'Data Matrix' for GSVA
setwd("M:/Richard B/TCGA_data/Pancancer/gsea_correct_format")
df <- read.csv("total_rna_seq_EL_expression_data.csv", header = T, row.names = 1)
data_matrix <- as.matrix(df)

# Import Gene Set (GMT Text Format, save an Excel doc as .txt)
setwd("M:/Richard B/Analysis/2019/May 2019/gsea_th17_signatures/Input")
th17_geneset_gmt <- getGmt("th17_genesets_gmt.txt")

# Run the test
## PLAGE, SSGSEA, GSVA
setwd("M:/Richard B/Analysis/2019/May 2019/gsea_th17_signatures/Output")

gsva_output <- gsva(data_matrix, th17_geneset_gmt, method = "gsva")
ssgsea_output <- gsva(data_matrix, th17_geneset_gmt, method = "ssgsea")
plage_output <- gsva(data_matrix, th17_geneset_gmt, method = "plage")

write.csv(gsva_output, "gsva_output.csv", row.names = T)
write.csv(ssgsea_output, "ssgsea_output.csv", row.names = T)
write.csv(plage_output, "plage_output.csv", row.names = T)

# Merge With Mutation Data
setwd("M:/Richard B/Analysis/2019/April 2019/Mutation_TCGA_Pancancer_Analysis/Subset")
mt_data <- read.csv("Mut_tidy_vs_other.csv", header = T)

gsva_transposed <- as.data.frame(t(gsva_output))
gsva_transposed <- rownames_to_column(gsva_transposed, var = "Patient.ID")  
gsva_merged <- merge(gsva_transposed, mt_data)  
  
ssgsea_transposed <- as.data.frame(t(ssgsea_output))
ssgsea_transposed <- rownames_to_column(ssgsea_transposed, var = "Patient.ID")  
ssgsea_merged <- merge(ssgsea_transposed, mt_data)    
  
plage_transposed <- as.data.frame(t(plage_output))
plage_transposed <- rownames_to_column(plage_transposed, var = "Patient.ID")  
plage_merged <- merge(plage_transposed, mt_data)    
  
setwd("M:/Richard B/Analysis/2019/May 2019/gsea_th17_signatures/Output")
write.csv(gsva_merged, "gsva_output_merged.csv", row.names = F)  
write.csv(ssgsea_merged, "ssgsea_output_merged.csv", row.names = F)  
write.csv(plage_merged, "plage_output_merged.csv", row.names = F)  



### Stats Comparison and Plotting ###

# Subsetting Based on Mutational Status
gsva_merged$Patient.ID <- NULL
ssgsea_merged$Patient.ID <- NULL
plage_merged$Patient.ID <- NULL

gsva_wt_data <- droplevels(subset(gsva_merged, Mutation_Status == "Other_Mut"))
gsva_mt_data <- droplevels(subset(gsva_merged, Mutation_Status == "Mut"))

ssgsea_wt_data <- droplevels(subset(ssgsea_merged, Mutation_Status == "Other_Mut"))
ssgsea_mt_data <- droplevels(subset(ssgsea_merged, Mutation_Status == "Mut"))

plage_wt_data <- droplevels(subset(plage_merged, Mutation_Status == "Other_Mut"))
plage_mt_data <- droplevels(subset(plage_merged, Mutation_Status == "Mut"))

# Subsetting Based on Comparison

gsva_th17_standard_wt <- gsva_wt_data$Th17_Standard
gsva_th17_standard_mt <- gsva_mt_data$Th17_Standard
gsva_th17_enhanced_wt <- gsva_wt_data$Th17_Enhanced
gsva_th17_enhanced_mt <- gsva_mt_data$Th17_Enhanced
gsva_th17_codeset_wt <- gsva_wt_data$Th17_CodeSet
gsva_th17_codeset_mt <- gsva_mt_data$Th17_CodeSet

ssgsea_th17_standard_wt <- ssgsea_wt_data$Th17_Standard
ssgsea_th17_standard_mt <- ssgsea_mt_data$Th17_Standard
ssgsea_th17_enhanced_wt <- ssgsea_wt_data$Th17_Enhanced
ssgsea_th17_enhanced_mt <- ssgsea_mt_data$Th17_Enhanced
ssgsea_th17_codeset_wt <- ssgsea_wt_data$Th17_CodeSet
ssgsea_th17_codeset_mt <- ssgsea_mt_data$Th17_CodeSet

plage_th17_standard_wt <- plage_wt_data$Th17_Standard
plage_th17_standard_mt <- plage_mt_data$Th17_Standard
plage_th17_enhanced_wt <- plage_wt_data$Th17_Enhanced
plage_th17_enhanced_mt <- plage_mt_data$Th17_Enhanced
plage_th17_codeset_wt <- plage_wt_data$Th17_CodeSet
plage_th17_codeset_mt <- plage_mt_data$Th17_CodeSet

# Wilcox Tests
gsva_statresult_th17_standard <- wilcox.test(gsva_th17_standard_wt, gsva_th17_standard_mt)
gsva_statresult_th17_enhanced <- wilcox.test(gsva_th17_enhanced_wt, gsva_th17_enhanced_mt)
gsva_statresult_th17_codeset <- wilcox.test(gsva_th17_codeset_wt, gsva_th17_codeset_mt)
gsva_statresult_th17_standard
gsva_statresult_th17_enhanced
gsva_statresult_th17_codeset

ssgsea_statresult_th17_standard <- wilcox.test(ssgsea_th17_standard_wt, ssgsea_th17_standard_mt)
ssgsea_statresult_th17_enhanced <- wilcox.test(ssgsea_th17_enhanced_wt, ssgsea_th17_enhanced_mt)
ssgsea_statresult_th17_codeset <- wilcox.test(ssgsea_th17_codeset_wt, ssgsea_th17_codeset_mt)
ssgsea_statresult_th17_standard
ssgsea_statresult_th17_enhanced
ssgsea_statresult_th17_codeset

plage_statresult_th17_standard <- wilcox.test(plage_th17_standard_wt, plage_th17_standard_mt)
plage_statresult_th17_enhanced <- wilcox.test(plage_th17_enhanced_wt, plage_th17_enhanced_mt)
plage_statresult_th17_codeset <- wilcox.test(plage_th17_codeset_wt, plage_th17_codeset_mt)
plage_statresult_th17_standard
plage_statresult_th17_enhanced
plage_statresult_th17_codeset

pVal_1_gsva_th17_standard <- wilcox.test(gsva_th17_standard_wt, gsva_th17_standard_mt)$p.value
pVal_2_gsva_th17_standard <- format(round(pVal_1_gsva_th17_standard, 4), nsmall = 4)
pVal_3_gsva_th17_standard <- paste("p = ", pVal_2_gsva_th17_standard, sep = "")

pVal_1_gsva_th17_enhanced <- wilcox.test(gsva_th17_enhanced_wt, gsva_th17_enhanced_mt)$p.value
pVal_2_gsva_th17_enhanced <- format(round(pVal_1_gsva_th17_enhanced, 4), nsmall = 4)
pVal_3_gsva_th17_enhanced <- paste("p = ", pVal_2_gsva_th17_enhanced, sep = "")

pVal_1_gsva_th17_codeset <- wilcox.test(gsva_th17_codeset_wt, gsva_th17_codeset_mt)$p.value
pVal_2_gsva_th17_codeset <- format(round(pVal_1_gsva_th17_codeset, 4), nsmall = 4)
pVal_3_gsva_th17_codeset <- paste("p = ", pVal_2_gsva_th17_codeset, sep = "")

pVal_1_ssgsea_th17_standard <- wilcox.test(ssgsea_th17_standard_wt, ssgsea_th17_standard_mt)$p.value
pVal_2_ssgsea_th17_standard <- format(round(pVal_1_ssgsea_th17_standard, 4), nsmall = 4)
pVal_3_ssgsea_th17_standard <- paste("p = ", pVal_2_ssgsea_th17_standard, sep = "")

pVal_1_ssgsea_th17_enhanced <- wilcox.test(ssgsea_th17_enhanced_wt, ssgsea_th17_enhanced_mt)$p.value
pVal_2_ssgsea_th17_enhanced <- format(round(pVal_1_ssgsea_th17_enhanced, 4), nsmall = 4)
pVal_3_ssgsea_th17_enhanced <- paste("p = ", pVal_2_ssgsea_th17_enhanced, sep = "")

pVal_1_ssgsea_th17_codeset <- wilcox.test(ssgsea_th17_codeset_wt, ssgsea_th17_codeset_mt)$p.value
pVal_2_ssgsea_th17_codeset <- format(round(pVal_1_ssgsea_th17_codeset, 4), nsmall = 4)
pVal_3_ssgsea_th17_codeset <- paste("p = ", pVal_2_ssgsea_th17_codeset, sep = "")

pVal_1_plage_th17_standard <- wilcox.test(plage_th17_standard_wt, plage_th17_standard_mt)$p.value
pVal_2_plage_th17_standard <- format(round(pVal_1_plage_th17_standard, 4), nsmall = 4)
pVal_3_plage_th17_standard <- paste("p = ", pVal_2_plage_th17_standard, sep = "")

pVal_1_plage_th17_enhanced <- wilcox.test(plage_th17_enhanced_wt, plage_th17_enhanced_mt)$p.value
pVal_2_plage_th17_enhanced <- format(round(pVal_1_plage_th17_enhanced, 4), nsmall = 4)
pVal_3_plage_th17_enhanced <- paste("p = ", pVal_2_plage_th17_enhanced, sep = "")

pVal_1_plage_th17_codeset <- wilcox.test(plage_th17_codeset_wt, plage_th17_codeset_mt)$p.value
pVal_2_plage_th17_codeset <- format(round(pVal_1_plage_th17_codeset, 4), nsmall = 4)
pVal_3_plage_th17_codeset <- paste("p = ", pVal_2_plage_th17_codeset, sep = "")

# Pre-Plotting Cleaning Up and factoring
gsva_merged2 <- gsva_merged
ssgsea_merged2 <- ssgsea_merged
plage_merged2 <- plage_merged

gsva_merged2$Mutation_Status <- gsub(pattern = "Other_Mut", replacement = "WT", x = gsva_merged2$Mutation_Status)
gsva_merged2$Mutation_Status <- gsub(pattern = "Mut", replacement = "MT", x = gsva_merged2$Mutation_Status)

ssgsea_merged2$Mutation_Status <- gsub(pattern = "Other_Mut", replacement = "WT", x = ssgsea_merged2$Mutation_Status)
ssgsea_merged2$Mutation_Status <- gsub(pattern = "Mut", replacement = "MT", x = ssgsea_merged2$Mutation_Status)

plage_merged2$Mutation_Status <- gsub(pattern = "Other_Mut", replacement = "WT", x = plage_merged2$Mutation_Status)
plage_merged2$Mutation_Status <- gsub(pattern = "Mut", replacement = "MT", x = plage_merged2$Mutation_Status)

gsva_merged2$Mutation_Status <- factor(gsva_merged2$Mutation_Status, levels = c("WT", "MT"))
ssgsea_merged2$Mutation_Status <- factor(ssgsea_merged2$Mutation_Status, levels = c("WT", "MT"))
plage_merged2$Mutation_Status <- factor(plage_merged2$Mutation_Status, levels = c("WT", "MT"))

# Violin Plot GSVA Th17 Standard
setwd("M:/Richard B/Analysis/2019/May 2019/gsea_th17_signatures/Output")
cbcols <- c("WT" = "#0000FF",
            "MT" = "#FF0000")
my_comparisons <- list(c("WT", "MT"))

cairo_pdf("./th17_standard_gsva_violin.pdf")
violin_1 <- ggplot(gsva_merged2, aes(x = Mutation_Status, y = Th17_Standard)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "GSVA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("GSVA Enrichment of Th17 Standard Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_gsva_th17_standard, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()
# Button to delete the violin file: unlink("gsva_violin.pdf")

# Violin Plot GSVA Th17 Enhanced
cairo_pdf("./th17_enhanced_gsva_violin.pdf")
violin_1 <- ggplot(gsva_merged2, aes(x = Mutation_Status, y = Th17_Enhanced)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
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
  annotate("text", label = pVal_3_gsva_th17_enhanced, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# Violin Plot GSVA Th17 CodeSet
cairo_pdf("./th17_codeset_gsva_violin.pdf")
violin_1 <- ggplot(gsva_merged2, aes(x = Mutation_Status, y = Th17_CodeSet)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "GSVA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("GSVA Enrichment of Th17 CodeSet Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_gsva_th17_codeset, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# Violin Plot SSGSEA Th17 Standard
cairo_pdf("./th17_standard_ssgsea_violin.pdf")
violin_1 <- ggplot(ssgsea_merged2, aes(x = Mutation_Status, y = Th17_Standard)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "SSGSEA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("SSGSEA Enrichment of Th17 Standard Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_ssgsea_th17_standard, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# Violin Plot SSGSEA Th17 Enhanced
cairo_pdf("./th17_enhanced_ssgsea_violin.pdf")
violin_1 <- ggplot(ssgsea_merged2, aes(x = Mutation_Status, y = Th17_Enhanced)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
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
  annotate("text", label = pVal_3_ssgsea_th17_enhanced, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# Violin Plot SSGSEA Th17 CodeSet
cairo_pdf("./th17_codeset_ssgsea_violin.pdf")
violin_1 <- ggplot(ssgsea_merged2, aes(x = Mutation_Status, y = Th17_CodeSet)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "SSGSEA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("SSGSEA Enrichment of Th17 CodeSet Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_ssgsea_th17_codeset, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# Violin Plot PLAGE Th17 Standard
cairo_pdf("./th17_standard_plage_violin.pdf")
violin_1 <- ggplot(plage_merged2, aes(x = Mutation_Status, y = Th17_Standard)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "PLAGE Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("PLAGE Enrichment of Th17 Standard Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_plage_th17_standard, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# Violin Plot PLAGE Th17 Enhanced
cairo_pdf("./th17_enhanced_plage_violin.pdf")
violin_1 <- ggplot(plage_merged2, aes(x = Mutation_Status, y = Th17_Enhanced)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "PLAGE Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("PLAGE Enrichment of Th17 Enhanced Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_plage_th17_enhanced, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# Violin Plot PLAGE Th17 CodeSet
cairo_pdf("./th17_codeset_plage_violin.pdf")
violin_1 <- ggplot(plage_merged2, aes(x = Mutation_Status, y = Th17_CodeSet)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(Mutation_Status, fill = Mutation_Status),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "PLAGE Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("PLAGE Enrichment of Th17 CodeSet Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_plage_th17_codeset, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

