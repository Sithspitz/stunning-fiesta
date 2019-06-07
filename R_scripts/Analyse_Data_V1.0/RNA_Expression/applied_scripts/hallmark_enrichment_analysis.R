#### ----------------------------------------------------------------------------------- ####
############################### Hallmark Enrichment Analysis ################################
#### ----------------------------------------------------------------------------------- ####

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
setwd("H:/My Documents/BEAR Datashare/broad_institute_gene_signatures_gsea/Broad Institute Genesets/Hallmark")
geneset_gmt <- getGmt("hallmark_gmt_code.txt")

# Run the test
setwd("H:/My Documents/BEAR Datashare/broad_institute_gene_signatures_gsea/GSEA_comparisons/HALLMARK/Output/")

ssgsea_output <- gsva(data_matrix, geneset_gmt, method = "ssgsea")
write.csv(ssgsea_output, "hallmark_output.csv", row.names = T)

# Merge With Mutation Data
## Ordered by MT status
setwd("M:/Richard B/Analysis/2019/April 2019/Mutation_TCGA_Pancancer_Analysis/Subset")
mt_data <- read.csv("KRAS_STK11_plus_dbl_tidy_vs_other.csv", header = T)

ssgsea_transposed <- as.data.frame(t(ssgsea_output))
ssgsea_transposed <- rownames_to_column(ssgsea_transposed, var = "Patient.ID")  
ssgsea_merged <- merge(ssgsea_transposed, mt_data)    
ordered_df <- arrange(ssgsea_merged, as.character(Mutation_Status))

setwd("H:/My Documents/BEAR Datashare/broad_institute_gene_signatures_gsea/GSEA_comparisons/HALLMARK/Output/")
write.csv(ordered_df, "hallmark_output_merged.csv", row.names = F)  


### Stats Comparison and Plotting ###

# Screening exercise to see enrichment before plotting #

# Restructuring
## See mt_position_notes for col numbers to call during comparisons
gsea_df <- ordered_df
gsea_df$Mutation_Status <- NULL
gsea_df$Patient.ID <- NULL
gsea_trans <- t(gsea_df)
write.csv(gsea_trans, "hallmark_output_temp.csv", row.names = T)
df <- read.csv("hallmark_output_temp.csv", row.names = 1, header = T)
df2 <- as.matrix(df)
data_names <- df
file.remove("hallmark_output_temp.csv")

mt_position_notes <- c("Double MT are from col 1:36, KRAS MT are from col 37:154, Other MT are from 155:469, STK11 MT are from 470:506, Total MT are from cl c(1:154,470:506)")

# Wilcox test KRAS MT vs Other MT
kras_output <- vector()
for (i in 1:nrow(df2)) {
  print(rownames(data_names[i, ]))
  wil_tested <- wilcox.test(df2[i,37:154], df2[i,155:469])
  if(wil_tested$p.value <=0.05) {
    kras_output[length(kras_output)+1] <- rownames(data_names[i, ])
    kras_output[length(kras_output)+1] <- wil_tested$p.value
  }
}

kras_sig <- as.data.frame(kras_output)
setwd("H:/My Documents/BEAR Datashare/broad_institute_gene_signatures_gsea/GSEA_comparisons/HALLMARK/Output/")
write.csv(kras_sig, "kras_mt_enriched_pathways_code.csv", row.names = F)

# Wilcox test STK11 MT vs Other MT
stk11_output <- vector()
for (i in 1:nrow(df2)) {
  print(rownames(data_names[i, ]))
  wil_tested <- wilcox.test(df2[i,470:506], df2[i,155:469])
  if(wil_tested$p.value <=0.05) {
    stk11_output[length(stk11_output)+1] <- rownames(data_names[i, ])
    stk11_output[length(stk11_output)+1] <- wil_tested$p.value
  }
}

stk11_sig <- as.data.frame(stk11_output)
setwd("H:/My Documents/BEAR Datashare/broad_institute_gene_signatures_gsea/GSEA_comparisons/HALLMARK/Output/")
write.csv(stk11_sig, "stk11_mt_enriched_pathways_code.csv", row.names = F)

# Wilcox test Double MT vs Other MT
double_output <- vector()
for (i in 1:nrow(df2)) {
  print(rownames(data_names[i, ]))
  wil_tested <- wilcox.test(df2[i,1:36], df2[i,155:469])
  if(wil_tested$p.value <=0.05) {
    double_output[length(double_output)+1] <- rownames(data_names[i, ])
    double_output[length(double_output)+1] <- wil_tested$p.value
  }
}

double_sig <- as.data.frame(double_output)
setwd("H:/My Documents/BEAR Datashare/broad_institute_gene_signatures_gsea/GSEA_comparisons/HALLMARK/Output/")
write.csv(double_sig, "double_mt_enriched_pathways_code.csv", row.names = F)

# Wilcox test Total MT vs Other MT
total_output <- vector()
for (i in 1:nrow(df2)) {
  print(rownames(data_names[i, ]))
  wil_tested <- wilcox.test(df2[i,c(1:154,470:506)], df2[i,155:469])
  if(wil_tested$p.value <=0.05) {
    total_output[length(total_output)+1] <- rownames(data_names[i, ])
    total_output[length(total_output)+1] <- wil_tested$p.value
  }
}

total_sig <- as.data.frame(total_output)
setwd("H:/My Documents/BEAR Datashare/broad_institute_gene_signatures_gsea/GSEA_comparisons/HALLMARK/Output/")
write.csv(total_sig, "total_mt_enriched_pathways_code.csv", row.names = F)

### Then look through the data and decode the 'GS' GENESET numbers and work out what to plot ###

