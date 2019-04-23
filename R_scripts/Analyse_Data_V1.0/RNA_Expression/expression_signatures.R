### Lung Adenocarcinoma PanCancer Atlas RNA Expression ###
## Should be adopted dependent on the analysis in question ##
source("./R_scripts/Functions/functions.R")
source("./R_scripts/Functions/GeneSets.R")
setwd("~/DataShare/TCGA_RNA_Analysis/Input/")
rna <- read.csv("luad_tcga_pancancer_rna.csv", header = T)



# Label Immune Gene Signatures
rna$aDC <- ifelse(
  (rna$Hugo_Symbol %in% aDC_genes), T, F)
rna$B_Cells <- ifelse(
  (rna$Hugo_Symbol %in% BCells_genes), T, F)
rna$CD8 <- ifelse(
  (rna$Hugo_Symbol %in% CD8_genes), T, F)
rna$Cytotoxic <- ifelse(
  (rna$Hugo_Symbol %in% Cytotoxic_genes), T, F)
rna$DC <- ifelse(
  (rna$Hugo_Symbol %in% DC_genes), T, F)
rna$Eosinophils <- ifelse(
  (rna$Hugo_Symbol %in% Eosinophils_genes), T, F)
rna$iDC <- ifelse(
  (rna$Hugo_Symbol %in% iDC_genes), T, F)
rna$Macrophages <- ifelse(
  (rna$Hugo_Symbol %in% Macrophages_genes), T, F)
rna$Mast_Cells <- ifelse(
  (rna$Hugo_Symbol %in% MastCells_genes), T, F)
rna$Neutrophils <- ifelse(
  (rna$Hugo_Symbol %in% Neutrophils_genes), T, F)
rna$NK <- ifelse(
  (rna$Hugo_Symbol %in% NK_genes), T, F)
rna$CD56_Bright_NK <- ifelse(
  (rna$Hugo_Symbol %in% NKCD56bright_genes), T, F)
rna$CD56_Dim_NK <- ifelse(
  (rna$Hugo_Symbol %in% NKCD56dim_genes), T, F)
rna$Tcm <- ifelse(
  (rna$Hugo_Symbol %in% Tcm_genes), T, F)
rna$Tem <- ifelse(
  (rna$Hugo_Symbol %in% Tem_genes), T, F)
rna$Tfh <- ifelse(
  (rna$Hugo_Symbol %in% Tfh_genes), T, F)
rna$Tgd <- ifelse(
  (rna$Hugo_Symbol %in% Tgd_genes), T, F)
rna$Th1 <- ifelse(
  (rna$Hugo_Symbol %in% Th1_genes), T, F)
rna$Th17 <- ifelse(
  (rna$Hugo_Symbol %in% Th17_genes), T, F)
rna$Th2 <- ifelse(
  (rna$Hugo_Symbol %in% Th2_genes), T, F)
rna$Tregs <- ifelse(
  (rna$Hugo_Symbol %in% Tregs_genes), T, F)

# Calculate mean scores for given genesets
aDCdata <- CalGeneMeanRNASeqZScore(rna, aDC)
colnames(aDCdata)[colnames(aDCdata)=="value"] <- "aDC_Mean"
Bcelldata <- CalGeneMeanRNASeqZScore(rna, B_Cells)
colnames(Bcelldata)[colnames(Bcelldata)=="value"] <- "B_Cell_Mean"
CD8data <- CalGeneMeanRNASeqZScore(rna, CD8)
colnames(CD8data)[colnames(CD8data)=="value"] <- "CD8_Mean"
Cytotoxicdata <- CalGeneMeanRNASeqZScore(rna, Cytotoxic)
colnames(Cytotoxicdata)[colnames(Cytotoxicdata)=="value"] <- "Cytotoxic_Mean"
DCdata <- CalGeneMeanRNASeqZScore(rna, DC)
colnames(DCdata)[colnames(DCdata)=="value"] <- "DC_Mean"
Eosinophilsdata <- CalGeneMeanRNASeqZScore(rna, Eosinophils)
colnames(Eosinophilsdata)[colnames(Eosinophilsdata)=="value"] <- "Eosinophil_Mean"
iDCdata <- CalGeneMeanRNASeqZScore(rna, iDC)
colnames(iDCdata)[colnames(iDCdata)=="value"] <- "iDC_Mean"
Macrophagesdata <- CalGeneMeanRNASeqZScore(rna, Macrophages)
colnames(Macrophagesdata)[colnames(Macrophagesdata)=="value"] <- "Macrophage_Mean"
Mastcelldata <- CalGeneMeanRNASeqZScore(rna, Mast_Cells)
colnames(Mastcelldata)[colnames(Mastcelldata)=="value"] <- "Mast_Cell_Mean"
Neutrophilsdata <- CalGeneMeanRNASeqZScore(rna, Neutrophils)
colnames(Neutrophilsdata)[colnames(Neutrophilsdata)=="value"] <- "Neutrophil_Mean"
NKdata <- CalGeneMeanRNASeqZScore(rna, NK)
colnames(NKdata)[colnames(NKdata)=="value"] <- "NK_Mean"
NKCD56brightdata <- CalGeneMeanRNASeqZScore(rna, CD56_Bright_NK)
colnames(NKCD56brightdata)[colnames(NKCD56brightdata)=="value"] <- "CD56_Bright_NK_Mean"
NKCD56dimdata <- CalGeneMeanRNASeqZScore(rna, CD56_Dim_NK)
colnames(NKCD56dimdata)[colnames(NKCD56dimdata)=="value"] <- "CD56_Dim_NK_Mean"
Tcmdata <- CalGeneMeanRNASeqZScore(rna, Tcm)
colnames(Tcmdata)[colnames(Tcmdata)=="value"] <- "Tcm_Mean"
Temdata <- CalGeneMeanRNASeqZScore(rna, Tem)
colnames(Temdata)[colnames(Temdata)=="value"] <- "Tem_Mean"
Tfhdata <- CalGeneMeanRNASeqZScore(rna, Tfh)
colnames(Tfhdata)[colnames(Tfhdata)=="value"] <- "Tfh_Mean"
Tgddata <- CalGeneMeanRNASeqZScore(rna, Tgd)
colnames(Tgddata)[colnames(Tgddata)=="value"] <- "Tgd_Mean"
Th1data <- CalGeneMeanRNASeqZScore(rna, Th1)
colnames(Th1data)[colnames(Th1data)=="value"] <- "Th1_Mean"
Th17data <- CalGeneMeanRNASeqZScore(rna, Th17)
colnames(Th17data)[colnames(Th17data)=="value"] <- "Th17_Mean"
Th2data <- CalGeneMeanRNASeqZScore(rna, Th2)
colnames(Th2data)[colnames(Th2data)=="value"] <- "Th2_Mean"
Tregsdata <- CalGeneMeanRNASeqZScore(rna, Tregs)
colnames(Tregsdata)[colnames(Tregsdata)=="value"] <- "Treg_Mean"

# Merge
m1 <- merge(aDCdata, Bcelldata, by = "Patient.ID")
m2 <- merge(m1, CD8data, by = "Patient.ID")
m3 <- merge(m2, Cytotoxicdata, by = "Patient.ID")
m4 <- merge(m3, DCdata, by = "Patient.ID")
m5 <- merge(m4, Eosinophilsdata, by = "Patient.ID")
m6 <- merge(m5, iDCdata, by = "Patient.ID")
m7 <- merge(m6, Macrophagesdata, by = "Patient.ID")
m8 <- merge(m7, Mastcelldata, by = "Patient.ID")
m9 <- merge(m8, Neutrophilsdata, by = "Patient.ID")
m10 <- merge(m9, NKdata, by = "Patient.ID")
m11 <- merge(m10, NKCD56brightdata, by = "Patient.ID")
m12 <- merge(m11, NKCD56dimdata, by = "Patient.ID")
m13 <- merge(m12, Tcmdata, by = "Patient.ID")
m14 <- merge(m13, Temdata, by = "Patient.ID")
m15 <- merge(m14, Tfhdata, by = "Patient.ID")
m16 <- merge(m15, Tgddata, by = "Patient.ID")
m17 <- merge(m16, Th1data, by = "Patient.ID")
m18 <- merge(m17, Th17data, by = "Patient.ID")
m19 <- merge(m18, Th2data, by = "Patient.ID")
m20 <- merge(m19, Tregsdata, by = "Patient.ID")
rm(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, 
   m12, m13, m14, m15, m16, m17, m18, m19)

# Heatmap Organisation
heat1 <- m20 %>% gather(-contains("Patient.ID"), key = "Geneset", value = "Zscore")
Bound2 <- acast(heat1, Patient.ID ~ Geneset, value.var = "Zscore")
mxCell <- t(Bound2)
xx <- mxCell[rowSums(!is.na(mxCell))!=0, colSums(!is.na(mxCell))!=0] 

# Heatmap Clustering 
hr <- hclust(as.dist(1-cor((xx), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.12))
clusterRows <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterRows[mycl]
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

setwd("~/DataShare/TCGA_RNA_Analysis/Output/")
cairo_pdf("./cytotoxic_immune_signature.pdf")
heatmap.2(xx,
          Rowv = F,
          Colv = as.dendrogram(hr),
          dendrogram = "column",
          #dendrogram = "col",
          scale = "row",
          col = my_palette,
          density.info = "none",
          trace = "none",
          ColSideColors = myClusterSideBar,
          key = F,
          labCol = F,
          margins = c(13,12))
dev.off()


# Output with the correct clusters
test <- heat1
test1 <- spread(test, key = "Geneset", value = "Zscore")
test1$cluster <- mycl
write.csv(test1, "luad_pancancer_total_immune_gene_clusters.csv", row.names = F)



# Not currently used: Doing a Linear Discriminant Analysis (LDA)
## dlyr note, if want to positively select then: 
### e.g. attrich <- expsig1 %>% dplyr:: select(contains(".ID"))
# att2 <- expsig1 %>% dplyr:: select(-contains(".ID"))
# att1 <- att2 %>% dplyr:: select(-contains("Subtype"))

