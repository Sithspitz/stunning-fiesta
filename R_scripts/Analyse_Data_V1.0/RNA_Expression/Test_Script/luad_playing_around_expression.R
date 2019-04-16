### Lung Adenocarcinoma PanCancer Atlas RNA Expression Playing Around ###
source("./R_scripts/Functions/functions.R")
source("./R_scripts/Functions/GeneSets.R")
setwd("~/DataShare/TCGA_RNA_Analysis/Testing/")
rna <- read.csv("luad_tcga_pancancer_rna.csv", header = T)




# Label Gene Signatures
Class2genes <- c("HLA-DOB", "HLA-DQA1", "HLA-DQA2", "HLA-DRB5", "HLA-DPB1",
                 "HLA-DMB", "HLA-DRA", "HLA-DMA", "HLA-DOA")

rna$Class2 <- ifelse(
  (rna$Hugo_Symbol %in% Class2genes), T, F)

rna$Th17 <- ifelse(
  (rna$Hugo_Symbol %in% Th17_genes), T, F)

class2datatemp <- CalGeneMeanRNASeqZScore(rna, Class2)
colnames(class2datatemp)[colnames(class2datatemp)=="value"] <- "Class_2_mean"

th17datatemp <- CalGeneMeanRNASeqZScore(rna, Th17)
colnames(th17datatemp)[colnames(th17datatemp)=="value"] <- "Th17_mean"

# Merge
expsig1 <- merge(th17datatemp, class2datatemp, by = "Patient.ID")

# Basic Heatmap
heat1 <- expsig1 %>% gather(-contains("Patient.ID"), key = "Geneset", value = "Zscore")
Bound2 <- acast(heat1, Patient.ID ~ Geneset, value.var = "Zscore")
mxCell <- t(Bound2)
xx <- mxCell[rowSums(!is.na(mxCell))!=0, colSums(!is.na(mxCell))!=0] 

hr <- hclust(as.dist(1-cor((xx), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.12))
clusterRows <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterRows[mycl]
my_palette <- colorRampPalette(c("blue", "black", "gold"))

heatmap.2(xx, 
          main = "Hierarchical Clustering of Gene Pathways",
          #Rowv = NA,
          Colv = as.dendrogram(hr),
          #dendrogram = "col",
          scale = "row",
          col = my_palette,
          density.info = "none",
          trace = "none",
          ColSideColors = myClusterSideBar,
          margins = c(12,12))

# Output with the correct clusters
test <- heat1
test1 <- spread(test, key = "Geneset", value = "Zscore")
test1$cluster <- mycl


# Doing a Linear Discriminant Analysis (LDA)
## dlyr note, if want to positively select then: 
### e.g. attrich <- expsig1 %>% dplyr:: select(contains(".ID"))
att2 <- expsig1 %>% dplyr:: select(-contains(".ID"))
att1 <- att2 %>% dplyr:: select(-contains("Subtype"))

att1[is.na(att1)] <- 0
pca <- prcomp(att1[,-1],
              center = TRUE,
              scale. = TRUE) 
tmp <- cor(att1[,-1])
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0
data.new <- att2[,!apply(tmp,2,function(x) any(abs(x) > 0.945))]



