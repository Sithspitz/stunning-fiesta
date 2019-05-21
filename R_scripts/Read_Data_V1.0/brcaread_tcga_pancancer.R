### DLBCL Read/Pre-proccessing ###
# Key commands
source("./R_scripts/Functions/functions.R")
setwd("M:/Richard B/TCGA_data/Breast Invasive Carcinoma/raw")



## RNA-Seq ##
# RSEM - RPKM EL
RNAEL1 <- read.csv("data_RNA_Seq_v2_expression_median.csv", header = T)

# Remove patients with NA
## Generally fix NA
RNAEL1[RNAEL1 == "#DIV/0!" | RNAEL1 == "NaN" | RNAEL1 == "#N/A" | RNAEL1 == "#VALUE!"] <- NA

# Keeps columns where NA is = 0
RNAEL2 <- RNAEL1[,colMeans(is.na(RNAEL1)) == 0]

# Convert all columns but Hugo_Symbol to numerics
Leave <- c("Hugo_Symbol")
selections <- !names(RNAEL2) %in% Leave
RNAEL2[,selections] <- as.numeric(as.character(unlist(RNAEL2[,selections])))
RNAEL3 <- droplevels(RNAEL2)

# Collate patients into a column called Patient.ID and values into 'RNASeqEL'
RNAEL4 <- RNAEL3 %>% gather(contains("TCGA"), key = "Patient.ID", value = "RNASeqEL")

# Factor the columns
NoFactor <- c("RNASeqEL")
RNAEL5 <- droplevels(nofactorthese(RNAEL4, NoFactor))

# Arrange into "DataType"
RNAEL6 <- RNAEL5 %>% gather(contains("RNASeqEL"), key = "DataType", value = "Value")

# Merge (rbind() here) and Write
setwd("M:/Richard B/TCGA_data/Breast Invasive Carcinoma/processed_csv")
write.csv(RNAEL6, "brca_tcga_pancancer_rna.csv", row.names = F)

