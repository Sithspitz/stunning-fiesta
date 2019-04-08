### Lung Adenocarcinoma PanCancer Atlas Read/Pre-proccessing ###
# Key commands
source("./R_scripts/Functions/functions.R")
setwd("L:/Richard B/TCGA_data/Pancancer/raw_csv")



## Clinical Read and Tidy ##

clinical_pat <- read.csv("data_clinical.csv", header = T)

# Remove uneeded columns
uneeded <- c("Sample.ID", "Cancer.Study", "Cancer.Type.Detailed", "Sample.Type", 
             "Ethnicity.Category", "Race.Category", "Subtype", 
             "Last.Alive.Less.Initial.Pathologic.Diagnosis.Date.Calculated.Day.Value",
             "Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date")
clin1 <- droplevels(clinical_pat[,!(names(clinical_pat) %in% uneeded)])

# Fix factors/numerics
## DontFactor vector is columns with numeric values - R may call the str() of these as integers
### The numericthese function converts any integers to numeric values
DontFactor <- c("Mutation.Count", "Fraction.Genome.Altered", "Diagnosis.Age", "Aneuploidy.Score")
clin2 <- nofactorthese(clin1, DontFactor)
clin3 <- numericthese(clin2, DontFactor)

# Change Patient IDs from - to . and write csv
clin3$Patient.ID <- gsub("-", "\\.", clin3$Patient.ID)
luad_tcga_pancancer_clinical <- clin3
setwd("L:/Richard B/R_WD/stunning-fiesta/Test_Output_WD")
write.csv(luad_tcga_pancancer_clinical, "luad_tcga_pancancer_clinical.csv", row.names = F)
setwd("L:/Richard B/TCGA_data/Pancancer/raw_csv")



## Copy Number Read and Tidy ##

# Create a blank df to allocate memory to read the CNA .csv
blankdf <- data.frame()[1:12000000, ]
CNA <- read.csv("data_CNA.csv", header = T)

# Collate patients into a column called "Patient.ID" and values into "CNA_NL" 
## Correct strcture of the TCGA Patient ID by changing - to . and removing the last three characters from each element
### Re-order columns
CNA1 <- CNA %>% gather(contains("TCGA"), key = "Patient.ID", value = "CNA_NL")
# ARCHIVE AS DONE THIS STEP ON EXCEL FOR THIS DATA: CNA1$Patient.ID <- gsub("-", "\\.", CNA1$Patient.ID)
# ARCHIVE AS DONE THIS STEP ON EXCEL FOR THIS DATA: CNA1$Patient.ID <- gsub(".{3}$", "", CNA1$Patient.ID)
CNA1 <- CNA1[,c(3,1,2,4)]

# Fix factors/numerics
DoFactor <- c("Patient.ID", "Entrez_Gene_Id")
CNA2 <- factorthese(CNA1, DoFactor)
NoFactor <- c("CNA_NL")
CNA3 <- nofactorthese(CNA2, NoFactor)

# Write csv
luad_tcga_pancancer_cna <- CNA3
setwd("L:/Richard B/R_WD/stunning-fiesta/Test_Output_WD")
write.csv(luad_tcga_pancancer_cna, "luad_tcga_pancancer_cna.csv", row.names = F)
rm(blankdf)


