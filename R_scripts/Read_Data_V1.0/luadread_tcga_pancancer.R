### Lung Adenocarcinoma PanCancer Atlas Read/Pre-proccessing ###
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

