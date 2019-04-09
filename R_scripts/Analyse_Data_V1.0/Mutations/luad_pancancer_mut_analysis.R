### Lung Adenocarcinoma PanCancer Atlas Mutation Analysis ###
source("./R_scripts/Functions/functions.R")
# CHANGE TO CORRECT WD DIR setwd()
mut <- read.csv("luad_tcga_pancancer_mut.csv", header = T)

# Fix factors/numerics
## WRITE IF NEEDED



# Count number of mutations #

# Seperate - subset anything that isn't a "Silent mutation"
Intermediate <- droplevels(subset(mut, Variant_Classification != "Silent"))
Intermediate$parameter <- as.factor(paste(Intermediate$Patient.ID, Intermediate$Variant_Type
                                          , sep = ","))

# Count indiviual types of MT/patient in TCGA
patmut1 <- data.frame(parameter = character(),
                      Number = double(),
                      stringsAsFactors = F)
c <- 1
for(i in levels(Intermediate$parameter)){
  print(i)
  work <- droplevels(subset(Intermediate, parameter == i))
  num <- nrow(work)
  patmut1[c, "parameter"] <- i
  patmut1[c, "Number"] <- num
  c <- c + 1
}


patmut2 <- separate(patmut1, parameter, c("Patient.ID", "Variant_Type"), sep = ",")
patmut3 <- spread(patmut2, key = "Variant_Type", value = "Number")

# Count total number of MTs
patmut2$Patient.ID <- as.factor(patmut2$Patient.ID)
total <- data.frame(Patient.ID = character(),
                    TotalMut = double(),
                    stringsAsFactors = F)
c <- 1
for(i in levels(patmut2$Patient.ID)){
  print(i)
  work <- droplevels(subset(patmut2, Patient.ID == i))
  Total <- sum(work$Number)
  total[c, "Patient.ID"] <- i
  total[c, "TotalMut"] <- Total
  c <- c + 1 
}

MutNum <- merge(patmut3, total, by = "Patient.ID")

# Convert NA to 0 and change column order
MutNum[is.na(MutNum)] <- 0

# Write CSV
## CHANGE SET WD()
write.csv(MutNum, "mut_num_luad_tcga_pancancer.csv", row.names = F)



# MUTATIONS WHICH RESULT IN A PROTEIN CHANGE?!?!


