### Lung Adenocarcinoma PanCancer Atlas Mutation Analysis ###
source("./R_scripts/Functions/functions.R")
setwd("~/DataShare/Mut_TCGA_Temp/Mutation_TCGA_Pancancer_Analysis/Input")
mut <- read.csv("luad_tcga_pancancer_mut.csv", header = T)



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
setwd("~/DataShare/Mut_TCGA_Temp/Mutation_TCGA_Pancancer_Analysis/Output")
write.csv(MutNum, "mut_num_luad_tcga_pancancer.csv", row.names = F)



# Most common variants in the TCGA #

# Firstly take the variants I am interested in
## Remove the levels I don't care about using droplevels()
Intermediate2 <- droplevels(subset(Intermediate, Variant_Classification != "RNA"))
Intermediate2 <- droplevels(subset(Intermediate2, Variant_Classification != "Intron"))
Intermediate2 <- droplevels(subset(Intermediate2, Variant_Classification != "3'UTR"))
Intermediate2 <- droplevels(subset(Intermediate2, Variant_Classification != "5'UTR"))
Intermediate2 <- droplevels(subset(Intermediate2, Variant_Classification != "5'Flank"))
Intermediate2 <- droplevels(subset(Intermediate2, Variant_Classification != "3'Flank"))
Intermediate2 <- droplevels(subset(Intermediate2, Variant_Classification != "Translation_Start_Site"))

Intermediate3 <- Intermediate2
Intermediate4 <- Intermediate2


Intermediate2$class <- as.factor(paste(Intermediate2$Patient.ID, Intermediate2$Variant_Classification, 
                                    sep = ","))
Intermediate3$protein_change <- as.factor(paste(Intermediate3$Patient.ID, Intermediate3$HGVSp_Short, 
                                       sep = ","))
Intermediate4$amino_acid_change <- as.factor(paste(Intermediate4$Patient.ID, Intermediate4$Amino_acids, 
                                                sep = ","))





# Count indivual types of alteration, protein change and amino acid change
var1 <- data.frame(parameter = character(),
                   Number = double(),
                   stringsAsFactors = F)
var2 <- data.frame(parameter = character(),
                   Number = double(),
                   stringsAsFactors = F)
var3 <- data.frame(parameter = character(),
                   Number = double(),
                   stringsAsFactors = F)

c <- 1
for(i in levels(Intermediate2$class)){
  print(i)
  work <- droplevels(subset(Intermediate2, class == i))
  num <- nrow(work)
  var1[c, "parameter"] <- i
  var1[c, "Number"] <- num
  c <- c + 1
}

c <- 1
for(i in levels(Intermediate3$protein_change)){
  print(i)
  work <- droplevels(subset(Intermediate3, protein_change == i))
  num <- nrow(work)
  var2[c, "parameter"] <- i
  var2[c, "Number"] <- num
  c <- c + 1
}

c <- 1
for(i in levels(Intermediate4$amino_acid_change)){
  print(i)
  work <- droplevels(subset(Intermediate4, amino_acid_change == i))
  num <- nrow(work)
  var3[c, "parameter"] <- i
  var3[c, "Number"] <- num
  c <- c + 1
}

# Seperate and spread variant class
## Remove NA
### Write CSV
var1_1 <- separate(var1, parameter, c("Patient.ID", "Variant_Classification"), sep = ",")
varclass <- spread(var1_1, key = "Variant_Classification", value = "Number")

varclass[is.na(varclass)] <- 0

setwd("~/DataShare/Mut_TCGA_Temp/Mutation_TCGA_Pancancer_Analysis/Output")
write.csv(varclass, "mut_class_luad_tcga_pancancer.csv", row.names = F)

# Seperate and spread protein changes
## Remove NA
### Write CSV - protein_change....csv is excel readable row-wise and protein_change...2.csv is not (column-wise), but same data
var2_1 <- separate(var2, parameter, c("Patient.ID", "Variant_Classification"), sep = ",")
proteinchange1 <- spread(var2_1, key = "Variant_Classification", value = "Number")

proteinchange1[is.na(proteinchange1)] <- 0

setwd("~/DataShare/Mut_TCGA_Temp/Mutation_TCGA_Pancancer_Analysis/Output")
write.csv(var2_1, "protein_change_luad_tcga_pancancer.csv", row.names = F)
write.csv(proteinchange1, "protein_change_luad_tcga_pancancer2.csv", row.names = F)

# Seperate and spread amino acid changes
## Remove NA
### Write CSV
var3_1 <- separate(var3, parameter, c("Patient.ID", "Variant_Classification"), sep = ",")
aminoacidchange1 <- spread(var3_1, key = "Variant_Classification", value = "Number")

aminoacidchange1[is.na(aminoacidchange1)] <- 0

setwd("~/DataShare/Mut_TCGA_Temp/Mutation_TCGA_Pancancer_Analysis/Output")
write.csv(aminoacidchange1, "amino_acid_change_luad_tcga_pancancer.csv", row.names = F)

