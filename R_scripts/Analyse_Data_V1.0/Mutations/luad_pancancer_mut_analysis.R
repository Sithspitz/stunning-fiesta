### Lung Adenocarcinoma PanCancer Atlas Mutation Analysis ###
source("./R_scripts/Functions/functions.R")
setwd("L:/Richard B/TCGA_data/Pancancer/processed_csv")
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
setwd("H:/My Documents/Analysis/2019/April 2019/temp/Mutation_TCGA_Pancancer_Analysis/Output")
write.csv(MutNum, "mut_num_luad_tcga_pancancer.csv", row.names = F)



# Most common variants in the TCGA #

# Firstly take the variants I am interested in
## Variant types as shown in 'levels' vector
levels <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
            "In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins",
            "Splice_Site", "Splice_Region")
var_subset <- subset(mut, Variant_Classification == levels, select = c(Patient.ID:PolyPhen))
var_subset$class <- as.factor(paste(var_subset$Patient.ID, var_subset$Variant_Classification, 
                                    sep = ","))
var_subset$protein_change <- as.factor(paste(var_subset$Patient.ID, var_subset$HGVSp_Short,
                                             sep = ","))
var_subset$amino_acid_change <- as.factor(paste(var_subset$Patient.ID, var_subset$Amino_acids,
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
for(i in levels(var_subset$class)){
  print(i)
  work <- droplevels(subset(var_subset, class == i))
  num <- nrow(work)
  var1[c, "parameter"] <- i
  var1[c, "Number"] <- num
  c <- c + 1
}

c <- 1
for(i in levels(var_subset$protein_change)){
  print(i)
  work <- droplevels(subset(var_subset, protein_change == i))
  num <- nrow(work)
  var2[c, "parameter"] <- i
  var2[c, "Number"] <- num
  c <- c + 1
}

c <- 1
for(i in levels(var_subset$amino_acid_change)){
  print(i)
  work <- droplevels(subset(var_subset, amino_acid_change == i))
  num <- nrow(work)
  var3[c, "parameter"] <- i
  var3[c, "Number"] <- num
  c <- c + 1
}

# Seperate and spread variant class
var1_1 <- separate(var1, parameter, c("Patient.ID", "Variant_Classification"), sep = ",")
var1_2 <- spread(var1_1, key = "Variant_Classification", value = "Number")

# Count total variant class 

######### NOT SURE ANYTHING BELOW HERE IS WORKING
######## OR IS EVEN NEEDED?!?!
var1_2$Patient.ID <- as.factor(var1_2$Patient.ID)
total1 <- data.frame(Patient.ID = character(),
                    TotalMut = double(),
                    stringsAsFactors = F)
c <- 1
for(i in levels(var1_2$Patient.ID)){
  print(i)
  work <- droplevels(subset(var1_2, Patient.ID == i))
  Total <- sum(work$Number)
  total1[c, "Patient.ID"] <- i
  total1[c, "TotalMut"] <- Total
  c <- c + 1 
}



# Similar to that:
## Can use amino acid change, so G/D







