### Lung Adenocarcinoma PanCancer Atlas Mutation Subsetting ###
source("./R_scripts/Functions/functions.R")
setwd("~/DataShare/TCGA_Mut_Analysis_Temporary/Mutation_TCGA_Pancancer_Analysis/Input")
mut <- read.csv("luad_tcga_pancancer_mut.csv", header = T)

# Seperate - subset anything that isn't a "Silent mutation"
Intermediate <- droplevels(subset(mut, Variant_Classification != "Silent"))



# Master Subsetting Loop Tidy - MT vs Other #

## Groups into which MT they have, or double MT vs others without MTs ##
### No silent MT included ###
#### Being used for STK11 and KRAS here, but can be changed ####
setwd("~/DataShare/TCGA_Mut_Analysis_Temporary/Mutation_TCGA_Pancancer_Analysis/Subset")
output <- data.frame(stringsAsFactors = F)

c <- 1
for(i in levels(Intermediate$Patient.ID)){
  print(paste0("working on: ", i))
  work <- droplevels(subset(Intermediate, Patient.ID == i))
  work$Hugo_Symbol <- as.factor(work$Hugo_Symbol)
  output[c, "Patient.ID"] <- i
  output[c, "Mutation_Status"] <- ifelse((length(is.na(levels(work$Hugo_Symbol)[levels(work$Hugo_Symbol) %in% "STK11"])) != 0 &
                                            length(is.na(levels(work$Hugo_Symbol)[levels(work$Hugo_Symbol) %in% "KRAS"])) != 0), "Double_Mut",
                                         ifelse((length(is.na(levels(work$Hugo_Symbol)[levels(work$Hugo_Symbol) %in% "STK11"])) == 0 &
                                                   length(is.na(levels(work$Hugo_Symbol)[levels(work$Hugo_Symbol) %in% "KRAS"])) != 0), "KRAS_Mut",
                                                ifelse((length(is.na(levels(work$Hugo_Symbol)[levels(work$Hugo_Symbol) %in% "STK11"])) != 0 &
                                                          length(is.na(levels(work$Hugo_Symbol)[levels(work$Hugo_Symbol) %in% "KRAS"])) == 0), "STK11_Mut", "Other_Mut")))
  c <- c + 1}

write.csv(output, "KRAS_STK11_plus_dbl_tidy_vs_other.csv", row.names = F)



# Subset KRAS, STK11 and KRAS/STK11 untidy #
MutChanges1 <- subset(Intermediate, Hugo_Symbol == "KRAS" | Hugo_Symbol == "STK11")
setwd("~/DataShare/TCGA_Mut_Analysis_Temporary/Mutation_TCGA_Pancancer_Analysis/Subset")
write.csv(MutChanges1, "KRAS_STK11_plus_dbl_untidy.csv" , row.names = F)

# Subset KRAS, STK11 and KRAS/STK11 tidy #

# tidy1 drops levels, and tidy2 removes 'true duplicates' from these levels 
## i.e. X2 TCGA.same with STK11 removed, whilst x2 TCGA.same with KRAS and STK11 not
tidy1 <- MutChanges1[, c("Patient.ID", "Hugo_Symbol")] %>% droplevels()
tidy2 <- tidy1[!duplicated(tidy1), ]

# The loop creates a df with whether STK11 single, KRAS single or double_mut, tidy
loop_output <- data.frame("Patient.ID" = character(),
                     "Mutation_status" = character(),
                     stringsAsFactors = F)
c <- 1
for(i in levels(tidy2$Patient.ID)){
  print(paste0("working on: ", i))
  work <- droplevels(subset(tidy2, Patient.ID == i))
  mut_ <- ifelse((nlevels(work$Hugo_Symbol) >= 2), "double_mut", as.character(levels(work$Hugo_Symbol)))
  loop_output[c, "Patient.ID"] <- i
  loop_output[c, "Mutation_status"] <- mut_
  c <- c + 1
}

write.csv(loop_output, "KRAS_STK11_plus_dbl_tidy.csv", row.names = F)



# Basic subset any other mutations am interested in #

# ATM below (for example)
ATM_sub1 <- subset(Intermediate, Hugo_Symbol == "ATM")
ATM_sub2 <- ATM_sub1[, c("Patient.ID", "Hugo_Symbol")] %>% droplevels()
ATM_sub3 <- ATM_sub2[!duplicated(ATM_sub2), ]
write.csv(ATM_sub3, "ATM_mut_tidy.csv", row.names = F)

