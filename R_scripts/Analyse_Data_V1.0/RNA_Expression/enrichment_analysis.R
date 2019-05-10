### Lung Adenocarcinoma PanCancer Atlas RNA Enrichment ###
## Should be adopted dependent on the analysis in question ##
mypackages <- c("GSEABase", "GSVA", "Biobase", "genefilter",
                "limma", "RColorBrewer", "GSVAdata")
# INSTALL reshape PACKAGE?!?
lapply(mypackages, library, character.only = T)
source("./R_scripts/Functions/functions.R")

setwd("~/DataShare/TCGA_RNA_Analysis/Input/")
rna <- read.csv("luad_tcga_pancancer_rna.csv", header = T)
mut <- read.csv("Mut_tidy_vs_other.csv", header = T)


# Merge and clean up Entrez Ids
m1 <- merge(rna, mut, by = "Patient.ID")
Intermediate2 <- m1
Intermediate2$Entrez_Gene_Id <- NULL 

# Droplevels RNASeqZScore and remove the datatype column
Intermediate3 <- droplevels(subset(Intermediate2, DataType == "RNASeqEL"))
Intermediate4 <- Intermediate3
Intermediate4$DataType <- NULL
Intermediate5 <- Intermediate4

# Subset just the Mut and the WT and drop the Mutation_Status column
my_Mut_subset <- droplevels(subset(Intermediate5, Mutation_Status == "Mut"))
my_WT_subset <- droplevels(subset(Intermediate5, Mutation_Status == "Other_Mut"))
my_Mut_subset2 <- my_Mut_subset
my_WT_subset2 <- my_WT_subset

my_Mut_subset2$Mutation_Status <- NULL
my_WT_subset2$Mutation_Status <- NULL


# Tidyr idea
my_Mut_subset3 <- my_Mut_subset2
dofactor <- c("Patient.ID", "Hugo_Symbol")
my_factored <- factorthese(my_Mut_subset3, dofactor)

my_extra <- my_factored
my_extra$extra <- as.factor(paste(my_extra$Patient.ID, my_extra$Hugo_Symbol, sep = ","))

dups <- my_extra
dropped <- dups[!duplicated(dups$extra), ]
dropped$extra <- NULL

t1 <- recast(dropped, Patient.ID ~ Hugo_Symbol, id.var = 1:2)
write.csv(t1, "t1_please.csv")








no_pid <- my_factored
no_pid$Patient.ID <- NULL

no_pid_spread <- spread(no_pid, Hugo_Symbol, Value)



my_factored %>% group_by(Hugo_Symbol) %>% mutate(grouped_id = row_number())
my_factored %>% spread(Patient.ID, Hugo_Symbol) %>% select(-grouped_id)

my_spread <- spread(my_factored, Hugo_Symbol, Value)

melted <- melt(my_Mut_subset3, id.vars = c("Patient.ID"))
t2 <- dcast(my_Mut_subset3, Patient.ID ~ Value, value.var = "Value")

write.csv(t1, "t1_please.csv")

t2 <- dcast(my_Mut_subset3, Patient.ID ~ Hugo_Symbol, value.var = "Value")

#### TRY AND GET THE DCAST TO WORK!!!!!! 

# t1 <- melt(my_Mut_subset3, id.vars =c("Patient.ID", "Hugo_Symbol"), 
  #         measure.vars = "Value")





#my_Mut_subset3 %>% spread(key = Hugo_Symbol, value = Value)

#my_Mut_subset4 <- rowid_to_column(my_Mut_subset3)

# TEST RESHAPE WITH MUT SUBSET
#md <- melt(my_Mut_subset2, id=(c("Patient.ID", "Hugo_Symbol")))



### I have managed to process the data to how I want it 
### But need to manipulate it futher to ensure that the input is correct for GSVA
### Is there anyway of doing this without losing information?
### Perhaps look at creating an ExpressionSet object to keep metadata?
### Need to transpose the data into the correct data.matrix in paper




