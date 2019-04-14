### Gene Sets ###
source("./R_scripts/Functions/functions.R")



# Import Pre-made gene sets #

# Tested with "Th17_geneset"
predir <- c("./R_scripts/Gene_Sets/")
setwd(predir)
temp <- list.files(pattern="*.csv")
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
         read.csv), envir = .GlobalEnv)

# Clean the geneset as characters
Th17_genes <- takegenelevels(Th17_geneset)

