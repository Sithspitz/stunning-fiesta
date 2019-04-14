### Gene Sets ###
source("./R_scripts/Functions/functions.R")



# Import Pre-made gene sets #
predir <- c("./R_scripts/Gene_Sets/")
setwd("./R_scripts/Gene_Sets/")
temp <- list.files(predir, pattern="*.csv")

