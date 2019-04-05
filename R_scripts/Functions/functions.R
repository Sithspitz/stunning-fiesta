### Functions ###

## Packages ###
# Install and load required packages
required <- c("tidyverse",
              "reshape2",
              "ggpubr",
              "survival",
              "data.table",
              "gplots",
              "devtools",
              "randomForest",
              "factoextra",
              "ggrepel",
              "cluster",
              "magrittr",
              "e1071",
              "ROCR",
              "cowplot",
              "car",
              "gridExtra")
for (lib in required)
{
  if (!require(lib, character.only = T))
  {
    install.packages(lib, dependencies = T)
    suppressMessages(library(lib, character.only = T, quietly = T))
  }
}

## General Functions ##
# Write CSV for data 
writeCsvD <- function(df){
  fn <- deparse(substitute(df))
  write.csv(df, file = paste0("L:/Richard B/TCGA_data/Pancancer/raw_csv",fn,".csv"), row.names = F)}

# Remove duplicate rows
deduplicate <- function(df){
  df1 <- df[!duplicated(df),]
  return(df1)
}

# Trim whitespace
trimWS <- function (x) gsub("^\\s+|\\s+$", "", x)

# Function to remove the addition of X. column at the start of a dataframe
Correct_Colnames <- function(df) {
  delete.columns <- grep("(^X$)|(^X\\.)(\\d+)($)", colnames(df), perl = T)
  if (length(delete.columns) > 0) {
    row.names(df) <- as.character(df[, grep("^X$", colnames(df))])
    df <- df[,-delete.columns]
    colnames(df) <- gsub("^X", "",  colnames(df))
  }
  return(df)
}

# Factor a list of given columns
factorthese <- function(df, somecolumns){
  Fctr <- names(df) %in% somecolumns
  df[,Fctr] <- lapply(df[,Fctr], function(column) as.factor(as.character(column)))
  return(df)
}

# Don't factor a list of given columns
nofactorthese <- function(df, somecolumns){
  Fctr <- !names(df) %in% somecolumns
  df[,Fctr] <- lapply(df[,Fctr], function(column) as.factor(as.character(column)))
  return(df)
}

# Don't character a list of given columns
nocharacterthese <- function(df, somecolumns){
  chrc <- !names(df) %in% somecolumns
  df[,chrc] <- lapply(df[,chrc], function(column) as.character(column))
  return(df)
}

# Numeric these
numericthese <- function(df, somecolumns){
  nmber <- names(df) %in% somecolumns
  df[,nmber] <- lapply(df[,nmber], function(column) as.numeric(as.character(column)))
  return(df)
}

## TCGA Specific Functions ##

# Create patient ID from sample ID
samptopat <- function(splist){
  splist <- as.character(splist)
  lslist <- vector(mode = "character",length = length(splist))
  samples <- regexpr("^TCGA",splist)
  spselect <- samples != -1
  lslist[!spselect] <- splist[!spselect]
  sp_pos_fix <- regexpr("TCGA[[:punct:]]{1}[[:alnum:]]{2}[[:punct:]]{1}[[:alnum:]]{4}",
                        splist[spselect],
                        perl = T)
  lslist[spselect] <- substr(splist[spselect],
                             sp_pos_fix,sp_pos_fix+attributes(sp_pos_fix)[[1]]-1)
  return(factor(lslist))
}

###### NEXT FUNCITON TO EXAMINE ON 08_04_19 IS SPLITTING AMINO ACID CHANGES #######



