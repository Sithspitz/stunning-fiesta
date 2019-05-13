### Gene Sets ###
source("./R_scripts/Functions/functions.R")



# Import Pre-made gene sets #

# All .csv files in the directory
predir <- c("./R_scripts/Gene_Sets/csv_format/")
setwd(predir)
temp <- list.files(pattern="*.csv")
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
         read.csv), envir = .GlobalEnv)

# Clean the immune genesets as characters
aDC_genes <- takegenelevels(aDC_geneset)
BCells_genes <- takegenelevels(BCells_geneset)
CD8_genes <- takegenelevels(CD8_geneset)
Cytotoxic_genes <- takegenelevels(Cytotoxic_geneset)
DC_genes <- takegenelevels(DC_geneset)
Eosinophils_genes <- takegenelevels(Eosinophils_geneset)
iDC_genes <- takegenelevels(iDC_geneset)
Macrophages_genes <- takegenelevels(Macrophages_geneset)
MastCells_genes <- takegenelevels(MastCells_geneset)
Neutrophils_genes <- takegenelevels(Neutrophils_geneset)
NK_genes <- takegenelevels(NK_geneset)
NKCD56bright_genes <- takegenelevels(NKCD56bright_geneset)
NKCD56dim_genes <- takegenelevels(NKCD56dim_geneset)
Tcm_genes <- takegenelevels(Tcm_geneset)
Tem_genes <- takegenelevels(Tem_geneset)
Tfh_genes <- takegenelevels(Tfh_geneset)
Tgd_genes <- takegenelevels(Tgd_geneset)
Th1_genes <- takegenelevels(Th1_geneset)
Th17_genes <- takegenelevels(Th17_geneset)
Th2_genes <- takegenelevels(Th2_geneset)
Tregs_genes <- takegenelevels(Tregs_geneset)

