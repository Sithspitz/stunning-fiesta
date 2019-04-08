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

# Splitting Amino acid changes
## Take last letter
lastvalue <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## Find position of AA
AApos <- function(splist){
  splist <- as.character(splist)
  lslist <- vector(mode = "character",length = length(splist))
  samples <- regexpr("[[:digit:]]",splist)
  spselect <- samples != -1
  lslist[!spselect] <- splist[!spselect]
  sp_pos_fix <- regexpr("[[:digit:]]{1,10}",
                        splist[spselect],
                        perl = T)
  lslist[spselect] <- substr(splist[spselect],
                             sp_pos_fix,sp_pos_fix+attributes(sp_pos_fix)[[1]]-1)
  return(factor(lslist))
}

# Remove NA from columns specified
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Clean Genesets
takegenelevels <- function(df){
  df1 <- df
  geneset <- as.character(df1[, 1])
  return(geneset)
}

# Calculate gene means
CalGeneMean <- function(df, genecol){
  genes <-    df[eval(substitute(genecol), df) == T, ]
  genes1 <- droplevels(subset(genes, DataType == "MicroarrayZScore"))
  GeneScore <- dcast(genes1, Hugo_Symbol ~ Patient.ID, value.var = "Value")
  GeneMean <- colMeans(GeneScore[-1], na.rm = T)
  GeneData <- melt(GeneMean)
  GeneData1 <- cbind(rownames(GeneData), data.frame(GeneData, row.names = NULL))
  colnames(GeneData1)[colnames(GeneData1) == "rownames(GeneData)"] <- "Patient.ID"
  GeneData2 <- GeneData1
  return(GeneData2)
}

## Creating Plots ##
# Add regression formula
lm_eqn <- function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

# Onco
ggOnco <- function(df){
  p <- ggplot(df,aes(x = OncoKRAS, y = Value)) +
    geom_boxplot(
      alpha = 0.5,
      width = 0.2)+ 
    geom_violin(aes(OncoKRAS, fill = OncoKRAS),
                scale = "width",
                alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "MutantKRAS",
         y = "Microarray Value"
    )+ 
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw()+
    theme(
      axis.text = element_text(size = 16))
  sig <- p + stat_compare_means(comparisons = my_comparisons,
                                label = "p.signif")
  return(sig)}

# Plot to make a violin plot across subtype for any Y attribute
ggAny <- function(df, column, YTitle){
  column <- eval(substitute(column),df, parent.frame())
  p <- ggplot(df, aes(x = Subtype, y = column)) +
    geom_boxplot(
      alpha = 0.5,
      width = 0.2)+ 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width",
                alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype",
         y = YTitle
    )+ 
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw()+
    theme(
      axis.text = element_text(size = 16)) + ylim(-2, 3)
  sig <- p + stat_compare_means(comparisons = my_comparisons,
                                label = "p.signif")
  return(sig)}

# Survival
ggsurv <- function(s, CI = "def", plot.cens = T, surv.col = "gg.def",
                   cens.col = "red", lty.est = 1, lty.ci = 2,
                   cens.shape = 3, back.white = F, xlab = "Time",
                   ylab = "Survival", main = "")
{
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) == T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = "def", plot.cens = T, surv.col = "gg.def",
                       cens.col = "red", lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = "Time",
                       ylab = "Survival", main = ""){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == "gg.def", "black", surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl <- if(CI == T | CI == "def") {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ("There are no censored observations")
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = "def", plot.cens = T, surv.col = "gg.def",
                       cens.col = "red", lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = "Time",
                       ylab = "Survival", main = "") {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), "="))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), "="))[1]
    gr.df <- vector("list", strata)
    ind <- vector("list", strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if(surv.col[1] != "gg.def"){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl <- pl + line
    
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop("Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata")
      }else if((length(surv.col) > 1 | surv.col == "gg.def")[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low,lty = group), col = surv.col)}
    } else {pl}
    
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ("There are no censored observations")
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}

# PC Score Plots
ggComponent <- function(df){
  p <- ggplot(df,aes(x = Subtype, y = Rank)) +
    geom_boxplot(
      alpha = 0.5,
      width = 0.2)+ 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width",
                alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype",
         y = "Rank of Component Score"
    )+ 
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw()+
    theme(axis.text = element_text(size = 16)) +
    theme(legend.direction = "horizontal", 
          legend.position = "top")
  sig <- p + stat_compare_means(comparisons = my_comparisons,
                                label = "p.signif")
  return(sig)}

# Plot Averages for bunch of factors
ggFactors <- function(df, YTitle){
  p <- ggplot(df,aes(x = Subtype, y = Rank)) +
    geom_boxplot(
      alpha = 0.5,
      width = 0.2)+ 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width",
                alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype",
         y = YTitle
    ) + 
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw() +
    theme(
      axis.text = element_text(size = 16)) +
    theme(legend.direction = "horizontal", 
          legend.position = "top")
  sig <- p + stat_compare_means(comparisons = my_comparisons,
                                label = "p.signif")
  return(sig)}

# Find the correlation between variable and PC
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}

# Calculate the contribution of variable to the PC
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}

