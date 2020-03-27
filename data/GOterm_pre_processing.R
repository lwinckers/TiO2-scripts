## ---------------------------
##
## Script name: data_pre_processing.R
##
## Purpose of script: pre-process GO-term genelists that will be used in: https://github.com/laurent2207/TiO2-scripts.
## GO-term genelists: Apoptopic process (GO:0006915), Inflammatory response (GO:0006954)
## GO-term genelists: Cellular response to DNA damage stimulus (GO:0006974), Cellular response to oxidative stress (GO:0034599)
##
## Author: Laurent Winckers, Martina Kutmon
##
## Date Created: 2020-03-23
##
## Session info:
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18362)
## Packages: biomaRt_2.42.0
##
## ---------------------------

##### Set up environment #####

### clear workspace and set string as factors to false
rm(list=ls())
options(stringsAsFactors = F)
gc()

### set working directroy
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

### install packages 
library(biomaRt)

### link function
source("./functions/GO_annotation.R")

### load in GO-term genelists, only GO.TERM and SYMBOL column 
path = paste0(getwd(), "/data-input")

ls_GO <- list.files(path = path, pattern = "GO00")

for (i in 1:length(ls_GO)){
  assign(ls_GO[i], read.table(paste0(path,"/", ls_GO[i]), header = T, sep = "\t")[c("GO.TERM", "SYMBOL")])
}

### combine them in one list to loop over them
ls_GO <- mget(ls(pattern = "GO00"))

### remove .txt at end of dataframe variable names
ls_GO <- setNames(ls_GO, gsub(pattern = ".txt", replacement = "", x = names(ls_GO)))

### map HGNC-symbols to Entrezgene IDs using the GO_annotation function
for (i in 1:length(ls_GO)){
  GO_annotation(fileName = paste0("ann_", names(ls_GO[i])), 
                data = ls_GO[[i]], 
                values = ls_GO[[i]]$SYMBOL,
                attributes = c('hgnc_symbol', 'entrezgene_id'),
                filter = "hgnc_symbol",
                path = paste0(getwd(), "/data-output/"))
}
