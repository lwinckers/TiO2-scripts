## ---------------------------
##
## Script name: data_pre_processing.R
##
## Purpose of script: pre-process pathway database-data that will be used in: https://github.com/laurent2207/TiO2-scripts.
## Pathway databases: WikiPathways-20190910-gmt-Homo_sapiens, KEGG-v7.0-entrezgene, Reactome-v7.0-entrez
##
## Author: Laurent Winckers, Martina Kutmon
##
## Date Created: 2020-03-23
##
## Session info:
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18362)
## Packages: qusage_2.20.0, data.table_1.12.8
##
## ---------------------------

##### Set up environment #####

### clear workspace and set string as factors to false
rm(list=ls())
options(stringsAsFactors = F)
gc()

### instal libraries
library(qusage)
library(data.table)

### set working directroy
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

### load gmt files and transpose the obtained lists to data frames
path = paste0(getwd(), "/data-input")
gmtFile <- list.files(path = path, pattern = ".gmt")
for (i in 1:length(gmtFile)){
  assign(gmtFile[i], read.gmt(paste0(path,"/",gmtFile[i])))
}

gmtFile <- mget(ls(pattern = ".gmt"))
gmtFile <- lapply(gmtFile, function(x){plyr::ldply(x,data.frame)})

### combine pathway database files and clean it
databases <- do.call(rbind, gmtFile)
colnames(databases) <- c("pathway", "entrezgene")
databases$entrezgene <- as.numeric(databases$entrezgene)

### save file of combined pathway databases
write.table(databases, "./data-output/pw_databases.txt", quote = F, row.names = F, sep = "\t")
