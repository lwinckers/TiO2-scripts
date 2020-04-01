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
library(plyr)

### set working directroy
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

### load gmt files and transpose the obtained lists to data frames
gmtFile <- qusage::read.gmt("./data-input/gmt_wp_Homo_sapiens.gmt")
gmtFile <- plyr::ldply(gmtFile, data.frame)

### Clean data frame
colnames(gmtFile) <- c("pathway", "entrezgene")
gmtFile$entrezgene <- as.numeric(gmtFile$entrezgene)

### save file of combined pathway databases
write.table(gmtFile, "./data-output/wp_database.txt", quote = F, row.names = F, sep = "\t")

