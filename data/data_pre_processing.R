## ---------------------------
##
## Script name: data_pre_processing.R
##
## Purpose of script: pre-process gene-expression data that will be used in: https://github.com/laurent2207/TiO2-scripts.
## Gene-expression data: GEO:GSE42069
##
## Author: Laurent Winckers, Martina Kutmon
##
## Date Created: 2020-03-23
##
## Session info:
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18362)
## Packages: biomaRt_2.42.0, dplyr_0.8.5
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
library(dplyr)

##### Gene-expression data pre-processing #####

### Load in gene-expression data, previously analysed using ArrayAnalysis (arrayanalysis.org)
# Caco2
caco2_low <- read.table("./data-input/TiO2_24hrs_10_Caco2.txt", sep = "\t", header = T)
caco2_high <- read.table("./data-input/TiO2_24hrs_100_Caco2.txt", sep = "\t", header = T)
# SAE
SAE_low <- read.table("./data-input/TiO2_24hrs_10_SAE.txt", sep = "\t", header = T)
SAE_high <- read.table("./data-input/TiO2_24hrs_100_SAE.txt", sep = "\t", header = T)
# THP1
THP1_low <- read.table("./data-input/TiO2_24hrs_10_THP1.txt", sep = "\t", header = T)
THP1_high <- read.table("./data-input/TiO2_24hrs_100_THP1.txt", sep = "\t", header = T)

### select necessary columns, remove columns that are not needed further down the process
# Caco2
caco2_low <- caco2_low[c(1,2,3,6)]
caco2_high <- caco2_high[c(1,2,3,6)]
# SAE
SAE_low <- SAE_low[c(1,2,3,6)]
SAE_high <- SAE_high[c(1,2,3,6)]
# THP1
THP1_low <- THP1_low[c(1,2,3,6)]
THP1_high <- THP1_high[c(1,2,3,6)]

### change column names adressing respective cell line and either 10 ug/ml, low (L) or 100 ug/ml, high (H) concentration
# Caco2
colnames(caco2_low)[c(2,3,4)] <- c("caco2_L_logFC", "caco2_L_FC", "caco2_L_pval")
colnames(caco2_high)[c(2,3,4)] <- c("caco2_H_logFC", "caco2_H_FC", "caco2_H_pval")
# SAE
colnames(SAE_low)[c(2,3,4)] <- c("SAE_L_logFC", "SAE_L_FC", "SAE_L_pval")
colnames(SAE_high)[c(2,3,4)] <- c("SAE_H_logFC", "SAE_H_FC", "SAE_H_pval")
# THP1
colnames(THP1_low)[c(2,3,4)] <- c("THP1_L_logFC", "THP1_L_FC", "THP1_L_pval")
colnames(THP1_high)[c(2,3,4)] <- c("THP1_H_logFC", "THP1_H_FC", "THP1_H_pval")

### select ensembl IDs from one of the datasets
## 11825
ids <- as.data.frame(caco2_low$ENSG_ID)
colnames(ids) <- "ENSG_ID"

### map Ensemble IDs to hgnc symbols
## 11791
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

genes <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'), 
  filters = 'ensembl_gene_id',
  values = ids$ENSG_ID,
  mart = ensembl
)

### remove NAs from gene file
## 11324
genes2 <- genes[!(is.na(genes$entrezgene_id)),]

### remove rows without HGNC-symbol
## 11298
genes3 <- genes2[-which(genes2$hgnc_symbol == ""),]

### merge datasets with annotated gene identifiers
colnames(genes3)[1] <- "ENSG_ID"
bound <- merge(genes3, caco2_low, by = "ENSG_ID")
bound <- merge(bound, caco2_high, by = "ENSG_ID")
bound <- merge(bound, SAE_low, by = "ENSG_ID")
bound <- merge(bound, SAE_high, by = "ENSG_ID")
bound <- merge(bound, THP1_low, by = "ENSG_ID")
bound <- merge(bound, THP1_high, by = "ENSG_ID")

### create seperate files for each cell line
#caco2 <- bound[,c(1:9)]
#SAE <- bound[,c(1:3,10:15)]
#THP1 <- bound[,c(1:3,16:21)]

### save annotated gene-expression files for each cell line and combined
write.table(bound, "./data-output/merged_TiO2.txt", sep = "\t", quote = F, row.names = F)
#write.table(caco2, "./data-output/merged_caco2.txt", sep = "\t", quote = F, row.names = F)
#write.table(SAE, "./data-output/merged_SAE.txt", sep = "\t", quote = F, row.names = F)
#write.table(THP1, ".//data-output/merged_THP1.txt", sep = "\t", quote = F, row.names = F)

### calculate score  for ranking of genes in GSEA analysis
# use specific rank-score calculation 
data_rnk <- bound
data_rnk <- data_rnk %>% mutate(caco2_L_rnk = caco2_L_FC * -log10(caco2_L_pval))
data_rnk <- data_rnk %>% mutate(caco2_H_rnk = caco2_H_FC * -log10(caco2_H_pval))
data_rnk <- data_rnk %>% mutate(SAE_L_rnk = SAE_L_FC * -log10(SAE_L_pval))
data_rnk <- data_rnk %>% mutate(SAE_H_rnk = SAE_H_FC * -log10(SAE_H_pval))
data_rnk <- data_rnk %>% mutate(THP1_L_rnk = THP1_L_FC * -log10(THP1_L_pval))
data_rnk <- data_rnk %>% mutate(THP1_H_rnk = THP1_H_FC * -log10(THP1_H_pval))

### clean ranked dataset and select columns that contain rank score and Entrezgene ID column
data_rnk <- data_rnk[c(1,2,3,(grep("_rnk", names(data_rnk))))]

### save file with ranked scores
write.table(data_rnk, "./data-output/rankscore_TiO2.txt", sep = "\t", quote = F, row.names = F)
