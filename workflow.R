## ---------------------------
##
## Script name: workflow.R
##
## Purpose of script: Takes one gene list for a process, finds relevant pathways and runs enrichment analysis (using transcriptomics data) to study how process is affected in dataset
##
## Author: Laurent Winckers, Martina Kutmon
##
## Date Created: 2020-03-24
##
## Session info:
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17763)
## Packages: clusterProfiler_3.14.3, plyr_1.8.6, biomaRt_2.42.0, dplyr_0.8.5, data.table_1.12.8, pheatmap_1.0.12,
## RColorBrewer_1.1-2, colorRamps_2.3 
##
## ---------------------------

# Step 1: Data import and pre-processing
# See data folder which holds three data pre-processing scripts
# For additional information see the README.md file

## ---------------------------

### set up environment
rm(list=ls())
options(stringsAsFactors = F)
gc()

### set working directroy
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

### load libraries
## step 2
library(clusterProfiler)
library(plyr)
library(biomaRt)

## step 3
library(clusterProfiler)
library(dplyr)

## step 4
library(data.table)
library(pheatmap)
library(colorRamps)
library(RColorBrewer)

### provide name of GO-term which you want to use
fileName <- "GO0006974"

## ---------------------------

# Step 2: GSEA 
# Run GSEA per comparison with selected gene set collection
# Select significant results in at least one comparison

### link to GSEA analysis function
source("./functions/GSEA.R")

### load gene-expression file with specific ranked score
data <- read.table("./data/data-output/rankscore_TiO2.txt", header = T, sep ="\t")

### load geneset
geneset <- read.table("./data/data-output/wp_database.txt", header = T, sep = "\t")

### perform GSEA analysis
GSEAanalysis(GENESET = geneset, fileName = paste0(fileName), data = data)

## ---------------------------

# Step 3: Pathway selection
# Read process gene list
# Read pathway gene set collections
# Run overrepresentation analysis
# Create gene set with all selected pathways

### load in gene GO-term genelist
goterm <- read.table(paste0("./data/data-output/ann_", fileName ,".txt"), header = T, sep ="\t")

### load in pathway databases file
databases <- read.table("./data/data-output/wp_database.txt", header = T, sep ="\t")


### perform enricher analysis
res_enr <- as.data.frame(enricher(gene = goterm$entrezgene_id, TERM2GENE = databases,
                              minGSSize = 10, maxGSSize = 500, 
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.05))

### save result
write.table(res_enr, paste0("./output/wpres_enricher_", fileName, ".txt"), quote = F, sep = "\t", row.names = F)

## ---------------------------

### select enriched pathways from combined pathway databases file
res_pw <- as.data.frame(databases[databases$pathway %in% res_enr$ID,])

### select only unique rows
res_pw <- unique(res_pw)

### remove NAs and empty values as they are pseudogenes, microRNAs or discontinued genes
res_pw <- res_pw[!is.na(res_pw$entrezgene),]
res_pw <- res_pw[-which(res_pw$entrezgene == ""),]

### save result
write.table(res_pw, paste0("./output/wppws_", fileName, ".txt"), quote = F, sep = "\t", row.names = F)

## ---------------------------

# Step 4: Create heatmap
# Columns = all comparisons
# Rows = selected pathways in gene set collection (only pathways which are significant in at least one comparison)
# NES is used for cell color
# Significance (maybe possible) as stars
# Cluster rows

### load in GSEA result files
files <- list.files(path = paste0(getwd(), "/output/GSEA"), pattern = paste0("wptest-", fileName))
for (i in 1:length(files)){
  assign(paste0("file", i), read.table(paste0(getwd(), "/output/GSEA/", files[i]), header = T, sep = "\t"))
}

## need to find nicer way for this
colnames(file1)[c(2:11)] <- paste(colnames(file1)[c(2:11)], "1", sep = "_")
colnames(file2)[c(2:11)] <- paste(colnames(file2)[c(2:11)], "2", sep = "_")
colnames(file3)[c(2:11)] <- paste(colnames(file3)[c(2:11)], "3", sep = "_")
colnames(file4)[c(2:11)] <- paste(colnames(file4)[c(2:11)], "4", sep = "_")
colnames(file5)[c(2:11)] <- paste(colnames(file5)[c(2:11)], "5", sep = "_")
colnames(file6)[c(2:11)] <- paste(colnames(file6)[c(2:11)], "6", sep = "_")

res_merge <- merge(file1, file2, by = "ID")
res_merge <- merge(res_merge, file3, by = "ID")
res_merge <- merge(res_merge, file4, by = "ID")
res_merge <- merge(res_merge, file5, by = "ID")
res_merge <- merge(res_merge, file6, by = "ID")

res_sig <- subset(res_merge, pvalue_1 < 0.01 | pvalue_2 < 0.01 | pvalue_3 < 0.01 | pvalue_4 < 0.01 | pvalue_5 < 0.01 | pvalue_6 < 0.01)
res_sig <- res_sig[c(1,(grep("enrichmentScore_", names(res_sig))), (grep("pvalue_", names(res_sig))))]

write.table(res_sig, paste0("./output/wpsigGSEA_", fileName, ".txt"), quote = F, sep = "\t", row.names = F)

### link to GSEA analysis function
source("./functions/heatmap.R")

res_sig < res_sig[c(1:6)]

heatmap(fileName = fileName, data = res_sig)

### session information
sessionInfo()

## ---------------------------

# Step 5: Network visualization
# Manual selection before
# Pathway - Gene view
# Pathway - Pathway view