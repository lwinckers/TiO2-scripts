## ---------------------------
##
## Script name: workflow.R
##
## Purpose of script: Takes one gene list for a process, finds relevant pathways and runs enrichment analysis (using transcriptomics data) to study how process is affected in dataset
##
## Author: Laurent Winckers, Martina Kutmon
##
## Date Created: xxxx-xx-xx
##
## Session info:
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17763)
## Packages: 
##
## ---------------------------

# Installation R packages
# install.packages("BiocManager")
## TODO: add all packages which need to be installed!
# BiocManager::install(c("rstudioapi","RCy3","clusterProfiler"))

## ---------------------------

# Step 1: Data import and pre-processing
# Read data file (describe how data needs to look like)
# Data pre-processing if needed
# Define number of experiments/comparisons


## ---------------------------

# Step 2: Pathway selection
# Read process gene list
# Read pathway gene set collections
# Run overrepresentation analysis
# Create gene set with all selected pathways

### set up environment
rm(list=ls())
options(stringsAsFactors = F)

### load libraries
library(clusterProfiler)
library(qusage)
library(plyr)

### load in gene GO-term genelist
fileName <- "ann_GO0006915.txt"

goterm <- read.table(paste0("./data-output/", fileName), header = T, sep ="\t")

### load in pathway databases file
databases <- read.table("./data-output/pw_databases.txt", header = T, sep ="\t")

### perform enricher analysis
res <- as.data.frame(enricher(gene = goterm, minGSSize = 10, TERM2GENE = databases, pvalueCutoff = 1))

### save result
resName <- "GO0006915"

write.table(res_pw , paste0("./data-output/res_enricher_", resName), quote = F, sep = "\t", row.names = F)

## ---------------------------

### select enriched pathways from combined pathway databases file
res_pw <- as.data.frame(databases[databases$pathway %in% res$ID,])

### select only unique rows
res_pw <- unique(res_pw)

### map entrezgene IDs to hgnc symbols
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

genes <- getBM(
  attributes = c('hgnc_symbol', 'entrezgene_id'), 
  filters = 'entrezgene_id',
  values = edge_table$entrezgene,
  mart = ensembl
)

### merge edge table with hgnc symbols
res_pw$hgnc_symbol <- genes$hgnc_symbol[match(res_pw$entrezgene, genes$entrezgene)]

### remove NAs and empty values as they are pseudogenes, microRNAs or discontinued genes
res_pw <- res_pw[!is.na(res_pw$hgnc_symbol),]
res_pw <- res_pw[-(res_pw$hgnc_symbol == ""),]

### save result
resName <- "GO0006915"

write.table(res_pw , paste0("./data-output/enriched_", resName), quote = F, sep = "\t", row.names = F)

## ---------------------------

# Step 3: GSEA 
# Run GSEA per comparison with selected gene set collection
# Select significant results in at least one comparison

### set up environment
rm(list=ls())
options(stringsAsFactors = F)

### load libraries
library(clusterProfiler)
library(dplyr)

### link to GSEA analysis function
source("./functions/GSEA.R")

### load gene-expression file with specific ranked score
data <- read.table("./data-output/ranked_TiO2.txt", header = T, sep ="\t")

### load geneset
fileName <- "enriched_"

geneset <- read.table(paste0("./data-output/", fileName), header = T, sep = "t")

### perform GSEA analysis
GSEAan(GENESET = geneset, fileName = "")

## ---------------------------

# Step 4: Create heatmap
# Columns = all comparisons
# Rows = selected pathways in gene set collection (only pathways which are significant in at least one comparison)
# NES is used for cell color
# Significance (maybe possible) as stars
# Cluster rows

## ---------------------------

# Step 5: Network visualization
# Manual selection before
# Pathway - Gene view
# Pathway - Pathway view