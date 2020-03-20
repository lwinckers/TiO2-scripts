## ---------------------------
##
## Script name: workflow.R
##
## Purpose of script: Takes one gene list for a process, finds relevant pathways and runs enrichment analysis (using transcriptomics data) to study how process is affected in dataset
## Dataset: 
## Pathway collections:
##
## Author: Laurent Winckers, Martina Kutmon
##
## Date Created: xxxx-xx-xx
##
## Session info:
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17763)
## Packages: readxl_1.3.1, rstudioapi_0.11, RColorBrewer_1.1-2, RCy3_2.6.3
##
## Cytoscape version 3.7.2
##
## ---------------------------

# Installation R packages
# install.packages("BiocManager")
## TODO: add all packages which need to be installed!
# BiocManager::install(c("rstudioapi","RCy3","clusterProfiler"))

## ---------------------------

# Step 1: Pathway selection
# Read process gene list
# Read pathway gene set collections
# Run overrepresentation analysis
# Create gene set with all selected pathways

## ---------------------------

# Step 2: Data import and pre-processing
# Read data file (describe how data needs to look like)
# Data pre-processing if needed
# Define number of experiments/comparisons

## ---------------------------

# Step 3: GSEA 
# Run GSEA per comparison with selected gene set collection
# Select significant results in at least one comparison

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