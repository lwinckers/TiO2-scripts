# Workflow

## How to run
This workflow can entirely be run in the [workflow.R script](https://github.com/laurent2207/TiO2-scripts/blob/master/workflow.R). The necessary data files are already pre-processed and can be found in the [data-output folder](https://github.com/laurent2207/TiO2-scripts/tree/master/data/data-output) in the [data folder](https://github.com/laurent2207/TiO2-scripts/tree/master/data). The [workflow.R script](https://github.com/laurent2207/TiO2-scripts/blob/master/workflow.R) loads in these pre-processed files and the script can be used from the onset.  

## Structure

### Main repository
* #### workflow.R script
Main workflow script. For details see below.

* #### Data folder
Folder which contains the raw data and three pre-processing scripts. Details for these scripts can be found below. 

* #### Output folder
Folder which is dedicated for the output of the [workflow.R script](https://github.com/laurent2207/TiO2-scripts/blob/master/workflow.R).

* #### Functions folder
Folder which contains the fucntions that will be used in the [workflow.R script](https://github.com/laurent2207/TiO2-scripts/blob/master/workflow.R).

## Scripts

### Pre-processing steps
* #### GO-term gene lists
Script that loads in GO-term genelists, and annotates the HGNC-symbols to entrezgene identifiers.

* #### Pathway GMT files
Script that loads in GMT files from WikiPathways, KEGG and Reactome and combines these three into one data-frame. This data-frame contains two columns; pathway and gene.

* #### Data pre-processing
Script that loads in mircoarray gene-expression data of six conditions, of which the raw files were pre-processed using [ArrayAnalysis](arrayanalysis.org), are combined into one data-frame. This data-frame contains fold change, log fold change and p-value for each conditions. Enseble identifiers, entrezgene identifiers and HGNC-symbols are annotated for all measured genes. 
Script also calculates scores which will be used to rank the genes for GSEA analysis. Scores are calculated according to the formula; signed Fold change * -log10(p-value).

### Workflow script
* #### Pathway selection
Selects pathways based on GO-term genelists. Enricher function from the [clusterProfiler package](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) is used to get significant enriched pathways from the three pathwaydatabses. 

* #### GSEA analysis
Performs GSEA analysis based on the significant enriched pathways, which are used as genesets, and the ranked genes. Genes are ranked accoring to the rankingscore which was calculated in the [data pre-processing script](https://github.com/laurent2207/TiO2-scripts/blob/master/data/data_pre_processing.R).

* #### Heatmap creation
Creates heatmap for each GO-term. Heatmap depicts the normalised enriched-score obtained fomr the GSEA analysis for pathways which are significant for at least one of the conditions. 
