# Inflammation network for titanium dioxide gene expression data analysis
This repository contains all the scripts used for the inflammation network analysis for titanium dioxide gene expression data.
The scripts are split up in various big tasks.

### Order of scripts
* [Common genes](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/common_genes.R) - All unique genes from both GO terms
* [Pathway selection](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/pathway_selection.r) - Selection of pathways based on the [unique genes](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/common_genes.R) from both GO terms
* [Network file creation](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/network_file_creation.r) - Creation of files used in for creation of networks in Cytoscape
* [GO terms network creation](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/GO-network-creation.R) - Creation of network based on the two GO terms
* [RCy3 - Cytoscape network creation](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/RCY3_network_creation.r) - Automated creation of network based on the selected pathways from [Pathway selection](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/pathway_selection.r) using the RCy3 package 
* [GSEA analysis in R](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/GSEA_analysis.R) - GSEA analysis in R based on [Cytoscape network](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/RCY3_network_creation.r) and [titanium dioxide gene expression data](https://github.com/laurent2207/TiO2_inflammation_network/tree/master/gene-expression-data)
* [Heatmap visualisation of GSEA analysis](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/Heatmap.R) - Visualisation of GSEA analysis in a heatmap

### Gene expression data
* Gene expression data was obtained from [GEO:GSE42069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42069) 
* Gene expression data was pre-processed and analysed with [ArrayAnalysis](http://www.arrayanalysis.org/)
* Cleaning and merging of gene expression data can be found in the [gene-expression-data folder](https://github.com/laurent2207/TiO2_inflammation_network/tree/master/gene-expression-data)

### Functions used
* [Install packages](https://github.com/laurent2207/TiO2_inflammation_network/blob/master/functions/autoInstallPackages.R) - Install packages which were not installed yet. For more information [click here!](https://github.com/laurent2207/autoInstallPackages)

### Authors
* Laurent Winckers - [GitHub](https://github.com/laurent2207)

### License
*No license yet*

### Acknowledgments
* Dr. Martina Summer-Kutmon - [GitHub](https://github.com/mkutmon)
