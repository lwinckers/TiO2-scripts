# Investigating the cell-specific toxicity response to titanium dioxide nanobelts
## How to run
The [ScriptTiO2.Rmd file](https://github.com/laurent2207/TiO2-scripts/blob/master/ScriptTiO2.Rmd) is seperated in three parts *i.e.* Part 1, Part 2 and Part 3 (described below). The necessary data files can be found in the [data folder](https://github.com/laurent2207/TiO2-scripts/tree/master/data). The orginal data, retrieved from [GEO, accession:GSE42069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42069) and analysed with [ArrayAnalysis](https://arrayanalysis.org), can be found in the [original_data subfolder](https://github.com/laurent2207/TiO2-scripts/tree/master/data/original_data) of the aforementioned [data folder](https://github.com/laurent2207/TiO2-scripts/tree/master/data). The files used in the Data pre-processing part (described below) can be found in the [data_pre_processing subfolder](https://github.com/laurent2207/TiO2-scripts/tree/master/data/data_pre_processing) of the aforementioned [data folder](https://github.com/laurent2207/TiO2-scripts/tree/master/data). The resulting files can be found in the [output folder](https://github.com/laurent2207/TiO2-scripts/tree/master/output). 
The functions that are used in the [ScriptTiO2.Rmd file](https://github.com/laurent2207/TiO2-scripts/blob/master/ScriptTiO2.Rmd) are found in the [functions folder](https://github.com/laurent2207/TiO2-scripts/tree/master/functions).

## Part descriptions

* #### Data pre-processing
Pre-processing of data files, which can be found in the [original_data subfolder](https://github.com/laurent2207/TiO2-scripts/tree/master/data/original_data). This part gene identifiers so that the files include Ensemble identifiers, HGNC symbols and Entrez Gene identifiers. Moreover, it merges the separate files into one file named TiO2-dataset.txt and saves this file in the [output folder](https://github.com/laurent2207/TiO2-scripts/tree/master/output).
Furthermore it creates depict signififcantly differentially expressed genes in volcano plots for all conditions, using the EnhanvedVolcano package in R (version 3.6.1) [[Blighe, Rana, and Lewis (2018)](https://github.com/kevinblighe/EnhancedVolcano)].

#### For the next three parts use the [ScriptTiO2.Rmd file](https://github.com/laurent2207/TiO2-scripts/blob/master/ScriptTiO2.Rmd)

* #### Part 1 - Find affected pathways
The first part is used to find affected pathways. This part uses overrepresentation analysis (ORA) based on gene expression data. This method is used to identify significantly enriched pathways from the WikiPathways database. 
In this part the enricher function of the clusterProfiler package in R [PMID:[22455463](https://pubmed.ncbi.nlm.nih.gov/22455463/)]. Almost all variables were kept standard except the Adjusted p-value cut-off and q-value cut-off were set at lower than 0.05 and minimal gene set and maximal gene set size were set to 10 and 300 respectively.
As mentioned above, this part makes use of the pathway models from the WikiPathways database (www.wikipathways.org). It is an open platform for the curation of biological pathways and hosts a pathway database with custom graphical pathway editing tools [PMID:[29136241](https://pubmed.ncbi.nlm.nih.gov/29136241/)]. 
For both this part and the second part, pathway-gene information was obtained in the Gene Matrix Transposed (GMT) file format from http://data.wikipathways.org/. The pathway dataset included the Curated Collection and the Reactome pathways.

* #### Part 2 - Find toxicity related pathways
The second part is used to find toxicity related pathways based on their genetic overlap with genes from the Gene Ontology (GO) terms. For this part the GO-terms were selected based on their relation to toxicologic processes i.e. responses in the human body. Genes were obtained from the Gene Ontology (GO) (AmiGO 2 version 2.5.12) categories “apoptotic process” (GO:0006915), “inflammatory response” (GO:0006954), “cellular response to DNA damage stimulus” (GO:0006974) and “response to oxidative stress” (GO:0006979) available at www.geneontology.org [PMID:[10802651](https://pubmed.ncbi.nlm.nih.gov/10802651/), [30395331](https://pubmed.ncbi.nlm.nih.gov/30395331/)]. 
The genes related to these four GO-terms are retrieved using the biomaRt package available in R (version 2.42.0) [PMID:[19617889](https://pubmed.ncbi.nlm.nih.gov/19617889/), [16082012](https://pubmed.ncbi.nlm.nih.gov/16082012/)]. ORA for comparison as an approach compared to genetich overlap was done similar as described above in part 1 and the WikiPathways database was also used as described above. 
 
* #### Part 3 - Study the effect on toxicity related pathways
The third part combines the results of the previous two parts. This part uses the results of the second part to filter those of the first part. This approach leads to toxicity related pathways which are used for further biological analysis. 

Cite the code: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5039276.svg)](https://doi.org/10.5281/zenodo.5039276)
