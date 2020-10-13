# Investigating the cell-specific toxicity response to titanium dioxide nanobelts
## How to run
The script is seperated in four parts i.e. Data pre-processing, Part 1, Part 2 and Part 3. The necessary data files can be found in the [data folder](https://github.com/laurent2207/TiO2-scripts/tree/master/data). The orginal data, retrieved from [GEO, accession:GSE42067](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42067) and analysed with [ArrayAnalysis](https://arrayanalysis.org), can be found in the [original_data subfolder]() of the aforementioned [data folder](https://github.com/laurent2207/TiO2-scripts/tree/master/data). The files used in the [Data pre-processing]() part of the script can be found in the [data_pre_processing subfolder]() of the aforementioned [data folder](https://github.com/laurent2207/TiO2-scripts/tree/master/data). The resulting files can be found in the [output folder](https://github.com/laurent2207/TiO2-scripts/tree/master/output). 
The functions that are used in the script are found in the [functions folder](https://github.com/laurent2207/TiO2-scripts/tree/master/functions).

## Part descriptions

* #### Data pre-processing


* #### Part 1 - Find affected pathways
The first module is used to find affected pathways. This module uses overrepresentation analysis (ORA) based on gene expression data. This method is used to identify significantly enriched pathways from the WikiPathways database. 
In this module the enricher function of the clusterProfiler package in R (version 3.14.3) [PMID:22455463].  Almost all variables were kept standard except the Adjusted p-value cut-off and q-value cut-off were set at lower than 0.05 and minimal gene set and maximal gene set size were set to 10 and 300 respectively.
As mentioned above, this module makes use of the pathway models from the WikiPathways database (www.wikipathways.org). It is an open platform for the curation of biological pathways and hosts a pathway database with custom graphical pathway editing tools [PMID:29136241]. 
For both this module and the second module, pathway-gene information was obtained in the Gene Matrix Transposed (GMT) file format from http://data.wikipathways.org/. The pathway dataset included the Curated Collection and the Reactome pathways.

* #### Part 2 - Find toxicity related pathways
The second module is used to find toxicity related pathways based on Gene Ontology (GO) terms and their respective genes. For this module the GO-terms were selected based on their relation to toxicologic processes i.e. responses in the human body. Genes were obtained from the Gene Ontology (GO) (AmiGO 2 version 2.5.12) categories “apoptotic process” (GO:0006915), “inflammatory response” (GO:0006954), “cellular response to DNA damage stimulus” (GO:0006974) and “cellular response to oxidative stress” (GO:0034599) available at www.geneontology.org [PMID:10802651][PMID:30395331]. 
The genes related to these four GO-terms are retrieved using the biomaRt package available in R (version 2.42.0) [PMID:19617889, 16082012]. ORA was done as described above and the WikiPathways database was also used as described above. 
 
* #### Part 3 - Study the effect on toxicity related pathways
The third module combines the results of the previous two modules. This module uses the results of the second module to filter those of the first module. This approach leads to toxicity related pathways which are used for further biological analysis. 

Cite the code: [![DOI](https://zenodo.org/badge/236458629.svg)](https://zenodo.org/badge/latestdoi/236458629)
