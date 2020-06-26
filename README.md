# Investigating the cell-specific toxicity response to titanium dioxide nanobelts
## How to run
This workflow is sperated in three different modules i.e. [module 1](https://github.com/laurent2207/TiO2-scripts/blob/master/Module1_Enrichment.Rmd), [module 2](https://github.com/laurent2207/TiO2-scripts/blob/master/Module2_GOterm_processing.Rmd) and [module 3](https://github.com/laurent2207/TiO2-scripts/blob/master/Module3_FilterResults.Rmd). The necessary data files can be found in the [data subfolder](https://github.com/laurent2207/TiO2-scripts/tree/master/data). They are already pre-processed and the resulting files can be found, after running the pre-processing script, in the [output folder](https://github.com/laurent2207/TiO2-scripts/tree/master/output). 
The functions that are used in the modules are found in the [functions subfolder](https://github.com/laurent2207/TiO2-scripts/tree/master/functions).

## Module descriptions

* #### Module 1 - Find affected pathways
The first module, to find affected pathways, is used to find affected pathways via overrepresentation analysis (ORA) based on gene expression data. Overrepresentation analysis is done to identify significantly enriched pathways from the WikiPathways database. The enricher function of the clusterProfiler package in R (version 3.14.3) [PMID:22455463] is used in this module. Adjusted p-value cut-off and q-value cut-off were set at lower than 0.05. Minimal gene set size was set to 10 and maximal gene set size was set to 300. All other variables were kept standard.
This module makes use of the pathway models database WikiPathways. WikiPathways (www.wikipathways.org) is an open platform for the curation of biological pathways and hosts a pathway database with custom graphical pathway editing tools [PMID:29136241]. 
For this module and the second module, pathway-gene information was obtained in the Gene Matrix Transposed (GMT) file format from http://data.wikipathways.org/, including the Curated Collection and the Reactome pathways .

* #### Module 2 - Find toxicity related pathways
The second module, to find toxicity related pathways, is used to find toxicity related pathways based on Gene Ontology (GO) terms and their respective genes. In this repository, GO-terms were selected based on relation to a toxicologic response in the human body. Genes were obtained from the Gene Ontology (GO) (AmiGO 2 version 2.5.12) categories “apoptotic process” (GO:0006915), “inflammatory response” (GO:0006954), “cellular response to DNA damage stimulus” (GO:0006974) and “cellular response to oxidative stress” (GO:0034599) available at www.geneontology.org [PMID:10802651][PMID:30395331]. Genes are retrieved using the biomaRt package available in R (version 2.42.0) [PMID:19617889, 16082012]. Genes are selected based on organism i.e. Homo sapiens and for all chromosomes. After retrieval genes are filtered based on evidence. Only genes with Gene Ontology evidence ND, NR, NAS, IEA, ISS, ISO, ISA, ISM, IGC and RCA are selected (http://geneontology.org/docs/guide-go-evidence-codes/). ORA was done as described above and WikiPathways was also used as described above. 
 
* #### Module 3 - Study the effect on toxicity related pathways
The third module, to study the effect on toxicity related pathways, combines the results of the previous two modules to come up with toxicity related pathways which are used for further biological analysis. The toxicity related pathways found in the second module were used to filter the affected pathways which were retrieved in the first module.

Cite the code: [![DOI](https://zenodo.org/badge/236458629.svg)](https://zenodo.org/badge/latestdoi/236458629)
