### Creation of network in Cytoscape, based on combined genes from GO terms
### Laurent Winckers, Maastricht University - BiGCaT
### 14-11-2019

# clean work space
rm(list=ls())

# set working directory to where script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# install packages 
source("./functions/autoInstallPackages.R")
using("RCy3", "data.table")

# check if cytoscape is open and check version
cytoscapePing()
cytoscapeVersionInfo()

# load GO term network files
edges <- read.table("./data-output/entrezgene/edges_GO_terms.txt", header = T, sep = "\t")
nodes <- read.table("./data-output/entrezgene/nodes_GO_terms.txt", header = T, sep = "\t")

edges <- as.data.frame(lapply(edges, as.character))
nodes <- as.data.frame(lapply(nodes, as.character))

# create network from the data frames and map column for entrezgene IDs
createNetworkFromDataFrames(nodes, edges, title = "GO_network", collection = "MyCollection")

# Add organic layout to the network for better visualization.

# load data set with gene expression values (logFC, p.value)
expr_data <- read.table("./gene-expression-data/data-output/merged_TiO2.txt", header = T, sep ="\t")

# load data into cytoscape
loadTableData(expr_data, data.key.colum = "entrezgene_id", table.key.column = "shared name")

# create heatmaps of data for respective nodes
setNodeCustomHeatMapChart(c("TiO2_24hrs_10_Caco2-logFC","TiO2_24hrs_10_SAE-logFC",
                            "TiO2_24hrs_10_THP1-logFC"), 
                          slot = 1, range = c(-0.26, 0.26))
setNodeCustomPosition(nodeAnchor = "N", graphicAnchor = "S", slot = 1)

setNodeCustomHeatMapChart(c("TiO2_24hrs_100_Caco2-logFC","TiO2_24hrs_100_SAE-logFC",
                            "TiO2_24hrs_100_THP1-logFC"), 
                          slot = 2, range = c(-0.26, 0.26))
setNodeCustomPosition(nodeAnchor = "S", graphicAnchor = "N", slot = 2)

# set node shape
lockNodeDimensions(TRUE)
setNodeShapeMapping(table.column = "type", 
                    table.column.values = c("Gene", "Process"),
                    shapes =  c("DIAMOND", "ELLIPSE"))

# set defaults
setNodeBorderColorDefault(new.color = "#999999")
setNodeBorderWidthDefault(new.width = 7)
setNodeColorDefault(new.color = "#FFFFFF")

# create subnetwork based on only significant genes in at least one of the datasets
createColumnFilter(filter.name = "sig genes1", column = "TiO2_24hrs_10_Caco2-P.Value", 0.05, "LESS_THAN")
createColumnFilter(filter.name = "sig genes2", column = "TiO2_24hrs_100_Caco2-P.Value", 0.05, "LESS_THAN")
createColumnFilter(filter.name = "sig genes3", column = "TiO2_24hrs_10_THP1-P.Value", 0.05, "LESS_THAN")
createColumnFilter(filter.name = "sig genes4", column = "TiO2_24hrs_100_THP1-P.Value", 0.05, "LESS_THAN")
createColumnFilter(filter.name = "sig genes5", column = "TiO2_24hrs_10_SAE-P.Value", 0.05, "LESS_THAN")
createColumnFilter(filter.name = "sig genes6", column = "TiO2_24hrs_100_SAE-P.Value", 0.05, "LESS_THAN")

sigexpr <- createCompositeFilter('combined filter', filter.list = c("sig genes1", "sig genes2",
                                                                    "sig genes3", "sig genes4",
                                                                    "sig genes5", "sig genes6"), 
                                                                      type = "ANY")

process <- createColumnFilter(filter.name = "Pathway filter", column = "type", "Process", "IS")

selectNodes(nodes = c(process$nodes, sigexpr$nodes), by.col = "shared name")

createSubnetwork(nodes = "selected", subnetwork.name = "subnetwork_GO_sig-genes")

# save networks
setCurrentNetwork(network = "GO_network")
fitContent()
png.file <- file.path("./data-output/images/GO_network.png")
exportImage(png.file, type = "png", resolution=600, zoom=500)  

setCurrentNetwork(network = "subnetwork_GO_sig-genes")
fitContent()
png.file <- file.path("./data-output/images./GO_sig-genes_network.png")
exportImage(png.file, type = "png", resolution=600, zoom=500)  

# save session
session.file <- file.path("./data-output/network/GO_network.cys")
saveSession(session.file)

# information about session
sessionInfo()
