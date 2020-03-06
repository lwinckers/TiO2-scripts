### Network creation in Cytoscape - using RCy3 package
### Laurent Winckers, Maastricht University - BiGCaT
### 14-11-2019

# clean work space
rm(list=ls())
options(stringsAsFactors = F)

# set working directory to where script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# install packages 
library(RCy3)

# check if cytoscape is open and check version
cytoscapePing()
cytoscapeVersionInfo()

# load in data of network
nodes <- as.data.frame(read.table("./data-output/nodes.txt", header = T, sep = "\t", stringsAsFactors = FALSE))
edges <- as.data.frame(read.table("./data-output/edges.txt", header = T, sep = "\t", stringsAsFactors = FALSE))

colnames(nodes)[1] <- c("id")
nodes <- nodes[-2]

colnames(edges)[c(1,2)] <- c("source", "target")
edges <- edges[-3]
edges$target <- as.character(edges$target)

# create network from the data frames and map column for entrezgene IDs
createNetworkFromDataFrames(nodes, edges, title = "MyNetwork", collection = "MyCollection")

# load data set with gene expression values (logFC, p.value)
expr_data <- read.table("./gene-expression-data/data-output/merged_TiO2.txt", header = T, sep ="\t")

# load data into cytoscape
loadTableData(expr_data, data.key.colum = "entrezgene_id", table.key.column = "shared name")

# check if tables are well merged
nodeTable <- getTableColumns(table = "node")
head(nodeTable)

# create heatmaps of data for respective nodes
setNodeCustomHeatMapChart(c("TiO2_24hrs_10_Caco2-logFC","TiO2_24hrs_10_SAE-logFC",
                            "TiO2_24hrs_10_THP1-logFC"), 
                          slot = 1, range = c(-0.26, 0.26))
setNodeCustomPosition(nodeAnchor = "N", graphicAnchor = "S", slot = 1)

setNodeCustomHeatMapChart(c("TiO2_24hrs_100_Caco2-logFC","TiO2_24hrs_100_SAE-logFC",
                            "TiO2_24hrs_100_THP1-logFC"), 
                          slot = 2, range = c(-0.26, 0.26))
setNodeCustomPosition(nodeAnchor = "S", graphicAnchor = "N", slot = 2)

# set node shape and label
lockNodeDimensions(TRUE)
setNodeShapeMapping(table.column = "type", 
                    table.column.values = c("Gene", "Pathway"),
                   shapes =  c("DIAMOND", "ELLIPSE"))
setNodeLabelMapping(table.column = "nodes")

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

process <- createColumnFilter(filter.name = "Pathway filter", column = "type", "Pathway", "IS")

selectNodes(nodes = c(process$nodes, sigexpr$nodes), by.col = "shared name")

createSubnetwork(nodes = "selected", subnetwork.name = "subnetwork-sig-genes")

# create community clusters
#glay.cmd = 'clusterAttribute = __glayCluster, createGroups = NO, network = current, restoreEdges = NO, selectedOnly = NO, showUI = NO, undirectedEdges = NO'
#commandsRun(glay.cmd)

# save networks
setCurrentNetwork(network = "MyNetwork")
fitContent()
png.file <- file.path("./data-output/images/pathway_network.png")
exportImage(png.file, type = "png", resolution=600, zoom=500)  

setCurrentNetwork(network = "subnetwork-sig-genes")
fitContent()
png.file <- file.path("./data-output/images/pathway_sig-genes_network.png")
exportImage(png.file, type = "png", resolution=600, zoom=500)  

# retrieve node table
setCurrentNetwork(network = "MyNetwork")
node_table <- getTableColumns(table = "node")

setCurrentNetwork(network = "subnetwork-sig-genes")
sig_node_table <- getTableColumns(table = "node")

write.table(node_table, "./data-output/network/cyto_node_table.txt", row.names = F, sep = "\t", quote = F)
write.table(sig_node_table, "./data-output/network/cyto_sig_node_table.txt", row.names = F, sep = "\t", quote = F)

# save Cytoscape session
session.file <- file.path("./data-output/network/pathway_network.cys")
saveSession(session.file)

# information about session
sessionInfo()
