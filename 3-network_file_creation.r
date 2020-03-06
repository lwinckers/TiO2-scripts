### Creation of files used for network creation in Cytoscape
### Laurent Winckers, Maastricht University - BiGCaT
### 14-11-2019

### qusage 2.18.0
### plyr 1.8.4
### biomaRt 2.40.5
### data.table 1.12.2
### limma 3.40.6

# clean work space
rm(list=ls())
options(stringsAsFactors = F)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# install packages 
library(qusage)
library(plyr)
library(biomaRt)
library(data.table)

##### GO TERM NETWORK #####

# load in GO term files and clean the files to make them according edge tables
path = paste0(getwd(), "/data-output")
ls_GO <- list.files(path = path, pattern = "^ann_GO")
for (i in 1:length(ls_GO)){
  assign(ls_GO[i], read.table(paste0(path,"/",ls_GO[i]), header = T, sep = "\t")[c("GO.TERM","entrezgene_id")])
  }
ls_GO <- mget(ls(pattern = "ann_GO"))
colnames <- c("source", "target") 
ls_GO <- lapply(ls_GO, setNames, colnames)

# combine edge tables of individual GO terms and make them overall character
names(ls_GO) <- NULL
edges_GO_terms <- do.call("rbind",ls_GO)
edges_GO_terms[] <- lapply(edges_GO_terms, as.character)

# create nodes file with type-annotation
nodes_GO_terms <- unique(stack(edges_GO_terms))[-2]
colnames(nodes_GO_terms) <- "id"
nodes_GO_terms$type <- "Gene"
nodes_GO_terms$type[nodes_GO_terms$id %like% "^GO:.*"] <- "Process"

# save files
write.table(edges_GO_terms, "./data-output/edges_GO_terms.txt", row.names = F, quote = F, sep = "\t")
write.table(nodes_GO_terms, "./data-output/nodes_GO_terms.txt", row.names = F, quote = F, sep = "\t")

##### PATHWAY NETWORK #####

# clean work space
rm(list=ls())
options(stringsAsFactors = F)

# load in pathway database gmt files
path = paste0(getwd(), "/data-input")
gmtFile <- list.files(path = path, pattern = ".gmt")
for (i in 1:length(gmtFile)){
  assign(gmtFile[i], read.gmt(paste0(path,"/",gmtFile[i])))
  }
gmtFile <- mget(ls(pattern = ".gmt"))
list2env(lapply(gmtFile, function(x){ldply(x,data.frame)}), envir = .GlobalEnv)

databases <- do.call(rbind, lapply(ls(pattern = ".gmt"), get))
colnames(databases) <- c("pathway", "entrezgene")
databases[] <- lapply(databases, gsub, pattern='[[:punct:]]', replacement=' ')

# load in genes file
genes <- read.table("./data-output/ann_genes.txt", sep = "\t", header = T)

# load selected pathways
selectedpws <- as.data.frame(read.table("./data-output/sigselpws.txt", header = T, sep = "\t"))
colnames(selectedpws)[1] <- "ID"

# load in selected pathways
#selectedpws1 <- as.data.frame(read.table("./data-output/results_ann_GO0006915.txt", header = T, sep = "\t")[,1])
#selectedpws2 <- as.data.frame(read.table("./data-output/results_ann_GO0006954.txt", header = T, sep = "\t")[,1])
#selectedpws3 <- as.data.frame(read.table("./data-output/results_ann_GO0006974.txt", header = T, sep = "\t")[,1])
#selectedpws4 <- as.data.frame(read.table("./data-output/results_ann_GO0034599.txt", header = T, sep = "\t")[,1])

#colnames(selectedpws1) <- "ID"
#colnames(selectedpws2) <- "ID"
#colnames(selectedpws3) <- "ID"
#colnames(selectedpws4) <- "ID"

#path = paste0(getwd(), "/data-output")
#res_GO <- list.files(path = path, pattern = "^results_ann_GO")
#for (i in 1:length(res_GO)){assign(res_GO[i], read.table(paste0(path,"/",res_GO[i]), header = T, sep = "\t"))}
#res_GO <- mget(ls(pattern = "^results_ann_GO"))

#colnames <- c("ID")
#list2env(lapply(res_GO, setNames, colnames), .GlobalEnv)

# edge table
edge_table <- as.data.frame(databases[databases$pathway %in% selectedpws$ID,])

#edges1 <- as.data.frame(databases[databases$pathway %in% selectedpws1$ID,])
#edges2 <- as.data.frame(databases[databases$pathway %in% selectedpws2$ID,])
#edges3 <- as.data.frame(databases[databases$pathway %in% selectedpws3$ID,])
#edges4 <- as.data.frame(databases[databases$pathway %in% selectedpws4$ID,])

#ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

#genes1 <- getBM(
#  attributes = c('hgnc_symbol', 'entrezgene_id'), 
#  filters = 'entrezgene_id',
#  values = edges1$entrezgene,
#  mart = ensembl
#)
#edges1$hgnc_symbol <- genes1$hgnc_symbol[match(edges1$entrezgene, genes1$entrezgene_id)]
#edges1 <- edges1[!is.na(edges1$hgnc_symbol),]
#edges1 <- edges1[-(edges1$hgnc_symbol == ""),]

#genes2 <- getBM(
#  attributes = c('hgnc_symbol', 'entrezgene_id'), 
#  filters = 'entrezgene_id',
#  values = edges2$entrezgene,
#  mart = ensembl
#)
#edges2$hgnc_symbol <- genes2$hgnc_symbol[match(edges2$entrezgene, genes2$entrezgene_id)]
#edges2 <- edges2[!is.na(edges2$hgnc_symbol),]
#edges2 <- edges2[-(edges2$hgnc_symbol == ""),]

#genes3 <- getBM(
#  attributes = c('hgnc_symbol', 'entrezgene_id'), 
#  filters = 'entrezgene_id',
#  values = edges3$entrezgene,
#  mart = ensembl
#)
#edges3$hgnc_symbol <- genes3$hgnc_symbol[match(edges3$entrezgene, genes3$entrezgene_id)]
#edges3 <- edges3[!is.na(edges3$hgnc_symbol),]
#edges3 <- edges3[-(edges3$hgnc_symbol == ""),]

#genes4 <- getBM(
#  attributes = c('hgnc_symbol', 'entrezgene_id'), 
#  filters = 'entrezgene_id',
#  values = edges4$entrezgene,
#  mart = ensembl
#)
#edges4$hgnc_symbol <- genes4$hgnc_symbol[match(edges4$entrezgene, genes4$entrezgene_id)]
#edges4 <- edges4[!is.na(edges4$hgnc_symbol),]
#edges4 <- edges4[-(edges4$hgnc_symbol == ""),]

#edges1 <- unique(edges1)
#edges2 <- unique(edges2)
#edges3 <- unique(edges3)
#edges4 <- unique(edges4)

#write.table(edges1 , "./data-output/edges-GO0006915-pw.txt", row.names=F, quote=F, sep="\t")
#write.table(edges2 , "./data-output/edges-GO0006954-pw.txt", row.names=F, quote=F, sep="\t")
#write.table(edges3 , "./data-output/edges-GO0006974-pw.txt", row.names=F, quote=F, sep="\t")
#write.table(edges4 , "./data-output/edges-GO0034599-pw.txt", row.names=F, quote=F, sep="\t")

# only unique rows
edge_table <- unique(edge_table)

# map entrezgene IDs to hgnc symbols
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

genes <- getBM(
  attributes = c('hgnc_symbol', 'entrezgene_id'), 
  filters = 'entrezgene_id',
  values = edge_table$entrezgene,
  mart = ensembl
)

# merge edge table with hgnc symbols
# remove NAs and empty values as they are pseudogenes, microRNAs or discontinued genes
edge_table$hgnc_symbol <- genes$hgnc_symbol[match(edge_table$entrezgene, genes$entrezgene)]
edge_table <- edge_table[!is.na(edge_table$hgnc_symbol),]
edge_table <- edge_table[-(edge_table$hgnc_symbol == ""),]

edge_entrez <- edge_table[-3]
edge_hgnc <- edge_table[-2]

colnames(edge_entrez)[c(1,2)] <- c("source", "target")
colnames(edge_hgnc)[c(1,2)] <- c("source", "target")

# save edge table
write.table(edge_table, "./data-output/sigseledges.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(edge_entrez, "./data-output/sigseledges_entrezgene.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(edge_hgnc, "./data-output/sigseledges_hgnc.txt", col.names = T, row.names = F, sep = "\t", quote = F)

# node table
# split edge table in two seperate data frames, rbind, unique, add type.
pathwayNodes <- as.data.frame(unique(edge_table$pathway))
colnames(pathwayNodes)[1] <- "nodes"

geneNodes_entrezgene <- as.data.frame(unique(edge_table$entrezgene))
colnames(geneNodes_entrezgene)[1] <- "nodes"

geneNodes_hgnc <- as.data.frame(unique(edge_table$hgnc_symbol))
colnames(geneNodes_hgnc)[1] <- "nodes"

nodes_entrezgene <- as.data.frame(rbind(pathwayNodes, geneNodes_entrezgene))
nodes_hgnc <- as.data.frame(rbind(pathwayNodes, geneNodes_hgnc))

both <- unique(edge_table[-1])
colnames(both) <- c("nodes", "nodes2")
pathwayNodes2 <- pathwayNodes
pathwayNodes2$nodes2 <- pathwayNodes$nodes
nodes <- as.data.frame(rbind(pathwayNodes2, both))

# add type
nodes$type <- "Gene"
nodes$type[nodes$nodes2 %in% genes$hgnc_symbol] <- "SelGene"
nodes$type[nodes$nodes %in% edge_table$pathway] <- "Pathway"

nodes_entrezgene$type <- "Gene"
nodes_entrezgene$type[nodes_entrezgene$nodes %in% genes$entrezgene] <- "SelGene"
nodes_entrezgene$type[nodes_entrezgene$nodes %in% edge_table$pathway] <- "Pathway"

nodes_hgnc$type <- "Gene"
nodes_hgnc$type[nodes_hgnc$nodes %in% genes$hgnc_symbol] <- "SelGene"
nodes_hgnc$type[nodes_hgnc$nodes %in% edge_table$pathway] <- "Pathway"

# save node table
write.table(nodes, "./data-output/sigselnodes.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(nodes_entrezgene, "./data-output/sigselnodes_entrezgene.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(nodes_hgnc, "./data-output/sigselnodes_hgnc.txt", col.names = T, row.names = F, sep = "\t", quote = F)

# information about this session
sessionInfo()


##### SIGNIFICANT PATHWAYS #####
sigPWs <- read.table("./data-output/GSEA/sig_pws.txt", header = T, sep = "\t")
sigPWs$unique.sigPWs.ID. <- as.character(sigPWs$unique.sigPWs.ID.)
edge_table <- as.data.frame(databases[databases$pathway %in% sigPWs[,1],])

write.table(edge_table, "./data-output/sig_edges.txt", col.names = T, row.names = F, sep = "\t", quote = F)
