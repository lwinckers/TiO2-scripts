 ### Creation of heatmap figure
### Laurent Winckers, Maastricht University - BiGCaT
### 14-11-2019

# clean work space
rm(list=ls())

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# install packages 
source("./functions/autoInstallPackages.R")
using("pheatmap", "colorRamps", "RColorBrewer")

# load necessary files
for (i in 1:length(list.files(path="./data-output/GSEA",pattern="GSEA"))){ 
          assign(list.files(path="./data-output/GSEA",pattern="GSEA")[i], 
          read.table(paste0(getwd(), "/data-output/GSEA/", 
                            list.files(path="./data-output/GSEA",pattern="GSEA")[i]), 
                     header = T, sep = "\t"))
}

##### PREPARATION DATA #####

# remove unnecessary columns
# change column names
gseaList <- mget(ls(pattern = 'GO0034599'))
gseaList <- lapply(gseaList, function(x){x[c(1,5,6)]}) 
for (i in 1:length(gseaList)){
  colnames(gseaList[[i]]) <- c("ID", paste0("NES_",i), paste0("pvalue_",i))
}
list2env(gseaList, envir = .GlobalEnv)

# merge data frames together
merged_GSEA <- Reduce(function(x, y) merge(x, y, all=TRUE), gseaList)
rownames(merged_GSEA) <- merged_GSEA[,1]
merged_GSEA <- merged_GSEA[-1]

##### CREATION HEATMAP #####

# load Cytoscape node table file and clean the file to create cluster file
#clusters <- read.table("./data-output/network/cyto_node_table.txt", 
#                         header = T, sep = "\t")
#clusters <- clusters[c(3,7,23)]
#clusters <- clusters[clusters$type == "Pathway",c(1,3)]
#colnames(clusters) <- c("ID", "cluster")

#GO_clusters <- data.frame(ID = c("GO:0006954", "GO:0006915", "GO:0006974", "GO:0034599"),
#                          cluster = 0)
#clusters <- rbind(clusters, GO_clusters)

#clusters <- clusters[order(clusters$cluster),]

# save file before setting first column as row names
#write.table(clusters, "./data-output/clusters.txt", row.names = F, sep = "\t", quote = F)

# set first column as row names
#rownames(clusters) <- clusters[,1]
#clusters <- clusters[-1]

# merge clusters with NES scores data frame, and order the rows based on the clusters
merged_GSEA <- merged_GSEA[-(grep("pvalue", names(merged_GSEA)))]
#merged_GSEA <- merge(merged_GSEA, clusters, by = "row.names")
#merged_GSEA$cluster <- as.numeric(as.character(merged_GSEA$cluster))
#merged_GSEA <- merged_GSEA[order(merged_GSEA$cluster),]
# <- merged_GSEA[-8]

# clean row names for heatmap
#merged_GSEA$ID <- gsub("%.*", "", merged_GSEA$ID)
#merged_GSEA$ID <- gsub("_", " ", merged_GSEA$ID)
#merged_GSEA$ID <- gsub("REACTOME ", "", merged_GSEA$ID)
#merged_GSEA$ID <- gsub("KEGG ", "", merged_GSEA$ID)
#merged_GSEA$ID <- toupper(merged_GSEA$ID)
#merged_GSEA$ID <- gsub(" CONTROL OF IMMUNE TOLERANCE BY VASOACTIVE INTESTINAL PEPTIDE",
#                       "CONTROL OF IMMUNE TOLERANCE BY VASOACTIVE INTESTINAL PEPTIDE",
#                       merged_GSEA$ID)
#
#rownames(merged_GSEA) <- merged_GSEA[,1]
#merged_GSEA <- merged_GSEA[-1]

# change order of columns
#merged <- merged[c(1,3,5,2,4,6)]

# create color palette
#newCols <- colorRampPalette(grDevices::rainbow(length(unique(clusters$cluster))))
#mycolors <- newCols(length(unique(clusters$cluster)))
#names(mycolors) <- unique(clusters$cluster)
#mycolors <- list(cluster = mycolors)

# create labels for rows of heatmap
labels_row <- rownames(merged_GSEA)

# create heatmap
pheatmap(merged_GSEA, cluster_cols = F, cluster_rows = F, 
         #color = col(10000),
         fontsize_row = 5, na_col = "#DDDDDD", scale = "column",
         #legend_breaks = c(),
         #legend_labels = c(),
         labels_row = labels_row,
         labels_col = c("Caco2 10µg/ml",
                        "Caco2 100µg/ml",
                        "SAE 10µg/ml",
                        "SAE 100µg/ml",
                        "THP1 10µg/ml",
                        "THP1 100µg/ml"),
         angle_col = 90,
         gaps_col = c(2,4),
         fontsize_col = 8,
         #annotation_row = clusters,
         #annotation_colors = mycolors,
         annotation_names_row = F,
         annotation_legend = T,
         #gaps_row = c(),
         cellheight = 10, cellwidth = 20,
         filename = "./data-output/images/GO0034599_heatmap.pdf")

# information about session
sessionInfo()
