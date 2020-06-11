filteredPlot <- function(x, cell, conc) {
  fileName = paste("ORA", x, cell, conc, sep = "_")
  data <- read.table(paste0("output/", fileName, ".txt"), header = T, sep = "\t", quote = "")
  data <- data[data$pvalue < 0.05,]
  
  data <- data[data$ID %in% sigORA$ID,]
  write.table(data, paste0("output/", fileName, "_filtered.txt"), quote = F, sep = "\t", row.names = F)
  
  s <- strsplit(data$geneID, split = "/")
  edges <- data.frame(WP = rep(data$ID, sapply(s, length)), Gene = unlist(s))
  nodes <- unique(data.frame(ID = c(edges[,1], edges[,2])))
  nodes$type <- ifelse(nodes$ID %like% "^WP" , 'Pathway', 'Gene')
  
  net <- graph_from_data_frame(d = edges, vertices = nodes, directed = F) 
  V(net)$color <- ifelse(V(net)$type == "Pathway", "orange", "lightblue")
  
  png(paste0("output/filteredPlot_",cell,"_",x,"_",conc,".png"), width = 2000, height = 2000)
  plot(net, edge.width = 5, vertex.size = 6, rescale = T, vertex.label = NA, xlim = c(-1.2,1.2), asp = 0)
  dev.off()
}