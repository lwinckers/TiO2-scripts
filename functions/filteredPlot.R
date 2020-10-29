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
  #V(net)$color <- ifelse(V(net)$type == "Pathway", "orange", "lightblue")
  
  # to create colors for each GO-term the pathway is affiliated with we need to create a matrix which depicts this affiliation
  # additionally we need to point out which ones are genes
  m <- sigORA
  m <- m[m$ID %in% data$ID,]
  row.names(m) <- m[,1]
  m <- m[-1]
  m$pid5 <- 0
  m <- as.data.frame(sapply(m, as.integer))
  
  k <- matrix(0, ncol = 4, nrow = (nrow(nodes)-nrow(data)))
  k <- data.frame(k)
  k$pid5 <- 1
  colnames(k) <- c("pid1", "pid2", "pid3" , "pid4", "pid5")
  
  mm <- rbind(m,k)
  am <- data.matrix(mm)
  
  values <- lapply(seq_len(nrow(am)), function(i) am[i,])
  
  # create colorblind friendly color palette
  pal <- c("#FFC20A", "#006CD1", "#D35FB7", "#D41159", "#D3D3D3")
  
  # increase font size
  V(net)$label.cex = 1.5
  # set label color to black
  V(net)$label.color = "#000000"
  
  # SVG
  svg(paste0("output/filteredPlot_",cell,"_",x,"_",conc,".svg"), width = 20, height = 20)
  plot(net, vertex.shape="pie", vertex.pie=values, vertex.pie.color=list(pal),
       edge.width = 5, vertex.size = ifelse(nodes$type == "Pathway", 6, 3), rescale = T, 
       vertex.label = ifelse(nodes$type == "Pathway", nodes$ID, NA), 
       vertex.frame.color="lightgray",
       xlim = c(-1.2,1.2), asp = 0)
  dev.off()
  
  # png
#  png(paste0("output/filteredPlot_",cell,"_",x,"_",conc,".png"), width = 2000, height = 2000)
#  plot(net, vertex.shape="pie", vertex.pie=values, vertex.pie.color=list(pal),
#       edge.width = 5, vertex.size = ifelse(nodes$type == "Pathway", 6, 3), rescale = T, 
#       vertex.label = ifelse(nodes$type == "Pathway", nodes$ID, NA), 
#       vertex.frame.color="lightgray",
#       xlim = c(-1.2,1.2), asp = 0)
#  dev.off()
}