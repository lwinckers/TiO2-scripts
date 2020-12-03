### ORA analysis function 
### Laurent Winckers, Maastricht University - Department of Bioinformatics BiGCaT
### 2020-03-23

enrichment <- function(wpid2gene,wpid2name, data, comparisons, fccutoff) {
  for (i in comparisons) {
    buffer <- data %>% dplyr::select(starts_with("entrez") | starts_with(i))
    
    # ORA
     up.genes <- buffer[buffer[,2] > fccutoff & buffer[,4] < 0.05, 1] 
     dn.genes <- buffer[buffer[,2] < -fccutoff & buffer[,4] < 0.05, 1]

     ewp.up <- clusterProfiler::enricher(
      up.genes,
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 10, maxGSSize = 300,
      TERM2GENE = wpid2gene,
      TERM2NAME = wpid2name)
    
    write.table(ewp.up, file = paste("output/ORApw_up_", i, ".txt", sep=""), sep="\t", quote = F, row.names = F)
  
    ewp.down <- clusterProfiler::enricher(
      dn.genes,
      pvalueCutoff = 1, 
      qvalueCutoff = 1,
      minGSSize = 10, maxGSSize = 300,
      TERM2GENE = wpid2gene,
      TERM2NAME = wpid2name)
    
    write.table(ewp.down, file = paste("output/ORApw_down_", i, ".txt", sep=""), sep="\t", quote = F, row.names = F)
    
    #dp.up <- dotplot(ewp.up, showCategory= 20, color = "p.adjust")
    cp.up <- cnetplot(ewp.up, showCategory=20, node_label="category")  
    #dp.down <- dotplot(ewp.down, showCategory= 20, color = "p.adjust")
    cp.down <- cnetplot(ewp.down, showCategory=20, node_label="category") 
    #figure <- ggarrange(dp.up, cp.up, dp.down, cp.down, labels = c("A", "B","C", "D"), ncol = 2, nrow = 2)
    figure <- ggarrange(cp.up, cp.down, labels = c("A", "B"), nrow = 2)
    
    # png
    #ggexport(figure, filename = paste("output/ORA-",i,".png",sep=""), width = 2000, height = 2000)
    
    # svg
    ggsave(file=paste("output/ORA-",i,".svg",sep=""), plot=figure, width=18, height=18, limitsize = FALSE)
  }
}
