### GSEA analysis function 
### Laurent Winckers, Maastricht University - Department of Bioinformatics BiGCaT
### 2020-03-23

enrichment <- function(wpid2gene,wpid2name, data, comparisons, fccutoff) {
  for (i in comparisons) {
    buffer <- data %>% dplyr::select(starts_with("entrez") | starts_with(i))
    
    # GSEA
    dat <- as.numeric(buffer[,5])
    names(dat) <- as.character(buffer[,1])
    dat <- sort(dat, decreasing = T)
    
    fcdat <- as.numeric(buffer[,3])
    names(fcdat) <- as.character(buffer[,1])
    
    res <- clusterProfiler::GSEA(dat, pvalueCutoff = 3,
                minGSSize = 10, maxGSSize = 300,
                TERM2GENE = wpid2gene,
                TERM2NAME = wpid2name,
                nPerm = 30000)
    write.table(res, file = paste("output/GSEA_", i, ".txt", sep=""), sep="\t", quote = F, row.names = F)
    
    dp <- dotplot(res, showCategory= 30, x = "NES")
    ep <- emapplot(res, showCategory = 30)
    cp <- cnetplot(res, showCategory = 30, node_label="category", foldChange = fcdat)    
    figure <- ggarrange(ggarrange(dp, ep, ncol = 2, labels = c("A","B")), cp, labels = c("C"), nrow = 2)
    ggexport(figure, filename = paste("output/GSEA-",i,".png",sep=""), width = 2000, height = 1800)
    
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
    
    write.table(ewp.up, file = paste("output/ORA_up_", i, ".txt", sep=""), sep="\t", quote = F, row.names = F)
  
    ewp.down <- clusterProfiler::enricher(
      dn.genes,
      pvalueCutoff = 1, 
      qvalueCutoff = 1,
      minGSSize = 10, maxGSSize = 300,
      TERM2GENE = wpid2gene,
      TERM2NAME = wpid2name)
    
    write.table(ewp.down, file = paste("output/ORA_down_", i, ".txt", sep=""), sep="\t", quote = F, row.names = F)
    
    #dp.up <- dotplot(ewp.up, showCategory= 20, color = "p.adjust")
    cp.up <- cnetplot(ewp.up, showCategory=20, node_label="category")  
    #dp.down <- dotplot(ewp.down, showCategory= 20, color = "p.adjust")
    cp.down <- cnetplot(ewp.down, showCategory=20, node_label="category") 
    #figure <- ggarrange(dp.up, cp.up, dp.down, cp.down, labels = c("A", "B","C", "D"), ncol = 2, nrow = 2)
    figure <- ggarrange(cp.up, cp.down, labels = c("A", "B"), nrow = 2)
    ggexport(figure, filename = paste("output/ORA-",i,".png",sep=""), width = 2000, height = 2000)
  }
}
