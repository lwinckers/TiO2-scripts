### Heatmap creation function
### Laurent Winckers, Maastricht University - Department of Bioinformatics BiGCaT
### 2020-04-01

heatmap <- function(fileName, data){
  rownames(data) <- data[,1]
  data <- data[-1]
  labels_row <- rownames(data)
  
  pheatmap(data, cluster_cols = F, cluster_rows = T, 
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
         filename = paste0("./output/images/", fileName, "_heatmap.pdf"))
  }