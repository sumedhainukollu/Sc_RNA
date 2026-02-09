library(Seurat)
library(dplyr)

sc <- readRDS("results/sc_umap.rds")

markers <- FindAllMarkers(sc,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

# Top markers per cluster
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

write.csv(markers, "results/markers/all_cluster_markers.csv",
          row.names = FALSE)

DoHeatmap(sc, features = top_markers$gene) +
  ggtitle("Top Marker Genes Heatmap")
