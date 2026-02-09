library(Seurat)

sc <- readRDS("results/sc_pca.rds")

# Neighbors & clustering
sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc, resolution = 0.5)

# UMAP
sc <- RunUMAP(sc, dims = 1:10)

DimPlot(sc, reduction = "umap", label = TRUE) +
  ggtitle("UMAP - scRNA clusters")

saveRDS(sc, "results/sc_umap.rds")
