library(Seurat)

sc <- readRDS("results/sc_umap.rds")

# Canonical markers
celltype_markers <- list(
  CD8_T = c("CD3D", "CD8A"),
  Treg = c("FOXP3", "IL2RA"),
  NK = c("NKG7", "GNLY"),
  B_cell = c("MS4A1", "CD79A"),
  Monocyte = c("LYZ", "S100A8"),
  Endothelial = c("PECAM1", "VWF"),
  Oligodendrocyte = c("MBP", "MOG"),
  Cycling = c("MKI67", "TOP2A")
)

DotPlot(sc, features = celltype_markers) +
  RotatedAxis()

# Manual annotation (example)
new.cluster.ids <- c(
  "CD8_T_Cell", "Treg_Cell", "NK_Cell",
  "Monocyte", "Endothelial_Cell",
  "Oligodendrocyte_like", "Cycling_Cell"
)

names(new.cluster.ids) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids)

DimPlot(sc, reduction = "umap",
        label = TRUE, repel = TRUE) +
  ggtitle("UMAP - Annotated Cell Types")

saveRDS(sc, "results/sc_annotated.rds")
