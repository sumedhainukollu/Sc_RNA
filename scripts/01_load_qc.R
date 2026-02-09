library(Seurat)
library(ggplot2)
library(patchwork)

# Load 10X data
data_dir <- "data/raw/Glioma_10X/"
counts <- Read10X(data.dir = data_dir)

sc <- CreateSeuratObject(
  counts = counts,
  project = "Glioma_scRNA",
  min.cells = 3,
  min.features = 200
)

# Mitochondrial percentage
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# QC plots
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtering
sc <- subset(sc,
             subset = nFeature_RNA > 200 &
                      nFeature_RNA < 6000 &
                      percent.mt < 10)

saveRDS(sc, "results/sc_qc_filtered.rds")
